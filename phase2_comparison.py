"""
Phase 2: Three-tier matching of v1_EquCab3_Master against v2.
"""
import logging
import pandas as pd
from normalization import normalize_address_id, normalize_chr, make_pos_key, classify_variant
from config import MAX_COORD_SHIFT_BP

log = logging.getLogger(__name__)


def compare_v1_v2(v1_master: pd.DataFrame, df_v2: pd.DataFrame) -> pd.DataFrame:
    """Perform three-tier matching and return the final master comparison table.

    Tiers:
        1. strict     — same IlmnID + Name + EquCab3 position
        2. id_shift   — same Name + same position, different IlmnID
        3. positional — same Chr:Pos, different IDs and Names
        4. v1_only    — in v1 but not matched in v2 (control probes)
        5. v2_only    — new markers in v2

    Returns:
        DataFrame with all markers from both manifests, with match_tier and
        in_v1 / in_v2 boolean columns.
    """
    log.info("=== Phase 2: Three-Tier Matching ===")

    # --- Prepare v2 ---
    df_v2 = df_v2.copy()
    for col in ["AddressA_ID", "AddressB_ID"]:
        if col in df_v2.columns:
            df_v2[col] = normalize_address_id(df_v2[col])
    df_v2["variant_type"] = classify_variant(df_v2["SNP"])
    df_v2["chr_v2_norm"] = normalize_chr(df_v2["Chr"])
    df_v2["pos_key_v2"] = make_pos_key(df_v2["chr_v2_norm"], df_v2["MapInfo"])

    v1_ids = set(v1_master["IlmnID"])
    v2_ids = set(df_v2["IlmnID"])
    log.info("v1 unique IlmnIDs: %d | v2 unique IlmnIDs: %d", len(v1_ids), len(v2_ids))
    log.info("IlmnIDs in both: %d | v1-only: %d | v2-only: %d",
             len(v1_ids & v2_ids), len(v1_ids - v2_ids), len(v2_ids - v1_ids))

    v1_names = set(v1_master["Name"])
    v2_names = set(df_v2["Name"])
    log.info("Names in both: %d | v1-only: %d | v2-only: %d",
             len(v1_names & v2_names), len(v1_names - v2_names), len(v2_names - v1_names))

    # --- Build match_tier map {IlmnID → tier} ---
    match_tier_v1 = {}   # IlmnID from v1 perspective
    match_tier_v2 = {}   # IlmnID from v2 perspective

    # --- Tier 1: Strict Match ---
    tier1 = pd.merge(
        v1_master[["IlmnID", "Name", "pos_key_ec3"]],
        df_v2[["IlmnID", "Name", "pos_key_v2"]],
        on=["IlmnID", "Name"],
        how="inner",
    )
    if MAX_COORD_SHIFT_BP == 0:
        tier1_matched = tier1[tier1["pos_key_ec3"] == tier1["pos_key_v2"]]
    else:
        # Parse positions and apply numeric tolerance
        tier1["pos_v1"] = tier1["pos_key_ec3"].str.split(":").str[1].astype(float)
        tier1["pos_v2"] = tier1["pos_key_v2"].str.split(":").str[1].astype(float)
        tier1_matched = tier1[
            (tier1["pos_key_ec3"].str.split(":").str[0] == tier1["pos_key_v2"].str.split(":").str[0])
            & ((tier1["pos_v1"] - tier1["pos_v2"]).abs() <= MAX_COORD_SHIFT_BP)
        ]

    strict_ids = set(tier1_matched["IlmnID"])
    for iid in strict_ids:
        match_tier_v1[iid] = "strict"
        match_tier_v2[iid] = "strict"
    log.info("Tier 1 (Strict): %d markers", len(strict_ids))

    # --- Tier 2: ID Shift (same Name + same position, different IlmnID) ---
    v1_unmatched = v1_master[~v1_master["IlmnID"].isin(strict_ids)]
    v2_unmatched = df_v2[~df_v2["IlmnID"].isin(strict_ids)]

    tier2 = pd.merge(
        v1_unmatched[["IlmnID", "Name", "pos_key_ec3"]].rename(columns={"IlmnID": "IlmnID_v1"}),
        v2_unmatched[["IlmnID", "Name", "pos_key_v2"]].rename(columns={"IlmnID": "IlmnID_v2"}),
        on="Name",
        how="inner",
    )
    tier2 = tier2[
        (tier2["IlmnID_v1"] != tier2["IlmnID_v2"])
        & (tier2["pos_key_ec3"] == tier2["pos_key_v2"])
    ].copy()
    # Store the paired IlmnID for later reference
    for _, row in tier2.iterrows():
        match_tier_v1[row["IlmnID_v1"]] = "id_shift"
        match_tier_v2[row["IlmnID_v2"]] = "id_shift"
    log.info("Tier 2 (ID Shift): %d pairs", len(tier2))

    # --- Tier 3: Positional Only ---
    id_shift_v1_names = set(tier2["Name"])
    v1_still = v1_unmatched[~v1_unmatched["Name"].isin(id_shift_v1_names)]
    v2_still = v2_unmatched[~v2_unmatched["Name"].isin(id_shift_v1_names)]

    # Exclude unplaced markers (Un:-1) from positional matching
    v1_placed = v1_still[v1_still["pos_key_ec3"] != "Un:-1"]
    v2_placed = v2_still[v2_still["pos_key_v2"] != "Un:-1"]

    tier3 = pd.merge(
        v1_placed[["IlmnID", "Name", "pos_key_ec3"]].rename(columns={"IlmnID": "IlmnID_v1", "Name": "Name_v1"}),
        v2_placed[["IlmnID", "Name", "pos_key_v2"]].rename(columns={"IlmnID": "IlmnID_v2", "Name": "Name_v2"}),
        left_on="pos_key_ec3",
        right_on="pos_key_v2",
        how="inner",
    )

    # Flag multi-match (many-to-many Cartesian products at same position)
    tier3["multi_match_v1"] = tier3.duplicated(subset=["IlmnID_v1"], keep=False)
    tier3["multi_match_v2"] = tier3.duplicated(subset=["IlmnID_v2"], keep=False)
    tier3["multi_match"] = tier3["multi_match_v1"] | tier3["multi_match_v2"]

    for _, row in tier3.iterrows():
        match_tier_v1[row["IlmnID_v1"]] = "positional"
        match_tier_v2[row["IlmnID_v2"]] = "positional"
    n_v1_pos = tier3["IlmnID_v1"].nunique()
    n_v2_pos = tier3["IlmnID_v2"].nunique()
    n_unique_pos = len(set(tier3["IlmnID_v1"]) | set(tier3["IlmnID_v2"]))
    log.info("Tier 3 (Positional): %d pairs → %d unique IlmnIDs (%d v1-side, %d v2-side); multi-match pairs: %d",
             len(tier3), n_unique_pos, n_v1_pos, n_v2_pos, tier3["multi_match"].sum())

    # --- Build outer join master table ---
    log.info("Building outer-join master comparison table...")

    # Add suffixes to avoid column collisions
    v1_for_merge = v1_master.copy()
    v2_for_merge = df_v2.copy()

    # Rename overlapping columns that aren't join keys
    v1_rename = {c: f"{c}_v1" for c in v1_for_merge.columns
                 if c in v2_for_merge.columns and c != "IlmnID"}
    v2_rename = {c: f"{c}_v2" for c in v2_for_merge.columns
                 if c in v1_for_merge.columns and c != "IlmnID"}
    v1_for_merge = v1_for_merge.rename(columns=v1_rename)
    v2_for_merge = v2_for_merge.rename(columns=v2_rename)

    master = pd.merge(v1_for_merge, v2_for_merge, on="IlmnID", how="outer")

    # --- Boolean membership flags ---
    master["in_v1"] = master["IlmnID"].isin(v1_ids)
    master["in_v2"] = master["IlmnID"].isin(v2_ids)

    # --- Assign match_tier ---
    def _assign_tier(row):
        iid = row["IlmnID"]
        # If the same IlmnID exists in both manifests, it can only be strict or coord_mismatch.
        # A Tier 3 cross-match may have added this IlmnID to match_tier_v1/v2 as "positional"
        # (because its v1 position coincidentally matched a *different* v2 IlmnID), but that
        # cross-match should not override the IlmnID-level classification.
        if row["in_v1"] and row["in_v2"]:
            t = match_tier_v1.get(iid) or match_tier_v2.get(iid)
            return t if t == "strict" else "coord_mismatch"
        if iid in match_tier_v1:
            return match_tier_v1[iid]
        if iid in match_tier_v2:
            return match_tier_v2[iid]
        if row["in_v1"] and not row["in_v2"]:
            return "v1_only"
        if row["in_v2"] and not row["in_v1"]:
            return "v2_only"
        return "unmatched"

    master["match_tier"] = master.apply(_assign_tier, axis=1)

    # Log how Tier 3 cross-match candidates were finally classified
    _tier3_ids = (
        {iid for iid, t in match_tier_v1.items() if t == "positional"} |
        {iid for iid, t in match_tier_v2.items() if t == "positional"}
    )
    _tier3_master = master.loc[master["IlmnID"].isin(_tier3_ids), "match_tier"]
    _n_pos = int((_tier3_master == "positional").sum())
    _n_coord = int((_tier3_master == "coord_mismatch").sum())
    log.info("Tier 3 final classification: %d positional (v2-only IlmnIDs at shared Chr:Pos), "
             "%d reclassified to coord_mismatch (same IlmnID exists in both manifests)",
             _n_pos, _n_coord)

    # --- Attach multi_match flag ---
    multi_ids = set(tier3.loc[tier3["multi_match"], "IlmnID_v1"]) | \
                set(tier3.loc[tier3["multi_match"], "IlmnID_v2"])
    master["multi_match"] = master["IlmnID"].isin(multi_ids)

    # --- Coordinate match flag (for common markers) ---
    pos_key_v1_col = "pos_key_ec3" if "pos_key_ec3" in master.columns else None
    pos_key_v2_col = "pos_key_v2" if "pos_key_v2" in master.columns else None
    if pos_key_v1_col and pos_key_v2_col:
        master["coord_match"] = (
            master[pos_key_v1_col].notna()
            & master[pos_key_v2_col].notna()
            & (master[pos_key_v1_col] == master[pos_key_v2_col])
        )
    else:
        master["coord_match"] = False

    # --- Summary ---
    tier_counts = master["match_tier"].value_counts()
    log.info("Match tier breakdown:\n%s", tier_counts.to_string())

    return master, tier2, tier3
