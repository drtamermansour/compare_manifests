"""
Phase 1: Integrate v1 manifest with v1.remapped to produce v1_EquCab3_Master.
"""
import logging
import pandas as pd
from normalization import (
    normalize_address_id, normalize_chr, make_pos_key,
    flag_ncbi_scaffolds, classify_variant,
)
from config import EQUCAB3_BUILD, MAPQ_MIN_PROBE, MAPQ_MIN_GENOMIC

log = logging.getLogger(__name__)


def build_v1_master(df_v1: pd.DataFrame, df_rem: pd.DataFrame) -> pd.DataFrame:
    """Join v1 + v1.remapped and apply normalizations and quality flags.

    Returns the v1_EquCab3_Master dataframe.
    """
    log.info("=== Phase 1: Building v1_EquCab3_Master ===")

    # --- Normalize Address IDs ---
    for col in ["AddressA_ID", "AddressB_ID"]:
        if col in df_v1.columns:
            df_v1[col] = normalize_address_id(df_v1[col])
        if col in df_rem.columns:
            df_rem[col] = normalize_address_id(df_rem[col])

    # --- Classify variant types in v1 ---
    df_v1["variant_type"] = classify_variant(df_v1["SNP"])

    # --- Deduplicate join keys (defensive) ---
    dup_v1 = df_v1.duplicated(subset=["IlmnID", "Name"])
    if dup_v1.any():
        log.warning("v1: %d duplicate IlmnID+Name rows — keeping first", dup_v1.sum())
        df_v1 = df_v1.drop_duplicates(subset=["IlmnID", "Name"], keep="first")

    dup_rem = df_rem.duplicated(subset=["IlmnID", "Name"])
    if dup_rem.any():
        log.warning("v1.remapped: %d duplicate IlmnID+Name rows — keeping first", dup_rem.sum())
        df_rem = df_rem.drop_duplicates(subset=["IlmnID", "Name"], keep="first")

    # --- Select remapped columns to bring in ---
    remapped_cols = [
        "IlmnID", "Name",
        "Chr_equCab3", "MapInfo_equCab3", "Strand_equCab3",
        "Ref_equCab3", "Alt_equCab3",
        "MAPQ_TopGenomicSeq", "MAPQ_Probe",
    ]
    # Only keep columns that actually exist in the remapped file
    remapped_cols = [c for c in remapped_cols if c in df_rem.columns]

    # --- Left join v1 + v1.remapped ---
    try:
        master = pd.merge(
            df_v1,
            df_rem[remapped_cols],
            on=["IlmnID", "Name"],
            how="left",
            validate="1:1",
        )
    except pd.errors.MergeError as e:
        log.error("Merge validation failed: %s. Falling back to non-validated merge.", e)
        master = pd.merge(df_v1, df_rem[remapped_cols], on=["IlmnID", "Name"], how="left")

    log.info("After join: %d rows in v1_master", len(master))

    # --- is_remapped flag ---
    master["is_remapped"] = master["Chr_equCab3"].notna()
    log.info("Markers with EquCab3 remapping: %d / %d", master["is_remapped"].sum(), len(master))

    # --- EC3 conflict detection for markers already in EquCab3 in v1 ---
    is_ec3_native = master["GenomeBuild"] == EQUCAB3_BUILD
    master["is_ec3_native"] = is_ec3_native
    n_ec3_native = is_ec3_native.sum()
    log.info("v1 markers with native EquCab3 coords (GenomeBuild='%s'): %d", EQUCAB3_BUILD, n_ec3_native)

    if "Chr_equCab3" in master.columns and "MapInfo_equCab3" in master.columns:
        chr_v1_norm = normalize_chr(master["Chr"])
        chr_rem_norm = normalize_chr(master["Chr_equCab3"].fillna(""))
        pos_v1 = pd.to_numeric(master["MapInfo"], errors="coerce")
        pos_rem = pd.to_numeric(master["MapInfo_equCab3"], errors="coerce")

        coord_conflict = (
            is_ec3_native
            & master["is_remapped"]
            & (
                (chr_v1_norm != chr_rem_norm)
                | (pos_v1 != pos_rem)
            )
        )
        master["coord_conflict_ec3"] = coord_conflict
        n_conflicts = coord_conflict.sum()
        log.info("EC3 coordinate conflicts (v1 native vs v1.remapped): %d", n_conflicts)
    else:
        master["coord_conflict_ec3"] = False

    # --- Normalized chromosome columns ---
    master["chr_ec3_norm"] = normalize_chr(master.get("Chr_equCab3", pd.Series("", index=master.index)).fillna(""))
    master["flag_chr_unresolvable"] = flag_ncbi_scaffolds(master["chr_ec3_norm"])

    # --- Positional key (EquCab3) ---
    master["pos_key_ec3"] = make_pos_key(master["chr_ec3_norm"], master.get("MapInfo_equCab3", pd.Series(0, index=master.index)))

    # --- MAPQ flags ---
    master["mapq_probe_int"] = pd.to_numeric(master.get("MAPQ_Probe", 0), errors="coerce").fillna(0).astype(int)
    master["mapq_genomic_int"] = pd.to_numeric(master.get("MAPQ_TopGenomicSeq", 0), errors="coerce").fillna(0).astype(int)

    master["flag_failed_probe_map"]  = master["mapq_probe_int"] == 0
    master["flag_low_mapq_probe"]    = master["mapq_probe_int"] < MAPQ_MIN_PROBE
    master["flag_low_mapq_genomic"]  = master["mapq_genomic_int"] < MAPQ_MIN_GENOMIC
    master["flag_unreliable_position"] = master["mapq_genomic_int"] == 0

    log.info("MAPQ_Probe = 0:           %d (%.1f%%)", master["flag_failed_probe_map"].sum(),
             100 * master["flag_failed_probe_map"].mean())
    log.info("MAPQ_TopGenomicSeq = 0:   %d (unreliable position)", master["flag_unreliable_position"].sum())

    return master
