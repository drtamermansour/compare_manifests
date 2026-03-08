"""
Output writers: CSV reports, summary statistics, and PLINK files.
"""
import logging
import pandas as pd
from datetime import datetime
from config import (
    MASTER_TABLE_CSV, DISCREPANCY_REPORT_CSV, SUMMARY_TXT,
    PLINK_FLIP_TXT, PLINK_AMBIG_TXT, PLINK_STRAND_CSV,
    EXPECTED_V1_ROWS, EXPECTED_REMAPPED_ROWS, EXPECTED_V2_ROWS,
)

log = logging.getLogger(__name__)

# Columns to include in the master comparison table (order matters for readability)
MASTER_COLS_PRIORITY = [
    "IlmnID", "Name_v1", "Name_v2",
    "in_v1", "in_v2", "is_remapped", "is_ec3_native", "match_tier", "coord_match", "multi_match",
    "Chr_equCab3", "MapInfo_equCab3", "Strand_equCab3",
    "chr_v2_norm", "pos_key_v2",
    "coord_conflict_ec3",
    "MAPQ_Probe", "MAPQ_TopGenomicSeq",
    "flag_failed_probe_map", "flag_low_mapq_probe",
    "flag_low_mapq_genomic", "flag_unreliable_position", "flag_chr_unresolvable",
    "sequence_match", "AlleleA_seq_match", "AlleleB_seq_match",
    "strand_flip_detected", "allele_complement_check",
    "is_ambiguous_snp", "is_indel", "indel_strand_check_skipped",
    "variant_type",
]


def write_master_table(master: pd.DataFrame) -> None:
    log.info("Writing master comparison table: %s", MASTER_TABLE_CSV)
    # Write priority columns first, then append remaining columns
    present_priority = [c for c in MASTER_COLS_PRIORITY if c in master.columns]
    remaining = [c for c in master.columns if c not in present_priority]
    ordered = present_priority + remaining
    master[ordered].to_csv(MASTER_TABLE_CSV, index=False)
    log.info("Master table: %d rows, %d columns", len(master), len(ordered))


def write_discrepancy_report(master: pd.DataFrame) -> None:
    log.info("Writing discrepancy report: %s", DISCREPANCY_REPORT_CSV)

    # Severity ordering: compute a severity score for sorting
    severity_flags = [
        "flag_unreliable_position",
        "coord_conflict_ec3",
        "flag_chr_unresolvable",
    ]
    bool_flags = severity_flags + [
        "flag_low_mapq_genomic",
        "flag_failed_probe_map",
        "flag_low_mapq_probe",
    ]

    # Boolean columns only (exclude string columns from "any flag" check)
    flag_cols = [c for c in bool_flags if c in master.columns]
    seq_mismatch = master.get("sequence_match", pd.Series(True, index=master.index)).eq(False).fillna(False)
    strand_flip  = master.get("strand_flip_detected", pd.Series("na", index=master.index)) == "flip"
    id_shift     = master.get("match_tier", pd.Series("", index=master.index)) == "id_shift"
    multi_match  = master.get("multi_match", pd.Series(False, index=master.index))

    any_flag = (
        master[flag_cols].any(axis=1)
        | seq_mismatch
        | strand_flip
        | id_shift
        | multi_match
    ) if flag_cols else (seq_mismatch | strand_flip | id_shift | multi_match)

    disc = master[any_flag].copy()

    # Severity score for sorting (higher = more severe)
    def _severity(row):
        score = 0
        if row.get("flag_unreliable_position"):   score += 100
        if row.get("coord_conflict_ec3"):          score += 80
        if row.get("sequence_match") is False:      score += 60
        if row.get("strand_flip_detected") == "flip": score += 40
        if row.get("match_tier") == "id_shift":    score += 20
        if row.get("multi_match"):                 score += 10
        return score

    disc["_severity"] = disc.apply(_severity, axis=1)
    disc = disc.sort_values("_severity", ascending=False).drop(columns=["_severity"])

    present_priority = [c for c in MASTER_COLS_PRIORITY if c in disc.columns]
    remaining = [c for c in disc.columns if c not in present_priority]
    disc[present_priority + remaining].to_csv(DISCREPANCY_REPORT_CSV, index=False)
    log.info("Discrepancy report: %d rows", len(disc))


def write_plink_files(master: pd.DataFrame) -> None:
    log.info("Writing PLINK files...")
    common = master["in_v1"] & master["in_v2"]

    # Name column may be suffixed to Name_v1 after the outer join
    name_col = "Name" if "Name" in master.columns else (
        "Name_v1" if "Name_v1" in master.columns else "Name_v2"
    )

    # plink_flip.txt — confirmed flips, excluding ambiguous SNPs and indels
    flip_mask = (
        common
        & (master.get("strand_flip_detected", "na") == "flip")
        & (~master.get("is_ambiguous_snp", pd.Series(False, index=master.index)))
        & (~master.get("is_indel", pd.Series(False, index=master.index)))
    )
    flip_names = master.loc[flip_mask, name_col].dropna()
    flip_names.to_csv(PLINK_FLIP_TXT, index=False, header=False)
    log.info("plink_flip.txt: %d markers", len(flip_names))

    # plink_ambiguous_exclude.txt — A/T and C/G SNPs
    ambig_mask = common & master.get("is_ambiguous_snp", pd.Series(False, index=master.index))
    ambig_names = master.loc[ambig_mask, name_col].dropna()
    ambig_names.to_csv(PLINK_AMBIG_TXT, index=False, header=False)
    log.info("plink_ambiguous_exclude.txt: %d markers", len(ambig_names))

    # plink_strand_report.csv — full strand action table
    strand_cols = [name_col, "IlmnID"]
    for c in ["IlmnStrand_v1", "IlmnStrand_v2", "Strand_equCab3", "RefStrand_v2",
              "strand_flip_detected", "is_ambiguous_snp", "is_indel",
              "allele_complement_check"]:
        if c in master.columns:
            strand_cols.append(c)

    def _plink_action(row):
        if row.get("is_indel"):                          return "indel_skip"
        if row.get("is_ambiguous_snp"):                  return "exclude_ambiguous"
        flip = row.get("strand_flip_detected", "na")
        if flip == "flip":                               return "flip"
        if flip == "consistent":                         return "keep"
        return "na"

    strand_df = master.loc[common, strand_cols].copy()
    strand_df["plink_action"] = master.loc[common].apply(_plink_action, axis=1)
    strand_df.to_csv(PLINK_STRAND_CSV, index=False)
    log.info("plink_strand_report.csv: %d rows", len(strand_df))


def write_summary(master: pd.DataFrame, v1_master: pd.DataFrame) -> None:
    log.info("Writing summary statistics: %s", SUMMARY_TXT)
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    common = master["in_v1"] & master["in_v2"]
    tier_counts = master["match_tier"].value_counts().to_dict()

    def _get(col, condition=None):
        if col not in master.columns:
            return "N/A"
        s = master[col] if condition is None else master.loc[condition, col]
        # Boolean flag columns become float64 after outer join (NaN for v2-only rows).
        # s.eq(True) returns False for NaN, 0.0, False — correctly counts True values
        # without triggering FutureWarning from fillna downcasting.
        return int(s.eq(True).sum())

    def _v1(col):
        return int(v1_master[col].sum()) if col in v1_master.columns else "N/A"

    n_ec3_native = int(v1_master["is_ec3_native"].sum()) if "is_ec3_native" in v1_master.columns else 0
    n_ec2_native = len(v1_master) - n_ec3_native

    lines = [
        f"=== Equine80select Manifest Comparison Summary ===",
        f"Generated: {now}",
        f"",
        f"FILE STATISTICS",
        f"  v1 assay markers:               {EXPECTED_V1_ROWS:>10,}",
        f"  v1.remapped markers:            {EXPECTED_REMAPPED_ROWS:>10,}",
        f"  v2 assay markers:               {EXPECTED_V2_ROWS:>10,}",
        f"  New markers in v2:              {tier_counts.get('v2_only', 0):>10,}",
        f"",
        f"PHASE 1: v1 Integration",
        f"  EquCab2 coords (v1 native):     {n_ec2_native:>10,}",
        f"  EquCab3 coords (v1 native):     {n_ec3_native:>10,}",
        f"  EC3 coordinate conflicts:       {_get('coord_conflict_ec3'):>10}",
        f"  MAPQ_Probe = 0 (flagged only):  {_get('flag_failed_probe_map'):>10}",
        f"  MAPQ_TopGenomicSeq = 0:         {_get('flag_unreliable_position'):>10}",
        f"  Chr unresolvable (scaffolds):   {_get('flag_chr_unresolvable'):>10}",
        f"",
        f"PHASE 2: v1 vs v2 Matching",
        f"  Strict matches:                 {tier_counts.get('strict', 0):>10,}",
        f"  ID Shifts (probe renumbered):   {tier_counts.get('id_shift', 0):>10,}",
        f"  Positional-only matches:        {tier_counts.get('positional', 0):>10,}",
        f"  Coord mismatches (pos shifted): {tier_counts.get('coord_mismatch', 0):>10,}",
        f"  v1-only markers (control):      {tier_counts.get('v1_only', 0):>10,}",
        f"  v2-only markers (new):          {tier_counts.get('v2_only', 0):>10,}",
        f"  Multi-match (positional):       {int(master['multi_match'].sum()):>10,}",
        f"",
        f"PHASE 3: Sequence Audit (common markers: {int(common.sum()):,})",
        f"  Probe sequence identical:       {int(master.loc[common, 'sequence_match'].sum()) if 'sequence_match' in master.columns else 'N/A':>10}",
        f"  Probe sequence changed:         {int(master.loc[common, 'sequence_match'].eq(False).fillna(False).sum()) if 'sequence_match' in master.columns else 'N/A':>10}",
        f"  Strand flips (non-ambiguous):   {int((master.loc[common, 'strand_flip_detected'] == 'flip').sum()) if 'strand_flip_detected' in master.columns else 'N/A':>10}",
        f"  Ambiguous SNPs (A/T or C/G):    {int(master.loc[common, 'is_ambiguous_snp'].sum()) if 'is_ambiguous_snp' in master.columns else 'N/A':>10}",
        f"  Indels (strand check skipped):  {int(master.loc[common, 'is_indel'].sum()) if 'is_indel' in master.columns else 'N/A':>10}",
    ]

    text = "\n".join(lines)
    with open(SUMMARY_TXT, "w") as f:
        f.write(text + "\n")
    print("\n" + text + "\n")
