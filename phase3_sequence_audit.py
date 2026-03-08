"""
Phase 3: Probe sequence comparison, strand flip detection, allele audit.
"""
import logging
import pandas as pd
from normalization import (
    is_ambiguous_snp, check_allele_complement, classify_variant, COMPLEMENT
)

log = logging.getLogger(__name__)


def _detect_ilmn_strand_flip(strand_v1, strand_v2) -> str:
    """Compare IlmnStrand (TOP/BOT) between v1 and v2.
    Returns: 'consistent', 'flip', or 'na'
    """
    s1 = (strand_v1 or "").strip().upper() if isinstance(strand_v1, str) else ""
    s2 = (strand_v2 or "").strip().upper() if isinstance(strand_v2, str) else ""
    if s1 not in ("TOP", "BOT") or s2 not in ("TOP", "BOT"):
        return "na"
    return "consistent" if s1 == s2 else "flip"


def run_sequence_audit(master: pd.DataFrame) -> pd.DataFrame:
    """Add sequence and strand audit columns to the master dataframe.

    Expects columns suffixed _v1 and _v2 for overlapping fields.
    The Strand_equCab3 column (from v1.remapped, no suffix) is used for strand comparison.

    Returns the master dataframe with new audit columns.
    """
    log.info("=== Phase 3: Sequence & Allele Audit ===")

    # --- Resolve column names with possible suffixes ---
    def _col(name, prefer_v1=True):
        """Find column by name, trying _v1/_v2 suffixes if bare name absent."""
        if name in master.columns:
            return name
        suffix = "_v1" if prefer_v1 else "_v2"
        return name + suffix if (name + suffix) in master.columns else None

    allele_a_v1 = _col("AlleleA_ProbeSeq", prefer_v1=True)
    allele_b_v1 = _col("AlleleB_ProbeSeq", prefer_v1=True)
    allele_a_v2 = _col("AlleleA_ProbeSeq", prefer_v1=False)
    allele_b_v2 = _col("AlleleB_ProbeSeq", prefer_v1=False)
    snp_v1      = _col("SNP", prefer_v1=True)
    snp_v2      = _col("SNP", prefer_v1=False)
    ilmn_strand_v1 = "IlmnStrand_v1" if "IlmnStrand_v1" in master.columns else _col("IlmnStrand", prefer_v1=True)
    ilmn_strand_v2 = "IlmnStrand_v2" if "IlmnStrand_v2" in master.columns else _col("IlmnStrand", prefer_v1=False)

    # Only audit common markers (in both v1 and v2)
    common = master["in_v1"] & master["in_v2"]
    n_common = common.sum()
    log.info("Common markers for sequence audit: %d", n_common)

    # --- Probe Sequence Comparison ---
    if allele_a_v1 and allele_a_v2:
        a_match = _seq_eq(master[allele_a_v1], master[allele_a_v2])
    else:
        log.warning("AlleleA_ProbeSeq columns not found; skipping AlleleA comparison")
        a_match = pd.Series(True, index=master.index)

    if allele_b_v1 and allele_b_v2:
        b_match = _seq_eq(master[allele_b_v1], master[allele_b_v2])
    else:
        b_match = pd.Series(True, index=master.index)

    master["AlleleA_seq_match"] = a_match.astype("boolean")
    master["AlleleB_seq_match"] = b_match.astype("boolean")
    master["sequence_match"] = (a_match & b_match).astype("boolean")
    # Only meaningful for common markers; pd.NA requires nullable boolean dtype
    master.loc[~common, ["AlleleA_seq_match", "AlleleB_seq_match", "sequence_match"]] = pd.NA

    n_seq_mismatch = (~master.loc[common, "sequence_match"]).sum()
    log.info("Probe sequence identical: %d | Changed: %d",
             n_common - n_seq_mismatch, n_seq_mismatch)

    # --- Variant Type & Ambiguous SNP flags ---
    snp_col = snp_v1 or snp_v2
    if snp_col:
        master["is_ambiguous_snp"] = master[snp_col].apply(
            lambda x: is_ambiguous_snp(x) if isinstance(x, str) else False
        )
        variant_col = "variant_type" if "variant_type" in master.columns else None
        if variant_col is None and snp_col:
            master["variant_type"] = classify_variant(master[snp_col])
    else:
        master["is_ambiguous_snp"] = False

    # --- Strand Flip Detection ---
    if ilmn_strand_v1 and ilmn_strand_v2:
        master["strand_flip_detected"] = master.apply(
            lambda r: _detect_ilmn_strand_flip(r.get(ilmn_strand_v1), r.get(ilmn_strand_v2)),
            axis=1,
        )
        master.loc[~common, "strand_flip_detected"] = "na"
    else:
        log.warning("Strand columns not found; strand flip detection skipped")
        master["strand_flip_detected"] = "na"

    n_flips = (master.loc[common, "strand_flip_detected"] == "flip").sum()
    log.info("Strand flips detected (common markers): %d", n_flips)
    log.info("Ambiguous SNPs (A/T or C/G, common): %d",
             master.loc[common, "is_ambiguous_snp"].sum())

    # --- Allele Complement Validation ---
    if snp_v1 and snp_v2:
        master["allele_complement_check"] = master.apply(
            lambda r: (
                check_allele_complement(r.get(snp_v1, ""), r.get(snp_v2, ""))
                if (common[r.name] and r.get("strand_flip_detected") == "flip")
                else "not_assessed"
            ),
            axis=1,
        )
    else:
        master["allele_complement_check"] = "not_assessed"

    # --- Indel flag ---
    if snp_col:
        from config import INDEL_SNPS
        master["is_indel"] = master[snp_col].isin(INDEL_SNPS)
        master["indel_strand_check_skipped"] = master["is_indel"] & common
    else:
        master["is_indel"] = False
        master["indel_strand_check_skipped"] = False

    log.info("Indels (strand check skipped): %d",
             master.loc[common, "indel_strand_check_skipped"].sum())

    return master


def _seq_eq(s1: pd.Series, s2: pd.Series) -> pd.Series:
    """Case-insensitive string equality; NaN == NaN is True (both absent = match).

    Uses .where() to avoid in-place dtype-mismatch assignments that raise
    FutureWarning in pandas when mixing numpy bool and nullable boolean dtypes.
    """
    both_null = s1.isna() & s2.isna()
    either_null = s1.isna() | s2.isna()
    str_match = s1.str.upper() == s2.str.upper()   # NaN inputs produce NaN here
    # one null → False; both null → True; both present → string comparison result
    result = str_match.where(~either_null, other=False)
    result = result.where(~both_null, other=True)
    return result.astype("boolean")
