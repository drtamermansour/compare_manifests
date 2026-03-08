"""
Normalization utilities for manifest columns.
"""
import pandas as pd
from config import CHR_NORM_MAP, INDEL_SNPS, CONTROL_SNPS, AMBIGUOUS_ALLELE_PAIRS


def normalize_address_id(series: pd.Series) -> pd.Series:
    """Strip leading zeros from Illumina AddressA/B_ID strings.
    v1 uses zero-padded 10-digit strings; v2 uses unpadded integers.
    Canonical form: unpadded integer string (e.g. '12758484').
    """
    return series.str.lstrip("0").replace("", "0")


def normalize_chr(series: pd.Series) -> pd.Series:
    """Normalize chromosome name strings to a canonical form.
    - X_NC_009175.3  → X
    - '0' or ''      → Un
    - Un_NW_*        → kept as-is (will be flagged as Chr_unresolvable downstream)
    """
    normed = series.fillna("").str.strip()
    normed = normed.replace(CHR_NORM_MAP)
    # Blank values after replacement
    normed = normed.replace("", "Un")
    return normed


def make_pos_key(chr_series: pd.Series, pos_series: pd.Series) -> pd.Series:
    """Create a string positional key 'CHR:POS' for hash-based joins.
    Unplaced markers (pos=0, NaN, or chr=Un) get key 'Un:-1'.
    """
    chr_norm = normalize_chr(chr_series)
    pos_int = pd.to_numeric(pos_series, errors="coerce").fillna(0).astype(int)

    # Treat position=0 as unplaced
    unplaced = (chr_norm.isin(["Un", ""])) | (pos_int == 0)
    key = chr_norm + ":" + pos_int.astype(str)
    key[unplaced] = "Un:-1"
    return key


def flag_ncbi_scaffolds(chr_series: pd.Series) -> pd.Series:
    """Return boolean mask: True if chromosome is an NCBI scaffold (Un_NW_* etc.).
    These cannot be matched to v2's simplified unplaced representation.
    """
    return chr_series.str.startswith("Un_", na=False) | chr_series.str.startswith("NW_", na=False)


def classify_variant(snp_series: pd.Series) -> pd.Series:
    """Classify each SNP field into: snp, indel, control, or unknown."""
    def _classify(val):
        if not isinstance(val, str):
            return "unknown"
        val = val.strip()
        if val in INDEL_SNPS:
            return "indel"
        if val in CONTROL_SNPS:
            return "control"
        if val.startswith("[") and "/" in val and val.endswith("]"):
            return "snp"
        return "unknown"
    return snp_series.map(_classify)


def parse_alleles(snp_str: str):
    """Extract (allele1, allele2) from '[A/G]' format.
    Returns None for indels, controls, or unparseable values.
    """
    if not isinstance(snp_str, str):
        return None
    snp_str = snp_str.strip()
    if not (snp_str.startswith("[") and snp_str.endswith("]") and "/" in snp_str):
        return None
    inner = snp_str[1:-1]
    parts = inner.split("/")
    if len(parts) == 2 and all(len(p) == 1 for p in parts):
        return (parts[0].upper(), parts[1].upper())
    return None


def is_ambiguous_snp(snp_str: str) -> bool:
    """Return True if SNP is a palindromic (A/T or C/G) pair.
    These cannot be strand-resolved from alleles alone.
    """
    alleles = parse_alleles(snp_str)
    if alleles is None:
        return False
    return set(alleles) in AMBIGUOUS_ALLELE_PAIRS


COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}


def check_allele_complement(snp_v1: str, snp_v2: str) -> str:
    """Check relationship between v1 and v2 SNP alleles.
    Returns: 'same', 'complement', 'mismatch', or 'cannot_assess'
    """
    a1 = parse_alleles(snp_v1)
    a2 = parse_alleles(snp_v2)
    if a1 is None or a2 is None:
        return "cannot_assess"
    comp_a1 = (COMPLEMENT.get(a1[0], "?"), COMPLEMENT.get(a1[1], "?"))
    if a2 == a1 or a2 == (a1[1], a1[0]):
        return "same"
    if a2 == comp_a1 or a2 == (comp_a1[1], comp_a1[0]):
        return "complement"
    return "mismatch"
