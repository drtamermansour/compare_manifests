"""
Configuration for the Equine80select manifest comparison pipeline.
"""
import os

# --- Input File Paths ---
V1_PATH = "Equine80select.v1.csv"
V1_REMAPPED_PATH = "Equine80select.v1.remapped.csv"
V2_PATH = "Equine80select.v2.csv"

# --- Output Directory ---
OUTPUT_DIR = "output"
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs("logs", exist_ok=True)

# --- Output File Paths ---
MASTER_TABLE_CSV        = os.path.join(OUTPUT_DIR, "master_comparison_table.csv")
DISCREPANCY_REPORT_CSV  = os.path.join(OUTPUT_DIR, "discrepancy_report.csv")
SUMMARY_TXT             = os.path.join(OUTPUT_DIR, "summary_statistics.txt")
PLINK_FLIP_TXT          = os.path.join(OUTPUT_DIR, "plink_flip.txt")
PLINK_AMBIG_TXT         = os.path.join(OUTPUT_DIR, "plink_ambiguous_exclude.txt")
PLINK_STRAND_CSV        = os.path.join(OUTPUT_DIR, "plink_strand_report.csv")
VENN_PNG                = os.path.join(OUTPUT_DIR, "venn_overlap.png")
CHR_DIST_PNG            = os.path.join(OUTPUT_DIR, "chromosome_distribution.png")

# --- Thresholds ---
MAPQ_MIN_PROBE   = 20   # below this: flag as low-quality probe remapping
MAPQ_MIN_GENOMIC = 20   # below this: flag as low-quality genomic remapping

# Maximum allowed coordinate difference (bp) for positional matching.
# Set to 0 for exact match, or >0 to allow liftover wobble.
MAX_COORD_SHIFT_BP = 0

# --- Expected Row Counts (for validation assertions) ---
EXPECTED_V1_ROWS       = 81_974
EXPECTED_REMAPPED_ROWS = 81_974
EXPECTED_V2_ROWS       = 84_319

# --- Genome Build Constants ---
EQUCAB2_BUILD = "2"
EQUCAB3_BUILD = "3"

# --- Chromosome sort order for visualization ---
CHR_SORT_ORDER = [str(i) for i in range(1, 32)] + ["X", "Y", "Un", "Chr_unresolvable"]

# --- Chromosome name normalization map ---
CHR_NORM_MAP = {
    "X_NC_009175.3": "X",
    "0": "Un",
    "": "Un",
}

# --- IlmnStrand values that indicate control probes (not assay rows) ---
CONTROL_STRAND_VALUES = {"Green", "Red", "Blue", "Purple", "Black"}

# --- SNP/variant classification ---
INDEL_SNPS    = {"[D/I]", "[I/D]"}
CONTROL_SNPS  = {"DNP (High)", "DNP (Bgnd)"}

# --- Palindromic allele pairs (ambiguous for strand resolution) ---
AMBIGUOUS_ALLELE_PAIRS = [{"A", "T"}, {"C", "G"}]
