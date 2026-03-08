"""
Equine80select Manifest Comparison Pipeline
============================================
Reconciles three Illumina microarray manifest files:
  - Equine80select.v1.csv        (legacy, EquCab2 primary coords)
  - Equine80select.v1.remapped.csv (bridge file, EquCab3 remapping)
  - Equine80select.v2.csv        (new manifest, native EquCab3)

Outputs:
  output/master_comparison_table.csv
  output/discrepancy_report.csv
  output/summary_statistics.txt
  output/plink_flip.txt
  output/plink_ambiguous_exclude.txt
  output/plink_strand_report.csv
  output/venn_overlap.png
  output/chromosome_distribution.png

Usage:
  python manifest_compare.py
"""
import logging
import sys
from datetime import datetime
from pathlib import Path

# ── Logging setup ────────────────────────────────────────────────────────────
log_path = Path("logs") / f"pipeline_{datetime.now():%Y%m%d_%H%M%S}.log"
log_path.parent.mkdir(exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(module)s: %(message)s",
    handlers=[
        logging.FileHandler(log_path),
        logging.StreamHandler(sys.stdout),
    ],
)
log = logging.getLogger(__name__)


def main():
    log.info("=" * 60)
    log.info("Equine80select Manifest Comparison Pipeline")
    log.info("=" * 60)

    # ── Phase 1: Load & Integrate ─────────────────────────────────────────────
    from io_utils import load_v1, load_remapped, load_v2
    from phase1_integration import build_v1_master

    df_v1  = load_v1()
    df_rem = load_remapped()
    df_v2  = load_v2()

    v1_master = build_v1_master(df_v1, df_rem)

    # ── Phase 2: Three-Tier Matching ──────────────────────────────────────────
    from phase2_comparison import compare_v1_v2

    master, tier2_id_shifts, tier3_positional = compare_v1_v2(v1_master, df_v2)

    # ── Phase 3: Sequence & Strand Audit ─────────────────────────────────────
    from phase3_sequence_audit import run_sequence_audit

    master = run_sequence_audit(master)

    # ── Outputs ───────────────────────────────────────────────────────────────
    from reporting import (
        write_master_table, write_discrepancy_report,
        write_plink_files, write_summary,
    )
    from visualization import plot_venn, plot_chromosome_distribution

    write_master_table(master)
    write_discrepancy_report(master)
    write_plink_files(master)
    write_summary(master, v1_master)
    plot_venn(master)
    plot_chromosome_distribution(master)

    log.info("=" * 60)
    log.info("Pipeline complete. Outputs written to: output/")
    log.info("Log saved to: %s", log_path)
    log.info("=" * 60)


if __name__ == "__main__":
    main()
