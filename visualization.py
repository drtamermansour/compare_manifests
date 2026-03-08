"""
Visualization: Venn diagram and chromosome distribution bar chart.
"""
import logging
import pandas as pd
from config import VENN_PNG, CHR_DIST_PNG, CHR_SORT_ORDER

try:
    import matplotlib
    matplotlib.use("Agg")  # non-interactive backend for server/script use
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker
    _MATPLOTLIB_AVAILABLE = True
except ImportError:
    _MATPLOTLIB_AVAILABLE = False

log = logging.getLogger(__name__)


def plot_venn(master: pd.DataFrame) -> None:
    """2-circle Venn diagram of IlmnID overlap between v1 and v2."""
    if not _MATPLOTLIB_AVAILABLE:
        log.warning("matplotlib not installed. Skipping Venn diagram. "
                    "Install with: pip install matplotlib matplotlib-venn")
        return
    try:
        from matplotlib_venn import venn2
    except ImportError:
        log.warning("matplotlib_venn not installed. Skipping Venn diagram. "
                    "Install with: pip install matplotlib-venn")
        return

    v1_ids = set(master.loc[master["in_v1"], "IlmnID"])
    v2_ids = set(master.loc[master["in_v2"], "IlmnID"])
    both   = v1_ids & v2_ids
    v1_only = v1_ids - both
    v2_only = v2_ids - both

    fig, ax = plt.subplots(figsize=(7, 5))
    venn2(
        subsets=(len(v1_only), len(v2_only), len(both)),
        set_labels=("v1 (Equine80select)", "v2 (Equine80select)"),
        ax=ax,
    )
    ax.set_title("IlmnID Overlap: v1 vs v2 Manifests", fontsize=13, pad=12)
    fig.tight_layout()
    fig.savefig(VENN_PNG, dpi=150)
    plt.close(fig)
    log.info("Venn diagram saved: %s", VENN_PNG)


def plot_chromosome_distribution(master: pd.DataFrame) -> None:
    """Grouped bar chart of marker counts per chromosome for v1 vs v2."""
    if not _MATPLOTLIB_AVAILABLE:
        log.warning("matplotlib not installed. Skipping chromosome distribution plot. "
                    "Install with: pip install matplotlib")
        return
    # Use normalized chromosome columns
    chr_v1_col = "chr_ec3_norm" if "chr_ec3_norm" in master.columns else None
    chr_v2_col = "chr_v2_norm" if "chr_v2_norm" in master.columns else None

    if not chr_v1_col or not chr_v2_col:
        log.warning("Normalized chromosome columns not found; skipping chromosome distribution plot.")
        return

    def _collapse_un(series: pd.Series) -> pd.Series:
        """Collapse Un_NW_* scaffold names into 'Un' for plotting.
        v1 uses NCBI scaffold names; v2 uses chr=0 → already normalized to 'Un'.
        Both are logically unplaced and should appear as a single 'Un' bar.
        """
        return series.where(~series.str.startswith("Un_", na=False), other="Un")

    v1_counts = _collapse_un(master.loc[master["in_v1"], chr_v1_col]).value_counts().rename("v1")
    v2_counts = _collapse_un(master.loc[master["in_v2"], chr_v2_col]).value_counts().rename("v2")

    df_plot = pd.DataFrame({"v1": v1_counts, "v2": v2_counts}).fillna(0).astype(int)

    # Sort chromosomes canonically
    present_chrs = list(df_plot.index)
    order = [c for c in CHR_SORT_ORDER if c in present_chrs]
    order += sorted(set(present_chrs) - set(order))  # any unexpected values at end
    df_plot = df_plot.reindex(order).fillna(0)

    x = range(len(df_plot))
    width = 0.4

    fig, ax = plt.subplots(figsize=(max(14, len(df_plot) * 0.5), 5))
    bars_v1 = ax.bar([i - width / 2 for i in x], df_plot["v1"], width=width, label="v1", color="#4C72B0", alpha=0.85)
    bars_v2 = ax.bar([i + width / 2 for i in x], df_plot["v2"], width=width, label="v2", color="#DD8452", alpha=0.85)

    ax.set_xticks(list(x))
    ax.set_xticklabels(df_plot.index, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Marker Count")
    ax.set_title("Marker Distribution by Chromosome: v1 vs v2 (EquCab3)", fontsize=12, pad=10)
    ax.legend()
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda val, _: f"{int(val):,}"))
    ax.grid(axis="y", linestyle="--", alpha=0.4)

    fig.tight_layout()
    fig.savefig(CHR_DIST_PNG, dpi=150)
    plt.close(fig)
    log.info("Chromosome distribution plot saved: %s", CHR_DIST_PNG)
