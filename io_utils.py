"""
File loading utilities for Illumina manifest CSV files.
"""
import logging
import pandas as pd
from config import (
    V1_PATH, V1_REMAPPED_PATH, V2_PATH,
    EXPECTED_V1_ROWS, EXPECTED_REMAPPED_ROWS, EXPECTED_V2_ROWS,
    CONTROL_STRAND_VALUES,
)

log = logging.getLogger(__name__)


def _find_controls_line(path: str, skip_metadata: bool) -> int:
    """Return nrows (assay data rows before [Controls]) by pre-scanning the file.

    Returns None if [Controls] is not found (read the whole file).

    IMPORTANT: Do NOT use pandas comment='[' to stop at [Controls].
    The comment parameter truncates every data row at the first '[' character,
    corrupting SNP values like [A/G] and every column that follows
    (GenomeBuild, Chr, MapInfo, AlleleA_ProbeSeq, etc.).

    File layout with skip_metadata=True:
        Lines 0-6  : 7 Illumina metadata lines  (skipped by skiprows=7)
        Line 7     : column header IlmnID,Name,... (used as header by pandas)
        Lines 8+   : assay data rows
        [Controls] : marks end of assay data

    pre_data_lines = 8 (7 skipped + 1 header), so:
        nrows = controls_absolute_lineno - 8
    """
    pre_data_lines = 8 if skip_metadata else 1  # skipped lines + 1 header
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for lineno, line in enumerate(fh):
            if lineno < pre_data_lines:
                continue
            if line.startswith("[Controls]"):
                return lineno - pre_data_lines  # number of assay data rows
    return None


def _load_illumina_manifest(path: str, skip_metadata: bool) -> pd.DataFrame:
    """Load an Illumina manifest CSV.

    Args:
        path: Path to the CSV file.
        skip_metadata: If True, skip the 7 Illumina metadata lines before the
                       column header (line 7 = IlmnID,Name,...). Uses nrows to
                       stop before [Controls] instead of comment='['.
                       Set False for v1.remapped which has no metadata preamble.
    """
    nrows = _find_controls_line(path, skip_metadata)
    if nrows is not None:
        log.info("%s: [Controls] found — reading %d assay rows", path, nrows)

    kwargs = dict(
        dtype=str,
        keep_default_na=False,
        na_values=[""],
        nrows=nrows,
    )
    if skip_metadata:
        kwargs["skiprows"] = 7  # skip lines 0-6 (metadata); line 7 becomes header
    df = pd.read_csv(path, **kwargs)

    # Strip whitespace from all string columns
    str_cols = df.select_dtypes("object").columns
    df[str_cols] = df[str_cols].apply(lambda c: c.str.strip())

    # Drop rows where IlmnID is null (artifact of comment parser or blank trailing lines)
    before = len(df)
    df = df.dropna(subset=["IlmnID"])
    if len(df) < before:
        log.warning("%s: dropped %d null-IlmnID rows", path, before - len(df))

    return df


def load_v1() -> pd.DataFrame:
    log.info("Loading v1: %s", V1_PATH)
    df = _load_illumina_manifest(V1_PATH, skip_metadata=True)
    _assert_rows(df, EXPECTED_V1_ROWS, V1_PATH)
    # Remove any control-probe artifacts that slipped through (IlmnStrand color names)
    mask = df["IlmnStrand"].isin(CONTROL_STRAND_VALUES)
    if mask.any():
        log.warning("v1: removing %d control-artifact rows (IlmnStrand color values)", mask.sum())
        df = df[~mask].reset_index(drop=True)
    log.info("v1 loaded: %d rows, %d columns", len(df), df.shape[1])
    return df


def load_remapped() -> pd.DataFrame:
    log.info("Loading v1.remapped: %s", V1_REMAPPED_PATH)
    df = _load_illumina_manifest(V1_REMAPPED_PATH, skip_metadata=False)
    _assert_rows(df, EXPECTED_REMAPPED_ROWS, V1_REMAPPED_PATH)
    log.info("v1.remapped loaded: %d rows, %d columns", len(df), df.shape[1])
    return df


def load_v2() -> pd.DataFrame:
    log.info("Loading v2: %s", V2_PATH)
    df = _load_illumina_manifest(V2_PATH, skip_metadata=True)
    _assert_rows(df, EXPECTED_V2_ROWS, V2_PATH)
    mask = df["IlmnStrand"].isin(CONTROL_STRAND_VALUES)
    if mask.any():
        log.warning("v2: removing %d control-artifact rows", mask.sum())
        df = df[~mask].reset_index(drop=True)
    log.info("v2 loaded: %d rows, %d columns", len(df), df.shape[1])
    return df


def _assert_rows(df: pd.DataFrame, expected: int, path: str) -> None:
    if len(df) != expected:
        log.warning(
            "%s: expected %d assay rows, got %d. Proceeding with caution.",
            path, expected, len(df)
        )
    else:
        log.info("%s: row count OK (%d)", path, expected)
