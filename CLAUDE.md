# Equine80select Manifest Comparison Pipeline

## Project Purpose

Reconciles three Illumina Equine80select microarray manifest files to migrate and compare markers from an older array version (v1) to a newer version (v2), using **EquCab3** as the unified coordinate system.

## How to Run

```bash
cd /home/tahmed/compare_manifests

# Install dependencies (first time only)
pip install pandas matplotlib matplotlib-venn

# Run the full pipeline
python manifest_compare.py
```

All outputs are written to `output/`. Timestamped logs are written to `logs/`.

---

## Input Files

| File | Rows | Cols | Notes |
|---|---|---|---|
| `Equine80select.v1.csv` | 81,974 | 19 | Legacy manifest. EquCab2 primary coords; ~6,022 markers already in EquCab3 (GenomeBuild='3'). 7 metadata rows before header + `[Controls]` section at tail. |
| `Equine80select.v1.remapped.csv` | 81,974 | 26 | Bridge file. All v1 markers with EquCab3 remapping. No metadata rows, no controls section. Adds: Chr_equCab3, MapInfo_equCab3, Strand_equCab3, Ref_equCab3, Alt_equCab3, MAPQ_TopGenomicSeq, MAPQ_Probe. |
| `Equine80select.v2.csv` | 84,319 | 21 | New manifest. Native EquCab3 coords. Manufactured 6/17/2025. 7 metadata rows before header + `[Controls]` tail. Adds: Exp_Clusters (=3), RefStrand (+/-). 2,345 new markers vs v1. |

### Critical Data Quirks

- **AddressA/B_ID format**: v1 and v1.remapped use zero-padded 10-digit strings (`0012758484`); v2 uses unpadded integers (`12758484`). Normalized to unpadded in all joins.
- **MAPQ_Probe = 0 for 51% of v1 markers**: Normal for probes near repetitive regions. NOT filtered — flagged only. The truly problematic subset is `MAPQ_TopGenomicSeq = 0` (~1,361 markers).
- **Chromosome naming mismatch**: v1.remapped uses full NCBI scaffold names (`Un_NW_019641767v1`) for unplaced markers; v2 uses `0` or blank. These are **unmatchable by position** and flagged as `flag_chr_unresolvable`.
- **Chromosome X**: v2 uses `X_NC_009175.3` for some X markers — normalized to `X`.
- **IlmnStrand color names**: Control probe artifacts (`Green`, `Red`, `Blue`, `Purple`, `Black`) in the IlmnStrand column of v1/v2. Filtered out after loading.
- **ID-shift markers**: 0 ID-shift markers. All 81,974 v1 IlmnIDs exist in v2 with the same IlmnID.
- **v1-only IlmnIDs**: 0 — ALL 81,974 v1 assay marker IlmnIDs exist in v2. There are no v1-only rows in the master table.
- **New IlmnIDs in v2**: 2,345 IlmnIDs are unique to v2 (set difference: v2_ids - v1_ids). Of these, 10 were positionally matched to v1 markers at the same Chr:Pos (labeled `positional`); the remaining 2,335 appear as `v2_only` tier (no match of any kind to v1).

---

## Module Structure

```
manifest_compare.py        # Entry point and logging setup
config.py                  # All paths, thresholds, and constants
io_utils.py                # File loading with metadata-skip and validation
normalization.py           # Address ID, chromosome, position, allele utilities
phase1_integration.py      # v1 + v1.remapped join → v1_EquCab3_Master
phase2_comparison.py       # Three-tier matching of v1_master vs v2
phase3_sequence_audit.py   # Probe sequence, strand flip, allele complement
reporting.py               # CSV/TXT/PLINK output writers
visualization.py           # Venn diagram + chromosome distribution chart
```

---

## Pipeline Phases

### Phase 1 — Data Integration (`phase1_integration.py`)

**Entry point**: `build_v1_master(df_v1, df_rem) → DataFrame`

1. Normalizes `AddressA_ID` and `AddressB_ID` (strip leading zeros).
2. Classifies variant types in v1 SNP column (`snp`, `indel`, `control`, `unknown`).
3. Deduplicates `IlmnID + Name` keys defensively before merging.
4. Left-joins v1 + v1.remapped on `[IlmnID, Name]` with `validate='1:1'`.
5. Detects EC3 coordinate conflicts: for the ~6,022 markers already in EquCab3 in v1 (`GenomeBuild='3'`), compares v1 `Chr`/`MapInfo` vs v1.remapped `Chr_equCab3`/`MapInfo_equCab3`. Priority: **v1.remapped always wins** as canonical EquCab3 coords.
6. Creates normalized chromosome column `chr_ec3_norm` and positional key `pos_key_ec3` (`CHR:POS` string).
7. Flags unresolvable NCBI scaffold chromosomes (`flag_chr_unresolvable`).
8. Computes MAPQ quality flags.

**Output columns added**:
- `is_remapped` — has EquCab3 coords in v1.remapped
- `coord_conflict_ec3` — v1 native EC3 coords differ from v1.remapped coords
- `chr_ec3_norm` — normalized chromosome name
- `pos_key_ec3` — `CHR:POS` positional key for joining
- `flag_chr_unresolvable` — NCBI scaffold (Un_NW_*) that can't match v2
- `flag_failed_probe_map` — MAPQ_Probe == 0
- `flag_low_mapq_probe` — MAPQ_Probe < 20 (threshold in config)
- `flag_low_mapq_genomic` — MAPQ_TopGenomicSeq < 20
- `flag_unreliable_position` — MAPQ_TopGenomicSeq == 0
- `variant_type` — snp / indel / control / unknown
- `is_ec3_native` — True for the ~6,022 markers whose `GenomeBuild == '3'` in v1 (had native EquCab3 coords before remapping). Stored as a column so `write_summary` can read it directly rather than re-evaluating `GenomeBuild == '3'`, which proved unreliable after the merge.

---

### Phase 2 — Three-Tier Matching (`phase2_comparison.py`)

**Entry point**: `compare_v1_v2(v1_master, df_v2) → (master_df, tier2_df, tier3_df)`

Markers are matched hierarchically — a marker is only evaluated for a lower tier if it failed all higher tiers.

| Tier | Label | Logic |
|---|---|---|
| 1 | `strict` | Same IlmnID + Name + EquCab3 position (`pos_key_ec3 == pos_key_v2`) |
| 2 | `id_shift` | Same Name + same position, different IlmnID |
| 3 | `positional` | A v2-only IlmnID (in_v1=False) placed at a Chr:Pos already occupied by a v1 marker with a **different** IlmnID. Both sides of the pair are different IlmnIDs from different manifests. The v2 row has `in_v1=False` because that specific IlmnID is absent from v1 — but the genomic **position** is shared with a v1 marker. `in_v1` is an IlmnID flag, not a position flag. **10 markers.** |
| — | `coord_mismatch` | IlmnID in both v1 and v2 (`in_v1=True AND in_v2=True`), but EquCab3 positions differ (v1.remapped liftover vs v2 native). **2,004 markers.** Includes ~648 NCBI-scaffold markers whose position in v1.remapped cannot be reconciled with v2's unplaced notation. |
| — | `v1_only` | In v1 but IlmnID absent from v2. In practice: **0** for Equine80select — all v1 assay IlmnIDs exist in v2. |
| — | `v2_only` | In v2 but IlmnID absent from v1 (new markers). **2,335 markers** (10 of the 2,345 unique-to-v2 IlmnIDs were positionally matched → labeled `positional`). |

Positional key for v2: normalized `Chr` + `MapInfo` (native EquCab3 in v2).
Unplaced markers (`Un:-1`) are excluded from tier-3 positional matching.
Many-to-many positional matches are flagged as `multi_match=True`.

**Output**: outer join of v1_master + df_v2 on `IlmnID`, with `in_v1`, `in_v2`, `match_tier`, `coord_match`, `multi_match` columns. Overlapping column names are suffixed `_v1` / `_v2`.

---

### Phase 3 — Sequence & Allele Audit (`phase3_sequence_audit.py`)

**Entry point**: `run_sequence_audit(master) → DataFrame`

Runs only on common markers (`in_v1 & in_v2`).

**Probe sequence comparison**: case-insensitive string equality on `AlleleA_ProbeSeq` and `AlleleB_ProbeSeq`. `NaN == NaN` is treated as a match (single-channel probes with no AlleleB).

**Strand flip detection**: compares `IlmnStrand_v1` vs `IlmnStrand_v2` (Illumina TOP/BOT convention, present in both manifests with identical semantics).
- `consistent` — same strand designation (both TOP or both BOT)
- `flip` — strand convention changed (TOP↔BOT) — indicates genuine strand-convention flip relevant for PLINK `--flip`
- `na` — one or both values missing or not TOP/BOT

**Why not `Strand_equCab3` vs `RefStrand`**: These measure different things. `Strand_equCab3` is the alignment strand of TopGenomicSeq from the liftover tool (+/-); `RefStrand` in v2 is Illumina's reference strand designation (+/-). With 100% identical probe sequences between v1 and v2, comparing these produced 40,949 spurious flips (~50%), which is biologically implausible. `IlmnStrand` TOP/BOT is the correct field for PLINK strand correction.

**Allele complement validation**: for flipped markers, checks whether v2 alleles are the complement of v1 alleles (e.g., `[A/G]` → `[T/C]`). Returns `same`, `complement`, `mismatch`, or `cannot_assess`.

**Special cases**:
- **Ambiguous SNPs** (A/T or C/G pairs): `is_ambiguous_snp=True` — strand cannot be inferred from alleles alone. Excluded from `plink_flip.txt`.
- **Indels** (`[D/I]`, `[I/D]`): `indel_strand_check_skipped=True` — no allele complement check.
- **Control SNPs** (`DNP (High)`, `DNP (Bgnd)`): classified in `variant_type` but not audited.

---

## Output Files

| File | Description |
|---|---|
| `output/master_comparison_table.csv` | Union of all markers (outer join). All flag columns included. Priority columns written first. |
| `output/discrepancy_report.csv` | Subset of master where any flag is True. Sorted by severity score (unreliable position > EC3 conflict > seq mismatch > strand flip > ID shift > multi-match). |
| `output/summary_statistics.txt` | Human-readable counts for all three phases. Also printed to stdout. |
| `output/plink_flip.txt` | Marker Names to pass to PLINK `--flip`. Confirmed flips only; excludes ambiguous SNPs and indels. |
| `output/plink_ambiguous_exclude.txt` | Marker Names for A/T and C/G SNPs. Recommended for exclusion from PLINK analyses. |
| `output/plink_strand_report.csv` | Full per-marker strand action table: `flip` / `keep` / `exclude_ambiguous` / `indel_skip` / `na`. |
| `output/venn_overlap.png` | 2-circle Venn diagram of IlmnID overlap (v1 vs v2). |
| `output/chromosome_distribution.png` | Grouped bar chart: marker counts per chromosome in v1 vs v2, sorted 1–31, X, Y, Un. |
| `logs/pipeline_YYYYMMDD_HHMMSS.log` | Full timestamped log of pipeline run. |

### Master Table Column Order (Priority Columns)
```
IlmnID, Name,
in_v1, in_v2, is_remapped, is_ec3_native, match_tier, coord_match, multi_match,
Chr_equCab3, MapInfo_equCab3, Strand_equCab3,
chr_v2_norm, pos_key_v2,
coord_conflict_ec3,
MAPQ_Probe, MAPQ_TopGenomicSeq,
flag_failed_probe_map, flag_low_mapq_probe,
flag_low_mapq_genomic, flag_unreliable_position, flag_chr_unresolvable,
sequence_match, AlleleA_seq_match, AlleleB_seq_match,
strand_flip_detected, allele_complement_check,
is_ambiguous_snp, is_indel, indel_strand_check_skipped,
variant_type,
[all remaining v1 and v2 original columns with _v1/_v2 suffixes]
```

---

## Configuration (`config.py`)

All tunable parameters live here. Edit this file to change behaviour without touching pipeline logic.

| Parameter | Default | Effect |
|---|---|---|
| `MAPQ_MIN_PROBE` | `20` | Threshold for `flag_low_mapq_probe` |
| `MAPQ_MIN_GENOMIC` | `20` | Threshold for `flag_low_mapq_genomic` |
| `MAX_COORD_SHIFT_BP` | `0` | Exact position match. Set >0 to allow liftover wobble (e.g., `5` for ±5 bp tolerance). |
| `EXPECTED_V1_ROWS` | `81974` | Row count assertion for v1.csv |
| `EXPECTED_REMAPPED_ROWS` | `81974` | Row count assertion for v1.remapped.csv |
| `EXPECTED_V2_ROWS` | `84319` | Row count assertion for v2.csv |
| `CHR_NORM_MAP` | `{X_NC_009175.3→X, 0→Un, ""→Un}` | Chromosome name normalization map |
| `CONTROL_STRAND_VALUES` | `{Green,Red,Blue,Purple,Black}` | IlmnStrand values that mark control probe artifact rows |
| `INDEL_SNPS` | `{[D/I],[I/D]}` | SNP values classified as indels |
| `CONTROL_SNPS` | `{DNP (High),DNP (Bgnd)}` | SNP values classified as control probes |
| `AMBIGUOUS_ALLELE_PAIRS` | `[{A,T},{C,G}]` | Palindromic pairs that cannot be strand-resolved |

---

## File Loading Logic

Illumina manifest CSVs have two structural hazards:
1. **Metadata preamble**: 7 lines of Illumina metadata before the column header row (v1 and v2 only).
2. **`[Controls]` section**: Appended after the assay rows. Starts with a line `[Controls]`.

Loading strategy in `io_utils.py`:

```python
# v1.csv and v2.csv
# Pre-scan finds [Controls] line → nrows stops reading there
nrows = _find_controls_line(path, skip_metadata=True)
pd.read_csv(path, skiprows=7, nrows=nrows, dtype=str, keep_default_na=False, na_values=[""])

# v1.remapped.csv — no metadata preamble, no [Controls] section
pd.read_csv(path, dtype=str, keep_default_na=False, na_values=[""])
```

**Why NOT `comment='['`**: Although `comment='['` would conveniently skip the `[Controls]` line,
pandas truncates EVERY row at the first `[` character — not just lines that start with `[`. Since
the SNP column contains values like `[A/G]`, every assay row would be truncated at that field,
leaving `GenomeBuild`, `Chr`, `MapInfo`, `AlleleA_ProbeSeq`, and all subsequent columns as `NaN`.
The fix is to pre-scan the file with `_find_controls_line()` and pass `nrows` instead.

File line layout (v1.csv and v2.csv):
```
Lines 0-6  → 7 Illumina metadata lines        (skiprows=7 skips these)
Line  7    → IlmnID,Name,...  (column header)  (pandas uses this as header)
Lines 8+   → assay data rows                   (nrows limits to these)
[Controls] → section marker                    (excluded by nrows)
23 rows    → control probe data                (excluded by nrows)
```

All columns loaded as `str` initially; numeric casting happens in phase modules.

---

## Normalization Functions (`normalization.py`)

| Function | Purpose |
|---|---|
| `normalize_address_id(series)` | Strip leading zeros from AddressA/B_ID: `'0012758484'` → `'12758484'` |
| `normalize_chr(series)` | Map `X_NC_009175.3→X`, `0/""→Un`; NCBI scaffold names kept as-is |
| `make_pos_key(chr, pos)` | Build `'CHR:POS'` string key; unplaced markers → `'Un:-1'` |
| `flag_ncbi_scaffolds(chr)` | Boolean mask: True for `Un_NW_*` or `NW_*` scaffold names |
| `classify_variant(snp)` | Returns `'snp'`, `'indel'`, `'control'`, or `'unknown'` |
| `parse_alleles(snp_str)` | Extract `(A1, A2)` from `'[A/G]'`; returns `None` for indels/controls |
| `is_ambiguous_snp(snp_str)` | True for A/T or C/G pairs |
| `check_allele_complement(snp_v1, snp_v2)` | Returns `'same'`, `'complement'`, `'mismatch'`, or `'cannot_assess'` |

---

## Adding a New Array Version (v3+)

To extend the pipeline for a future `Equine80select.v3.csv`:

1. **`config.py`**: Add `V3_PATH`, `EXPECTED_V3_ROWS`, and output paths for v3 comparison results.
2. **`io_utils.py`**: Add `load_v3()` following the same pattern as `load_v2()`. Check whether v3 has metadata rows and a `[Controls]` tail.
3. **`phase2_comparison.py`**: The `compare_v1_v2()` function is generic enough — duplicate and rename it to `compare_v2_v3(v2_master, df_v3)`. The v2 manifest would serve as the new "master" baseline.
4. **`manifest_compare.py`**: Add the new phase call after the existing Phase 2.
5. **Update `EXPECTED_*_ROWS`** constants to match the new file's actual row counts.

If v3 introduces new columns (like v2 added `RefStrand`), add them to `MASTER_COLS_PRIORITY` in `reporting.py` for clean output ordering.

---

## Adding a New Remapping File

If a new liftover remapping file is generated (e.g., `Equine80select.v2.remapped.csv`):

1. Add its path to `config.py`.
2. Add `load_v2_remapped()` to `io_utils.py`.
3. In `phase1_integration.py`, update `remapped_cols` to match the new file's column names (confirm whether it uses `Chr_equCab3` or a different naming convention).

---

## Known Limitations

- **NCBI scaffold markers**: ~648 markers in v1.remapped have `Un_NW_*` chromosome names that cannot be positionally matched against v2's `0`/blank unplaced representation. These are flagged `flag_chr_unresolvable=True` and will appear as `coord_mismatch` in Phase 2 (IlmnID exists in both, but the scaffold position in v1.remapped is incomparable to v2's unplaced notation).
- **Ambiguous SNPs in strand analysis**: A/T and C/G palindromic SNPs (`is_ambiguous_snp=True`) cannot be strand-resolved from alleles alone. They are excluded from `plink_flip.txt` and listed in `plink_ambiguous_exclude.txt`.
- **MAPQ_Probe = 0**: 42,162 markers (~51%) have this value. This is expected behavior for probes in repetitive regions — the genomic position (from `MAPQ_TopGenomicSeq`) is still reliable. These are flagged but not filtered.
- **Indel strand check**: `[D/I]` and `[I/D]` markers cannot have their strand validated via allele complement logic. Flagged as `indel_strand_check_skipped`.
- **Coordinate shift tolerance**: `MAX_COORD_SHIFT_BP=0` by default (exact match). If future liftover produces small positional shifts (±1–5 bp), increase this threshold in `config.py`.

---

## Verification Checklist

After running the pipeline, verify:

1. Logs show row count assertions pass: `v1: row count OK (81974)`, `v2: row count OK (84319)`.
2. `summary_statistics.txt` shows `New markers in v2 = 2,335` (v2_only tier; 10 additional v2-only IlmnIDs are labeled `positional`).
3. `summary_statistics.txt` shows `ID Shifts = 0` (all v1 IlmnIDs exist in v2 with the same IlmnID).
4. `master_comparison_table.csv` has exactly **84,319 rows** (= all v2 IlmnIDs; no v1-only rows exist since all v1 IlmnIDs are present in v2).
5. `discrepancy_report.csv` contains only rows where at least one flag is True.
6. `plink_flip.txt` entries are NOT in `plink_ambiguous_exclude.txt` (these lists must be disjoint).
7. `venn_overlap.png` overlap (in_v1 & in_v2 = 81,974) = `strict` (79,970) + `coord_mismatch` (2,004) = 81,974. All `positional` rows have `in_v1=False` (they are v2-only IlmnIDs at positions shared with a different v1 IlmnID).

---

## Dependencies

```
pandas >= 1.5
numpy
matplotlib
matplotlib-venn      # optional; Venn diagram is skipped gracefully if absent
```

Python standard library: `logging`, `datetime`, `pathlib`, `os`.
