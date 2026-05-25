# Atlas rebuild plan — control-only, from raw on disk

**Date:** 2026-05-17

## Source of truth

The **per-deposit decision** for what to include is in:

```
data/atlas_rebuild_target_list.tsv
```

Each row is one GSE with a `decision` column:

| Decision | Action | Source |
|---|---|---|
| `INCLUDE_EXPLICIT` | Load only the GSMs in the `gsms` column | xlsx (2) lists specific GSMs |
| `INCLUDE_ALL_SAMPLES` | Load every GSM from the deposit's raw dir | xlsx (2) Samples cell is exactly "all samples" |
| `SKIP_QUALIFIED` | Skip — deferred | xlsx (2) has qualified text ("All samples (VEH ONLY)" etc.) |
| `SKIP_DESCRIPTIVE` | Skip — deferred | xlsx (2) has paper-specific descriptive text |
| `SKIP_EMPTY` | Skip — deferred | xlsx (2) Samples cell is empty |
| `SKIP_NO_ANNOTATION` | Skip — deferred | Not in xlsx (2) (added after last annotation pass) |

**This plan does NOT enumerate which GSEs are which.** Edit the TSV if a decision needs to change; the rebuild script picks up whatever the current TSV says.

## Workflow

### 1. Resolve raw dirs + loader types (~30 min)

Extend `atlas_rebuild_target_list.tsv` with two derived columns:

- `raw_dir` — locate the deposit's raw files. Check, in order: `data/raw/<slug>/`, `data/_blocked/<slug>/`, `data/_pending_r_interop/<slug>/`, `data/_pending_fastq_reprocess/<slug>/`.
- `loader_type` — one of the format families in `docs/loader_format_reference.md` (`10x_mtx`, `dge_text`, `csv`, `alevin`, `rds`, `h5ad`, `cellranger_h5`, `streaming_csv`, `series_level_mtx`).

Save as `data/rebuild_config.tsv`. Only INCLUDE rows need these columns populated; SKIP rows can leave them blank.

### 2. Archive current atlases (~5 min)

```
data/processed/      → data/_archive/processed_pre_rebuild_2026-05-17/
data/processed_v2/   → data/_archive/processed_v2_pre_rebuild_2026-05-17/
data/manifest.tsv    → data/_archive/manifest_pre_rebuild.tsv
data/manifest_v2.tsv → data/_archive/manifest_v2_pre_rebuild.tsv
```

`data/raw/` and `data/reference/` stay untouched.

### 3. Write `scripts/rebuild_atlas.py` (~1-2 hours code)

For each INCLUDE row in `data/rebuild_config.tsv`:

1. Dispatch to the loader template per `loader_type` (see `docs/loader_format_reference.md`).
2. Apply GSM filter: keep cells whose `obs['gsm']` is in `gsm_filter` (or all if `INCLUDE_ALL_SAMPLES`).
3. Apply uniform QC: `n_counts ≥ 500`, `n_genes ≥ 200`, `pct_mito ≤ 20%` (MT detection via HNOCA's 13 Ensembl IDs in `data/reference/hnoca_var_canonical.tsv`).
4. Project to HNOCA-compatible schema (reuse `scripts/migrate_to_hnoca_schema.py`'s gene mapping + obs harmonization + `sanitize_for_h5()`).
5. Set `is_control=True` on every cell (the filter already excludes non-controls).
6. Write `data/processed/<slug>.h5ad` and append to a new `data/manifest.tsv`.

### 4. Dry-run on 3 deposits (~15 min)

Pick one row per major loader_type to validate end-to-end before scaling. Verify cell counts match the GSM list, schema is HNOCA-compatible, `ad.concat()` works across any pair.

### 5. Full rebuild (~3-6 hours runtime)

Run over all INCLUDE rows in background with progress monitor.

### 6. Validate (~30 min)

- Total cell / GSM / GSE counts match TSV expectations.
- 3 random concat tests across deposits.
- Per-deposit pct_mito sanity (median 2–8% for healthy organoids).
- Manifest row count = number of INCLUDE rows that loaded successfully.

### 7. Update docs (~15 min)

- `docs/atlas_overview.md`: refresh headline numbers; drop v2 references.
- `docs/atlas_compilation_notes.md`: append rebuild section.
- `docs/data_sources.md`: update totals.

## QC policy (uniform across rebuild)

```
n_counts ≥ 500
n_genes  ≥ 200
pct_mito ≤ 20%
```

Smart-seq2 deposits (`gse195692`, `gse185052`) additionally get `layers['counts_lengthnorm']` for cross-tech comparability.

## What this plan does NOT do

- Re-download any raw data (everything on disk).
- Touch `data/raw/` or `data/reference/`.
- Build the integrated / scPoli atlas (downstream step after rebuild).
- Train any DL models.
- Decide what to do with the 34 SKIP rows — those wait on `data/brain_organoid (2).xlsx` updates from the colleague.

## To add a deferred deposit back later

1. Open `data/brain_organoid (2).xlsx`, find the row, fill in the Samples cell with either explicit GSMs or `"all samples"`.
2. Re-run the script that produces `atlas_rebuild_target_list.tsv` (the same one that built the current TSV).
3. Run the rebuild script with `--slugs <slug>` (or just for the changed rows).
4. Append to the manifest.

## Pointers

| What | Where |
|---|---|
| Per-deposit decisions | `data/atlas_rebuild_target_list.tsv` |
| Authoritative control-GSM annotation | `data/brain_organoid (2).xlsx`, sheet "Paper Datasets", col M |
| Canonical gene set (36,842 HGNC ↔ Ensembl + MT flag + gene length) | `data/reference/hnoca_var_canonical.tsv` |
| Format-family loader recipes | `docs/loader_format_reference.md` |
| Existing converter templates | `/tmp/conv_*.py`, `/tmp/rds_to_mtx.R`, `/tmp/run_a[23].py`, `scripts/migrate_to_hnoca_schema.py` |
