# Handoff — Building the control-only Neural Organoid Atlas

For an LLM continuing the **raw data → control-group atlas** build. This is the clean
plate: only the scripts that take downloaded GEO/array datasets and produce the
integrated control-only atlas. (Everything from the integration/benchmark/poster phase
has been moved to `scripts/_archive/` — see bottom.)

## Goal
A control-only (no disease/mutant cells) human neural-organoid atlas: per-deposit raw
extraction → uniform QC → HGNC-canonical gene space → concatenation → finalized per-GSM
annotations → scVI-ready preprocessed file.

## Environment
- Python: `/opt/homebrew/Caskroom/miniforge/base/bin/python3.13` (miniforge base; scvi 1.4.1, torch 2.10 MPS, anndata 0.12, scanpy 1.12).
- Some control-fix steps need R (`fix_gse297594_control.py` reads Seurat `.rds`).

## Source-of-truth / config files
| File | Role |
|---|---|
| `data/rebuild_config.tsv` | **driver** — per deposit: `accession, slug, decision, gsms` (control whitelist), `raw_dir`, `loader_type` (13 formats), `is_smartseq2`, `organoid_type`, `multi_lineage`, `samples_cell` |
| `data/raw/<slug>/` | downloaded raw data (per `raw_dir`) |
| `data/reference/hnoca_var_canonical.tsv` | 36,842-gene **HGNC canonical** target space (+ ensembl/length/mt/HVG flags) |
| `data/brain_organoid_GSMannotations.xlsx` (sheet `GSM Annotations`) | finalized per-GSM annotation source of truth (schema: `docs/annotation_schema.md`) |
| `data/manifest.tsv` | output registry (per-deposit path, n_cells, control counts, status) |

## Clean-plate scripts & order
All under `scripts/`. Run from repo root with the miniforge python.

1. **`rebuild_atlas.py`** — per-deposit raw extraction. Dispatches by `loader_type`
   (`10x_mtx_per_gsm`, `cellranger_h5`, `series_level_mtx`, `csv`, `dge_text`,
   `tar_archive`, `alevin`, `streaming_csv`, `rds`, `h5ad`, …) → filter to control GSM
   whitelist → QC (`n_counts≥500, n_genes≥200, pct_mito≤20`) → map genes to HGNC
   canonical → HNOCA obs schema → `data/processed/<slug>.h5ad` + manifest.
   `--dry-run` to validate; `--slugs a b` to subset.
2. **`migrate_to_hnoca_schema.py`** — standardize var (HGNC + flags) + obs harmonization
   (`bio_sample/tech_sample/batch/individual`); raw counts unchanged.
3. **`reuse_from_archive.py`** — reuse SAFE_REUSE deposits from `data/_archive/` (skip
   re-download), re-filter to whitelist, re-attach canonical var.
4. **`fix_gse297594_control.py`**, **`fix_gse296775_control.py`** — pooled + cell-hashed
   (HTO) deposits whose raw MTX mixes control+mutant: demultiplex `condition==WT` from the
   Seurat `.rds` meta.data, then the same QC/HGNC/schema. (Contamination control —
   critical: the naive MTX load silently included mutant cells.) [`fix_gse239542.py`
   exists but is gitignored.]
5. **`concatenate_atlas.py`** — project each deposit onto the 36,842 canonical genes
   (sparse `X@P`; dup symbols summed, missing zero-filled), write temps with identical
   var + unique obs_names, `concat_on_disk` → **`data/atlas_v5_full.h5ad`** (memory-safe).
6. **`inject_finalized_annotations.py`** (+ **`recover_per_cell_annotations.py`** for
   barcode-encoded samples; **`patch_obs_annotations.py`** for obs-only refresh) — join the
   8 finalized fields per cell by GSM, set `annotation_level` (`gsm` vs `deposit`).
   **`membership_diff.py`** decides per-deposit reload-vs-obs-patch.
7. **`preprocess_atlas.py`** — final prep: HVG (3000, batch-aware) + `layers['counts']` +
   lognorm `X` → **`data/atlas_v5_preprocessed.h5ad`** (scVI-ready).
- **`regen_manifest.py`** — rebuild `manifest.tsv` from `data/processed/`.
- **`fetch_hnoca.sh`** — fetch HNOCA reference assets (canonical genes / comparison).

## Data flow
```
data/raw/<slug>/ ─rebuild_atlas(load+QC+HGNC)→ data/processed/<slug>.h5ad
                  (+reuse_from_archive, +fix_*_control demux)
   └─concatenate_atlas(project→36,842, concat_on_disk)→ data/atlas_v5_full.h5ad
   └─inject_finalized_annotations(per-GSM)→ (annotated, control-only)
   └─preprocess_atlas(HVG+counts+lognorm)→ data/atlas_v5_preprocessed.h5ad
```

## Current state (as of 2026-05-25)
- **`data/atlas_v5_full.h5ad`** = 4,079,890 cells × 36,842 genes, **118 deposits**, ~500
  GSMs, control-only, all 8 finalized annotation fields injected. X = raw counts.
- **`data/atlas_v5_preprocessed.h5ad`** = 4,016,176 × 3000 HVG (cell_ranger, batch-aware),
  `layers['counts']` raw, `X` lognorm.
- Atlas metadata: 40 organoid types, 45 protocols, 24% multi-lineage, age 0–276 d (median 70).

## Gotchas / decisions (don't re-discover these)
- **Pooled+hashed deposits** can't be split by genotype from the MTX — must demux from the
  Seurat object (`fix_*_control.py`). GSE297594 was contaminated before this fix.
- **Counts integrity**: some archived deposits had normalized data in the counts slot
  (broken pipeline); verify X is raw integer counts (NB models need true counts).
- **`multi_lineage` is mixed-dtype** in source (`'0'/'1'/'False'/'True'/'No'`) — normalize
  to bool: `{'1','True','Yes'}→multi`.
- **HGNC canonical (36,842)** is the gene space; atlas var_names are HGNC symbols, but Braun
  reference is Ensembl — bridge via `hnoca_var_canonical.tsv`.
- QC thresholds and HVG choice (3000) are tunable in `rebuild_atlas.py` / `preprocess_atlas.py`
  and materially affect downstream integration.

## What's archived (`scripts/_archive/`, the next phase, not needed to build the atlas)
Integration + benchmark + poster work: scVI/scANVI/scPoli training (`train_*`,
`scpoli_pilot`), Braun label transfer (`braun_*`), benchmarks Q1–Q3 (`benchmark_*`),
scIB metrics (`scib_*`), figures (`make_*umap`, `poster_panels`, `figure2_*`,
`fig1*`), stratified mixing, `calibration_ece`, and the shared `atlas_common.py` /
`_provenance.py` utilities they use. See `docs/benchmark_pipeline_review.md` for that phase.

## Pointers
- Plan docs: `docs/rebuild_plan_2026-05-17.md`, `docs/rebuild_plan_2026-05-23.md` (pooled/hashed), `docs/handoff_2026-05-20.md`.
- Annotation schema: `docs/annotation_schema.md`.
