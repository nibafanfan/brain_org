# brain_organoid

Benchmark of multi-lineage human brain organoid protocols, extending the
Human Neural Organoid Cell Atlas (HNOCA) framework to multi-lineage
differentiation methods published in 2024–2025.

See [`docs/proposal.md`](docs/proposal.md) for the full project proposal,
[`docs/rebuild_plan_2026-05-23.md`](docs/rebuild_plan_2026-05-23.md) for the
current build plan, and [`docs/handoff_2026-05-20.md`](docs/handoff_2026-05-20.md)
for the latest state.

## Pipeline

A control-only atlas is assembled per-deposit, harmonized to a canonical gene
space, annotated per cell, integrated with scVI, then labelled. Run order and
entry points (`scripts/`):

| Step | Script | What it does |
|------|--------|--------------|
| 1. Build | `rebuild_atlas.py` (+ `migrate_to_hnoca_schema.py`) | Per-deposit loaders → uniform QC (`n_counts≥500, n_genes≥200, pct_mito≤20`) → project onto canonical **36,842-gene HGNC** space → HNOCA-style `obs` schema → `data/processed/<slug>.h5ad` + manifest |
| 1b. Reuse | `reuse_from_archive.py` | Reuse previously-validated deposits, filtered to a finalized GSM whitelist (`rebuild_config.tsv`) |
| 1c. Contamination | `fix_gse297594_control.py`, `fix_gse296775_control.py` | Pooled/cell-hashed deposits: demux genotype/condition from the Seurat object (raw pooled MTX can't be split) — extract true controls only |
| 2. Concatenate | `concatenate_atlas.py` | Memory-safe `concat_on_disk` onto the canonical genes → `data/atlas_v5_full.h5ad` |
| 3. Annotate | `inject_finalized_annotations.py`, `recover_per_cell_annotations.py`, `patch_obs_annotations.py`, `membership_diff.py` | Inject the 8 finalized annotation fields per cell, keyed by GSM, with an `annotation_level` flag (`gsm` = authoritative per-GSM, `deposit` = coarser fallback) |
| 4. Preprocess | `preprocess_atlas.py` | Drop zero-count cells, prune micro-batches, normalize+log1p, **`cell_ranger` HVG** (n=3000, `batch_key=bio_sample`); raw counts kept in `layers['counts']` |
| 5. Integrate | `train_scvi.py` | scVI (NB, `n_latent=30`, `batch_key=tech_sample`, ~15 epochs, Apple MPS) → latent + model |
| 6. Evaluate | `eval_integration.py`, `braun_label_transfer.py` | Integration metrics; Braun-2023 label transfer via scANVI |

Two batch keys are used deliberately: **`bio_sample`** (per-organoid) for HVG
selection, **`tech_sample`** (sequencing library — the real technical batch) for
scVI integration.

## Current state & outcomes (atlas v5)

**→ Full methods, outcomes, and pilot validation: [`docs/v5_outcomes_and_validation.md`](docs/v5_outcomes_and_validation.md)**

- **Atlas:** 4.08M cells · 118 deposits · ~500 GSMs · canonical 36,842-gene space.
  Excludes HNOCA-overlapping datasets; pooled/hashed deposits use HTO-demuxed
  controls.
- **Annotations:** all 8 finalized fields per cell, **100% value-accurate** vs the
  source sheet; **94.5% at `gsm` granularity** (`annotation_level` flag marks the
  5.5% deposit-level fallback). GSM coverage 512/653 (rest are Tier-2 unbuilt or
  pooled/unsplit).
- **HVG:** `cell_ranger`, 3000 genes, `batch_key=bio_sample` — chosen because
  `pearson_residuals` stalled and `seurat_v3` hit a singular-matrix loess at this
  batch granularity.
- **scVI:** NB, `n_latent=30`, `batch_key=tech_sample`. **Epoch experiment (15 vs
  100):** reconstruction converged by epoch 15 (1124.1 → flat); the extra ELBO
  gain to epoch 100 is *pure KL shrinkage* and validation reconstruction actually
  drifts up — so 15 epochs is adequate, **not undertrained**. Longer training
  moved batch/bio mixing only 1.82→1.91; residual batch structure is **structural**
  (505 libraries confounded with biology), addressed by scANVI / coarser batch key,
  not more epochs.
- **Validation:** preprocessed matrix passes integrity checks (integer counts max
  55,140, log-norm `X` max 8.85, 0 degenerate cells/genes, clean monotonic ELBO).
  Cell-type markers are biologically coherent atlas-wide (SOX2 39%, VIM 73%,
  DCX 50%, MKI67 9.5%, RBFOX3 13%, GFAP 7.3%) and **per-deposit pan-markers appear
  in all 118 deposits → no cross-deposit gene mis-alignment**. One outlier
  (`gse290048_pineal`) is biologically coherent pineal tissue, flagged as a
  watch-item.
- **Label transfer (Braun 2023 → scANVI) pilot:** ~99% neural lineage (Radial glia
  45%, Neuron 34%, Glioblast 12%, Neuroblast 7%), CellClass mean confidence 0.943;
  Region softer (~48% forebrain, conf 0.55) as expected. Watch items: immune/
  vascular at 0% in the 200k pilot (resolve at full scale), and the atlas
  `cell_type` cross-check column is mostly `unknown`.

The ~800 GB data tree and annotation workbooks are gitignored.

## Layout

```
configs/                   # YAML/JSON configs for integration & benchmarking runs
data/
  raw/                     # untouched downloads (gitignored)
  processed/               # harmonized AnnData (gitignored)
  external/                # third-party data (gitignored)
  reference/               # primary brain reference atlases (gitignored)
docs/                      # proposal, design notes
notebooks/                 # exploratory analyses
scripts/                   # CLI entrypoints, data fetchers
src/brain_organoid/
  integration/             # scVI / scANVI / SAE training & embeddings
  benchmarks/              # mapping accuracy, alignment, coverage metrics
  models/                  # model definitions
  utils/                   # I/O, QC, plotting helpers
tests/
results/                   # output tables, embeddings (gitignored)
figures/                   # plots (gitignored)
```

## Setup

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
```

## Data

Primary brain references and HNOCA artifacts are downloaded by scripts under
`scripts/`. URLs and DOIs are listed in [`docs/data_sources.md`](docs/data_sources.md).

## For reviewers

Feedback wanted on the **research plan** and **model training**. **Please read
[`docs/v5_outcomes_and_validation.md`](docs/v5_outcomes_and_validation.md) first**
— it has the as-run methods, measured outcomes, and pilot validation (data
integrity + per-deposit cell-type marker checks), including evidence that
directly addresses common a-priori concerns:

1. **scVI epochs** — we **ran 15 vs 100 epochs**: reconstruction loss is flat
   (1124.1 → 1124.9), the extra ELBO gain is pure KL shrinkage, and validation
   reconstruction drifts up. So 15 epochs is adequate, not undertrained; residual
   batch structure is structural (505 libraries ≈ biology), not an optimization
   gap. Is our reasoning sound, and is scANVI / a coarser batch key the right next
   lever (vs more epochs)?
2. **HVG strategy** — `cell_ranger`, n=3000, `batch_key=bio_sample` (after
   `pearson_residuals` stalled and `seurat_v3` hit singular-matrix loess). The
   `MIN_BATCH=50` prune was a **no-op** here (smallest batch 162 cells) — keep or
   drop the guard? Better flavor for ~4M cells / 505 batches?
3. **Contamination control** — demuxing pooled/cell-hashed deposits from the
   Seurat object, not raw counts (`fix_gse297594_control.py`). We agree it should
   become a first-class `seurat_hto_control` loader in `rebuild_atlas.py` dispatch
   — review the approach and the generalization.
4. **Annotation provenance** — per-cell vs deposit-level (`annotation_level`)
   fallback (94.5% gsm-level) and its impact on label transfer / benchmarking.
5. **Reproducibility/portability** — absolute paths + hard-coded interpreter:
   what to parameterize?

Please flag methodological risks, statistical concerns, simpler/standard
alternatives, and concrete next steps. Entry points:
`docs/v5_outcomes_and_validation.md` · `docs/rebuild_plan_2026-05-23.md` →
`scripts/rebuild_atlas.py` → `scripts/preprocess_atlas.py` → `scripts/train_scvi.py`.

## Reproducibility via shared config

Use `configs/atlas_v5.yaml` and optional CLI overrides (`--config`, `--root`, `--in`, `--out`, `--out-tag`).

### Reproduce figures

```bash
python scripts/eval_integration.py --config configs/atlas_v5.yaml --out-tag v5_full
```

### Reproduce benchmarks

```bash
python scripts/preprocess_atlas.py --config configs/atlas_v5.yaml --out-tag v5_full
python scripts/train_scvi.py --config configs/atlas_v5.yaml --out-tag v5_full
python scripts/braun_label_transfer.py --config configs/atlas_v5.yaml --out-tag v5_full
```
