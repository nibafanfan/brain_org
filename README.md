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

## Current state (atlas v5)

4.08M cells · 118 deposits · ~500 GSMs · all 8 finalized annotation fields per
cell (94.5% at `gsm` granularity) · scVI latent trained (clean monotonic ELBO).
Excludes HNOCA-overlapping datasets and uses HTO-demultiplexed controls for
pooled deposits. The ~800 GB data tree and annotation workbooks are gitignored.

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

Feedback wanted on the **research plan** and **model training**, specifically:

1. **HVG strategy** — `cell_ranger` flavor, n=3000, `batch_key=bio_sample`:
   appropriate at this scale and batch granularity? (We moved off
   `pearson_residuals`/`seurat_v3`, which stalled or hit singular-matrix loess
   on degenerate per-batch variance.)
2. **scVI setup** — `n_latent=30`, ~15 epochs, NB likelihood, and the
   `bio_sample`-for-HVG vs `tech_sample`-for-scVI batch-key split.
3. **Contamination control** — demultiplexing pooled/cell-hashed deposits from
   the processed Seurat object rather than raw counts (`fix_gse297594_control.py`
   as the template). Is this robust / generalizable?
4. **Annotation provenance** — the per-cell vs deposit-level (`annotation_level`)
   fallback and its implications for label transfer and benchmarking.
5. **Reproducibility/portability** — scripts currently use absolute paths and a
   hard-coded interpreter; what should be parameterized?

Please flag methodological risks, statistical concerns, simpler/more standard
alternatives, and concrete next steps. Entry points: `docs/rebuild_plan_2026-05-23.md`
→ `scripts/rebuild_atlas.py` → `scripts/preprocess_atlas.py` → `scripts/train_scvi.py`.
