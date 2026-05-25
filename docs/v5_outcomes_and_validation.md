# Atlas v5 — methods, outcomes, and validation

**Date:** 2026-05-25 · supersedes earlier training notes for the v5 build.
Companion to `rebuild_plan_2026-05-23.md` (build plan) and `handoff_2026-05-20.md`.

This document records the **as-run** methods, the **measured outcomes**, and the
**pilot validation checks** (data integrity + cell-type sanity) for atlas v5, so
a reviewer can assess the actual pipeline rather than defaults.

---

## 1. Build summary (what's in the atlas)

- **`data/atlas_v5_full.h5ad`** — 4,079,890 cells × 36,842 canonical HGNC genes,
  118 deposits, ~500 GSMs, 508 `bio_sample` / 512 `tech_sample`.
- Control-only; **excludes HNOCA-overlapping datasets**; pooled/cell-hashed
  deposits use **HTO-demultiplexed controls** (not raw pooled MTX).
- **Annotations:** all 8 finalized fields (`cell_type_origin`, `age_days`,
  `organoid_type`, `protocol`, `unguided`, `multi_lineage`, `vascularized`,
  `slice`) injected **per cell, keyed by GSM**, plus `gsm` and an
  `annotation_level` provenance flag.
  - **100% value-accurate** vs the finalized `GSM Annotations` sheet on the
    3.82M gsm-level cells (every field; the only non-matches were sheet blanks
    rendered as `unknown`).
  - **94.5% of cells `annotation_level=gsm`** (authoritative per-GSM); 5.5%
    `deposit` (series-level deposits where per-cell GSM isn't resolvable).

### GSM coverage vs the finalized sheet (653 GSMs)
- **512 present** in the atlas.
- **105 absent** — accession not built yet (26 Tier-2 deposits: raw-build /
  FASTQ-reprocess / GEO re-download).
- **36 absent** — accession built but pooled/unsplit (per-GSM split needs raw
  reprocessing): GSE98201, GSE252522, GSE197887, GSE243015, GSE181518,
  GSE251684, GSE117512, GSE122342, GSE227640, GSE219245.

---

## 2. Preprocessing (as run)

`preprocess_atlas.py`, on the full atlas:

1. **Drop zero-count cells** — 63,714 cells (1.56%) had 0 counts in the
   canonical gene space (genes lost during projection). Removed before HVG.
   → 4,016,176 cells.
2. **Prune micro-batches** `< 50` cells — **no-op on this data** (smallest
   `bio_sample` = 162 cells; 0 batches / 0 cells removed). Kept as a guard.
3. **Normalize** `normalize_total(1e4)` + `log1p` → `X`; raw counts preserved in
   `layers['counts']`.
4. **HVG = `cell_ranger` flavor**, `n_top_genes=3000`, `batch_key=bio_sample`.

**HVG flavor rationale (why not the defaults):**
- `pearson_residuals` (the earlier directive) **stalled** — single-threaded
  per-batch residuals over 508 batches ran >1 h with no progress.
- `seurat_v3` **crashed** — singular-matrix (`reciprocal condition number
  3.7e-15`) in the per-batch loess on a degenerate-variance batch.
- `cell_ranger` uses **binned dispersion lookup** (no loess), so it cannot throw
  singular-matrix errors and is fast. Batch granularity is `bio_sample`
  (per-organoid); selection is robust (saturation analysis: 2–4k genes all valid).

Output `data/processed/atlas_v5_preprocessed.h5ad` — 4,016,176 × 3000.

---

## 3. scVI integration (as run) + the epoch experiment

`train_scvi.py`: `n_layers=2`, `n_latent=30`, `gene_likelihood='nb'`,
`layer='counts'`, **`batch_key='tech_sample'`** (505 sequencing libraries — the
real technical batch), Apple MPS.

Two batch keys are used deliberately: **`bio_sample`** (per-organoid) for HVG;
**`tech_sample`** (library) for scVI integration.

### Epoch experiment — 15 vs 100 (addresses "undertraining at ~4M cells")
We ran both and compared. **15 epochs is adequate; more epochs do not improve the
data fit.**

| metric | epoch 15 | epoch 100 |
|---|---|---|
| reconstruction_loss_train | 1124.13 | **1124.95** (flat / slightly worse) |
| kl_local_train | 76.54 | 39.98 |
| elbo_train | 1200.67 | 1164.92 |
| validation_loss | — | bottoms ≈ epoch 4 (1119), **drifts up to 1126.75** |

**Interpretation:** the entire epoch-15→100 ELBO improvement comes from the **KL
term shrinking** (latent tightening toward the prior), while **reconstruction was
already converged at epoch 15**. Validation reconstruction is best early and
degrades — so 100 epochs is mildly **over-regularized**, not underfit. (Early
stopping, patience 10, never fired because it monitors ELBO, which KL keeps
pulling down.)

### Integration outcome (why we don't chase epochs)
100ep (batch 512, ~784k steps, 2.2 h) vs 15ep (batch 2048, ~29k steps, 7 min):
batch mixing improved only **~10% relative**, **biology dropped ~7%**, net
batch/bio ratio **1.82 → 1.91**. The residual batch structure is **structural,
not an optimization gap**: `tech_sample` has 505 categories **confounded with real
biology** (each dataset is its own protocol / age / region). Global kNN
same-neighbor-fraction conflates batch with biology.

**Real levers (not epochs):** coarser `batch_key`, **scANVI** (semi-supervised
with transferred labels), or **cell-type-stratified** iLISI/kBET. The definitive
mixing litmus test is deferred to after Braun-2023 label transfer (so mixing can
be measured *within* cell type). See `rebuild_plan_2026-05-23.md`.

Final outputs: `data/scvi_model_v5_full/`, `data/scvi_latent_v5_full.h5ad`
(4,016,176 × 30; latent stored in `X` — copy to `obsm['X_scvi']` before
`sc.pp.neighbors`).

---

## 4. Pilot validation — data integrity (preprocessed matrix)

All checks on `atlas_v5_preprocessed.h5ad` passed:

| Check | Result |
|---|---|
| File loads, shape | (4,016,176, 3000), `layers=['counts']` |
| `counts` dtype / range | float32, min 0, **max 55,140**, **all integer**, no NaN/Inf |
| `X` is log-norm, distinct | max **8.85**, min 0, `X != counts` |
| Degenerate cells/genes | **0** zero-count cells, **0** zero-count genes |
| obs/var indices | unique; `tech_sample` 505 batches, no NaN, min 162 cells/batch |
| Training-loss signal | ELBO descends monotonically, **no NaN/Inf** |

The clean monotonic ELBO under an NB likelihood is strong evidence the counts
layer is genuine raw integers (the prior "normalized-not-counts" loader bug would
manifest as NaN/Inf or a non-decreasing loss).

---

## 5. Pilot validation — cell-type marker sanity

### Atlas-wide (% cells > 0, raw counts)
| marker | cell type | % cells > 0 |
|---|---|---|
| SOX2 | neural progenitor / radial glia | 38.7% |
| VIM | radial glia / progenitor | ~73% |
| DCX | immature neuron | 50.4% |
| MKI67 | proliferating | 9.5% |
| RBFOX3 (NeuN) | mature neuron | 13.0% |
| GFAP | astrocyte | 7.3% |

This composition — abundant progenitors (SOX2/VIM) and immature neurons (DCX), a
modest cycling fraction (MKI67), fewer mature neurons/astrocytes — is exactly
what's expected for a neurogenic organoid atlas.

### Per-deposit marker check (catches gene mis-alignment across 118 deposits)
Per-deposit % cells > 0, min / median / max across deposits:

| marker | min | median | max | deposits at 0% |
|---|---|---|---|---|
| SOX2 | 0.0 | 35.7 | 90.6 | 4 |
| DCX | 0.2 | 55.3 | 95.7 | **0** |
| MKI67 | 0.5 | 7.0 | 40.3 | **0** |
| GFAP | 0.0 | 1.3 | 72.8 | 2 |
| RBFOX3 | 0.1 | 9.9 | 50.9 | **0** |
| VIM | 1.9 | 73.2 | 99.9 | **0** |

**Conclusion: cross-deposit gene mis-alignment is ruled out.** Pan-markers
(DCX, MKI67, RBFOX3, VIM) are detected in **all 118 deposits**. A mis-aligned
gene column would blink to 0% in chunks of deposits — it does not. The few zeros
are biology (GFAP absent in young organoids; SOX2 absent in mature/non-cortical).

### Outlier flagged for the integration review: `gse290048_pineal`
Every *cortical* marker is near-zero — but this is **biology, not corruption**:
- It's pineal tissue (photoreceptor-like), not cortical: **PDC 35% (mean 5.70),
  RORB 49%, OTX2 11%, AANAT 12%, NEUROD1 12%** — pineal/photoreceptor genes map
  and express correctly (3 independent markers landing right ≠ chance).
- Normal depth (median 1,796 genes / 3,394 counts, 10x 3′ v3).
- **Caveat:** its top-expressed genes are unusual (`ENSG…` IDs, lncRNAs,
  `TNNI2`, `FOXD4`; **no ribosomal/mito dominance**) — suggests upstream
  processing differences (ribo/mito stripped and/or lncRNA-rich reference), not
  mis-mapping. **Watch item:** if it forms a lone outlier cluster in the scVI
  UMAP, review or exclude.

---

## 6. Known limitations / next steps

- **Tier-2 coverage:** 141 finalized GSMs not yet in the atlas (105 unbuilt
  accessions; 36 in pooled/unsplit deposits needing per-GSM raw splitting).
- **`seurat_hto_control` loader:** pooled/cell-hashed control extraction is
  currently standalone (`fix_gse297594_control.py`, `fix_gse296775_control.py`),
  **not** in `rebuild_atlas.py`'s loader dispatch — a future rebuild could
  accidentally route a hashed study through the pooled-MTX loader. Formalizing
  this loader is the top correctness item (already noted in the build plan).
- **Integration metric:** adopt cell-type-stratified mixing (post label transfer)
  rather than global kNN, which confounds batch with biology here.
- **Reproducibility:** scripts use absolute paths and a hard-coded interpreter —
  parameterize before external reuse.
