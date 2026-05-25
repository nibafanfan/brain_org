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
Head-to-head, **same-dataset/same-batch neighbor fraction** (kNN mixing proxy;
lower = better mixing for batch keys, higher = better conservation for biology):

| key (n categories) | old (2048 / 15ep, ~29k steps) | new (512 / 100ep, ~784k steps) | random baseline | Δ |
|---|---|---|---|---|
| tech_sample (505) | 0.304 | 0.269 | 0.0054 | −12% |
| bio_sample (501) | 0.305 | 0.270 | 0.0054 | −11% |
| dataset_slug (113) | 0.464 | 0.418 | 0.017 | −10% |
| organoid_type (33, biology) | 0.553 | 0.513 | 0.156 | −7% |

**Biology/batch ratio** (integration quality): old 0.553/0.304 = **1.82** → new
0.513/0.269 = **1.91** (~5% better). So **27× more compute (2.2 h vs 7 min) bought
~10% relative batch-mixing gain while biology dropped ~7%** — the latent didn't
restructure, everything just contracted slightly. Underfitting was a *minor*
factor; the 7-min run was already near the achievable mixing for this setup.

The residual batch structure is **structural, not an optimization gap**:
`tech_sample` has 505 categories **confounded with real biology** (each dataset
carries its own protocol / age / region), and the same-neighbor metric conflates
the two. scVI can't (and shouldn't) erase "same dataset = partly same biology."
**Keep the converged 100-epoch model** and move on; real levers for mixing are a
coarser `batch_key`, scANVI, or per-cell-type iLISI/kBET — not more epochs.
(An independent Codex methods review reached the same conclusion.)

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

## 6. Label transfer pilot (Braun 2023 → scANVI)

`braun_label_transfer.py` pilot on a 200k-cell subsample, ~7.4 min end-to-end.
**CellClass** via scANVI surgery on the scVI latent; **Region** via kNN.

### CellClass (scANVI) — biologically sane
| Class | % | |
|---|---|---|
| Radial glia | 44.7% | progenitors |
| Neuron | 34.0% | |
| Glioblast | 12.3% | glia |
| Neuroblast | 7.1% | |
| Neuronal IPC | 0.9% | |
| Fibroblast | 0.8% | |
| Immune / Vascular / Oligo / Erythrocyte | 0.0% | ⚠ see caveats |

**~99% neural lineage** (radial glia + neuron + glioblast + neuroblast + IPC) —
exactly what brain organoids should be. **Mean confidence 0.943; 89% of cells
> 0.8.** Transfer is confident and plausible.

### Region (kNN) — plausible, appropriately softer
Telencephalon 31.6% + Forebrain 16.1% ≈ **48% forebrain** (cortical protocols are
most common), Midbrain 20.7%, Diencephalon 15.1%, Medulla 12.1%. **Mean confidence
0.550** — much lower than CellClass, as expected: organoid regional patterning is
imprecise, so regional identity transfers noisily (honest, not broken).

### Caveats to resolve at full scale
1. **Immune & Vascular = 0.0%** — investigated (see §6b). **Not biology and not
   missing features**: the cell types are present (~1–3% in the atlas) and their
   markers were already in the transfer panel. The cause is **class imbalance +
   pilot subsampling + argmax** (Immune/Vascular are only 0.5–0.7% of the Braun
   reference too). The proper test is the **full, un-subsampled run**; if still
   suppressed, use class-balanced sampling / soft-probability thresholds.
2. **Validation crosstab was weak** — the atlas `cell_type` column is almost all
   `unknown`, giving no signal. Before the full run, confirm the **finalized
   per-cell annotations** (`cell_type_origin`, `organoid_type`, …,
   `annotation_level`, `gsm`) are present/populated in `atlas_v5_full.h5ad` so
   transfers can be cross-checked against real organoid labels.

### Pre-full-run checklist
- Confirm finalized annotation fields are populated (not `unknown`) for the
  cross-check; the injection field-map and `annotation_level` semantics apply.
- Full-atlas **rare-class audit** (immune/vascular/oligo presence).
- Benchmark **stratified by `annotation_level`** (`gsm` vs `deposit`) so coarse
  fallback cells don't bias results.

## 6b. Rare-lineage representation — HVG depth vs forced marker panel

The pilot's "0% Immune/Vascular" prompted an HVG sensitivity check (3k vs 5k) and
a marker-prevalence audit. **Key finding: those lineages exist; the features
didn't.**

### Marker expression prevalence (300k-cell sample, % cells > 0)
- microglia: P2RY12 **1.0%**, TMEM119 0.64%, CSF1R 0.75%, C1QB 0.17%, CX3CR1 0.09%
- myeloid/immune: AIF1 **1.35%**, PTPRC 0.79%
- endothelial: CLDN5 **2.95%**, KDR 0.76%, FLT1 0.67%, PECAM1 0.38%, VWF 0.23%

So there are **thousands of real microglia (~1%) and endothelial (~1–3%) cells**.

### 3k vs 5k HVG (cell_ranger, `batch_key=bio_sample`)
- 5k is a near-**superset** of 3k (2997/3000 retained) + 2003 new genes → low-risk.
- 5k **adds endothelial PECAM1, KDR**; 3k already had CLDN5, FLT1.
- 5k adds **zero microglia-specific markers** — CX3CR1/P2RY12/TMEM119/C1Q/CSF1R are
  expressed but too **low-variance** to be HVGs at *either* depth.
- Oligodendrocyte (OLIG1/2, SOX10, MBP, PLP1, PDGFRA) and mural (PDGFRB, RGS5,
  ACTA2) programs are already fully covered at 3k.

### Conclusion — it's class imbalance, not features (revised after testing)
We tested the forced-panel hypothesis and it was **wrong for the transfer**:
- The transfer's panel comes from **Braun's** HVGs (∩ atlas), and **42 of 47**
  curated markers were *already* selected — forcing the other 5 changed nothing
  for Immune/Vascular (still 0%) and shifted the composition (Glioblast 12%→34%),
  i.e. run-to-run instability on the 200k subsample.
- The real cause: **the rare classes are rare in the Braun reference too** —
  Immune **0.5%** (8,102 cells), Vascular **0.7%** (11,741), Oligo 0.4%,
  Erythrocyte 0.5%. In `--pilot` mode the reference is subsampled to 100k, leaving
  only ~500 Immune / ~700 Vascular training cells, and scANVI's **argmax**
  prediction essentially never selects a 0.5%-prior class.

So the 0% is **class imbalance + pilot subsampling + argmax**, not missing
features. Real levers (in priority):
1. **Full run, no subsampling** — all 1.67M Braun cells give the rare classes
   8–12k training examples each (vs ~500 in the pilot). This is the proper test
   and has **not yet been run**.
2. If still suppressed at full scale: **class-balanced scANVI sampling / class
   weights**, or **soft-probability thresholds** for rare classes (assign Immune
   if P(Immune) exceeds a per-class threshold, not just argmax), or an independent
   **marker-score** assignment for immune/endothelial cells.

The `RARE_PANEL` forcing in `braun_label_transfer.py` is kept only as cheap
insurance (guarantees the markers are in the panel even when a pilot subsamples
Braun); it is **not** the fix. 5k is still adopted as the low-risk preprocessing
baseline (near-superset of 3k, + endothelial PECAM1/KDR), but HVG depth was never
the rare-lineage bottleneck.

**Artifacts:** 3k = `atlas_v5_preprocessed.h5ad`; 5k =
`atlas_v5_5k_preprocessed.h5ad` (+ `*_hvg.tsv` each). Prior trained 3k model/latent
preserved untouched.

## 7. Known limitations / next steps

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
