# Brain-Organoid Atlas → Primary-Reference Benchmark: Pipeline, Outputs, Challenges & Gaps

**Status:** working end-to-end as of 2026-05-25. Written for a code reviewer (e.g. Codex)
to critique the implementation and propose improvements. It is deliberately **honest and
self-critical**: every known weakness, confound, and shortcut is called out. Where a claim
is shaky, it says so.

**TL;DR of the science:** we embedded a 4.08M-cell brain-organoid atlas with scVI, transferred
primary fetal-brain labels (Braun 2023) onto it via scArches+kNN, and answered three benchmark
questions about multi-lineage vs single-lineage protocols. The science conclusions are
*repertoire/correspondence-based* and hold up; the *integration quality* and several *methods
choices* are the weak points and are where review is most valuable.

---

## 1. Environment & invocation

- Python: `/opt/homebrew/Caskroom/miniforge/base/bin/python3.13` (miniforge base env).
- Libs: `scvi-tools 1.4.1`, `torch 2.10.0` (MPS backend), `anndata 0.12.10`, `scanpy 1.12`,
  `scikit-learn 1.8.0`. Hardware: Apple M2 Ultra, 192 GB RAM, MPS GPU.
- All scripts are standalone, run from repo root `/Users/eg/brain_organoid`, write logs to
  `data/logs/<name>.log`. No config system; **paths and hyperparameters are hardcoded as module
  constants or argparse defaults** (review item — see §7).

---

## 2. Data inputs

| Input | Shape | Notes |
|---|---|---|
| `data/processed/atlas_v5_preprocessed.h5ad` | 4,016,176 × 3000 | HVG-subset, `layers['counts']`=raw, `X`=lognorm, `batch_key=tech_sample` (505). scVI training input. |
| `data/atlas_v5_full.h5ad` | 4,079,890 × 36,842 | full-gene atlas, **`X`=raw counts** (CSR), var_names = **HGNC symbols**. Query source for label transfer + benchmarks. |
| `data/raw/braun_2023/braun_all.h5ad` | 1,665,937 × 58,226 | primary fetal brain reference, **`X`=raw counts**, var_names = **Ensembl**. Labels: `CellClass`(12), `Region`(10), `Subregion`(18), `Age`(PCW), `donor_id`(26). |
| `data/reference/hnoca_var_canonical.tsv` | 36,842 rows | gene bridge: `hgnc_symbol` ↔ `ensembl` (atlas var_names == hgnc_symbol exactly). |

Note the two atlas cell counts differ (4,016,176 preprocessed vs 4,079,890 full); the scVI
latent is a subset of the full atlas. They share `obs_names` (verified 100% overlap), so joins
are by barcode.

---

## 3. Pipeline stages

### Stage A — scVI embedding  (`scripts/train_scvi.py`, 126 ln)
- Trains scVI on `atlas_v5_preprocessed` (`layer=counts`, `batch_key=tech_sample`, `n_layers=2`,
  `n_latent=30`, `gene_likelihood=nb`).
- Final model used: `--batch-size 512 --max-epochs 100 --early-stopping` →
  `data/scvi_model_v5_full/`, `data/scvi_latent_v5_full.h5ad` (4,016,176 × 30, latent stored in `X`).
- Edits made this session: MPS-fallback env var set **before** torch import; `EpochLogger`
  callback for clean per-epoch ELBO; `--batch-size`/`--num-workers` args.
- **Known issues:**
  - `--early-stopping` monitors `elbo_validation`, which keeps creeping down via the KL term
    long after reconstruction plateaus → **early stopping never actually fired** (ran the full
    100 epochs). Needs a `min_delta` / different monitor.
  - `--num-workers > 0` is unsafe on macOS (`spawn` pickles the in-memory dataset per worker →
    OOM); defaulted to 0. Documented in the arg help.
  - Latent saved in `X`, not `obsm['X_scvi']`; every downstream consumer re-copies `X→obsm`. Minor but error-prone.

### Stage B — integration eval  (`scripts/eval_integration.py`, 74 ln)
- 200k-cell subsample, kNN same-neighbor fraction for `dataset_slug` (batch) and `organoid_type`
  (biology) vs random baseline, plus a 2-panel UMAP → `data/scvi_umap_eval_v5.png`.
- **This is a crude integration metric** (see §7).

### Stage C — Braun reference model  (`scripts/braun_label_transfer.py`, 282 ln)
- Stage 1 of this script trains scVI→scANVI on Braun: shared Ensembl genes (34,746) → 2000 HVG
  (`flavor='seurat'`) **+ a forced `RARE_PANEL` of ~6 lineage markers** (→ 2006 genes total),
  `batch_key=donor_id`, `label_key=CellClass`. Saves `data/braun_scanvi_full/` (the **healthy**,
  reusable reference: 88% self-accuracy on Braun).
- Stage 2 (query surgery + `model.predict()`) in this script is **DEPRECATED** — its scANVI-head
  prediction collapses on the organoid query (see §6.4). Superseded by Stage D.
- `RARE_PANEL` and HVG-forcing were added mid-session (by the user) to address rare-type
  detection; it does **not** fix the head collapse, but it usefully guarantees the rare markers
  are in the panel for Stage F gating.

### Stage D — production label transfer  (`scripts/braun_transfer_finalize.py`, 114 ln) ✅
- Reuses `braun_scanvi_full`; reindexes Braun to the model's stored 2006 genes (read from
  `model.pt`); streams the full 4M-cell query (memory-safe, see §6.3); scArches surgery (15 ep);
  then **kNN-on-the-joint-latent** for both `CellClass` and `Region` (+ confidence).
- Output: `data/braun_transfer_full_knn.h5ad` (4,079,890 × 30 latent + `CellClass_pred/conf`,
  `Region_pred/conf`, and all 54 atlas obs columns).
- **This is the canonical transfer output.** `data/braun_transfer_full.h5ad` is the broken
  scANVI-head run and should be deleted (see §7).

### Stage E — stratified batch mixing  (`stratified_mixing.py`, `stratified_mixing_region.py`)
- Joins `CellClass_pred` (+`Region_pred`) onto the scVI latent; within each cell type (and
  type×region) builds a kNN among only those cells and measures same-`dataset_slug` fraction vs
  the within-group baseline → `data/stratified_mixing{,_region}.tsv`.

### Stage F — benchmark Q1/Q2/Q3
- `benchmark_q1_coverage.py` + `benchmark_q1_sizecontrol.py`: per-deposit cell-type/region
  presence; multi vs single; equal-detection-power robustness check.
- `benchmark_q2_correspondence.py`: pseudobulk (mean raw counts→CP10K→log1p, 2006 genes)
  organoid-class × Braun-class correlation → `data/q2_correspondence_corr.tsv`.
- `benchmark_q2_markergate.py`: marker-score gating of microglia/endothelium/oligo, transfer
  miss-rate, correspondence of gated cells to Braun.
- `benchmark_q3_gaps.py`: reference-coverage (fraction of each Braun cell's kNN that are organoid)
  by CellClass/Region/Age, multi vs single → `data/q3_coverage.tsv`.
- Diagnostics (kept for provenance): `diag_scanvi_collapse.py`, `braun_query_fix_test.py`.

---

## 4. Outputs catalog

| File | What | Trust |
|---|---|---|
| `data/scvi_model_v5_full/`, `data/scvi_latent_v5_full.h5ad` | scVI atlas embedding | ✅ |
| `data/braun_scanvi_full/` | Braun reference scANVI (2006 genes) | ✅ reusable |
| `data/braun_transfer_full_knn.h5ad` | **canonical** organoid labels (CellClass/Region + conf + latent) | ✅ |
| `data/braun_transfer_full.h5ad` | scANVI-head run, 100% one class | ❌ delete |
| `data/braun_transfer_pilot.h5ad` | 200k pilot | sanity only |
| `data/stratified_mixing{,_region}.tsv` | within-type batch mixing | ✅ (crude metric) |
| `data/q2_correspondence_corr.tsv` | organoid×Braun pseudobulk corr | ⚠ partial circularity |
| `data/q3_coverage.tsv` | per-Braun-cell organoid coverage | ⚠ read relative only |
| `data/scvi_umap_eval_v5.png` | integration UMAP | qualitative |

---

## 5. Results summary (for context)

- **Integration:** global same-dataset kNN ratio 17.5×; stratified by CellClass ~13–16×;
  by CellClass×Region ~8–13× for major neural groups. → real residual batch effect (~10×),
  partly tamable, after controlling for type+region.
- **Q1:** multi-lineage doesn't broaden repertoire (NS) but shifts composition to support
  lineages (microglia +19pp deposit-presence, glia +21, oligo +10, endothelium +5/multi-only);
  single-lineage enriched for neurogenic trajectory (Neuroblast −27). Size-controlled.
- **Q2:** correspondence to primary graded — core neural 0.86–0.93; microglia 0.82 (clean);
  endothelium/oligo ~0.5 (immature).
- **Q3:** multi-lineage closes microglia (+0.09) & oligo (+0.19) coverage; vascular (0.011) and
  posterior CNS (cerebellum/hindbrain) gaps persist; no maturation-age gradient.

---

## 6. Challenges encountered (war stories — useful for reviewers)

**6.1 "Training too fast" → batch-size/step-count tradeoff.** Initial run used `batch_size=2048`
(GPU-efficient) but that's ~16× fewer gradient steps than scVI's default-128 over the same epochs.
Retraining at 512/100ep (~784k steps, 2.2h) over 2048/15ep (~29k steps, 7min) improved batch
mixing only ~10% relative and *dropped* biology ~7%; biology/batch ratio barely moved (1.82→1.91).
**Conclusion: residual batch structure is structural (505 batches confounded with biology), not
underfitting.** Don't chase epochs.

**6.2 Gene-namespace mismatch.** Atlas var_names are HGNC symbols; Braun is Ensembl. Naive overlap
was 57/3000. Bridged via `hnoca_var_canonical.tsv` (hgnc↔ensembl); 34,746-gene intersection.

**6.3 Memory.** Full query matrix is ~92 GB in RAM (4M × 36,842 @ 7.7% density); the first
finalize draft did `query.layers['counts']=query.X.copy()` → ~175 GB → OOM risk. Fixed by
**streaming `chunked_X` and keeping only the ~2006 model genes** before materializing. Braun loads
fine (~7 GB @ 0.9% density).

**6.4 scANVI classifier-head collapse (the big one).** The scANVI `model.predict()` on the
organoid query collapsed to **100% of a single class** (Radial glia in one run, Neuron in another —
the class is *unstable*). Diagnosis (`diag_scanvi_collapse.py`): the classifier predicts the **Braun
reference at 88% accuracy** with correct proportions — so it is **not** a training/imbalance
collapse; it's an **out-of-distribution-query** failure. `braun_query_fix_test.py` confirmed
reducing surgery epochs (40→10→0) does not fix it. **Fix: bypass the head, transfer labels by kNN
on the joint latent** (the latent generalizes; the classifier head does not). This is why Stage D
uses kNN, not `predict()`.

**6.5 Run-provenance confusion.** Two full runs existed concurrently (a 2000-gene run I launched
and a 2006-gene `RARE_PANEL` run started mid-session). Timestamps were inconsistent
(`braun_transfer_full.h5ad` predates its own model dir), proving they came from different runs.
Lesson: **out-tag/version every artifact and stamp the git SHA + params into the .h5ad `.uns`**
(not done — review item).

**6.6 Per-cell marker gating is ambient-contaminated.** `score_genes` thresholds caught 7–12% of
cells as "oligo" (MBP/PLP1 are high-expressors that leak into ambient RNA) and endothelial gating
was polluted by `A2M`/`FLT1`/`KDR`. Counts are unreliable; only the **correspondence asymmetry**
(microglia 0.82 ≫ endo/oligo ~0.5) is robust across marker sets.

---

## 7. Methodological limitations & gaps (read critically)

1. **Integration metric is crude.** Everything uses a kNN same-neighbor fraction. We never ran
   established metrics (**iLISI / kBET / scIB-metrics / kNN-graph connectivity**). The same-fraction
   conflates batch and biology and has an awkward baseline. *Replace with scib-metrics.*
2. **Q2 Part 1 circularity.** Organoid cells were labeled by latent proximity to Braun, then we
   correlate organoid-class vs Braun-class *expression*. Expression and latent aren't independent,
   so diagonal dominance is partly guaranteed. The marker-gated Part 2 is the de-circularized check,
   but it has its own (gating) problems. *A cleaner design: held-out marker-defined ground truth, or
   correlation on genes explicitly excluded from the scVI/scANVI feature set.*
3. **Q3 absolute coverage is uninterpretable.** Coverage is ~0.05 everywhere (vs 0.5 neutral)
   because organoid and primary occupy offset latent regions — part biology, part residual scArches
   query-batch separation. We only read *relative* differences. *Needs a calibrated null (e.g.
   reference-vs-reference coverage) to separate "gap" from "global domain shift".*
4. **No replicates / CIs.** Single random seed for most subsamples; no bootstrap CIs on any reported
   number (Q1 presence rates, correspondences, coverage). Small-n populations (Neural crest n=47,
   Vascular n=596 in Braun samples) are noisy and reported without error bars.
5. **kNN label transfer ignores class imbalance.** Rare Braun classes (Immune/Vascular 0.5–0.7%)
   are systematically under-called → the "~0% Immune" global result. No distance weighting,
   class-balancing, or abstention/OOD threshold (cells with no good match still get a label).
6. **scANVI is trained then thrown away.** We spend the compute to train the scANVI head, then
   bypass it. Either fix the head for OOD (class-balanced sampling, calibration, or scArches with
   `unfrozen`/different freeze flags) or drop scANVI and train plain scVI for the reference.
7. **Massive code duplication.** The gene bridge, `same_frac`, `pseudobulk`, chunked query reader,
   and Braun-reindex logic are re-implemented across ~6 scripts. *Extract a shared `organoid_atlas/`
   module* (io, gene-mapping, metrics). High-value refactor.
8. **No config / provenance.** Hardcoded paths, magic numbers (N_SUB, CAP, thresholds), no git-SHA
   or param stamping into outputs, no manifest. Reproducibility relies on reading each script.
9. **`multi_lineage` is a mixed-dtype mess** (`'0'/'1'/'False'/'True'/'No'` strings) hand-mapped in
   every script. Should be normalized once upstream in the atlas build.
10. **Alignment by positional/`reindex` assumptions.** Several joins assume `obs_names` align or use
    `.reindex(...).values`; no assertions guard against silent misalignment. *Add explicit
    index-equality asserts.*
11. **Pilot vs full divergence.** The pilot used a different (100k, 2000-gene) reference than the
    full run; the pilot "validated" a pipeline that then behaved differently at scale (the head
    collapse only showed at full scale). Pilots should mirror the full config more closely.
12. **No held-out evaluation of the transfer.** We never split Braun into train/test to quantify
    transfer accuracy independent of the self-prediction number.

---

## 8. Concrete asks for the reviewer (Codex)

Highest-value first:

1. **Replace the integration metric** in `eval_integration.py` / `stratified_mixing*.py` with
   `scib-metrics` (iLISI, kBET, cLISI, kNN-graph connectivity); keep the stratified variant.
2. **De-duplicate** into a shared module (gene bridge, chunked reader, `same_frac`, pseudobulk,
   Braun-reindex). Flag any subtle differences between the current copies (e.g. normalization).
3. **Harden the kNN label transfer** (`braun_transfer_finalize.py`): distance-weighted kNN,
   class-balanced reference sampling for rare types, and an **abstention threshold** so off-manifold
   organoid cells aren't force-labeled. Assess whether this recovers Immune/Vascular.
4. **Investigate / fix or remove the scANVI head** (`braun_label_transfer.py` Stage 2): is the OOD
   collapse fixable via class-balanced `n_samples_per_label`, calibration, or different scArches
   freeze flags? If not, drop scANVI for plain scVI + kNN and save the compute.
5. **Add provenance**: stamp git SHA + params + input file hashes into `.uns` of every output; add a
   `--out-tag`-versioned manifest. Delete `data/braun_transfer_full.h5ad` (broken).
6. **Add CIs / seeds**: bootstrap CIs on Q1 presence rates, Q2 correspondences, Q3 coverage; expose a
   `--seed` everywhere; flag small-n populations.
7. **Q3 null calibration**: implement reference-vs-reference (and organoid-vs-organoid) coverage as
   a baseline so "gap" is separated from the global organoid↔primary domain offset.
8. **Q2 circularity**: re-run correspondence on a gene set held out of the scVI/scANVI features.
9. **Early-stopping fix** in `train_scvi.py`: `min_delta` or monitor a metric that actually plateaus.
10. **Replace per-cell marker gating** (`benchmark_q2_markergate.py`) with **cluster-level**
    annotation (Leiden on the latent, score clusters) to escape ambient-RNA false positives.

---

## 9. Reproduction order

```bash
P=/opt/homebrew/Caskroom/miniforge/base/bin/python3.13
# A. scVI embedding (~50 min, MPS)
$P scripts/train_scvi.py --max-epochs 100 --batch-size 512 --early-stopping --out-tag v5_full --accelerator mps
# B. integration eval
$P scripts/eval_integration.py
# C. Braun reference scANVI (~part of the full transfer script; reused thereafter)
#    (braun_label_transfer.py Stage 1 builds data/braun_scanvi_full)
# D. canonical label transfer (~50 min): surgery + kNN-on-latent
$P scripts/braun_transfer_finalize.py
# E. stratified mixing
$P scripts/stratified_mixing.py ; $P scripts/stratified_mixing_region.py
# F. benchmarks
$P scripts/benchmark_q1_coverage.py ; $P scripts/benchmark_q1_sizecontrol.py
$P scripts/benchmark_q2_correspondence.py ; $P scripts/benchmark_q2_markergate.py
$P scripts/benchmark_q3_gaps.py
```

Diagnostics (optional, explain the design choices): `diag_scanvi_collapse.py`,
`braun_query_fix_test.py`.

---

## 10. What to trust vs distrust (one-paragraph honest summary)

**Trust:** the scVI embedding exists and is reproducible; the Braun reference scANVI is healthy
(88% self-accuracy); the kNN-on-latent transfer gives biologically coherent labels (validated
against `cell_type_original`); the Q1 composition shift and Q2 correspondence *asymmetry* are robust
to reasonable perturbations. **Distrust / caveat:** absolute integration quality (residual ~10×
within-type batch), any absolute coverage number in Q3 (domain offset), rare-type *counts* from
per-cell gating (ambient), and any single-run number without a CI. The science story
(*multi-lineage adds primary-faithful microglia, partially oligo, but not vasculature or posterior
identity*) rests on relative/repertoire comparisons that survive these caveats; the methods scaffold
underneath needs the hardening listed in §8.

---

## 11. Addendum — changes since commit `aac5cb4` (post first review)

This section responds to the first code review. **Review item #1 (calibrated kNN
transfer) is implemented; items #2–#5 are still open.**

### #1 Calibrated kNN transfer — DONE  (`scripts/braun_transfer_calibrated.py`)
Replaces the argmax / unweighted / unbalanced kNN in `braun_transfer_finalize.py`.
Reuses the **existing** joint latents (query `X` from `braun_transfer_full_knn.h5ad`
+ Braun reference latent from `braun_scanvi_full`) — **no surgery re-run** (~3 min).
Adds:
- **distance-weighted** votes (`w = 1/(dist+eps)`),
- **class-balanced reference sampling** (cap 6000/class → 61,744 ref cells) as prior
  correction for rare classes,
- **abstention**: max posterior < `TAU=0.4` → `CellClass_cal='Unknown'`,
- **OOD diagnostic**: query mean-kNN-distance vs the reference's own 95th-percentile
  in-distribution distance.

Output: `data/braun_transfer_full_calibrated.h5ad` (adds `CellClass_cal`,
`CellClass_cal_conf`, `abstain`, `ood`).

**Results:**
- Held-out Braun accuracy 0.965 (weighted) ≈ 0.964 (unweighted); rare-class recall
  (Immune/Vascular/Oligo) = **1.0 both** → the classifier separates rare types fine
  *when the reference is balanced*; the original under-calling was **reference
  imbalance**, not a classifier defect.
- Query distribution rebalanced: Radial glia 38→28%, Neuron 31→25% (were over-called);
  **Neuronal IPC 3.5→7.3%, Oligo 0.6→1.2% (both ~2×)**, Glioblast/Neuroblast/Neural
  crest up; abstain 3.6%.
- **Truly-rare types stay rare** even with balanced ref + perfect Braun recall:
  Immune 0.24→0.35%, Vascular 0.03→0.05%. ⇒ microglia/endothelium are **biologically
  scarce in organoids** (~0.2–0.4%, consistent with marker gating), not a transfer
  artifact. The Q1–Q3 rare-type gap is real.
- **OOD = 78.3%** — 78% of organoid cells are farther from Braun than 95% of Braun
  cells are from each other. This **quantifies the organoid↔primary domain offset**
  (= Q3's ~0.05 coverage, §5/§7.3) and is the strongest evidence yet that transferred
  labels are a **"nearest fetal correlate," not ground truth**. Recommend foregrounding
  this in any writeup.

**Caveat on the held-out check:** train/test were both drawn from the *balanced* subset,
so it demonstrates classifier capability, not the magnitude of the imbalance fix; the
real evidence of the fix is the query-distribution shift above.

### #2 scIB metrics — DONE  (`scripts/scib_metrics_eval.py`)
Replaces the crude same-neighbor ratio with the `scib-metrics` Benchmarker panel,
computed on the scVI latent (300k subsample), batch=`dataset_slug`,
label=`CellClass_cal` (abstained dropped). **Installed in an isolated venv
(`/Users/eg/.venvs/scib`, scib-metrics 0.5.9 + jax) to keep jax/numpy churn OUT of the
scvi base env** (the base env has known numpy-ABI fragility). Output `data/scib_metrics.tsv`.

Results (raw, single embedding → per-metric values interpretable; aggregates are
summaries): **Batch** — iLISI **0.015**, kBET 0.185, PCR 0.004, graph connectivity 0.683;
**Bio** — cLISI **0.962**, silhouette label 0.502, isolated labels 0.629, NMI/ARI 0.293/0.154.
→ confirms **poor batch integration + good biology preservation**; the earlier same-neighbor
ratios (~10–17×) were directionally correct. This is now the integration yardstick (the
same-neighbor scripts are demoted to supplementary).

### Single calibrated-label rerun — DONE  (`scripts/benchmark_rerun_calibrated.py`)
Recomputes Q1 + Q2 on `CellClass_cal` vs old `CellClass_pred` in one pass.
- **Q1 robust**: multi-lineage still adds microglia (+12pp), oligo (+16), neural crest (+22),
  endothelium (+5). **Correction:** Glioblast's +23 multi-enrichment was an argmax artifact
  → +4 under calibration. Neuroblast single-enrichment robust (−19).
- **Q2 robust**: graded correspondence preserved (core neural 0.88–0.93 > microglia 0.76 /
  endothelium 0.76); **diagonal dominance 12/12 both old and new**. Minor downshifts only.
- **Q3 unchanged**: aggregates by Braun labels on the organoid side → organoid-label
  calibration does not affect it.

**Net:** the benchmark conclusions survive both the calibrated labels and the proper metrics;
the only material change is correcting the overstated Glioblast enrichment.

### #3 De-circularize Q2 — DONE  (`scripts/benchmark_q2_heldout.py`)
Recomputes organoid↔Braun pseudobulk correspondence on **genes held out of the 2006
transfer feature set** (32,740 shared genes available; 4,000 sampled), with bootstrap
CIs. Result: **diagonal dominance 11/12 on independent genes**, self-corr 0.84–0.96,
tight CIs; `circular_inflation` is *negative* everywhere (held-out ≥ on-feature) → the
correspondence is **not** a circular artifact. Caveat: random held-out genes are
housekeeping-dominated (uniformly high corr) so the readout is diagonal dominance, not
magnitude; a held-out-*variable*-gene variant would be sharper (follow-up). Only failure:
Vascular (n=107) → Neural crest by 0.018 (the already-weak rarest class). Output
`data/q2_heldout_correspondence.tsv` (+ `.provenance.json`).

### #5-lite Provenance stamping — DONE  (`scripts/_provenance.py`)
Shared helper: `stamp(adata, script, params)` → `.uns['provenance']` (git SHA, script,
params, timestamp, python + lib versions); `write_sidecar(tsv, ...)` → `<tsv>.provenance.json`.
Wired into `braun_transfer_calibrated.py` (.uns) and `benchmark_q2_heldout.py` (sidecar).
Full #5 (central config + de-dup of gene-bridge/chunked-reader/metrics into a module) still open.

### #4 Cluster-level marker gating — DONE  (`scripts/benchmark_q2_clustergate.py`)
Replaces fragile per-cell `score_genes` thresholds: Leiden-cluster the query on the scVI
latent (400k subsample, res=2), score marker panels PER CLUSTER (ambient-robust), gate
whole clusters, report prevalence + transferred-label agreement (CellClass_cal = Braun
cross-ref) + confidence + OOD fraction. Output `data/q2_clustergate.tsv` (+ provenance).
Results — far cleaner than per-cell, and **corrected a misattribution**:
- **Microglia**: 1 coherent cluster, 0.18%, marker-score 1.14, **98% agree with transferred
  Immune**, conf 0.98, **90% multi-lineage** → robustly real.
- **Endothelium**: **no cluster above gate** → not a resolvable population (scattered/ambient/
  immature), consistent with Q1–Q3.
- **"Oligo"-marker cluster**: 1.71%, but **98% transferred as Neural crest, not Oligo**
  (PLP1/MBP also mark neural-crest/peripheral glia) → per-cell oligo gating was a
  misattribution; true mature CNS oligodendrocytes don't form a distinct cluster either.
- OOD-stratified: these rare support clusters are *less* OOD (microglia 28%, neural-crest 21%)
  than the global 78% — support lineages sit closer to the Braun manifold than the immature
  neural mass.

### τ-sweep (review #1 refinement) — DONE  (in `braun_transfer_calibrated.py`)
abstain% by τ: 0.2→0.0, 0.3→0.5, **0.4→3.6**, 0.5→9.3, 0.6→21.9. **τ=0.4 is the knee**
(steep rise above it); rare classes retained at all τ. Also now saves `CellClass_cal_argmax`
(pre-abstention) so τ can be re-analyzed without recompute. Provenance bug fixed
(`str()`-coerce lib versions; torch's version object broke `.uns` serialization).

### #4 robustness sweep — DONE  (`scripts/benchmark_q2_clustergate_sweep.py`)
Resolutions {1.0,1.5,2.0} × gate {0.25,0.30,0.35} × 2 seeds, + per-cluster **Braun
centroid-correlation** (top class + margin, on the 2006 genes — independent of the
kNN-transfer label). Output `data/q2_clustergate_sweep.tsv`. Conclusions are stable:
- **Microglia**: same ~0.18% single cluster every run, centroid → **Immune, margin 0.41–0.46**
  (decisive). Transfer label and expression-centroid agree.
- **Endothelium**: at/below the detection limit — seed0 finds **no cluster** at any res;
  seed1 finds a ~70-cell (0.02%) cluster that centroids → Vascular (margin ~0.17). Not a
  reliably resolvable population.
- **"Oligo"-marker cluster**: centroid → **Neural crest** in all 6 runs (small margin 0.03–0.13,
  since neural-crest/oligo share PLP1/MBP programs) → confirms the misattribution; no clean
  mature-CNS-oligo population.

### Baseline-PCA scIB — DONE  (`scripts/scib_metrics_baseline.py`)
Added unintegrated lognorm-PCA(50) as a comparator embedding in the Benchmarker (vs scVI
latent), pinned seed. Output `data/scib_metrics_baseline.tsv`. **Sobering result: scVI ≈
unintegrated PCA on every axis** — iLISI 0.0157 (PCA) vs 0.0152 (scVI), kBET 0.190 vs 0.189,
cLISI 0.961 vs 0.962, Total 0.444 vs 0.446. Conclusion (narrowed wording):
**no measurable integration gain vs PCA under current scIB panel**. Two readings (both consistent with prior findings):
(a) the 505-batch structure is biological/irremovable, so neither method mixes it — and neither
should; (b) if better cross-dataset mixing is genuinely needed, the lever is label-aware
integration (scANVI/scPoli), not this scVI run. Either way, the pipeline should **not claim
strong integration credit** from scVI on these metrics.

### Null-calibrated OOD — DONE  (`scripts/benchmark_ood_nullcalibrated.py`)
Built an in-distribution null (Braun test→train kNN-dist; p95=1.300, p99=1.581) instead of the
ref-self threshold. **Global query OOD: 75.0% (>p95), 42.1% (>p99)** — confirms the prior
self-method 78.3% (robust to calibration). **Per-class OOD** (each class vs its own Braun null):
every class majority-OOD (68–90%), and the **rare support lineages are *most* offset** —
Vascular 90%, Oligo 87%, Immune 79% vs dominant neural ~68–70%. I.e. even where organoids make
the right cell type, it is transcriptomically distinct from the primary counterpart (worst for
support lineages). Output `data/ood_nullcalibrated.tsv`. (NB: differs from clustergate's global
28%/21% OOD because that asks "near *any* primary cell" while this asks "near *same-class*
primary cells".)

### Still open
- **#5** full provenance/config + shared-module refactor — partial (stamping done; central
  config + de-dup of gene-bridge/chunked-reader/metrics still pending).
- Queued: held-out-*variable*-gene Q2; calibration reliability curve/ECE.

### Related reference
- `docs/annotation_schema.md` — schema of the (gitignored) annotation workbooks
  (`brain_organoid(_GSMannotations).xlsx`): sheets, columns/vocabularies, xlsx→obs field map,
  `annotation_level` (gsm vs deposit) semantics, and known data-quality quirks.

---

## 12. Methods-final robustness tables (collated from existing outputs)

### A. Cluster-gating robustness — microglia (the one robustly recovered rare lineage)
Source: `data/q2_clustergate_sweep.tsv` (recovery/n/centroid across the grid) +
`data/q2_clustergate.tsv` (transfer agreement / OOD / conf, canonical run seed0·res2).
Gate thresholds {0.25, 0.30, 0.35} give **identical** results (microglia cluster-mean
score ≈1.14 ≫ all gates), so threshold is collapsed below.

| seed | res | recovered | n cells (%) | centroid→Braun (margin) |
|---|---|---|---|---|
| 0 | 1.0 | Y | 735 (0.184%) | Immune (0.458) |
| 0 | 1.5 | Y | 735 (0.184%) | Immune (0.454) |
| 0 | 2.0 | Y | 736 (0.184%) | Immune (0.451) |
| 1 | 1.0 | Y | 749 (0.187%) | Immune (0.410) |
| 1 | 1.5 | Y | 755 (0.189%) | Immune (0.409) |
| 1 | 2.0 | Y | 774 (0.194%) | Immune (0.409) |

Canonical run (seed0·res2): **transfer agreement 98% Immune**, mean conf **0.984**,
**OOD 28%** (global ref-self threshold), **90% multi-lineage**.
Other lineages across the grid: **endothelium** recovered only at seed1 (~69–73 cells,
0.02%, →Vascular) and **not at seed0** → not robust; **"oligo"-marker** cluster recovered
in all runs but centroid → **Neural crest** (margin 0.03–0.13), not Oligo.

### B. OOD threshold sensitivity (null-calibrated)
Source: `data/ood_nullcalibrated.tsv` + provenance sidecar. Null = Braun test→train
mean kNN-distance; global thresholds p95=1.300, p99=1.581.

**Global query OOD:** **p95 = 75.0%**, **p99 = 42.1%** (prior ref-self-p95 method: 78.3%).

**Per-class OOD (vs own-class null):**

| CellClass | n_query | p95 thr | OOD% @p95 | OOD% @p99 |
|---|---|---|---|---|
| Vascular | 238 | 1.243 | 89.5 | *(pending)* |
| Oligo | 6,143 | 1.423 | 86.6 | *(pending)* |
| Immune | 1,691 | 1.346 | 79.2 | *(pending)* |
| Fibroblast | 23,165 | 1.427 | 77.0 | *(pending)* |
| Neuroblast | 60,604 | 1.322 | 75.3 | *(pending)* |
| Neuronal IPC | 36,202 | 1.391 | 73.0 | *(pending)* |
| Radial glia | 141,432 | 1.337 | 70.1 | *(pending)* |
| Glioblast | 66,638 | 1.552 | 69.0 | *(pending)* |
| Neuron | 126,206 | 1.280 | 68.2 | *(pending)* |
| Erythrocyte | 17 | 0.710 | 100.0 | *(n too small)* |

Rare support lineages (Vascular/Oligo/Immune) are the **most** per-class OOD — even where
organoids make the right cell type, it's transcriptomically distinct from primary.
*Per-class p99 column is now code-ready in `benchmark_ood_nullcalibrated.py` and will populate
on the next OOD run (not recomputed here to avoid the ~85-min job).*

---

## 13. Decision request — integration claim & comparator (GO / NO-GO)

**Finding:** under the current scIB panel there is **no measurable integration gain vs PCA**
(iLISI 0.0157 vs 0.0152, kBET 0.190 vs 0.189, Total 0.444 vs 0.446). Q1–Q3 are therefore kept
on repertoire/correspondence footing, which does not depend on tight cross-dataset mixing.

**Decision needed (reviewer):**
- **Option (a) — Freeze.** Keep scVI as latent/denoiser, **drop integration-gain claims**, ship
  on the repertoire/correspondence + OOD framing. No further compute.
- **Option (b) — Scoped comparator.** Run **one** label-aware comparator embedding
  (scANVI or scPoli) into the same Benchmarker, judged against predefined acceptance thresholds:

  | Criterion | GO threshold |
  |---|---|
  | ΔiLISI vs scVI | ≥ **+0.05** |
  | ΔkBET vs scVI | ≥ **+0.05** |
  | cLISI drop (biology) | ≤ **0.02** (no collapse) |
  | Global OOD @p95 | **≤ 75%** (not worse) |
  | Mean transfer confidence | not worse than current |

  **GO** only if iLISI **and** kBET both clear +0.05 **and** all guardrails hold; otherwise
  **NO-GO → revert to (a)**.

Thresholds are a proposal — adjust as needed. Remaining open items unchanged: full #5 (config +
shared-module de-dup); queued held-out-*variable*-gene Q2 and calibration reliability/ECE.

### Outcome — comparator RAN (Codex GO), result NO-GO → freeze scVI
Ran scANVI (label-aware, semi-supervised by `CellClass_cal`, 20 ep, from `scvi_model_v5_full`)
→ `scripts/train_scanvi_comparator.py`, evaluated vs PCA + scVI in one Benchmarker
→ `scripts/scib_metrics_comparator.py`, `data/scib_metrics_comparator.tsv`.

| Metric | X_pca | X_scvi | X_scanvi |
|---|---|---|---|
| iLISI (batch) | 0.0157 | 0.0152 | 0.0171 |
| kBET | 0.190 | 0.189 | **0.254** |
| Graph connectivity | 0.707 | 0.717 | 0.783 |
| cLISI (bio) | 0.961 | 0.962 | 0.990 |
| NMI / ARI | 0.33/0.19 | 0.33/0.19 | 0.45/0.27 |
| Total | 0.444 | 0.446 | 0.484 |

vs scVI: **ΔiLISI = +0.0019 (need ≥+0.05) ❌**, ΔkBET = +0.0651 ✅, cLISI drop = −0.028 ✅.
**EMBEDDING-LEVEL VERDICT: NO-GO → freeze scVI.** (OOD/confidence guardrails not run — embedding
gate already failed; no re-transfer.)

Interpretation: scANVI improved kBET / graph-connectivity / bio-conservation but **did not move
the batch-mixing bar (iLISI flat)** — reinforcing that the 505-batch structure is
biological/irremovable even for label-aware integration. The bio-conservation gains
(cLISI/NMI/ARI) are **partly circular** (scANVI is trained on the same `CellClass_cal` labels
those metrics score against), so they are not independent evidence. **Decision: freeze scVI as
latent/denoiser; do not claim integration gain vs PCA under the current scIB panel.**

**DECISION LOCKED (reviewer-confirmed, 2026-05-25).** Keep NO-GO; the pre-registered
**conjunctive** gate (ΔiLISI **AND** ΔkBET) was **not** relaxed to OR post hoc (doing so would
weaken credibility). **iLISI is the primary batch-mixing metric in this decision gate; the kBET
improvement is acknowledged but insufficient without iLISI movement.** scVI is retained as the
latent/denoiser backbone; all results are framed on **repertoire / correspondence + OOD-aware
transfer**, not integration-gain claims. No further integration experiments planned.
Next (compute idle until specified): full **#5** (config + shared-module de-dup),
**held-out-*variable*-gene Q2**, **calibration reliability/ECE**.
