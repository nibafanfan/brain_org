# Method review — atlas rebuild + integration pipeline

Date: 2026-05-25
Scope: review of current v5 plan/training choices as implemented in `scripts/preprocess_atlas.py`, `scripts/train_scvi.py`, and pooled-control handling template `scripts/fix_gse297594_control.py`.

## 1) HVG strategy (`cell_ranger`, n=3000, `batch_key='bio_sample'`)

### What is strong already
- Using `cell_ranger` flavor is pragmatic and robust at very large scale, and avoids the known local `seurat_v3` dependency/ABI issues noted in the docs.
- Using a biological replicate key for HVG (`bio_sample`) while reserving technical key for integration is conceptually reasonable: it favors genes recurring across biological units rather than just sequencing pools.
- Zero-count-cell removal before HVG is correct and should remain mandatory.

### Risks / concerns
- `bio_sample` cardinality is very high at this scale (~hundreds of samples) and heterogeneous in size; with `MIN_BATCH=50` you reduce instability, but you also remove small biological units entirely from feature selection.
- Dropping micro-batches before HVG can bias selected genes toward large protocols/datasets and away from rare but valid programs.
- `cell_ranger` HVG on log-normalized data can overweight high-mean housekeeping-like variance patterns in ultra-large mixtures. Not fatal, but worth sensitivity checks.
- Fixed 3000 HVGs may be slightly tight for broad cross-protocol neurodevelopment mixtures (neurons + glia + vascular + immune); 4k–6k often stabilizes cross-lineage resolution.

### Recommended adjustments
1. Keep current pipeline as baseline, but run a **feature-sensitivity panel**:
   - 3k (`cell_ranger`, current)
   - 5k (`cell_ranger`)
   - optionally Pearson-residual HVG (if environment can support it stably)
2. Replace hard removal of micro-batches with one of:
   - compute HVG on all cells but cap per-batch contribution via per-batch downsampling, or
   - retain micro-batches and aggregate HVG by prevalence rank (how many batches nominate each gene).
3. Report HVG provenance metrics each run: per-lineage marker retention, per-dataset contribution, and overlap/Jaccard with previous run.

## 2) scVI setup (latent dim, epochs, likelihood, batch-key split)

### Current choices
- `gene_likelihood='nb'`
- `n_latent=30`
- old run: batch_size=2048, ~15 epochs (~29k steps)
- retrain: batch_size=512, 100 epochs (~784k steps)
- HVG keyed on `bio_sample`, scVI batch keyed on `tech_sample`

### Assessment
- **Batch-key split is defensible**: biological key for HVG, technical key for generative correction is a standard and often beneficial separation.
- `n_latent=30` is a sensible default for multi-million-cell integration.
- `NB` is reasonable if QC and chemistry heterogeneity are already moderated.

### Risks / concerns
- Undertraining was a plausible hypothesis initially, but a 100-epoch retrain showed only modest gains in mixing and a similar contraction of biology, so epoch count is not the dominant bottleneck for this dataset.
- No explicit train/val diagnostics are written out for model selection; decisions may be based on runtime rather than convergence.
- `NB` vs `ZINB`: usually NB is fine for UMI data, but this should be justified empirically at least once (reconstruction/transfer metrics).
- Single latent dimensionality may miss nuanced regional lineages; 30 may be fine, but 20/50 should be checked.


### New empirical update (post-review retrain)
- Head-to-head comparison (old 2048/15ep vs new 512/100ep) shows ~10–12% lower same-neighbor batch metrics (`tech_sample`, `bio_sample`, `dataset_slug`) and ~7% lower `organoid_type` neighborhood purity.
- The biology-to-batch ratio improved only modestly (1.82 -> 1.91, ~5% relative), indicating the latent mostly contracted rather than fundamentally restructured.
- Practical conclusion: additional epochs alone are a diminishing-returns lever; the 100-epoch run is a reasonable converged model to keep, but further improvement likely requires metric and model/key changes rather than longer training.

### Recommended adjustments
1. Keep the converged 100-epoch model as the production default for downstream transfer unless/until a different objective function is adopted.
2. Add a **small hyperparameter sweep** on a stratified 200k–500k subset:
   - latent dims: 20, 30, 50
   - likelihood: NB (baseline), optional ZINB check
   - epochs: with early stopping enabled
3. Promote best config to full run and use convergence criteria rather than fixed epoch budget.
4. Save training diagnostics artifacts (per-epoch ELBO/loss, runtime, batch mixing metrics, label transfer metrics) in versioned files.
5. Consider covariates if available (`assay_sc`, chemistry, donor/individual) to reduce over-correction or residual technical structure.
6. If stronger mixing is still required, prioritize structural levers over more epochs: coarser batch definitions, semi-supervised scANVI with existing annotations, and scIB-style metrics (e.g., iLISI/kBET within cell types) that reduce biology-confounded interpretation.

## 3) Pooled/cell-hashed deposits (`fix_gse297594_control.py` pattern)

### What is strong
- Correctly recognizes that raw pooled MTX cannot infer genotype/control in hashed designs.
- Uses Seurat `meta.data` (`condition`, singlet filters) as source of truth for control extraction.
- Re-runs standard QC + canonical projection after extraction, preserving atlas consistency.

### Risks / concerns
- Reliance on processed Seurat metadata introduces a trust boundary: upstream demux/QC decisions may vary across studies.
- Current singlet logic allows all cells when `HTO_classification.global` is absent; this can admit ambiguous droplets in some datasets.
- Potential barcode reorder/mismatch risks are mitigated by reindexing, but should still be explicitly audited in output logs.
- Template is not yet generalized in main dispatcher; operational risk of accidental fallback to pooled-MTX loader remains.

### Recommended adjustments
1. Implement formal loader type (`seurat_hto_control`) in `rebuild_atlas.py` dispatch and require explicit config entry for flagged deposits.
2. Add hard validation checks per deposit:
   - all retained cells have allowed control labels,
   - singlet/demux status present or explicitly waived with reason,
   - cell count delta vs raw object recorded.
3. Emit a per-deposit contamination audit report (WT/MUT counts before/after, excluded classes, unknown labels).
4. Prefer storing both raw demux fields and normalized atlas fields in obs for traceability.

## 4) Deposit-level vs per-cell annotation fallback (`annotation_level`)

### Assessment
- Per-cell annotations should be primary whenever available; deposit-level fallback is useful but can blur heterogeneity.
- For mixed-type deposits (e.g., CO + TA in one slug), deposit-level labels are statistically unsafe for benchmarking and transfer evaluation.

### Risks / implications
- Label transfer benchmarking may appear better/worse depending on how many cells inherit coarse deposit-level labels.
- Cross-protocol comparisons can be confounded if one cohort is mostly per-cell and another mostly deposit-level fallback.

### Recommended adjustments
1. Treat `annotation_level` as a first-class covariate in every benchmark table/plot.
2. Run all transfer metrics in three strata:
   - per-cell annotated only,
   - deposit-level fallback only,
   - combined.
3. Add minimum confidence flags for fallback-derived labels and avoid using them as hard truth in supervised evaluation.

## 4b) New Braun/scANVI pilot readout (post-100-epoch decision)

### Observed pilot behavior
- Pilot transfer runtime was fast (~7.4 minutes end-to-end) with biologically plausible class composition and high class confidence (mean ~0.943; large majority >0.8).
- CellClass predictions were neural-dominant (~99% combined radial glia / neuron / glioblast / neuroblast / neuronal IPC), which is directionally consistent with a control-only brain organoid atlas.
- Region transfer was intentionally softer (mean confidence ~0.55), with plausible forebrain-heavy plus midbrain/diencephalon/medulla distribution; this pattern is expected because regional patterning in organoids is less discrete than major lineage identity.

### Caveats to carry into the full run
1. Immune/vascular classes at ~0% in the pilot are a key uncertainty for Benchmark Q2 (microglia/endothelium correspondence). A 200k pilot can miss rare classes; absence/presence should be decided only on the full atlas run.
2. Validation against existing organoid labels was weak because many deposits retain coarse `cell_type='unknown'` at the harmonized layer; transfer validation should rely on finalized per-GSM fields where available.

### Required pre-full-run check (annotation availability)
- Confirm finalized annotation fields are present and populated in the full atlas object before benchmarking transfer outputs. The finalized per-cell fields expected from the GSM reconciliation pipeline are:
  - `cell_type_origin`, `age_days`, `organoid_type`, `protocol`, `unguided`, `multi_lineage`, `vascularized`, `slice`, plus provenance `annotation_level` and per-cell `gsm`.
- Perform benchmark stratification by `annotation_level` (`gsm` vs `deposit`) so transfer quality claims are not inflated/deflated by coarse fallback labels.

## 5) Reproducibility / portability (absolute paths, hard-coded interpreter)

### Current issue
- Several scripts hard-code local paths/interpreters (e.g., `/Users/eg/...`, fixed Rscript/Python paths).

### Risks
- Non-portable outside one machine.
- Hard to run in CI/HPC/containers.
- Higher risk of silent path drift and partial rebuild mismatches.

### Recommended adjustments
1. Centralize configuration:
   - `--project-root`, `--input`, `--output`, `--config` CLI flags with sensible defaults.
   - environment-variable fallback (e.g., `BRAIN_ORG_ROOT`).
2. Remove interpreter path assumptions from scripts; rely on active environment and document required versions.
3. Add a lockfile/environment spec (`environment.yml` or `pyproject` + pinned versions) and a reproducible runbook.
4. Version every major artifact with a manifest hash (config + code commit + input manifest + timestamp).

## Methodological risks to flag now
- Residual batch structure appears largely structural/confounded with biology (protocol/age/region composition), not mainly an optimization-length issue.
- Feature-selection bias from micro-batch pruning in highly imbalanced datasets.
- Demux-dependent contamination control is robust only if metadata quality and singlet policies are explicitly audited.
- Mixed annotation granularity can confound downstream benchmarks if not stratified.

## Simple / standard alternatives
- Keep `cell_ranger` HVG + NB scVI as baseline, but add one standardized sensitivity matrix (HVG count × latent dim × epoch/early-stop) and pick by objective metrics.
- Use fully per-cell benchmarking subsets as primary headline results; use fallback-labeled cells only in secondary analyses.
- Migrate scripts to parameterized CLIs and one orchestration entrypoint (Makefile/Snakemake) before the next full rebuild.

## Concrete next steps (recommended order)
1. Parameterize paths/interpreters in preprocess/train/fix scripts.
2. Add `seurat_hto_control` loader path in core dispatcher with validation/audit outputs.
3. Run 200k–500k pilot sweeps (HVG 3k vs 5k; latent 20/30/50; early stopping).
4. Keep/use the 100-epoch converged model for Braun 2023 label transfer and proceed with downstream analyses.
5. If integration quality is revisited, evaluate coarser batch keys or scANVI + scIB metrics rather than extending epochs further.
6. Run Braun transfer on the full atlas and explicitly audit whether rare immune/vascular classes emerge at full scale.
7. Re-run label transfer benchmarks with stratification by `annotation_level` and with coverage reports for finalized per-GSM fields.
8. Freeze a release bundle: code commit, configs, manifests, metrics, and artifact checksums.


## Appendix — git reproducibility / checkout note

- Commit IDs shown in agent logs may come from an ephemeral local branch and may not exist on your machine unless that exact commit was pushed to the remote.
- If `git checkout <sha>` fails with `pathspec ... did not match`, first fetch remote heads/PR refs and then checkout a remote-tracking ref that contains the change (for example `origin/<branch>` or `origin/pr/<id>`), or cherry-pick the commit that *does* exist in your repo history.
- For this repository snapshot, the method-review update is present at commit `87bd618` on the current branch history.
