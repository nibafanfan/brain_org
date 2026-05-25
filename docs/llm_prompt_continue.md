# Prompt: continue the brain organoid atlas project

Paste this prompt into a new LLM session to pick up the project. Adjust the
"Today's task" section to whatever the next priority is.

---

## Project context

You are working on `/Users/eg/brain_organoid` — a brain-organoid scRNA-seq
benchmark project that **extends the HNOCA atlas (He, Dony, Fleck et al.
2024, *Nature*)** to a broader collection of brain-organoid protocols, with
a particular focus on **multi-lineage protocols published in 2024–2026** and
their correspondence to primary developing human brain (Braun et al. 2023,
*Science*).

The three benchmark questions driving every modeling decision are:

1. Do multi-lineage protocols (vascularized + microglia + multi-region) show
   **improved cell-type coverage** vs single-lineage HNOCA protocols?
2. Do the **added cell types** (endothelium, microglia, pericytes, regional
   neurons) transcriptomically correspond to their primary brain counterparts?
3. What **cell-type gaps remain** even in multi-lineage systems — i.e. what
   should future protocols target?

Healthy/control samples drive the integration; disease-model deposits are
kept with WT/isogenic-correction controls tagged per-sample.

## Current state (2026-05-14)

- **137 deposits in atlas** (`data/manifest.tsv`)
- **9,507,920 cells**, 6,305,852 control cells (66.3%), 992 samples (684
  control)
- Per-deposit AnnData files in `data/processed/<slug>.h5ad`
- Reference atlases downloaded: Braun 2023 (1.67M cells,
  `data/raw/braun_2023/braun_all.h5ad`) and HNOCA
  (`data/raw/hnoca_2024/hnoca_cleanedmeta.h5ad`)
- 9 deposits in `data/_pending_fastq_reprocess/` (need CellRanger from SRA)
- 6 deposits in `data/_blocked/` (4 permanent rejects + 2 deferred)
- 187-row candidate xlsx at `data/brain_organoid_partial_annotate.xlsx`
  with 28 manually-annotated rows; auto-annotated draft for the rest is at
  `data/brain_organoid_annotated_for_colleague.xlsx`

## Required reading before you touch anything

In this order:

1. **`docs/atlas_overview.md`** — what's in the atlas, schema, caveats, how
   to load. Start here.
2. **`docs/atlas_compilation_notes.md`** — narrative of how it was built;
   format-specific handlers, cleanup decisions, deposit-level quirks. Read
   §3 (schema), §4.1 (filtering decisions), §4.2 (control vs disease), §4.6
   (R-interop pipeline). This is the authoritative cleanup record.
3. **`docs/data_sources.md`** — per-paper accession map. Use this to find
   which paper a deposit came from.
4. **`docs/proposal.md`** — original project framing.
5. **`data/manifest.tsv`** — authoritative deposit index. Every cleanup or
   integration step should keep this in sync.

## Canonical `obs` schema (every `.h5ad` already has these)

| Column | Notes |
|---|---|
| `sample_id`, `gsm`, `is_control` | sample-level (HTO pools flagged as caveat) |
| `n_counts`, `n_genes`, `pct_mito` | per-cell QC |
| `organoid_type`, `accession`, `dataset_slug`, `dataset_filter` | identity |

`X` is always raw integer UMI counts (CSR sparse). `var_names` is
heterogeneous (mix of Ensembl IDs and HGNC symbols) — **gene namespace
harmonization has NOT been done yet** and is a prerequisite for integration.

## How HNOCA harmonized (we should mostly follow this)

HNOCA's atlas (`hnoca_cleanedmeta.h5ad`) is the template for what our final
integrated atlas should look like. Key elements:

- **CELLxGENE-style obs** — 48 columns with `_original` versions preserved
  (`cell_type` + `cell_type_original`, etc.). Sample hierarchy
  (`bio_sample`, `tech_sample`, `batch`, `individual`).
- **Gene namespace**: `var.index` is HGNC symbols; `var['ensembl']` is
  parallel Ensembl IDs; `var['gene_length']` is bp.
- **Layers**: `X` is raw counts; `layers['counts_lengthnorm']` is
  length-normalized for cross-tech (Smart-seq2 vs 10x) comparability.
- **Integration**: `obsm['X_scpoli']` is the scPoli latent embedding (NOT
  scVI). `obsm['X_umap_scpoli']` is UMAP from that.
- **Annotation**: `annot_level_1` → `annot_level_4_rev2` hierarchical
  (Snapseed); plus `annot_region_rev2`, `annot_ntt_rev2`.
- **HVG selection**: cross-batch (`var['highly_variable_nbatches']`
  records how many of 36 batches each HVG is variable in).

## Today's task

**Edit this section before pasting.** Pick one of:

### A. Finish cleanup
- Streaming-CSV loader for `gse180122_alsftd` batch 1 (~100k cells, 1.8 GB
  CSV that OOMs on naive `pd.read_csv`). Use `pd.read_csv(chunksize=)` with
  `scipy.sparse.lil_matrix` accumulation. See `data/_blocked/gse180122_alsftd/`.
- Chunked-sparse loader for `gse171344` (36 HTO-multiplexed Drop-seq DGE
  files, ~1.4 GB total). See `data/_blocked/gse171344/`.
- Spot-check the 159 auto-annotated rows in
  `data/brain_organoid_annotated_for_colleague.xlsx` against paper Methods
  sections — calibrated accuracy is ~66% overall, with Multi-Lineage at
  89% and Cell Type Origin at 78%.

### B. Schema migration to HNOCA-compatible
- Write `scripts/migrate_to_hnoca_schema.py`. For each `data/processed/*.h5ad`:
  - Add `_original` versions of existing obs columns
  - Add `bio_sample`, `tech_sample`, `batch`, `individual` columns
    (mostly = `gsm` for single-sample deposits)
  - Add `cell_type` (initially empty / `"unknown"`, populated after Snapseed)
- Gene namespace harmonization: build a mapping from Ensembl ↔ HGNC using
  HNOCA's `var` as ground truth. Rewrite `var.index` to HGNC symbols,
  populate `var['ensembl']` from the Ensembl ID. Drop genes not in HNOCA's
  set (or with `var['highly_variable_nbatches'] < 3` if using HVG filter).
- Add `var['gene_length']` from Ensembl (use `pyensembl` or download GTF).
- Add `layers['counts_lengthnorm']` for Smart-seq2 deposits.

### C. Integration via scPoli (transfer learning from HNOCA)
- Use scPoli's reference-mapping mode rather than retraining from scratch.
  scPoli model + reference embedding live in
  `data/raw/hnoca_2024/` (check the theislab GitHub for the exact files).
- For each of our 137 deposits: project onto HNOCA's latent space, save
  `obsm['X_scpoli']` per deposit and concat into a single atlas h5ad.
- Compute neighbor graph + UMAP in the integrated space.
- **Hardware**: Mac local CPU is fine for projection; training is not
  needed if using the published HNOCA scPoli model.

### D. Benchmark analysis (the headline)
Once integration is done:
- **Q1 (coverage)**: per-protocol, count primary brain cell types that have
  ≥N matched organoid cells (N=10 or 50). Compare single-lineage HNOCA
  protocols vs the multi-lineage subset (deposits flagged
  `Multi-Lineage=Yes` in the annotation xlsx). Output: a table and bar
  chart of coverage by protocol category.
- **Q2 (correspondence)**: for each added cell type (endothelium,
  microglia, pericytes, regional neurons), measure transcriptomic
  similarity to primary counterpart in Braun. Use kNN purity, marker-gene
  correlation, embedding distance.
- **Q3 (gaps)**: invert Q1 — list primary cell types with zero/near-zero
  representation in *any* protocol. This is the punch list for future
  protocol development.

### E. DL analysis (optional, if Q1-Q3 outputs need a deeper model)
- Train a small MLP or VAE on the integrated embedding to predict primary
  brain cell-type from organoid cell features. Use as a "mapping
  confidence" classifier.
- Or train an SAE (sparse autoencoder) on the integrated counts and look at
  per-protocol activation patterns of the learned features.

## Working conventions

- **Per-tool-call focus**. Don't bundle "download + convert + update docs"
  into one call. Each tool call does one logical step.
- **Update `data/manifest.tsv`** when a deposit gets added/changed/removed.
- **Keep `docs/atlas_compilation_notes.md` in sync** for narrative changes.
- **Update `docs/atlas_overview.md`** headline numbers (deposit count, cell
  count, control %) at the end of any session that changes them.
- **Healthy controls only** for atlas integration; disease samples kept in
  the AnnData but filtered out at integration time via `is_control`.
- **Don't re-download what's on disk.** Check `data/raw/<slug>/` before
  fetching from GEO.
- **scPoli, not scVI**, for HNOCA-compatibility — unless there's a specific
  reason to retrain.

## Known caveats (skim before you commit to a method)

- **Gene namespace heterogeneity is real and unfixed.** Some deposits use
  Ensembl IDs in `var_names`, others use HGNC symbols. Harmonize before
  integrating.
- **HTO-multiplexed deposits** (`gse189535`, `gse281452_iMG`,
  `gse296775_strada`, `gse297594_mecp2`) have `is_control=True` flagged
  at the sample level but the cells inside are a mixed cohort. Don't trust
  the cell-level control labels until you re-demultiplex the HTO tags.
- **`gse282644_hiv` is broken** (12 cells, orientation error in source
  CSV). Skip it.
- **Auto-annotated xlsx has ~66% strict match** with manual annotations.
  Multi-Lineage and Cell Line are the most reliable fields (~80-89%);
  Age, Protocol, Notes need most review.

## Output expectations

Every session should produce, at minimum:

1. **Code or AnnData artifact** matching the canonical schema.
2. **`data/manifest.tsv` row** added/updated for any new/changed deposit.
3. **Short narrative** appended to `docs/atlas_compilation_notes.md` if a
   processing decision was made.
4. **Update** `docs/atlas_overview.md` headline counts if they changed.

## When to ask the user vs decide yourself

**Decide yourself**: code changes, format-specific loader choices, internal
file layout, function naming, library choices.

**Ask the user before**:
- Modifying or deleting any `.h5ad` already in the atlas
- Re-running a >1-hour compute job
- Pulling >10 GB of new data
- Changing the canonical obs schema
- Making methodological choices that affect the benchmark answers (e.g.,
  cell-type coverage threshold, integration method swap, control flag
  redefinition)

---

## Starter actions for the agent

When you first read this prompt:

1. `cat data/manifest.tsv | head -5` and `wc -l data/manifest.tsv` to confirm
   current deposit count
2. `head -50 docs/atlas_compilation_notes.md` to learn the established
   pipeline patterns
3. `ls data/_pending_fastq_reprocess/ data/_blocked/` to see open queues
4. Confirm with user which of A/B/C/D/E is today's task, then proceed
