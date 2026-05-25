# Brain organoid scRNA-seq atlas — current state

## Two parallel atlases on disk

| Location | Schema | Use for |
|---|---|---|
| `data/processed/` (v1) | ad-hoc obs (sample_id, gsm, is_control, …); var is mix of Ensembl IDs and HGNC symbols | per-deposit work, before-migration analyses, rollback |
| `data/processed_v2/` (v2, **HNOCA-compatible**) | CELLxGENE-style obs with `_original` preservation; var is HGNC symbols in HNOCA's 36,842-gene set; has `gene_length`, `counts_lengthnorm` layer for Smart-seq2 | **integration with HNOCA**, scPoli transfer learning, atlas-level analyses |

### `data/processed/` (v1) — 137 deposits

**137** `.h5ad` files (updated 2026-05-14). **Each file is one deposit** (one paper's scRNA-seq dataset), pre-loaded as an [AnnData](https://anndata.readthedocs.io) object with:

- `X` — raw integer UMI counts (cells × genes, CSR sparse)
- `obs` — per-cell metadata with a canonical schema (see below)
- `var` — per-gene info (typically Ensembl IDs as `var_names`, gene symbols in `var['gene_symbol']`)

Each `.h5ad` is **independently loadable** — they were built so you can `ad.read_h5ad(...)` any one, or `ad.concat([...])` a set of them, without surprises.

The authoritative index is `data/manifest.tsv`: one row per deposit with `slug`, `accession`, `path`, `n_cells`, `n_genes`, `n_samples`, `n_control_samples`, `n_control_cells`, `organoid_type`, `filter`, `status`.

## Headline numbers

| Metric | Value |
|---|---:|
| Deposits | **137** |
| Total cells | **10,905,385** |
| Control cells | 7,641,866 (70.1%) |
| Total samples | 1,032 |
| Control samples | 722 |
| Largest deposit | `gse171344` — 1,404,709 cells (HTO Drop-seq, 36 pools) |
| Smallest useful deposit | `gse195692` — 72 cells (Smart-seq2, tiny pilot) |

## Coverage

Brain organoid scRNA-seq spanning **2017–2026** publications. Top organoid types by deposit count:

- **cerebral** (54), **cortical** (8), **midbrain** (4), **cerebellar** (3)
- specialized: **hypothalamic, pineal, nucleus basalis, ganglionic eminence, thalamocortical, hippocampal, brainstem, spinal, telencephalic**
- multi-lineage / multi-modal: **vascularized, microglia-containing, oligodendrocyte-bearing, neural-vascular, cortical-assembloids, multi-region, multi-protocol**

Includes the original 36 datasets from HNOCA (He et al. 2024), plus ~100 deposits from 2024–2026 that postdate HNOCA. Disease-model deposits are kept with healthy/WT controls tagged per-sample.

## Canonical `obs` schema (every h5ad has these)

| Column | Meaning |
|---|---|
| `sample_id` | per-deposit sample label (often GSM title or condition name) |
| `gsm` | NCBI GEO sample ID, or `'unknown'`/`'pool'` for aggregated cases |
| `is_control` | bool — True for healthy/wild-type/vehicle samples (sample-level, not always cell-level) |
| `n_counts`, `n_genes`, `pct_mito` | per-cell QC |
| `organoid_type` | coarse label (e.g. `cerebral`, `midbrain`, `cortical_assembloid`) — best-guess from paper title/Methods |
| `accession` | GEO accession (e.g. `GSE197887`) |
| `dataset_slug` | matches filename + manifest row |
| `dataset_filter` | `"authors-filtered"` (deposit was already cell-called) or `"min_counts=500, min_genes=200"` (we applied threshold filter) |

Deposit-specific extra columns are present when the GEO Sample Characteristics carried useful metadata: `genotype`, `diagnosis`, `treatment`, `cell_line`, `protocol_age_days`, etc. — not standardized; whatever the GEO record carried.

## Loading example

```python
import anndata as ad
import pandas as pd

manifest = pd.read_csv('data/manifest.tsv', sep='\t')

# Load one
a = ad.read_h5ad(manifest.query('slug == "gse271116_pd_midbrain"').path.iloc[0])

# Healthy controls only across multiple deposits
healthy_only = a[a.obs['is_control']].copy()

# Concat a bunch (canonical schema, so this works)
adatas = [ad.read_h5ad(p) for p in manifest.query('organoid_type == "cerebral"').path]
combined = ad.concat(adatas, axis=0, join='outer', merge='same')
```

## Reference atlases (for benchmarking)

| | Location | Cells |
|---|---|---:|
| **Braun et al. 2023** (primary first-trimester developing human brain) | `data/raw/braun_2023/braun_all.h5ad` | 1.67M |
| **HNOCA 2024** (He et al. — integrated organoid atlas, scPoli-integrated) | `data/raw/hnoca_2024/hnoca_cleanedmeta.h5ad` | 1.77M |

These are the references the organoid deposits get mapped against for the benchmark; not part of the 137-deposit count.

## Caveats worth knowing before integrating

1. **Gene namespace is heterogeneous.** Most deposits use Ensembl IDs in `var_names`, some use HGNC symbols; gene_symbol column normalizes display but not the index. **Atlas-level Ensembl harmonization hasn't been done yet** — recommend using Braun's var_names as canonical and dropping ultra-rare genes during integration.

2. **Cell calling varies.** ~20 deposits were "raw droplets only" deposits where we applied `n_counts ≥ 500 AND n_genes ≥ 200` (a CellRanger-knee approximation). The other ~117 used the authors' own cell calls. The `filter` column in manifest distinguishes them. The simple-threshold ones include more borderline cells.

3. **`is_control` is per-sample, not per-cell, for HTO-multiplexed pools.** A few deposits (`gse189535`, `gse281452_iMG`, `gse296775_strada`, `gse297594_mecp2`) demultiplex Dex+Vehicle conditions within the same GSM. The whole sample is flagged `is_control=True` but the cells inside are a mixed cohort. Don't treat their cell-level control labels as ground truth without re-demultiplexing the HTO tags.

4. **Disease state isn't a top-level obs column yet.** Right now it's implicit in `is_control` + deposit-specific genotype/diagnosis columns. We recommend adding an explicit `disease_state` (`healthy` / `<disease_name>` / `unknown`) at integration time.

5. **`gse282644_hiv` is broken** (12 cells; CSV orientation error in the source deposit). Skip it.

## What's not in the atlas

- **9 deposits** in `data/_pending_fastq_reprocess/` — ship only normalized/CPM values. To use, pull raw FASTQ from SRA and run CellRanger. Listed in `docs/data_sources.md`.
- **`data/_blocked/` (6 dirs)** — the only place that holds *genuinely* unconverted deposits. Each has a `REASON.md`:
  - **4 permanent rejects**: `gse97882` (ATAC peaks), `gse135634` (bedGraph), `gse167208` (bulk RNA-seq), `gse271118_spatial` (MERFISH). Two of these (`gse97882`, `gse271118_spatial`) have **sibling accessions that DID make it into the atlas** — the paper is represented, just via a different GEO ID (`gse98201_mge_cortical`, `gse271116_pd_midbrain`).
  - **2 deferred** (need a smarter loader): `gse171344` (HTO-multiplexed Drop-seq, pandas OOMs on the 36 DGE files); `gse180122_alsftd` (batch 2 is in atlas already; batch 1 needs streaming CSV).

> **Note on raw data:** previously-blocked deposits that were *recovered* keep their raw download in `data/raw/<slug>/`, not in `_blocked/`. `_blocked/` only holds the 6 dirs above — anything else with raw data is in `data/raw/`.

## `data/processed_v2/` (v2 — HNOCA-compatible)

**137 deposits** migrated to be `ad.concat`-compatible with HNOCA — full coverage as of 2026-05-15. (First pass migrated 129; the remaining 8 had gene-namespace edge cases — Ensembl version suffixes, build prefixes, numeric placeholders, transposed axes — all fixed; details in `atlas_compilation_notes.md` §5.5.)

What's different from v1:

```python
# var.index = HGNC symbols matching HNOCA's 36,842-gene set
# var columns: ensembl, gene_symbol, gene_length, mt, highly_variable*, hv_nbatches
# obs columns add HNOCA-style harmonization with _original preservation:
#   cell_type + cell_type_original, cell_line + cell_line_original,
#   disease + disease_original, assay_sc + assay_sc_original,
#   bio_sample, tech_sample, batch, individual
# layers['counts_lengthnorm'] present for Smart-seq2 deposits only
```

Index: `data/manifest_v2.tsv` (129 rows). Load with the same `ad.read_h5ad(path)` pattern as v1.

## Where to find more detail

- **Per-deposit narrative** including format quirks, cleanup decisions, and known issues: `docs/atlas_compilation_notes.md`
- **Per-paper accession reference** (paper → GEO accession → which file we used): `docs/data_sources.md`
- **What got rejected and why** (with sibling-accession recoveries called out): `data/triage_rejected.md`
