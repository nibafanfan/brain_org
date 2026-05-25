# Loader format reference

Concrete recipe per source-data format encountered in the brain organoid atlas. Use this when picking a `loader_type` for a new deposit, or when extending the rebuild script.

---

## Standard 10x trio — the canonical case

Three files per sample:

```
GSM####_barcodes.tsv.gz   ← cell barcodes (one per line, 16-bp)
GSM####_features.tsv.gz   ← gene IDs (col 1 = Ensembl ID, col 2 = symbol, col 3 = type)
GSM####_matrix.mtx.gz     ← Matrix Market sparse: genes × cells, integer counts
```

Loader:

```python
from scipy.io import mmread
import anndata as ad, pandas as pd
from scipy.sparse import csr_matrix

X = mmread(mtx_path).T.tocsr().astype('int32')  # transpose to cells × genes
barcodes = pd.read_csv(bc_path, sep='\t', header=None)[0].tolist()
features = pd.read_csv(feat_path, sep='\t', header=None)
adata = ad.AnnData(
    X=X,
    obs=pd.DataFrame({'sample_id': label, 'gsm': gsm}, index=barcodes),
    var=pd.DataFrame({'gene_symbol': features[1].values}, index=features[0].values))
```

For deposits with **per-GSM trios** (each GSM has its own file set), loop and `ad.concat(ads, join='outer')` on the cell axis.

Examples: `cai2025`, `gse276558`, `gse242329`, most modern deposits.

---

## Non-standard format families

### A. Drop-seq DGE text (`*.dge.txt.gz`)

Genes as rows, cells as columns, tab-delimited, integer counts.

```python
df = pd.read_csv(dge, sep='\t', index_col=0)   # genes are rows
X = csr_matrix(df.values.T.astype('int32'))    # transpose to cells × genes
adata = ad.AnnData(X=X, obs={'gsm': gsm}, var=df.index)
```

When files are large (>30 MB compressed each), use chunked streaming with `chunksize=` to avoid OOM. Apply CellRanger-knee filter (`n_counts ≥ 500, n_genes ≥ 200`) at the per-file boundary before concat.

Examples: `gse86153`, `gse132105`, `gse171344` (36 streamed files).

### B. Single CSV per sample (gene × cell or cell × gene)

Orientation depends on the deposit. Auto-detect by shape:

```python
df = pd.read_csv(path, index_col=0)
if df.shape[0] > df.shape[1]:
    X = csr_matrix(df.values.T.astype('int32'))  # genes were rows
    cell_ids = df.columns; gene_ids = df.index
else:
    X = csr_matrix(df.values.astype('int32'))
    cell_ids = df.index; gene_ids = df.columns
```

Examples: `gse224346`, `gse146878 (SC_CTRL3, SC_FXS3)`.

### C. Streaming for OOM-sized dense CSVs (>1 GB compressed)

```python
chunks = pd.read_csv(path, chunksize=2000, index_col=0)  # rows = genes
sparse_blocks = [csr_matrix(c.values.astype('int32')) for c in chunks]
X = vstack(sparse_blocks).T.tocsr()  # transpose to cells × genes
```

Memory peak ~6 GB instead of OOM on the ~14 GB dense matrix.

Example: `gse180122_alsftd` batch 1 (1.8 GB compressed, 24,583 genes × 105,071 cells).

### D. Alevin/salmon output — **axes reversed from CellRanger**

```
GSM####_quants_mat.mtx.gz       ← CELLS × genes (note: opposite of 10x!)
GSM####_quants_mat_rows.txt.gz  ← cell barcodes
GSM####_quants_mat_cols.txt.gz  ← gene IDs with version suffix (ENSG00000123456.5)
```

Do **not** transpose. Strip the version suffix from gene IDs:

```python
import re
_ENSG = re.compile(r'ENSG\d{11}')
def strip_version(g):
    m = _ENSG.search(str(g))
    return m.group(0) if m else None
```

Example: `gse253230_ube3a` (rebuilt from Alevin output after detecting the v1 conversion had swapped axes).

### E. Seurat RDS files (R-interop pipeline)

Use the `Rscript`-based pipeline at `/tmp/rds_to_mtx.R`:

```python
import subprocess
subprocess.run(['Rscript', '/tmp/rds_to_mtx.R', rds_path, out_dir])
# R script: readRDS → extract @assays$RNA@counts → write mtx trio + metadata.tsv
# Then load the resulting mtx trio via standard 10x loader (above).
```

**Double-gzipped RDS gotcha:** some deposits ship `*.rds.gz` where the inner `.rds` is *itself* gzipped. The Python wrapper handles this by doing one manual `gunzip` step before invoking the R script (since `readRDS(gzfile(...))` chokes on double-compression).

R 4.4.1 + Seurat 5.3.0 + Matrix installed locally.

Examples: `gse277968_4protocols`, `gse231546_arid`, `gse163952`, `gse219245`, `gse306010`, `gse310490`, `gse325956`, `gse271116_pd_midbrain`.

### F. Pre-built h5ad

```python
adata = ad.read_h5ad(h5ad_path)
```

After read, sanitize: coerce object/category obs columns to str with NaN→'' for HDF5 vlen-string compatibility (`sanitize_for_h5()`).

**Axis-swap detection:** if `obs_names` look like Ensembl IDs and `var_names` look like barcodes, the original loader transposed it. Rebuild from raw or do `adata = adata.T`.

**Misnamed extensions:** some files have `.h5ad` but are actually CellRanger v2 `.h5` (with `/<reference_name>/` group containing `barcodes`, `data`, `indices`, `indptr`, `shape`, `genes`, `gene_names`). Use `h5py` directly:

```python
with h5py.File(path) as h:
    g = h['unknown']  # or whatever the root key is
    X = csc_matrix((g['data'][:], g['indices'][:], g['indptr'][:]),
                   shape=tuple(g['shape'][:])).T.tocsr().astype('int32')
    barcodes = [b.decode() for b in g['barcodes'][:]]
    genes = [b.decode() for b in g['genes'][:]]
```

Example: `gse195510` was a v2 h5 misnamed as h5ad.

### G. CellRanger h5 (`filtered_feature_bc_matrix.h5`)

```python
import scanpy as sc
adata = sc.read_10x_h5(path)
```

Used for many CellRanger 3+ outputs.

### H. Custom format quirks (per-deposit special handling)

- **Numeric placeholders** (`var_names = '1','2',...,'36601'`): the upload lost gene names entirely. Borrow var order from another deposit using the same 10x reference. For 36,601 genes the donor is `gse171263` (10x GRCh38-2020-A). Example: `gse165577`.
- **Ensembl with build prefix** (`hg19_ENSG...`, `GRCh38_ENSG...`, `GRCh38_2020A_v5_ENSG...`): strip via regex `_ENSG_RE.search(name)`. Examples: `gse165975`, `gse190815`, `gse208710`.
- **Versioned Ensembl** (`ENSG00000123456.5`): strip `.<digits>` via the same regex. Examples: `gse183627_kmt2d`, `gse295097_radiation`, `gse301348_h5n1`.
- **HTO-multiplexed pools**: don't attempt cell-level demux (HTO tag matrix isn't deposited). Mark the whole pool as `is_control=True` only at the sample level with a `hto_pool=True` caveat flag. Example: `gse171344`.
- **Series-level merged matrix** (instead of per-GSM trios): `GSE####_barcodes/genes/matrix.gz` at series root. Treat as one large sample; cell-line / condition info has to come from barcode prefixes in `obs_names`. Example: `gse134049`.

### I. Permanently rejected formats

These can't drive a per-cell atlas:

| Format | What it is | Example |
|---|---|---|
| Pseudobulk TSV | STAR/salmon/DESeq2 length-scaled counts collapsed to per-sample columns | `gse288165` |
| xlsx with TPM/RPKM | CLC Workbench bulk RNA-seq, one per-sample workbook | `gse167208` |
| ATAC peaks `.bed.gz` | open-chromatin region calls, not RNA | `gse97882` |
| MERFISH/Xenium HDF5 | spatial transcriptomics with ~140 probe panel | `gse271118_spatial` |
| `*.fastq.bedGraph.gz` | coverage tracks from FASTQ alignment, no count matrix | `gse135634` |
| FASTQ-only SRA deposit (no processed counts) | needs `cellranger count` to be usable | `gse280341_msm`, many in `_pending_fastq_reprocess/` |

---

## Filter step (after loading, before write)

For the **control-group atlas rebuild**, do the GSM filter on the loaded AnnData:

```python
allowed_gsms = set(config_row['gsm_filter'].split(','))  # or skip filter if 'ALL'
adata = adata[adata.obs['gsm'].astype(str).str.upper().isin(allowed_gsms)].copy()
```

Then mark `adata.obs['is_control'] = True` (the filter already excludes non-controls).

## Uniform QC

After GSM filter, apply consistent QC:

```python
counts = np.asarray(adata.X.sum(axis=1)).ravel()
ngenes = np.asarray((adata.X > 0).sum(axis=1)).ravel()
# MT mask: use HNOCA canonical for Ensembl namespace, str.startswith('MT-') for symbols
mt_mask = build_mt_mask(adata)  # see scripts/migrate_to_hnoca_schema.py
mt_counts = np.asarray(adata.X[:, mt_mask].sum(axis=1)).ravel()
pct_mito = mt_counts / np.maximum(1, counts) * 100

keep = (counts >= 500) & (ngenes >= 200) & (pct_mito <= 20)
adata = adata[keep].copy()
adata.obs['n_counts'] = counts[keep]
adata.obs['n_genes'] = ngenes[keep]
adata.obs['pct_mito'] = pct_mito[keep]
```

## HNOCA-compatible schema (rebuild output)

Every rebuilt h5ad should have:

- `X` = raw integer counts (CSR sparse).
- `var.index` = HGNC symbols (map via `data/reference/hnoca_var_canonical.tsv`).
- `var` columns: `ensembl`, `gene_symbol`, `gene_length`, `mt` (bool), `highly_variable`, `highly_variable_nbatches`.
- `obs` columns: `sample_id`, `gsm`, `is_control=True`, `n_counts`, `n_genes`, `pct_mito`, `organoid_type`, `accession`, `dataset_slug`, `dataset_filter`, plus HNOCA harmonization (`cell_type='unknown'`, `cell_type_original`, `cell_line`, `cell_line_original`, `bio_sample`, `tech_sample`, `batch`, `individual`, `disease='healthy'`, `assay_sc`).
- `layers['counts_lengthnorm']` for Smart-seq2 deposits only.

See `scripts/migrate_to_hnoca_schema.py` for the canonical implementation; the rebuild script should fold this in.

---

## Existing templates to reuse

| Script | Format handled |
|---|---|
| `/tmp/conv_rds.py` + `/tmp/rds_to_mtx.R` | Seurat RDS (E) |
| `/tmp/run_a2.py` | Streaming CSV (C) |
| `/tmp/run_a3.py` | Per-file HTO DGE (A + filter) |
| `/tmp/convert_std.py` | Standard 10x mtx h5 (top + G) |
| `/tmp/conv_amb.py`, `/tmp/conv_amb2.py` | Ambiguous / mixed-format handlers (H) |
| `scripts/migrate_to_hnoca_schema.py` | Gene namespace + HNOCA obs (post-load schema) |
