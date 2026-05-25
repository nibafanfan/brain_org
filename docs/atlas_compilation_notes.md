# Brain organoid atlas — compilation notes for collaborators

**Status snapshot (2026-05-15, after Phase A cleanup + Phase B schema migration + gene-namespace patches):** **137 deposits**, **10,905,385 cells**, 7,641,866 control cells (70.1%), 1,032 samples (722 control). All 137 deposits now also present in v2 schema (`data/processed_v2/`). Plus 6 deposits in `_blocked/` (4 permanent rejects + 2 deferred) and 9 in `_pending_fastq_reprocess/`.

> **Cell-count correction (2026-05-15):** the original v1 conversion of `gse253230_ube3a` had inflated the cell count to 151,948 because the loader treated per-sample gene-row blocks as cells. Rebuild from raw Alevin output gave the correct count of 39,646 cells × 61,487 Ensembl IDs (post-version-strip), across 9 samples (`6week` + `11week_sub1..8`).

> **Recent additions since prior snapshot (121 → 137, +16):**
> - Cleanup recoveries: gse180122_alsftd batch 1+2 combined (121,480 cells; previously batch 2 only at 16k), gse171344 HTO Drop-seq (1.4M cells across 36 multiplexed pools).
> - Sibling-accession recoveries: gse98201_mge_cortical (paper at GSE97882 was ATAC-only), gse271116_pd_midbrain (paper at GSE271118 was MERFISH-only).
> - Auto-annotation batch: 7 deposits recovered from the "unaccounted" rejected list (gse86153, gse146878, gse220085, gse224346, gse306010, gse310490, gse325956).
> - Plus 7 recoveries from re-audited `_blocked/` deposits.
> - Manifest integrity fix on 2026-05-14: dropped duplicate `gse137877` row, broken `gse282644_hiv` (12-cell deposit), and 2 corrupt `gse171344` placeholder rows from prior partial conversions.

This is a working note for anyone joining the project mid-stream — it explains what was actually done to compile, clean, and annotate the atlas so you don't have to reverse-engineer it from code. Where the data has structural quirks you should know about *before* using it, those are flagged. The master record is `data/manifest.tsv`; this document is the narrative behind it.

---

## 1. Source material

- **Candidate list:** `data/brain_organoid.xlsx` — 189 paper rows curated by the project lead. Each row is a paper title + one or two GEO/SRA hyperlinks. The "Dataset" cells render as the text "GEO Accession viewer" but the real accession is in the hyperlink target URL.
- **Already-published organoid atlas:** HNOCA (He, Dony, Fleck et al. 2024, *Nature*) — 36 datasets, ~1.7M cells. Listed in the "HNOCA Paper Datasets" sheet; flagged as **skip** in our pipeline since they're already in HNOCA's published integration.
- **Primary brain reference:** Braun et al. 2023, *Science* — 1.67M first-trimester developing human brain cells. Downloaded to `data/raw/braun_2023/braun_all.h5ad` (11 GB).
- **HNOCA reference:** `hnoca_cleanedmeta.h5ad` (17 GB) downloaded to `data/raw/hnoca_2024/`. Carries Snapseed cell-type annotations.

## 2. Triage process

Each candidate row went through three filters:

1. **Format triage** (GEO/ENA API queries) — is this human single-cell? Is there a per-cell matrix or only pseudobulk? What's the file format?
2. **Scope triage** — brain organoid? Reject spinal/retinal/gut/cardiac/non-human. (Disease-model deposits are kept; we filter to controls at the sample level.)
3. **Compilation triage** — try to convert to AnnData; if format defeats automation, defer to a quarantine bucket.

Triage outputs are persisted as `data/triage_pass1.tsv`, `data/triage_pass2.tsv`, `data/triage_accepted.md`, `data/triage_rejected.md`. The rejected list went through a manual re-audit (user spot-checks) which recovered 26 of 38 originally-rejected deposits.

## 3. Per-deposit AnnData schema

Every h5ad in `data/processed/` follows the same `obs` schema, so you can `ad.concat()` them without surprises:

| Column | Type | Meaning |
|---|---|---|
| `sample_id` | str | Per-deposit sample label (often the GSM title or per-condition name) |
| `gsm` | str | NCBI GEO sample ID (e.g. `GSM8501109`). `'unknown'` or `'pool'` for aggregated/non-mappable cases. |
| `is_control` | bool | True if this cell's sample is a healthy/wild-type/vehicle control. Per-sample, not per-cell — see §4 for caveats. |
| `n_counts` | int | Total UMI/read count per cell |
| `n_genes` | int | Number of expressed genes per cell |
| `pct_mito` | float | Percent mitochondrial transcripts (MT-* genes); NaN if no mito genes detected |
| `organoid_type` | str | e.g. `cerebral`, `midbrain`, `cortical_assembloid`, `hypothalamic` (best guess from paper title/Methods) |
| `accession` | str | GEO accession (e.g. `GSE197887`) |
| `dataset_slug` | str | Internal slug matching the h5ad filename and `data/manifest.tsv` row |
| `dataset_filter` | str | Either `"authors-filtered"` (deposit was already cell-called) or `"min_counts=500, min_genes=200"` (we filtered raw droplets) — see §4 |

Many deposits also have deposit-specific obs columns when the GEO characteristics surfaced useful metadata: `genotype`, `diagnosis`, `treatment`, `cell_line`, `tissue`, `condition`, `protocol_age_days`, `condition_label`. These are not standardized — they're whatever the GSM characteristics field carried.

`X` is always raw integer UMI counts in CSR sparse format (one exception: GSE304918 had only log-normalized data and was quarantined to `_pending_fastq_reprocess/`; the remaining 112 are all raw).

`var.index` is whatever the original deposit used — Ensembl IDs (`ENSG…`) for most, gene symbols (`MALAT1`) for older ones, or a mix. `var['gene_symbol']` carries the human-readable symbol when the deposit provided it; otherwise it mirrors `var_names`. **Gene-namespace harmonization to a single Ensembl reference is a downstream integration step we have not yet done.**

## 4. Cleanup decisions made — pay attention before integrating

### 4.1 Empty-droplet filtering (5 deposits, ~27M nominal cells → ~500K real cells)

Several deposits uploaded the 10x **raw** `feature_bc_matrix` (every droplet ever recorded, mostly ambient RNA) instead of the **filtered** version (real cells called by CellRanger). For these we applied `n_counts ≥ 500 AND n_genes ≥ 200` as a CellRanger-knee approximation. The flag `filter_applied = "min_counts=500, min_genes=200"` distinguishes them in the manifest.

| Deposit | Pre-filter | Post-filter | Kept |
|---|---:|---:|---:|
| gse260532_kcnj2 | 6,794,880 | 7,098 | 0.10% |
| gse251684_striato_nigral | 3,983,974 | 30,145 | 0.76% |
| gse320222_cbp | 11,383,859 | 94,074 | 0.83% |
| gse231319_polaroid | 5,472,341 | 71,221 | 1.30% |
| gse286235_nbm | 5,438,521 | 33,661 | 0.62% |
| gse300486_gaucher | 9,108,846 | 253,228 | 2.78% |
| gse312664_microglia | 4,792,450 | 155,142 | 3.24% |

Total filtered: **20 deposits** (~20% of the atlas). The other **89 deposits** were `authors-filtered` (we trusted the deposit's cell-calling).

The simple thresholds are conservative. For higher-quality cell calling, swap in `emptyDrops` (R DropletUtils) or recompute against the published cell barcodes (which most papers don't expose). Treat the `min500c/200g` deposits' borderline cells with extra skepticism.

### 4.2 Disease-vs-control sample annotation

The project decision is "healthy/control samples only" for atlas integration, but disease-model deposits are kept because they carry healthy WT/isogenic-corrected controls. We mapped controls per-GSM from the GEO Sample Characteristics field. The detection logic looks at:

1. **Genotype** field — `WT`, `Wildtype`, `Parental`, `Isogenic`, `Euploid`, `Ctrl`, `Healthy control`, `+/+`, gene-corrected variants (`iCtrl`, `GC Miro1 R272Q`).
2. **Disease state / Diagnosis / Condition** — `Healthy control`, `Typical`, `Unaffected` vs `Case`, `Patient`, `Affected`, named disease.
3. **Treatment** — `None`, `DMSO`, `Vehicle`, `Untreated`, `Mock`, `Sham`, `Normoxia`, `Baseline`, `NoCort` vs anything else.
4. **Sample title keywords** — fallback heuristic for opaque GSM IDs (XHI002, KH003 etc.).

When the deposit was a **multiplexed pool** (hashtag-demultiplexed half-control / half-treated cells per GSM), we couldn't separate cell-level control vs disease without re-demultiplexing the HTOs. For those we flag the whole sample `is_control=True` with the understanding that it's actually a mixed cohort. **Affected deposits:** `gse296775_strada`, `gse297594_mecp2`, `gse281452_iMG`, `gse189535` (Dex+Vehicle pools), `gse285126_fndc4` (also blocked separately).

After initial conversion, several deposits showed `is_control = 0/N samples` because the keyword detector missed deposit-specific labels (`Ctrl`, `iCtrl`, `Parental`, `GC`, `NoCort`, `apoe3wt`, `nbh`, `NT*`, `Neg*`, `FCT`, `NBH`, `EpiC`, `WTS42`, etc.). Each was fixed by inspecting the GEO characteristics directly and overriding the GSM→control map. Across 3 fix passes we resolved **0 zero-control deposits** in the final atlas.

### 4.3 Format-specific handlers we built

The 121 conversions used multiple loader paths. If you re-process or extend, you'll likely need these for new deposits:

| Format | Example deposits |
|---|---|
| 10x mtx trio (standard `_barcodes/features/matrix`) | majority |
| 10x mtx trio with **dash separator** (`-barcodes.tsv.gz`) | GSE280812, GSE286054, GSE252522 |
| 10x mtx with **dot separator** (`.barcodes.tsv.gz`) | GSE164089, GSE241743 |
| 10x mtx with **long suffix** (`_human_filtered_feature_bc_matrix.mtx.gz`) | GSE273907 |
| 10x mtx nested in **per-sample subdirectory** | GSE304516 (`<gsm>/<label>/<label>_barcodes.tsv.gz`), GSE197887 (per-sample tar.gz extract) |
| 10x mtx with custom `_quants_mat_cols/_rows.txt.gz` naming | GSE253230 |
| 10x h5 (filtered_feature_bc_matrix.h5) | many |
| 10x h5 raw (need empty-droplet filter) | GSE286235, GSE198927 |
| Old 10x h5 (`filtered_gene_bc_matrices_h5.h5`) | GSE108571 |
| CellBender output h5 | GSE281452 |
| Per-sample zip with inner 10x trio | GSE281622 |
| Per-sample tar.gz with inner cellranger_outs | GSE208438, GSE183627 |
| Drop-seq DGE text (`*.dge.txt.gz`) | GSE207749, GSE208418, GSE237855, GSE195692 |
| Gene × cell CSV per GSM | GSE113089, GSE192405, GSE196423 |
| `.csv.bz2` (bzip2 not gzip) | GSE184878 |
| Aggregated cells-in-columns CSV (one big file) | GSE180122 (OOM'd, blocked) |
| h5ad gzipped per GSM | GSE239542 |
| h5ad pre-built whole-deposit | GSE304918 (normalized, quarantined), GSE195510 |
| Parse Biosciences split-seq with broken meta mapping | GSE285126 (blocked) |
| MERFISH/Xenium spatial transcriptomics panel | GSE271116/271118 (blocked, not scRNA-seq) |
| Seurat RDS only | GSE277968, GSE219245, GSE163952, GSE231546 — recovered via `Rscript → mtx → AnnData` pipeline (see §4.6) |

| Re-audited deferrals from new-68 (orientation-aware reads, GSE prefix variants, mixed format) | GSE187877, GSE183903, GSE243015, GSE165577, GSE132105 |

### 4.4 Quirks worth knowing per-deposit

- **`gse197887_typical_autism`** (10x, 10 samples): sample IDs are opaque codes (`XHI002` … `XHI061`). Control/case is determined by `diagnosis` field in each GSM. 5/10 are control.
- **`gse283473_ws`** (6 samples, Williams syndrome): all 6 GSMs labeled `disease state: Healthy control` so all marked is_control=True. The WS patient samples are elsewhere or never deposited.
- **`gse237133_midbrain`**: 3 GSMs — `Ctrl` (healthy), `PD Miro1 R272Q` (Parkinson patient), `iCtrl/GC Miro1 R272Q` (isogenic gene-corrected, which counts as a control). The fix uses isogenic correction logic; **2/3 samples are control**.
- **`gse180122_alsftd`**: only 1/2 batches converted (the other batch OOM'd on the 3 GB CSV and is blocked). The 16,410 cells are from C9 batch 2; ALS/FTD patient lines (`CS30/31/32/33`) co-exist with healthy controls (`EpiC*`, `WTS42*`) in the same matrix.
- **`gse282644_hiv`**: only **12 cells** — likely a transcript-count CSV that was oriented wrong on auto-detect. Effectively unusable; do not include in integration.
- **`gse195692`**: 72 cells across 12 samples. Smart-seq2-like format; per-sample cell counts are tiny.
- **`gse189535`** (16 samples, Glucocorticoid): each "sample" is a hashtagged pool of Dex+Vehicle conditions. All 16 marked is_control=True (≥half the cells are vehicle controls) but the labeling is sample-level, not cell-level. Downstream demultiplexing required for clean cell-level control vs treatment.
- **`gse195510`**: deposit ships a single h5ad named `…_COUNTS_FILTERED.h5ad`. Initial auto-loader produced an empty AnnData; cleanup pending — sits in `_blocked/`.
- **`gse286054_meth`**: titled "methamphetamine-induced neuroinflammation". `NT*` = no-treatment (control), `MA*` = methamphetamine. 8/14 are NT controls.
- **`gse303735_gwi`**: Gulf War Illness study. `vehicle_treated` = control, `toxicant_treated` = exposed. 4/8 controls.

### 4.5 Cell-count outliers in the manifest

- **`gse183627_kmt2d`** = 517,309 cells (6.2% of the atlas). Originally a 17 GB RAW.tar with mixed FASTQ-mtx-multiome content; we filtered to only scRNA-seq tarballs and Cell Ranger trios. Skewed cell-count is real.
- **`gse273907_oxygen`** = 278,137 cells (oxygen-challenge panel, 7 samples). Filtered_feature_bc_matrix per sample.
- **`gse300486_gaucher`** = 253,228 cells (Gaucher disease + spinal organoids, 16 samples). Filtered from 9.1M raw droplets.
- **`gse287254_3dprinted_vasc`** = 388,584 cells (perfused vascularized cerebral organoids, 8 samples).
- **`gse263652_cerebellar`** = 214,075 cells (medulloblastoma cerebellar organoids — but ALL 16 GSMs are labeled control here; verify before using as "healthy cerebellar reference").
- **`gse245719_thalcortx`** = 142,210 cells (thalamocortical assembloids, 22q11.2 deletion). Listed as 31/31 control samples — likely a pool labeling issue similar to STRADA; verify cell-level genotype before drawing conclusions about WT vs 22q11.2.

### 4.6 R-interop pipeline for Seurat RDS deposits (added 2026-05-13)

We have an `Rscript`-based pipeline at `/tmp/rds_to_mtx.R` that reads a Seurat `.rds` (or `.rds.gz` that's actually double-compressed — see below), extracts the raw `counts` slot, and writes a 10x-style mtx trio (`barcodes.tsv`, `features.tsv`, `matrix.mtx`) plus a `metadata.tsv` for `obj@meta.data`. Python then wraps the mtx into AnnData using the same `emit()` schema as the rest of the atlas.

**Gotcha:** several deposits ship `*.rds.gz` files where the inner `.rds` is *itself* gzipped (double-compressed). The Python wrapper handles this by doing one `gunzip` step manually before invoking the R script, since `readRDS(gzfile(...))` chokes on the double-compression.

**Deposits recovered via this pipeline (2026-05-13):**
- **GSE277968_4protocols** (69 per-GSM RDS files + 1 combined Seurat object): 91,260 cells, 4 protocols × 4 cell lines × replicates. This is the most architecturally important comparison study we had been blocked on.
- **GSE231546_arid** (multiome counts.tsv + RDS): 54,480 cells via the `cnts.tsv.gz` direct read.
- **GSE163952** (HSV1-infected cerebral organoids): 34,632 cells.
- **GSE219245**: 4,742 cells.

For future RDS deposits, the pipeline at `/tmp/conv_rds.py` is the template. R 4.4.1 + Seurat + Matrix are already installed locally.

### 4.6.5 Streaming-CSV loader for OOM-sized deposits (added 2026-05-14)

Two deposits that previously OOM'd on naive `pd.read_csv` were recovered with chunked streaming → CSR sparse construction:

- **gse180122_alsftd** batch 1 (1.8 GB compressed, 24,583 genes × 105,071 cells dense). Streamed in row-chunks of 2,000 genes via `pd.read_csv(chunksize=2000)`; each chunk converted to `csr_matrix` and `sp.vstack`'d. Final shape after concat with batch 2: **121,480 cells × 24,582 genes**. Memory peak ~6 GB instead of the ~14 GB OOM of dense parsing.
- **gse171344** HTO-multiplexed Drop-seq (36 DGE files, ~1.4 GB compressed total, up to 100k pre-filter droplets per file). Per-file chunked stream + `n_counts ≥ 500 AND n_genes ≥ 200` cell-call filter, then sparse concat across files. Final shape: **1,404,709 cells × 52,584 genes** — single largest deposit in the atlas. ~6.7h runtime serial.

The pattern (genes-as-rows chunked stream + sparse vstack + cell-call filter at the per-file boundary) is reusable for any future deposit shipping a single dense CSV that exceeds working memory. Script template lives at `/tmp/run_a2.py` / `/tmp/run_a3.py`.

**Manifest integrity fix (same date)**: Audit pass found 4 corrupt rows (column-shift errors, duplicates, malformed). Audit logic now in `data/manifest.tsv.bak_pre_a1` history; backup retained.

## 5. Quarantine — deposits NOT in the atlas

Each quarantined directory has a `REASON.md` next to the raw data.

### `data/_pending_fastq_reprocess/` (2)
Deposit ships only normalized counts. Need to pull raw FASTQ from SRA and re-run Cell Ranger.
- **gse304918_thyroid**: log-normalized AnnData only. SRA `PRJNA1303443`, 64.8 GB FASTQ, ~8 hr Cell Ranger.
- **gse280341_msm**: R RDS file with normalized Seurat counts.

### `data/_pending_r_interop/` (0)
**Empty as of 2026-05-13.** All four RDS-only deposits (`gse277968_4protocols`, `gse219245`, `gse163952`, `gse231546`) recovered via the new `Rscript → mtx → AnnData` pipeline (§4.5).

### `data/_blocked/` (~24 remaining after re-audit)
Format defeats automated conversion. After the re-audit + R-interop passes, the remaining blocked deposits are mostly fundamental data-quality issues. Highlights by reason:

| Reason category | Deposits |
|---|---|
| Not scRNA-seq (ATAC peaks, FASTQ bedGraphs, spatial transcriptomics) | gse97882, gse135634, gse271118_spatial |
| Smart-seq2 / plate-based with normalized-only output | gse75140, gse124299 (3096 GSMs, CPM only), gse185052, gse124174 |
| Single aggregated gene×cell matrix with normalized values (TPM/SCT/log-norm) | gse131094, gse150903 |
| OOM on large dense CSV | gse180122_alsftd (24K genes × 150K cells) — partial 16k-cell h5ad still in atlas from one of two batches |
| Mixed/complex format we haven't reverse-engineered | gse137877, gse167208 (xlsx), gse163018 |
| Parse split-seq with broken matrix↔meta mapping | gse285126 |
| Drop-seq / micro-Cellman naming with mixed types | gse195692 |
| Processed/metadata only — no raw counts | gse145306 (CPM only), gse181290, gse233295 |

The 9 deposits recovered between the previous snapshot and now (gse187877, gse183903, gse243015, gse165577, gse132105, plus the 4 from R-interop) have been **removed from `_blocked/`** and are now in `data/processed/`.

If a still-quarantined deposit becomes critical to a benchmark question, the path to recovery is in its `REASON.md`.

## 5.5 HNOCA-compatible schema migration (added 2026-05-14)

The atlas now has a **v2 schema variant at `data/processed_v2/`** with 129 deposits migrated to be `ad.concat`-compatible with the published HNOCA atlas. The original v1 schema in `data/processed/` is retained for rollback.

### What changed in v2

**Gene namespace**: `var.index` is HGNC symbols matching HNOCA's 36,842-gene set. Per-deposit, var was reindexed:
- Deposits with Ensembl-ID `var.index` → mapped via HNOCA's `var['ensembl']` (built into `data/reference/hnoca_var_canonical.tsv`)
- Deposits with symbol `var.index` → passed through, kept only symbols present in HNOCA

Each migrated h5ad has these new `var` columns: `ensembl`, `gene_symbol`, `gene_length` (bp), `mt` (bool), `highly_variable` (3,000 HVGs), `highly_variable_rank`, `highly_variable_nbatches`.

**Obs schema** — HNOCA-style harmonization with `_original` preservation:

| Column added | Source |
|---|---|
| `cell_type` + `cell_type_original` | existing `cell_type` (or `"unknown"`) preserved as `_original`; `cell_type` set to `"unknown"` pending Snapseed |
| `cell_line` + `cell_line_original` | existing `cell_line` preserved; HNOCA naming convention |
| `disease` + `disease_original` | inferred from `genotype`/`diagnosis` if present; defaults to `"healthy"` |
| `assay_sc` + `assay_sc_original` | `"10x 3' v3"` default; `"Smart-seq2"` for plate-based deposits |
| `bio_sample`, `tech_sample`, `batch`, `individual` | sample hierarchy (mostly = `gsm`, slug, slug, cell_line respectively) |
| `assay_differentiation`, `organoid_age_days`, `publication`, `sample_source`, `state_exact` | misc HNOCA columns |
| `organ` + `organ_original`, `organism` + `organism_original`, `sex` + `_original`, `ethnicity` + `_original`, `suspension_type` + `_original` | filled with sane defaults |

**Layers**: for Smart-seq2 deposits (`gse185052`, `gse195692`) only, added `layers['counts_lengthnorm'] = X / gene_length * mean(gene_length)` for cross-tech comparability.

### Migration outcome

- **137 / 137 deposits successfully migrated** (100% after 2026-05-15 patches; original first-pass was 129/137)
- Mean gene retention: **~80%** (range ~30% – 100%)
- Best retention: Ensembl-namespace deposits that already used GENCODE-compatible IDs (eg `gse197887_typical_autism` at 97.0%, `gse324211_xist` at 100% — 36,842/36,842).
- Worst retention: deposits whose source used legacy or atypical symbol versions (eg gse171344 at 36.4%).

**Originally 8 deposits failed migration, all "0 genes mapped". All recovered on 2026-05-15.** Five distinct failure modes were diagnosed and fixed:

| Slug | Failure mode | Example var_name | Fix |
|---|---|---|---|
| `gse183627_kmt2d` | Ensembl with version suffix | `ENSG00000000003.15` | regex extract `ENSG\d{11}` from any decorated string (added to script's `normalize_ensembl()`) |
| `gse295097_radiation` | Ensembl with version suffix | `ENSG00000243485.5` | same |
| `gse301348_h5n1` | Ensembl with version suffix | `ENSG00000000003.15` | same |
| `gse165975` | Ensembl with reference-build prefix | `hg19_ENSG00000223972` | regex extract handles prefix |
| `gse190815` | Ensembl with CellRanger ref prefix | `GRCh38_2020A_v5_ENSG00000243485` | regex extract handles prefix |
| `gse208710` | Ensembl with `GRCh38_` prefix | `GRCh38_ENSG00000243485` | regex extract handles prefix |
| `gse165577` | Numeric placeholders 1..36601 (gene names lost during prior conversion) | `'1','2',...,'36601'` | borrow var ordering from another deposit using same 10x reference (`gse171263`) and rewrite v1 h5ad |
| `gse253230_ube3a` | Transposed deposit (axes swapped); obs inflated with per-sample gene blocks | obs: `gse253230_ube3a__6week__ENSG00000228794.11`, var: 24-char codes | rebuild v1 from raw Alevin output (`data/raw/gse253230/*_quants_mat.{mtx,cols,rows}.gz`); rows=cells, cols=genes; strip version suffix; concat 9 samples |

Two additional bugs surfaced and were fixed in the script during the recovery:
- **HDF5 string write error on mixed-type var/obs columns** (NaN + str). Fixed by `sanitize_for_h5()` coercing all object/category columns to str with NaN→''.
- **Duplicate HGNC names from collapsed Ensembl mappings**: `var_names_make_unique()` was creating `GENE-1`, `GENE-2` suffixes that then failed canonical-var join. Fixed by dropping duplicates *before* uniquifying.

### Validation pass

After migration, every v2 h5ad was checked:
- **104/129 fully OK** — every required column present, gene set ⊂ HNOCA, X integer.
- **25/129 warnings** — X dtype is float32 instead of integer. Underlying values still integer-valued (CSV-loaded deposits store as float); scVI/scPoli will cast at load. Not a hard problem.
- **0 hard failures**.
- 3/3 random concat tests succeeded — confirms schemas are merge-compatible.

### Pointers

- Script: `scripts/migrate_to_hnoca_schema.py` (273 lines, command-line `--slugs` filter, `--dry-run` mode).
- Canonical gene table: `data/reference/hnoca_var_canonical.tsv` (36,842 × 7 cols, 2.2 MB).
- v2 manifest: `data/manifest_v2.tsv` (129 rows).
- v2 h5ads: `data/processed_v2/<slug>.h5ad` (~40 GB total).
- Migration log: `/tmp/migrate_log.json`.

## 6. Quality signals we DO have, and what to check before publishing

**Per-cell QC fields are present** in every h5ad: `n_counts`, `n_genes`, `pct_mito`. You can re-filter per-deposit at AnnData load time.

**Filter heterogeneity is explicit** in `manifest.tsv` (`filter` column). Treat `min500c/200g` deposits and `authors-filtered` deposits as potentially different cell populations; the simple-threshold ones likely retain more ambiguous / borderline-empty cells.

**Sample-level control flags** are reliable for deposits with one cell line per sample. For multiplexed pools (a few % of the atlas, listed in §4.2), the flag is a sample-level summary, not a cell-level truth.

**Gene namespace heterogeneity is real and unfixed.** Some deposits use Ensembl, some use HGNC symbols, some have both. Older deposits (pre-2019) may use deprecated gene symbols. Before integration, run a one-pass namespace harmonization (recommend: use the Braun atlas's `var_names` as the canonical Ensembl set, drop genes not present in ≥3 deposits).

**Donor-level metadata is sparse.** GEO Sample Characteristics frequently omit donor sex, age, ethnicity. Where present, it's in deposit-specific obs columns (`cell_line`, `donor_id`, `Sex`); standardize at integration time if needed.

**Disease state needs an explicit obs column** before publication. Right now disease state is implicit in the `is_control` flag plus deposit-specific columns. Recommend adding `disease_state` (str: `healthy`, `<disease_name>`, `unknown`) as a top-level obs column at integration.

**Mixed pools** (`gse189535`, `gse281452_iMG`, `gse296775_strada`, `gse297594_mecp2`): cell-level demultiplexing requires re-running on HTO/CMO tag matrices. Until then, do not treat their cell-level labels as ground truth.

**`gse282644_hiv` is broken** (12 cells, almost certainly a transpose error in the source CSV). Remove from integration.

## 7. What's compiled, in numbers

| Metric | Value |
|---|---:|
| Compiled deposits | **121** |
| Total cells | **8,690,390** |
| Control cells | 5,686,329 (65.4%) |
| Total samples | 916 |
| Control samples | 629 |
| Single-sample deposits | 17 |
| Multi-sample deposits with all-control | 46 |
| Multi-sample deposits with mixed control/disease | 58 |
| Deposits filtered with min500c/200g | 20 |
| Deposits using author-deposited filter | 98 |
| Cell-count buckets | 8 (<5K), 35 (5K–25K), 50 (25K–100K), 27 (100K–500K), 1 (>500K) |
| Quarantined deposits | ~24 in `_blocked/` + 2 in `_pending_fastq_reprocess/` |
| Reference atlases ready | Braun 2023 (1.67M cells) + HNOCA 2024 (1.77M cells) |

Top organoid types by deposit count: **cerebral** (54), **cortical** (8), **midbrain** (3), **cerebellar** (2), **striatal** (2), **forebrain** (2), **brain** (2), plus deposits tagged `cerebral_hsv`, `cerebral_arid`, `multi_protocol` and the long tail (hypothalamic, pineal, nucleus basalis, ganglionic eminence, hippocampal, thalamocortical, medullary spinal trigeminal, single rosette, multi-region, microglia/cortical assembloid, neural-vascular, etc.). The `organoid_type` column is a coarse first-pass label assigned from paper titles — re-annotate from Methods sections for any deposits that will be central to the analysis.

## 8. Recommended next steps before integration

1. **HNOCA dedup** — compute set intersection between HNOCA's 36 constituent dataset accessions and our 121. Drop or downweight overlaps so we're not double-counting.
2. **Gene namespace harmonization** — pick the Braun var_names as the canonical set, project all 121 deposits onto it, drop ultra-rare genes.
3. **Disease-state column** — add explicit `disease_state` obs column derived from `is_control`, `genotype`, and `diagnosis` per deposit.
4. **Spot-check the medium-confidence control flags** — particularly mixed-pool deposits (§4.2) and the 0-ctrl-deposits we fixed (§4.2 end).
5. **Drop `gse282644_hiv`** (12 cells, likely orientation error).
6. **Resolve `gse195510`** (h5ad reader failure during pass-1) — manual investigation, then re-add or formally block.
7. **Maximize recall before dedup** (per the project lead's instruction): work through the ~10 manual paper-only xlsx rows (find GEO links manually) and the remaining `_blocked/` deposits before deduplication.
8. **Pilot integration** — 3–5 representative deposits + Braun + HNOCA through scVI/scPoli to validate the schema before scaling to all 121.

The atlas is ready for integration work, with the caveats above documented per-deposit in `manifest.tsv` and `data_sources.md`. The `data/manifest.tsv` is the authoritative record; everything else (this note, `data_sources.md`, the `_blocked/REASON.md` files, `triage_rejected.md`) is supporting documentation.
