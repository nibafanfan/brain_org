#!/usr/bin/env python3
"""Rebuild gse297594 control-only from the two demultiplexed Seurat objects.

WHY: the raw GEO MTX for GSE297594 is pooled (pool1-4), and each pool mixes
isogenic-control (WT) and MeCP2-mutant (R133C/R255X) cells that were combined by
cell hashing (HTO). The pooled barcode matrix therefore CANNOT be split by
genotype. The old atlas build loaded the pooled MTX and marked all 164,189 cells
is_control=True -> it silently included thousands of Rett-mutant cells.

The genotype/control call exists ONLY in the Seurat objects' meta.data
(`condition` = WT/MUT, derived from `hash.ID`). This script:
  1. (R) reads each .rds, subsets `condition=="WT"` Singlets, writes a 10x trio
     + a per-cell metadata csv (orig.ident, genotype, library, final_annotation).
  2. (py) assembles cells x genes, runs the SAME apply_qc + map_to_hgnc +
     attach_canonical_var + migrate_obs + sanitize_for_h5 as rebuild_atlas.py,
     overrides bio_sample=orig.ident (per-organoid) / tech_sample=library (pool),
     and writes data/processed/gse297594_mecp2.h5ad.

Run with: /opt/homebrew/Caskroom/miniforge/base/bin/python (anndata 0.12 env).
R is the system Rscript 4.4.1 (Seurat available).
"""
import subprocess
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

ROOT = Path('/Users/eg/brain_organoid')
sys.path.insert(0, str(ROOT / 'scripts'))
from migrate_to_hnoca_schema import (  # noqa: E402
    load_canonical, map_to_hgnc, attach_canonical_var, migrate_obs, sanitize_for_h5,
)
from rebuild_atlas import apply_qc  # noqa: E402

SLUG = 'gse297594_mecp2'
ACCESSION = 'GSE297594'
RAW = ROOT / 'data' / 'raw' / 'gse297594'
WORK = ROOT / 'data' / '_pending_r_interop' / 'gse297594_wt'
OUT = ROOT / 'data' / 'processed' / f'{SLUG}.h5ad'
RSCRIPT = '/usr/local/bin/Rscript'

# pool/library number -> GEO GSM (pool1..4). For provenance only.
LIB_TO_GSM = {'1': 'GSM8995448', '2': 'GSM8995449', '3': 'GSM8995450', '4': 'GSM8995451'}

# source object -> (organoid_type, multi_lineage). CO = cortical organoid (single
# lineage); TA = telencephalic assembloid (cortical + MGE/ventral = multi-lineage).
SOURCE_META = {
    'CO': ('Cortical organoid', 'False'),
    'TA': ('Cortical-MGE telencephalic assembloid', 'True'),
}

t0 = time.time()
def log(m): print(f"[{time.time()-t0:7.1f}s] {m}", flush=True)


R_EXTRACT = r'''
suppressPackageStartupMessages({library(Seurat); library(Matrix)})
args <- commandArgs(trailingOnly=TRUE)
rds <- args[1]; src <- args[2]; outdir <- args[3]
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
obj <- readRDS(rds)
DefaultAssay(obj) <- "RNA"
obj <- tryCatch(JoinLayers(obj), error=function(e) obj)   # Seurat v5 multi-layer -> single
md <- obj@meta.data
keep <- rownames(md)[md$condition == "WT" &
                     (is.null(md$HTO_classification.global) | md$HTO_classification.global == "Singlet")]
cat(src, ": total", ncol(obj), "-> WT singlets", length(keep), "\n")
m <- tryCatch(GetAssayData(obj, assay="RNA", layer="counts"),
              error=function(e) GetAssayData(obj, assay="RNA", slot="counts"))
m <- m[, keep, drop=FALSE]
writeMM(m, file.path(outdir, "matrix.mtx"))
write.table(rownames(m), file.path(outdir, "features.tsv"),
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(colnames(m), file.path(outdir, "barcodes.tsv"),
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
mdk <- md[keep, , drop=FALSE]
out <- data.frame(
  cell           = keep,
  orig_ident     = as.character(mdk$orig.ident),
  genotype       = as.character(mdk$genotype),
  mutation       = as.character(mdk$mutation),
  condition      = as.character(mdk$condition),
  library        = as.character(mdk$library),
  hash_id        = as.character(mdk$hash.ID),
  final_annotation = as.character(mdk$final_annotation),
  stringsAsFactors = FALSE
)
write.csv(out, file.path(outdir, "meta.csv"), row.names=FALSE)
cat("OK\n")
'''


def extract_one(src, rds_name):
    outdir = WORK / src
    rds = RAW / rds_name
    if (outdir / 'matrix.mtx').exists() and (outdir / 'meta.csv').exists():
        log(f"{src}: extract cached at {outdir}")
        return outdir
    script_p = WORK / f'extract_{src}.R'
    WORK.mkdir(parents=True, exist_ok=True)
    script_p.write_text(R_EXTRACT)
    log(f"{src}: running R extraction on {rds_name} (reads ~{rds.stat().st_size/1e9:.1f} GB rds)...")
    r = subprocess.run([RSCRIPT, str(script_p), str(rds), src, str(outdir)],
                       capture_output=True, text=True, timeout=3600)
    print(r.stdout, flush=True)
    if r.returncode != 0 or 'OK' not in r.stdout:
        raise RuntimeError(f"R extraction failed for {src}:\n{r.stderr[-2000:]}")
    return outdir


def load_source(src, outdir):
    """Read the trio + meta -> cells x genes AnnData with raw obs from Seurat."""
    X = sio.mmread(str(outdir / 'matrix.mtx')).T.tocsr().astype('int32')  # genes x cells -> cells x genes
    barcodes = pd.read_csv(outdir / 'barcodes.tsv', header=None)[0].astype(str).tolist()
    features = pd.read_csv(outdir / 'features.tsv', header=None)[0].astype(str).tolist()
    meta = pd.read_csv(outdir / 'meta.csv').set_index('cell')
    meta = meta.reindex(barcodes)  # align to matrix column order
    assert X.shape == (len(barcodes), len(features)), (X.shape, len(barcodes), len(features))
    obs = pd.DataFrame(index=[f'{src}__{b}' for b in barcodes])
    obs['source'] = src
    obs['orig_ident'] = meta['orig_ident'].values
    obs['genotype'] = meta['genotype'].values
    obs['library'] = meta['library'].astype(str).values
    obs['gsm'] = obs['library'].map(LIB_TO_GSM).fillna('GSM_unknown').values
    obs['sample_id'] = meta['orig_ident'].astype(str).str.lower().values
    obs['final_annotation'] = meta['final_annotation'].astype(str).values
    var = pd.DataFrame(index=features)
    log(f"{src}: loaded {X.shape[0]} WT cells x {X.shape[1]} genes")
    return ad.AnnData(X=X, obs=obs, var=var)


def main():
    canonical, ensembl_to_hgnc = load_canonical()
    log(f'canonical: {canonical.shape[0]} genes')

    parts = []
    for src, rds_name in [('CO', 'GSE297594_CO_Seuratobj.rds'),
                          ('TA', 'GSE297594_TA_Seuratobj.rds')]:
        outdir = extract_one(src, rds_name)
        parts.append(load_source(src, outdir))
    a = ad.concat(parts, axis=0, join='outer', fill_value=0, merge='same')
    a.obs_names_make_unique()
    log(f"concatenated WT: {a.shape} | by source: {a.obs['source'].value_counts().to_dict()}")
    log(f"  by genotype: {a.obs['genotype'].value_counts().to_dict()}")

    # --- identical pipeline to rebuild_atlas.process_one ---
    a, qc = apply_qc(a, canonical, SLUG)
    log(f"QC: {qc['pre_qc']} -> {qc['post_qc']} cells ({qc['n_mt']} MT genes)")
    a, mp = map_to_hgnc(a, canonical, ensembl_to_hgnc)
    a = attach_canonical_var(a, canonical)
    log(f"gene map: {mp['in']} -> {mp['kept']} ({mp['ns']})")

    # per-cell organoid_type / multi_lineage by source object
    a.obs['organoid_type'] = a.obs['source'].map(lambda s: SOURCE_META[s][0]).astype(str)
    a.obs['multi_lineage'] = a.obs['source'].map(lambda s: SOURCE_META[s][1]).astype(str)
    a.obs['accession'] = ACCESSION
    a.obs['dataset_slug'] = SLUG
    a.obs['is_control'] = True  # condition=="WT" filter already applied in R
    a.obs['dataset_filter'] = 'seurat_hto_WT_2026-05-23'
    a = migrate_obs(a, SLUG)

    # Override sample hierarchy: bio_sample = per-organoid (orig_ident);
    # tech_sample = the multiplexed pool/library (the real technical batch).
    a.obs['bio_sample'] = a.obs['orig_ident'].astype(str)
    a.obs['tech_sample'] = (SLUG + '_lib' + a.obs['library'].astype(str)).astype(str)
    a.obs['batch'] = a.obs['tech_sample']
    # keep author cell types as the *_original (atlas-wide cell_type stays 'unknown')
    a.obs['cell_type_original'] = a.obs['final_annotation'].astype(str)
    a.obs['cell_type'] = 'unknown'

    a = sanitize_for_h5(a)
    a.write_h5ad(OUT, compression='gzip')
    log(f"WROTE {OUT}  shape={a.shape}")
    log(f"  is_control all True: {bool(a.obs['is_control'].all())}")
    log(f"  bio_sample (organoids): {a.obs['bio_sample'].nunique()} | tech_sample: {a.obs['tech_sample'].nunique()}")
    log(f"  organoid_type: {a.obs['organoid_type'].value_counts().to_dict()}")
    log("DONE")


if __name__ == '__main__':
    main()
