#!/usr/bin/env python3
"""Rebuild control-only atlas from raw data per docs/rebuild_plan_2026-05-17.md.

For each INCLUDE row in data/rebuild_config.tsv:
  1. Load via the appropriate loader_type
  2. Filter cells to allowed GSMs (or keep all if INCLUDE_ALL_SAMPLES)
  3. Apply uniform QC (n_counts>=500, n_genes>=200, pct_mito<=20)
  4. Map var to HGNC via HNOCA canonical (36,842 genes)
  5. Attach HNOCA-compatible obs schema
  6. Write data/processed/<slug>.h5ad and append to manifest

Usage:
  scripts/rebuild_atlas.py --dry-run        # validate without writing
  scripts/rebuild_atlas.py --slugs gse242329 gse218457  # subset
  scripts/rebuild_atlas.py                  # full run
"""
import argparse
import gzip
import json
import os
import re
import subprocess
import sys
import tarfile
import time
import traceback
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

ROOT = Path('/Users/eg/brain_organoid')
CONFIG = ROOT / 'data' / 'rebuild_config.tsv'
CANONICAL = ROOT / 'data' / 'reference' / 'hnoca_var_canonical.tsv'
OUT_DIR = ROOT / 'data' / 'processed'
MANIFEST = ROOT / 'data' / 'manifest.tsv'
LOG_PATH = ROOT / 'data' / 'rebuild_log.json'

# Reuse helpers from migrate_to_hnoca_schema
sys.path.insert(0, str(ROOT / 'scripts'))
from migrate_to_hnoca_schema import (  # noqa: E402
    load_canonical, normalize_ensembl, detect_namespace,
    map_to_hgnc, attach_canonical_var, migrate_obs,
    sanitize_for_h5, add_lengthnorm_layer, SMART_SEQ2,
)


# ---------- Loader registry ----------

def _read_trio(bc_path, feat_path, mtx_path, sample_id, gsm):
    """Standard 10x trio reader. Returns AnnData with obs[gsm], obs[sample_id]."""
    X = sio.mmread(str(mtx_path)).T.tocsr().astype('int32')  # transpose to cells x genes
    barcodes = pd.read_csv(bc_path, sep='\t', header=None)[0].astype(str).tolist()
    features = pd.read_csv(feat_path, sep='\t', header=None)
    if features.shape[1] >= 2:
        var_ids = features[0].astype(str).tolist()
        var_syms = features[1].astype(str).tolist()
        var_df = pd.DataFrame({'gene_symbol_orig': var_syms}, index=var_ids)
    else:
        var_ids = features[0].astype(str).tolist()
        var_df = pd.DataFrame(index=var_ids)
    if X.shape[0] != len(barcodes) and X.shape[0] == len(var_ids):
        X = X.T.tocsr()
    if X.shape[1] != len(var_ids):
        raise ValueError(f"shape mismatch: X={X.shape} barcodes={len(barcodes)} features={len(var_ids)}")
    obs = pd.DataFrame({'gsm': gsm, 'sample_id': sample_id}, index=[f'{sample_id}__{b}' for b in barcodes])
    return ad.AnnData(X=X, obs=obs, var=var_df)


def _extract_raw_tar(raw_dir):
    """If a GSE*_RAW.tar is present and per-GSM files aren't yet extracted, extract."""
    for f in os.listdir(raw_dir):
        if f.endswith('_RAW.tar'):
            tar_path = raw_dir / f
            # Check if any GSM* files already exist
            has_gsm = any(p.startswith('GSM') for p in os.listdir(raw_dir))
            if not has_gsm:
                with tarfile.open(tar_path) as t:
                    t.extractall(raw_dir)
            return


def load_10x_mtx_per_gsm(raw_dir, gsm_filter=None):
    """Per-GSM trios (barcodes/features/matrix). One AnnData per GSM, concat at end."""
    raw_dir = Path(raw_dir)
    _extract_raw_tar(raw_dir)
    files = os.listdir(raw_dir)
    # Group files by GSM
    gsm_files = {}
    for f in files:
        m = re.match(r'(GSM\d+)', f)
        if m:
            gsm_files.setdefault(m.group(1), []).append(f)
    parts = []
    for gsm, fs in sorted(gsm_files.items()):
        if gsm_filter and gsm not in gsm_filter:
            continue
        bc = next((f for f in fs if 'barcode' in f.lower() and ('.tsv' in f.lower() or f.endswith('.txt') or f.endswith('.txt.gz'))), None)
        ft = next((f for f in fs if ('feature' in f.lower() or 'gene' in f.lower()) and ('.tsv' in f.lower() or f.endswith('.txt') or f.endswith('.txt.gz')) and 'matrix' not in f.lower()), None)
        mx = (next((f for f in fs if 'matrix' in f.lower() and '.mtx' in f.lower()), None)
              or next((f for f in fs if '.mtx' in f.lower()), None))  # some GSMs name the mtx e.g. GSMxxxx_D07.mtx.gz
        if not (bc and ft and mx):
            print(f"    skip {gsm}: incomplete trio ({fs[:5]})", flush=True)
            continue
        try:
            sample_id = gsm.lower()
            a = _read_trio(raw_dir/bc, raw_dir/ft, raw_dir/mx, sample_id, gsm)
            parts.append(a)
        except Exception as e:
            print(f"    err loading {gsm}: {e}", flush=True)
    if not parts:
        raise RuntimeError("no GSM parts loaded")
    return ad.concat(parts, axis=0, join='outer', fill_value=0)


def _pick_cellranger_h5(fs):
    """Choose ONE .h5 per GSM. Prefer a filtered matrix; never a raw/unfiltered
    one (which contains empty droplets and, if loaded alongside filtered,
    duplicates every real cell)."""
    filt = [f for f in fs if 'filtered' in f.lower()]
    if filt:
        return sorted(filt)[0]
    nonraw = [f for f in fs if 'raw' not in f.lower() and 'unfiltered' not in f.lower()]
    if nonraw:
        return sorted(nonraw)[0]
    return sorted(fs)[0]


def load_cellranger_h5(raw_dir, gsm_filter=None):
    """CellRanger filtered_feature_bc_matrix.h5 per GSM."""
    import scanpy as sc
    raw_dir = Path(raw_dir)
    _extract_raw_tar(raw_dir)
    # Group .h5 by GSM; load only ONE per GSM (deposits often ship both a
    # filtered and a raw/unfiltered/cellbender matrix — loading both doubles cells).
    h5_by_gsm = {}
    for f in sorted(os.listdir(raw_dir)):
        if not f.endswith('.h5'):
            continue
        m = re.match(r'(GSM\d+)', f)
        if not m:
            continue
        h5_by_gsm.setdefault(m.group(1), []).append(f)
    parts = []
    for gsm, fs in sorted(h5_by_gsm.items()):
        if gsm_filter and gsm not in gsm_filter:
            continue
        f = _pick_cellranger_h5(fs)
        if len(fs) > 1:
            print(f"    {gsm}: {len(fs)} .h5 files, using {f}", flush=True)
        try:
            sample_id = gsm.lower()
            a = sc.read_10x_h5(str(raw_dir/f))
            a.obs['gsm'] = gsm
            a.obs['sample_id'] = sample_id
            a.obs_names = [f'{sample_id}__{b}' for b in a.obs_names]
            # Use Ensembl IDs (unique) as var_names; preserve symbols in gene_symbol_orig
            if 'gene_ids' in a.var.columns:
                a.var['gene_symbol_orig'] = a.var_names.astype(str)
                a.var_names = a.var['gene_ids'].astype(str)
            else:
                a.var_names_make_unique()
            parts.append(a)
        except Exception as e:
            print(f"    err loading {gsm}: {e}", flush=True)
    if not parts:
        raise RuntimeError("no h5 parts loaded")
    return ad.concat(parts, axis=0, join='outer', fill_value=0)


def load_csv(raw_dir, gsm_filter=None):
    """One CSV per sample (gene × cell or cell × gene). Auto-detect orientation."""
    raw_dir = Path(raw_dir)
    _extract_raw_tar(raw_dir)
    files = sorted(os.listdir(raw_dir))
    parts = []
    for f in files:
        if not (f.endswith('.csv') or f.endswith('.csv.gz') or f.endswith('.csv.bz2')):
            continue
        m = re.match(r'(GSM\d+)', f)
        if not m:
            continue
        gsm = m.group(1)
        if gsm_filter and gsm not in gsm_filter:
            continue
        try:
            df = pd.read_csv(raw_dir/f, index_col=0)
            # Heuristic: more rows than cols → genes are rows (transpose)
            if df.shape[0] > df.shape[1]:
                X = sp.csr_matrix(df.values.T.astype('int32'))
                cell_ids = df.columns.astype(str).tolist()
                gene_ids = df.index.astype(str).tolist()
            else:
                X = sp.csr_matrix(df.values.astype('int32'))
                cell_ids = df.index.astype(str).tolist()
                gene_ids = df.columns.astype(str).tolist()
            sample_id = gsm.lower()
            obs = pd.DataFrame({'gsm': gsm, 'sample_id': sample_id},
                              index=[f'{sample_id}__{c}' for c in cell_ids])
            var = pd.DataFrame(index=gene_ids)
            a = ad.AnnData(X=X, obs=obs, var=var)
            parts.append(a)
        except Exception as e:
            print(f"    err {f}: {e}", flush=True)
    if not parts:
        raise RuntimeError("no csv parts loaded")
    return ad.concat(parts, axis=0, join='outer', fill_value=0)


def load_dge_text(raw_dir, gsm_filter=None):
    """Drop-seq DGE text. Genes as rows, cells as cols, tab-delimited."""
    raw_dir = Path(raw_dir)
    _extract_raw_tar(raw_dir)
    files = sorted(os.listdir(raw_dir))
    parts = []
    for f in files:
        if not ('.dge.txt' in f.lower() or '.dge.tsv' in f.lower()):
            continue
        m = re.match(r'(GSM\d+)', f)
        if not m:
            continue
        gsm = m.group(1)
        if gsm_filter and gsm not in gsm_filter:
            continue
        try:
            df = pd.read_csv(raw_dir/f, sep='\t', index_col=0)
            X = sp.csr_matrix(df.values.T.astype('int32'))
            sample_id = gsm.lower()
            obs = pd.DataFrame({'gsm': gsm, 'sample_id': sample_id},
                              index=[f'{sample_id}__{c}' for c in df.columns.astype(str)])
            var = pd.DataFrame(index=df.index.astype(str))
            a = ad.AnnData(X=X, obs=obs, var=var)
            parts.append(a)
        except Exception as e:
            print(f"    err {f}: {e}", flush=True)
    if not parts:
        raise RuntimeError("no dge parts loaded")
    return ad.concat(parts, axis=0, join='outer', fill_value=0)


def load_series_level_mtx(raw_dir, gsm_filter=None):
    """Single series-level trio at the raw_dir root. All cells from one sample
    (cell-line/condition has to come from barcode prefixes if needed)."""
    raw_dir = Path(raw_dir)
    files = os.listdir(raw_dir)
    bc = next((f for f in files if 'barcode' in f.lower() and '.tsv' in f.lower() and not f.startswith('GSM')), None)
    ft = next((f for f in files if ('feature' in f.lower() or 'gene' in f.lower()) and '.tsv' in f.lower() and not f.startswith('GSM') and 'matrix' not in f.lower()), None)
    mx = next((f for f in files if 'matrix' in f.lower() and '.mtx' in f.lower() and not f.startswith('GSM')), None)
    if not (bc and ft and mx):
        raise RuntimeError(f"no series-level trio in {raw_dir}: {files[:5]}")
    slug = raw_dir.name
    a = _read_trio(raw_dir/bc, raw_dir/ft, raw_dir/mx, sample_id=slug, gsm=slug.upper())
    # For series-level, treat the slug as the sole GSM unless gsm_filter says otherwise
    if gsm_filter:
        # Series-level data has no per-cell GSM info — can't filter. Keep all.
        pass
    return a


def load_h5ad(raw_dir, gsm_filter=None):
    """Pre-built h5ad. Apply axis-swap detection. Supports .h5ad and .h5ad.gz."""
    raw_dir = Path(raw_dir)
    parts = []
    # Avoid double-reading when both X.h5ad and X.h5ad.gz are present (the .gz is
    # the original; the .h5ad is its decompression). Prefer the .h5ad; only fall
    # back to a .gz whose decompressed sibling is missing.
    all_files = os.listdir(raw_dir)
    h5ad_present = {f for f in all_files if f.endswith('.h5ad')}
    files = sorted(h5ad_present)
    files += sorted(f for f in all_files
                    if f.endswith('.h5ad.gz') and f[:-3] not in h5ad_present)
    if not files:
        raise RuntimeError(f"no .h5ad in {raw_dir}")
    for f in files:
        path = raw_dir / f
        if f.endswith('.h5ad.gz'):
            decompressed = raw_dir / f[:-3]
            if not decompressed.exists():
                with gzip.open(path, 'rb') as gz, open(decompressed, 'wb') as out:
                    out.write(gz.read())
            path = decompressed
        m = re.match(r'(GSM\d+)', f)
        gsm = m.group(1) if m else raw_dir.name.upper()
        if gsm_filter and gsm not in gsm_filter:
            continue
        a = ad.read_h5ad(path)
        # Prefer raw counts: some processed h5ads (e.g. Seurat exports) store
        # normalized values in X and raw integer counts in layers['counts'].
        # Using X directly would feed normalized data into QC/scVI as if counts.
        if 'counts' in a.layers:
            a.X = a.layers['counts'].copy()
        for _ln in list(a.layers.keys()):
            del a.layers[_ln]
        obs_sample = str(a.obs_names[0]) if a.n_obs > 0 else ''
        var_sample = str(a.var_names[0]) if a.n_vars > 0 else ''
        if re.search(r'ENSG\d', obs_sample) and len(var_sample) > 12 and var_sample.replace('-','').replace('_','').isalnum():
            a = a.T
        sample_id = gsm.lower() if m else raw_dir.name
        a.obs['gsm'] = gsm
        a.obs['sample_id'] = sample_id
        a.obs_names = [f'{sample_id}__{b}' for b in a.obs_names]
        parts.append(a)
    if not parts:
        raise RuntimeError(f"no h5ad parts kept after gsm_filter in {raw_dir}")
    if len(parts) == 1:
        return parts[0]
    return ad.concat(parts, axis=0, join='outer', fill_value=0)
    # Axis-swap check: if obs looks like genes (ENSG/short uppercase) and var looks like barcodes
    obs_sample = str(a.obs_names[0])
    var_sample = str(a.var_names[0])
    if re.search(r'ENSG\d', obs_sample) and len(var_sample) > 12 and var_sample.replace('-','').replace('_','').isalnum():
        a = a.T
    if 'gsm' not in a.obs.columns:
        a.obs['gsm'] = raw_dir.name.upper()
    if 'sample_id' not in a.obs.columns:
        a.obs['sample_id'] = raw_dir.name
    if gsm_filter:
        a = a[a.obs['gsm'].astype(str).str.upper().isin(gsm_filter)].copy()
    return a


def load_rds(raw_dir, gsm_filter=None):
    """Seurat RDS via Rscript subprocess. Handles double-gzipped RDS."""
    raw_dir = Path(raw_dir)
    parts = []
    # Avoid double-reading when both X.rds and X.rds.gz are present; prefer the
    # .rds, only fall back to a .gz whose decompressed sibling is missing.
    all_files = os.listdir(raw_dir)
    rds_present = {f for f in all_files if f.endswith('.rds')}
    rds_files = sorted(rds_present)
    rds_files += sorted(f for f in all_files
                        if f.endswith('.rds.gz') and f[:-3] not in rds_present)
    for f in rds_files:
        rds_in = raw_dir / f
        # Decompress double-gz if needed
        if f.endswith('.rds.gz'):
            rds_uncompressed = raw_dir / f[:-3]
            if not rds_uncompressed.exists():
                with gzip.open(rds_in, 'rb') as gz, open(rds_uncompressed, 'wb') as out:
                    out.write(gz.read())
            rds_in = rds_uncompressed
        out_dir = raw_dir / f'rds_extract_{f.replace(".rds.gz","").replace(".rds","")}'
        out_dir.mkdir(exist_ok=True)
        if not (out_dir / 'matrix.mtx').exists():
            script = f"""
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
obj <- readRDS("{rds_in}")
if (inherits(obj, "Seurat")) {{
  m <- GetAssayData(obj, slot="counts")
}} else if (inherits(obj, "dgCMatrix") || inherits(obj, "matrix") || inherits(obj, "dgTMatrix")) {{
  m <- obj
}} else {{ stop(paste("Unknown class:", paste(class(obj), collapse=","))) }}
writeMM(m, "{out_dir}/matrix.mtx")
write.table(rownames(m), "{out_dir}/features.tsv", sep="\\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(colnames(m), "{out_dir}/barcodes.tsv", sep="\\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
cat("OK\\n")
"""
            script_p = f'/tmp/rds_{raw_dir.name}_{f}.R'
            with open(script_p, 'w') as fp: fp.write(script)
            r = subprocess.run(['Rscript', script_p], capture_output=True, text=True, timeout=1800)
            if r.returncode != 0:
                print(f"    R err on {f}: {r.stderr[:300]}", flush=True)
                continue
        # Read back as 10x trio
        # GSM detection from filename
        m = re.match(r'(GSM\d+)', f)
        gsm = m.group(1) if m else raw_dir.name.upper()
        if gsm_filter and gsm not in gsm_filter:
            continue
        sample_id = gsm.lower() if m else f"{raw_dir.name}_{f.replace('.rds.gz','').replace('.rds','')}"
        try:
            a = _read_trio(out_dir/'barcodes.tsv', out_dir/'features.tsv', out_dir/'matrix.mtx', sample_id, gsm)
            parts.append(a)
        except Exception as e:
            print(f"    assembly err for {f}: {e}", flush=True)
    if not parts:
        raise RuntimeError(f"no rds parts loaded from {raw_dir}")
    return ad.concat(parts, axis=0, join='outer', fill_value=0)


def load_tar_archive(raw_dir, gsm_filter=None):
    """Nested tar within the deposit — extract then dispatch.

    Try in order: alevin output (if alevin subdir exists), 10x per-GSM, series-level.
    """
    raw_dir = Path(raw_dir)
    _extract_raw_tar(raw_dir)
    # After extraction, find any nested .tar.gz and extract
    for f in os.listdir(raw_dir):
        if f.endswith('.tar.gz') and not f.endswith('_RAW.tar'):
            with tarfile.open(raw_dir/f) as t:
                t.extractall(raw_dir)
    # Look for alevin structure (subdir with /alevin/ inside or quants_mat.mtx)
    has_alevin = any(
        (raw_dir/sub).is_dir() and
        ((raw_dir/sub/'alevin').exists() or (raw_dir/sub/'quants_mat.mtx').exists() or (raw_dir/sub/'quants_mat.mtx.gz').exists())
        for sub in os.listdir(raw_dir)
    )
    if has_alevin:
        return load_alevin(raw_dir, gsm_filter)
    # Try 10x per-GSM, then series-level
    try:
        return load_10x_mtx_per_gsm(raw_dir, gsm_filter)
    except RuntimeError:
        return load_series_level_mtx(raw_dir, gsm_filter)


def load_alevin(raw_dir, gsm_filter=None):
    """salmon-alevin: quants_mat.mtx + rows.txt + cols.txt per GSM. Cells × genes (no transpose)."""
    raw_dir = Path(raw_dir)
    parts = []
    # Search for alevin output dirs
    for sub in os.listdir(raw_dir):
        sub_p = raw_dir / sub
        if not sub_p.is_dir():
            continue
        alevin_dir = sub_p / 'alevin' if (sub_p / 'alevin').exists() else sub_p
        mtx = alevin_dir / 'quants_mat.mtx'
        if not mtx.exists():
            # Try gzipped variant
            mtx_gz = alevin_dir / 'quants_mat.mtx.gz'
            if mtx_gz.exists():
                with gzip.open(mtx_gz, 'rb') as gz, open(mtx, 'wb') as out:
                    out.write(gz.read())
            else:
                continue
        rows_p = alevin_dir / 'quants_mat_rows.txt'
        cols_p = alevin_dir / 'quants_mat_cols.txt'
        rows = pd.read_csv(rows_p, header=None)[0].astype(str).tolist()
        cols = pd.read_csv(cols_p, header=None)[0].astype(str).tolist()
        X = sio.mmread(str(mtx)).tocsr().astype('int32')
        if X.shape[0] != len(rows):
            X = X.T.tocsr()
        # GSM detection from parent dir
        m = re.search(r'(GSM\d+)', sub)
        gsm = m.group(1) if m else sub.upper()
        if gsm_filter and gsm not in gsm_filter:
            continue
        sample_id = gsm.lower()
        obs = pd.DataFrame({'gsm': gsm, 'sample_id': sample_id},
                          index=[f'{sample_id}__{r}' for r in rows])
        var = pd.DataFrame(index=cols)
        parts.append(ad.AnnData(X=X, obs=obs, var=var))
    if not parts:
        raise RuntimeError(f"no alevin parts in {raw_dir}")
    return ad.concat(parts, axis=0, join='outer', fill_value=0)


def load_streaming_csv(raw_dir, gsm_filter=None):
    """Large dense CSV (>1 GB) — chunked load."""
    raw_dir = Path(raw_dir)
    parts = []
    for f in sorted(os.listdir(raw_dir)):
        if not (f.endswith('.csv') or f.endswith('.csv.gz')):
            continue
        m = re.match(r'(GSM\d+)', f)
        if not m:
            continue
        gsm = m.group(1)
        if gsm_filter and gsm not in gsm_filter:
            continue
        # Stream-load: read in chunks, assemble sparse
        chunks = pd.read_csv(raw_dir/f, chunksize=2000, index_col=0)
        sparse_blocks = []
        col_names = None
        for c in chunks:
            sparse_blocks.append(sp.csr_matrix(c.values.astype('int32')))
            if col_names is None:
                col_names = c.columns.astype(str).tolist()
        # Vstack and transpose (genes were rows)
        X = sp.vstack(sparse_blocks).T.tocsr()
        gene_ids = [g for chunk_genes in [pd.read_csv(raw_dir/f, usecols=[0]).iloc[:,0].astype(str).tolist()] for g in chunk_genes]
        sample_id = gsm.lower()
        obs = pd.DataFrame({'gsm': gsm, 'sample_id': sample_id},
                          index=[f'{sample_id}__{c}' for c in col_names])
        var = pd.DataFrame(index=gene_ids)
        parts.append(ad.AnnData(X=X, obs=obs, var=var))
    if not parts:
        raise RuntimeError("no streaming_csv parts")
    return ad.concat(parts, axis=0, join='outer', fill_value=0)


# Loader dispatch table
LOADERS = {
    '10x_mtx_per_gsm': load_10x_mtx_per_gsm,
    '10x_mtx_positional': load_10x_mtx_per_gsm,  # same — only positional ordering of column names
    'cellranger_h5': load_cellranger_h5,
    'csv': load_csv,
    'csv_bz2': load_csv,  # gzip/bz2 transparent
    'dge_text': load_dge_text,
    'series_level_mtx': load_series_level_mtx,
    'h5ad': load_h5ad,
    'rds': load_rds,
    'tar_archive': load_tar_archive,
    'alevin': load_alevin,
    'streaming_csv': load_streaming_csv,
}


# ---------- QC ----------

def build_mt_mask(adata, canonical):
    """Detect MT genes using HNOCA canonical's mt column (Ensembl-based) or MT- prefix (HGNC-based)."""
    var_names = adata.var_names.astype(str)
    ns = detect_namespace(var_names)
    if ns == 'ensembl':
        mt_ensembls = set(canonical[canonical['mt'] == True]['ensembl'].dropna().astype(str).tolist())
        return np.array([normalize_ensembl(v) in mt_ensembls for v in var_names])
    else:
        return np.array([str(v).upper().startswith('MT-') for v in var_names])


def apply_qc(adata, canonical, slug):
    """Uniform QC: n_counts >= 500, n_genes >= 200, pct_mito <= 20."""
    if not sp.issparse(adata.X):
        adata.X = sp.csr_matrix(adata.X)
    counts = np.asarray(adata.X.sum(axis=1)).ravel()
    ngenes = np.asarray((adata.X > 0).sum(axis=1)).ravel()
    mt_mask = build_mt_mask(adata, canonical)
    n_mt = int(mt_mask.sum())
    if n_mt > 0:
        mt_counts = np.asarray(adata.X[:, mt_mask].sum(axis=1)).ravel()
        pct_mito = mt_counts / np.maximum(1, counts) * 100
    else:
        pct_mito = np.zeros(adata.n_obs, dtype=float)
    keep = (counts >= 500) & (ngenes >= 200) & (pct_mito <= 20)
    adata = adata[keep].copy()
    adata.obs['n_counts'] = counts[keep].astype(np.int32)
    adata.obs['n_genes'] = ngenes[keep].astype(np.int32)
    adata.obs['pct_mito'] = pct_mito[keep].astype(np.float32)
    return adata, {'n_mt': n_mt, 'pre_qc': len(counts), 'post_qc': adata.n_obs}


# ---------- Main pipeline ----------

def process_one(row, canonical, ensembl_to_hgnc, dry_run=False):
    slug = row['slug']
    if not slug or pd.isna(slug):
        return {'slug': None, 'accession': row.get('accession',''), 'status': 'skip_no_slug'}
    raw_dir = row.get('raw_dir', '')
    if not raw_dir or pd.isna(raw_dir):
        return {'slug': slug, 'status': 'skip_no_raw_dir'}
    raw_dir = ROOT / raw_dir
    if not raw_dir.exists():
        return {'slug': slug, 'status': 'skip_missing_raw_dir', 'raw_dir': str(raw_dir)}
    loader_type = row['loader_type']
    if loader_type in ('missing', 'custom'):
        return {'slug': slug, 'status': 'skip_unsupported_loader', 'loader': loader_type}
    if loader_type not in LOADERS:
        return {'slug': slug, 'status': 'skip_unknown_loader', 'loader': loader_type}

    # GSM filter
    gsm_filter = None
    if row['decision'] == 'INCLUDE_EXPLICIT':
        if isinstance(row.get('gsms'), str) and row['gsms']:
            gsm_filter = set(g.strip().upper() for g in row['gsms'].split(',') if g.strip())
    # INCLUDE_ALL_SAMPLES → no filter

    t0 = time.time()
    a = LOADERS[loader_type](raw_dir, gsm_filter)
    if a.n_obs == 0:
        return {'slug': slug, 'status': 'fail_zero_cells_post_load'}

    # Apply GSM filter at obs level (loader may have loaded everything)
    if gsm_filter:
        keep_mask = a.obs['gsm'].astype(str).str.upper().isin(gsm_filter).values
        a = a[keep_mask].copy()
        if a.n_obs == 0:
            return {'slug': slug, 'status': 'fail_zero_cells_post_gsm_filter'}

    n_loaded = a.n_obs
    # QC
    a, qc_stats = apply_qc(a, canonical, slug)
    if a.n_obs == 0:
        return {'slug': slug, 'status': 'fail_zero_cells_post_qc', **qc_stats}

    # Map to HGNC + attach canonical var
    a, map_stats = map_to_hgnc(a, canonical, ensembl_to_hgnc)
    if a is None:
        return {'slug': slug, 'status': 'fail_zero_genes_post_map', **qc_stats, **map_stats}
    a = attach_canonical_var(a, canonical)

    # Migrate obs to HNOCA schema
    a.obs['accession'] = row.get('accession', slug.upper())
    a.obs['dataset_slug'] = slug
    a.obs['organoid_type'] = str(row.get('organoid_type', '') or '')
    a.obs['multi_lineage'] = str(row.get('multi_lineage', '') or '')
    a.obs['is_control'] = True  # the filter already excluded non-controls
    a.obs['dataset_filter'] = f"rebuild_{row['decision']}_2026-05-17"
    a = migrate_obs(a, slug)

    # Lengthnorm layer for Smart-seq2
    if row.get('is_smartseq2', False) or slug in SMART_SEQ2:
        a = add_lengthnorm_layer(a)

    # Sanitize
    a = sanitize_for_h5(a)

    # Write
    if not dry_run:
        out_path = OUT_DIR / f'{slug}.h5ad'
        a.write_h5ad(out_path, compression='gzip')
        result_path = str(out_path)
    else:
        result_path = '(dry-run, not written)'

    return {
        'slug': slug,
        'status': 'ok',
        'loader': loader_type,
        'gsm_filter_size': len(gsm_filter) if gsm_filter else None,
        'cells_loaded': n_loaded,
        'cells_post_qc': qc_stats['post_qc'],
        'n_mt_genes': qc_stats['n_mt'],
        'genes_kept': map_stats['kept'],
        'genes_in': map_stats['in'],
        'final_shape': list(a.shape),
        'has_lengthnorm': 'counts_lengthnorm' in a.layers,
        'path': result_path,
        'seconds': round(time.time() - t0, 1),
    }


# ---------- Main driver ----------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--slugs', nargs='*', help='Limit to specific slugs (for dry-run)')
    ap.add_argument('--limit', type=int, default=0, help='Process at most N rows')
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    canonical, ensembl_to_hgnc = load_canonical()
    print(f'canonical: {canonical.shape[0]} genes; {len(ensembl_to_hgnc)} ENSG→HGNC', flush=True)

    config = pd.read_csv(CONFIG, sep='\t', dtype=str).fillna('')
    include_mask = config['decision'].str.startswith('INCLUDE')
    config = config[include_mask].copy()
    print(f'INCLUDE rows in config: {len(config)}', flush=True)

    if args.slugs:
        config = config[config['slug'].isin(args.slugs)]
        print(f'limited to slugs: {len(config)} rows', flush=True)
    if args.limit:
        config = config.head(args.limit)
        print(f'limited to first {args.limit} rows', flush=True)

    # Mark boolean
    config['is_smartseq2'] = config['is_smartseq2'].astype(str).str.lower() == 'true'

    results = []
    t_start = time.time()
    for i, row in config.iterrows():
        idx = list(config.index).index(i) + 1
        slug = row['slug'] or '(blank)'
        print(f'\n[{idx}/{len(config)}] {slug} ({row["decision"]}, loader={row["loader_type"]})', flush=True)
        try:
            r = process_one(row, canonical, ensembl_to_hgnc, dry_run=args.dry_run)
            print(f'  result: {r.get("status")} {r.get("final_shape", "")} {r.get("seconds", "")}s', flush=True)
        except Exception as e:
            r = {'slug': slug, 'status': 'err', 'error': str(e), 'tb': traceback.format_exc()}
            print(f'  ERR: {e}', flush=True)
            traceback.print_exc()
        results.append(r)

    # Build manifest from successful results. MERGE into any existing manifest so a
    # targeted --slugs run only updates the slugs it touched (never clobbers others).
    if not args.dry_run:
        cols = ['slug','accession','path','n_cells','n_genes','n_samples',
                'n_control_samples','n_control_cells','organoid_type','filter','status']
        existing = {}
        if MANIFEST.exists():
            prev = pd.read_csv(MANIFEST, sep='\t', dtype=str, keep_default_na=False)
            for _, pr in prev.iterrows():
                existing[pr['slug']] = {c: pr.get(c, '') for c in cols}
        for r in results:
            if r.get('status') == 'ok':
                existing[r['slug']] = {
                    'slug': r['slug'],
                    'accession': config[config['slug']==r['slug']].iloc[0].get('accession',''),
                    'path': r['path'],
                    'n_cells': r['final_shape'][0],
                    'n_genes': r['final_shape'][1],
                    'n_samples': '',
                    'n_control_samples': '',
                    'n_control_cells': r['final_shape'][0],
                    'organoid_type': config[config['slug']==r['slug']].iloc[0].get('organoid_type',''),
                    'filter': r.get('loader', ''),
                    'status': 'ok',
                }
        pd.DataFrame([existing[s] for s in sorted(existing)], columns=cols).to_csv(MANIFEST, sep='\t', index=False)
        print(f'\nWROTE manifest: {MANIFEST} ({len(existing)} rows, merged)', flush=True)

    # Log
    with open(LOG_PATH, 'w') as f:
        json.dump([{k: (list(v) if isinstance(v, tuple) else v) for k, v in r.items() if not isinstance(v, type)} for r in results], f, indent=2, default=str)

    n_ok = sum(1 for r in results if r.get('status') == 'ok')
    n_skip = sum(1 for r in results if r.get('status', '').startswith('skip'))
    n_fail = sum(1 for r in results if r.get('status', '') not in ('ok',) and not r.get('status', '').startswith('skip'))
    elapsed = time.time() - t_start
    print(f'\n===== SUMMARY =====', flush=True)
    print(f'Total: {len(results)}  ok: {n_ok}  skipped: {n_skip}  failed: {n_fail}', flush=True)
    print(f'Elapsed: {elapsed:.0f}s ({elapsed/60:.1f}m)', flush=True)
    print(f'Log: {LOG_PATH}', flush=True)


if __name__ == '__main__':
    main()
