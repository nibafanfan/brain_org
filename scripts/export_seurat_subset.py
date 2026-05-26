#!/usr/bin/env python3
"""Export a random 100k-cell subset of atlas_v5_full as 10x MTX for Seurat.

Run with the miniforge BASE env python:
  /opt/homebrew/Caskroom/miniforge/base/bin/python scripts/export_seurat_subset.py

Produces (in data/seurat_subset_100k/):
  matrix.mtx.gz       genes x cells raw counts (CellRanger v3 orientation)
  features.tsv.gz     ensembl  symbol  "Gene Expression"
  barcodes.tsv.gz     one barcode per line
  scvi_latent.csv     100k x 30 integrated embedding (reuse, no re-integration)
  metadata.csv        per-cell annotations for coloring clusters

Cells are sampled only from those present in BOTH the full atlas and the scVI
latent (zero-count cells were dropped before scVI), so every exported cell has a
latent vector. Order is identical across all 5 files.
"""
import gzip
import numpy as np
import scipy.io as sio
import scipy.sparse as sp
import anndata as ad
import h5py

SEED = 0
N = 100_000
ATLAS = "data/atlas_v5_full.h5ad"
LATENT = "data/scvi_latent_v5_full.h5ad"
OUTDIR = "data/seurat_subset_100k"
META_COLS = [
    "cell_type", "cell_type_origin", "gsm", "accession", "protocol",
    "organoid_type", "multi_lineage", "age_days", "cell_line", "n_genes",
    "n_counts", "pct_mito", "tech_sample", "bio_sample",
]

import os
os.makedirs(OUTDIR, exist_ok=True)

print("[1/6] reading barcodes from both files...")
with h5py.File(ATLAS, "r") as f:
    atlas_bc = f["obs"][f["obs"].attrs["_index"]][:].astype(str)
with h5py.File(LATENT, "r") as f:
    lat_bc = f["obs"][f["obs"].attrs["_index"]][:].astype(str)

lat_set = set(lat_bc.tolist())
eligible = np.flatnonzero(np.isin(atlas_bc, lat_bc))
print(f"    atlas={len(atlas_bc):,}  latent={len(lat_bc):,}  eligible(both)={len(eligible):,}")

rng = np.random.default_rng(SEED)
sel = np.sort(rng.choice(eligible, size=N, replace=False))
sel_bc = atlas_bc[sel]
print(f"    sampled {len(sel):,} cells (seed={SEED})")

print("[2/6] extracting raw counts (backed)...")
adata = ad.read_h5ad(ATLAS, backed="r")
sub = adata[sel].to_memory()
X = sub.X.tocsr()
X.eliminate_zeros()
print(f"    counts shape {X.shape}  nnz={X.nnz:,}  dtype={X.dtype}  max={X.max()}")
assert sub.obs_names.tolist() == sel_bc.tolist()

print("[3/6] writing matrix.mtx.gz (genes x cells)...")
genes_x_cells = X.T.tocsr()  # CellRanger MTX = features x barcodes
with gzip.open(f"{OUTDIR}/matrix.mtx.gz", "wb") as fh:
    sio.mmwrite(fh, genes_x_cells.astype(np.int32), field="integer")

print("[4/6] writing features.tsv.gz & barcodes.tsv.gz...")
ensembl = sub.var["ensembl"].astype(str).values
symbol = sub.var_names.astype(str).values
with gzip.open(f"{OUTDIR}/features.tsv.gz", "wt") as fh:
    for e, s in zip(ensembl, symbol):
        e = e if e and e != "nan" else s  # fall back to symbol as id if no ensembl
        fh.write(f"{e}\t{s}\tGene Expression\n")
with gzip.open(f"{OUTDIR}/barcodes.tsv.gz", "wt") as fh:
    for b in sel_bc:
        fh.write(f"{b}\n")

print("[5/6] writing scvi_latent.csv (aligned)...")
lat = ad.read_h5ad(LATENT, backed="r")
lat_pos = {b: i for i, b in enumerate(lat_bc)}
rows = np.array([lat_pos[b] for b in sel_bc])
Z = lat[rows].to_memory().X
Z = np.asarray(Z)
import csv
with open(f"{OUTDIR}/scvi_latent.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["barcode"] + [f"scvi_{i+1}" for i in range(Z.shape[1])])
    for b, row in zip(sel_bc, Z):
        w.writerow([b] + [f"{v:.6g}" for v in row])
print(f"    latent {Z.shape}")

print("[6/6] writing metadata.csv...")
cols = [c for c in META_COLS if c in sub.obs.columns]
meta = sub.obs[cols].copy()
meta.insert(0, "barcode", sel_bc)
meta.to_csv(f"{OUTDIR}/metadata.csv", index=False)
print(f"    metadata cols: {cols}")

print(f"\nDONE -> {OUTDIR}/")
print("Files:", os.listdir(OUTDIR))
