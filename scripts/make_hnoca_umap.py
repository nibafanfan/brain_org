#!/usr/bin/env python3
"""HNOCA-concept UMAP: layout from the scVI latent (expression-derived), colored
(a) by transferred cell-type/region labels and organoid_type (HNOCA-style atlas
panels), and (b) by canonical marker genes (gene-expression overlay on the SAME
coordinates). Demonstrates that layout = expression; coloring is an overlay choice.
"""
import time
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad, scanpy as sc
import scipy.sparse as sp
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

ROOT = Path('/Users/eg/brain_organoid')
N_SUB = 300_000
MARKERS = ['SOX2', 'DCX', 'RBFOX3', 'AQP4', 'P2RY12', 'CLDN5']   # RG, neuroblast, neuron, astro, microglia, endo
t0 = time.time()
def log(m): print(f"[{time.time()-t0:7.1f}s] {m}", flush=True)

lat = ad.read_h5ad(ROOT / 'data/scvi_latent_v5_full.h5ad')          # X = scVI latent
cal = ad.read_h5ad(ROOT / 'data/braun_transfer_full_calibrated.h5ad', backed='r')
for c in ['CellClass_cal', 'Region_pred']:
    lat.obs[c] = cal.obs[c].reindex(lat.obs_names).values
log(f"latent {lat.shape}; labels joined")

rng = np.random.default_rng(0)
idx = np.sort(rng.choice(lat.n_obs, min(N_SUB, lat.n_obs), replace=False))
a = lat[idx].copy()
a.obsm['X_scvi'] = np.asarray(a.X)

# marker expression for the same cells, from the full-gene atlas (markers aren't all HVG)
atlas = ad.read_h5ad(ROOT / 'data/atlas_v5_full.h5ad', backed='r')
pos = pd.Series(np.arange(atlas.n_obs), index=atlas.obs_names)
mk = [m for m in MARKERS if m in set(atlas.var_names)]
ge = atlas[pos.loc[a.obs_names].to_numpy()].to_memory()[:, mk]
X = ge.X.toarray() if sp.issparse(ge.X) else np.asarray(ge.X)
mexpr = ad.AnnData(X=X.astype('float32')); mexpr.var_names = mk
sc.pp.normalize_total(mexpr, target_sum=1e4); sc.pp.log1p(mexpr)
for g in mk:
    a.obs[g] = mexpr[:, g].X.toarray().ravel() if sp.issparse(mexpr.X) else np.asarray(mexpr[:, g].X).ravel()
log(f"markers loaded: {mk}")

sc.pp.neighbors(a, use_rep='X_scvi', n_neighbors=15, random_state=0)
sc.tl.umap(a, random_state=0)
log("UMAP computed")

# (a) HNOCA-style label panels
fig, ax = plt.subplots(1, 3, figsize=(24, 7))
sc.pl.umap(a, color='CellClass_cal', ax=ax[0], show=False, size=3,
           title='CellClass (Braun-transferred)', legend_loc='right margin')
sc.pl.umap(a, color='Region_pred', ax=ax[1], show=False, size=3,
           title='Region (Braun-transferred)', legend_loc='right margin')
sc.pl.umap(a, color='organoid_type', ax=ax[2], show=False, size=3,
           title='organoid_type (protocol)', legend_loc=None)
fig.suptitle('scVI-latent UMAP — colored by transferred labels (HNOCA-style)', fontsize=14)
fig.tight_layout(); fig.savefig(ROOT / 'data/scvi_umap_hnoca_celltype.png', dpi=110)
log("wrote data/scvi_umap_hnoca_celltype.png")

# (b) marker-gene overlay on the SAME layout
n = len(mk); ncol = 3; nrow = int(np.ceil(n / ncol))
fig, ax = plt.subplots(nrow, ncol, figsize=(7 * ncol, 6 * nrow))
ax = np.atleast_1d(ax).ravel()
for i, g in enumerate(mk):
    sc.pl.umap(a, color=g, ax=ax[i], show=False, size=3, cmap='viridis', title=g)
for j in range(n, len(ax)):
    ax[j].axis('off')
fig.suptitle('Same UMAP — colored by marker-gene expression (log1p CP10K)', fontsize=14)
fig.tight_layout(); fig.savefig(ROOT / 'data/scvi_umap_hnoca_markers.png', dpi=110)
log("wrote data/scvi_umap_hnoca_markers.png")
log("DONE")
