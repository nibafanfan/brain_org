#!/usr/bin/env python3
"""HNOCA-concept UMAP: layout from the scVI latent (expression-derived), colored
(a) by transferred cell-type/region labels and organoid_type (HNOCA-style atlas
panels), and (b) by canonical marker genes (gene-expression overlay on the SAME
coordinates). Demonstrates that layout = expression; coloring is an overlay choice.

Paths are CLI args (defaults below) for portability; outputs accept --out-tag and
get a provenance sidecar JSON.
"""
import argparse, sys, time
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad, scanpy as sc
import scipy.sparse as sp
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _provenance import write_sidecar

MARKERS = ['SOX2', 'DCX', 'RBFOX3', 'AQP4', 'P2RY12', 'CLDN5']   # RG, neuroblast, neuron, astro, microglia, endo
t0 = time.time()
def log(m): print(f"[{time.time()-t0:7.1f}s] {m}", flush=True)

ap = argparse.ArgumentParser()
ap.add_argument('--root', default='/Users/eg/brain_organoid')
ap.add_argument('--in-latent', default=None, help='scVI latent h5ad (default <root>/data/scvi_latent_v5_full.h5ad)')
ap.add_argument('--in-calibrated', default=None, help='calibrated transfer h5ad (labels)')
ap.add_argument('--in-atlas', default=None, help='full-gene atlas h5ad (marker expression)')
ap.add_argument('--out-dir', default=None, help='output dir (default <root>/data)')
ap.add_argument('--out-tag', default='', help='suffix for output filenames, e.g. v5')
ap.add_argument('--n-sub', type=int, default=300_000)
ap.add_argument('--legend-top-n', type=int, default=15, help='top-N organoid_type categories in legend (+Other)')
args = ap.parse_args()

ROOT = Path(args.root)
IN_LAT = Path(args.in_latent or ROOT / 'data/scvi_latent_v5_full.h5ad')
IN_CAL = Path(args.in_calibrated or ROOT / 'data/braun_transfer_full_calibrated.h5ad')
IN_ATLAS = Path(args.in_atlas or ROOT / 'data/atlas_v5_full.h5ad')
OUT = Path(args.out_dir or ROOT / 'data')
tag = (f"_{args.out_tag}" if args.out_tag else "")
prov = {'in_latent': str(IN_LAT), 'in_calibrated': str(IN_CAL), 'in_atlas': str(IN_ATLAS),
        'n_sub': args.n_sub, 'legend_top_n': args.legend_top_n, 'markers': MARKERS}

def top_n(series, n, other='Other'):
    keep = series.value_counts().index[:n]
    return series.where(series.isin(keep), other).astype('category')

lat = ad.read_h5ad(IN_LAT)
cal = ad.read_h5ad(IN_CAL, backed='r')
for c in ['CellClass_cal', 'Region_pred']:
    lat.obs[c] = cal.obs[c].reindex(lat.obs_names).values
log(f"latent {lat.shape}; labels joined")

rng = np.random.default_rng(0)
idx = np.sort(rng.choice(lat.n_obs, min(args.n_sub, lat.n_obs), replace=False))
a = lat[idx].copy()
a.obsm['X_scvi'] = np.asarray(a.X)

# marker expression for the same cells, from the full-gene atlas (markers aren't all HVG)
atlas = ad.read_h5ad(IN_ATLAS, backed='r')
missing = a.obs_names.difference(atlas.obs_names)
assert len(missing) == 0, f"{len(missing)} subsample cells absent from atlas — cannot align markers"
pos = pd.Series(np.arange(atlas.n_obs), index=atlas.obs_names)
mk = [m for m in MARKERS if m in set(atlas.var_names)]
ge = atlas[pos.loc[a.obs_names].to_numpy()].to_memory()[:, mk]
X = ge.X.toarray() if sp.issparse(ge.X) else np.asarray(ge.X)
mexpr = ad.AnnData(X=X.astype('float32')); mexpr.var_names = mk
sc.pp.normalize_total(mexpr, target_sum=1e4); sc.pp.log1p(mexpr)
for g in mk:
    a.obs[g] = mexpr[:, g].X.toarray().ravel() if sp.issparse(mexpr.X) else np.asarray(mexpr[:, g].X).ravel()
a.obs['organoid_type_topN'] = top_n(a.obs['organoid_type'].astype(str), args.legend_top_n)
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
sc.pl.umap(a, color='organoid_type_topN', ax=ax[2], show=False, size=3,
           title=f'organoid_type (top {args.legend_top_n}+Other)', legend_loc='right margin')
fig.suptitle('scVI-latent UMAP — colored by transferred labels (HNOCA-style)', fontsize=14)
fig.tight_layout()
p1 = OUT / f'scvi_umap_hnoca_celltype{tag}.png'
fig.savefig(p1, dpi=110); write_sidecar(p1, __file__, prov)
log(f"wrote {p1}")

# (b) marker-gene overlay on the SAME layout
n = len(mk); ncol = 3; nrow = int(np.ceil(n / ncol))
fig, ax = plt.subplots(nrow, ncol, figsize=(7 * ncol, 6 * nrow))
ax = np.atleast_1d(ax).ravel()
for i, g in enumerate(mk):
    sc.pl.umap(a, color=g, ax=ax[i], show=False, size=3, cmap='viridis', title=g)
for j in range(n, len(ax)):
    ax[j].axis('off')
fig.suptitle('Same UMAP — colored by marker-gene expression (log1p CP10K)', fontsize=14)
fig.tight_layout()
p2 = OUT / f'scvi_umap_hnoca_markers{tag}.png'
fig.savefig(p2, dpi=110); write_sidecar(p2, __file__, prov)
log(f"wrote {p2}")
log("DONE")
