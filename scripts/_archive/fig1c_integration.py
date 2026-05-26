#!/usr/bin/env python3
"""HNOCA Fig.1c-style integration view: the integrated UMAP colored by dataset and
protocol. Reuses the persisted umap_points coordinates (no UMAP recompute) and joins
dataset_slug from the latent obs.

HONEST FRAMING: this is an integration DIAGNOSTIC. Our scIB iLISI~0.015 shows datasets
mix only partially (residual batch structure), so expect dataset-structured regions
rather than HNOCA's fully-blended look.

  python scripts/fig1c_integration.py [--out-tag v5]
"""
import argparse, sys, time
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad, scanpy as sc
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parent))
from atlas_common import load_config

t0 = time.time()
def log(m): print(f"[{time.time()-t0:7.1f}s] {m}", flush=True)
ap = argparse.ArgumentParser()
ap.add_argument('--root', default=None)
ap.add_argument('--out-tag', default='v5')
ap.add_argument('--legend-top-n', type=int, default=15)
args = ap.parse_args()
cfg = load_config(root=args.root)
OUT = cfg.root / 'data/poster'
pts_path = OUT / f'umap_points_{args.out_tag}.tsv.gz'
assert pts_path.exists(), f"{pts_path} missing — run poster_panels.py first"

pts = pd.read_csv(pts_path, sep='\t', index_col=0)
lat = ad.read_h5ad(cfg.latent, backed='r')
ds = lat.obs['dataset_slug'].reindex(pts.index).astype(str)
log(f"{len(pts):,} cells; {ds.nunique()} datasets; protocols {pts['protocol'].nunique() if 'protocol' in pts else 'NA'}")

def topn(s, n, other='Other'):
    keep = s.value_counts().index[:n]
    return s.where(s.isin(keep), other).astype('category')

a = ad.AnnData(np.zeros((len(pts), 1), 'float32'))
a.obs_names = pts.index
a.obsm['X_umap'] = pts[['UMAP1', 'UMAP2']].to_numpy()
a.obs['dataset_slug'] = pd.Categorical(ds.values)          # 113 -> no legend
a.obs['protocol_topN'] = topn(pts['protocol'].astype(str), args.legend_top_n) if 'protocol' in pts else 'NA'

fig, ax = plt.subplots(1, 2, figsize=(26, 11))
sc.pl.umap(a, color='dataset_slug', ax=ax[0], show=False, size=6, frameon=False,
           legend_loc=None, title=f'by dataset ({ds.nunique()}) — integration diagnostic (want mixed)')
sc.pl.umap(a, color='protocol_topN', ax=ax[1], show=False, size=6, frameon=False,
           legend_loc='right margin', title=f'by protocol (top {args.legend_top_n})')
fig.suptitle('Integration view (scVI latent UMAP) — colored by batch/protocol', fontsize=16, y=1.0)
fig.tight_layout(rect=[0, 0, 1, 0.97])
for ext in ('png', 'pdf'):
    fig.savefig(OUT / f'fig1c_integration_{args.out_tag}.{ext}', dpi=150)
log(f"wrote fig1c_integration_{args.out_tag}.png/pdf -> {OUT}")
log("DONE")
