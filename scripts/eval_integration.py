#!/usr/bin/env python3
"""Validate the scVI embedding: did it mix batches while preserving biology?

Run with miniforge base python. Subsamples for tractable UMAP, builds a kNN graph
on X_scvi, and reports:
  - batch mixing: mean fraction of k-NN sharing the same dataset_slug (lower=better)
    vs the no-integration baseline expectation.
  - biology: mean fraction of k-NN sharing the same organoid_type (higher=preserved).
  - UMAP PNG colored by dataset_slug and organoid_type.
"""
import time
import numpy as np
import anndata as ad
import scanpy as sc
from brain_organoid.config import add_common_overrides, load_config, merge_cli_overrides, resolve_path, require_file, stamp_provenance

t0 = time.time()
def log(m): print(f"[{time.time()-t0:7.1f}s] {m}", flush=True)

ap = add_common_overrides(__import__('argparse').ArgumentParser())
ap.add_argument('--subsample', type=int, default=None)
ap.add_argument('--k', type=int, default=None)
args = ap.parse_args()
cfg = merge_cli_overrides(load_config(args.config), args)
lat_rel = args.in_override or cfg['paths']['scvi_latent_full']
LAT = require_file(resolve_path(cfg, lat_rel), 'scvi latent input')
N_SUB = args.subsample or cfg['defaults']['integration_eval']['subsample']
K = args.k or cfg['defaults']['integration_eval']['knn_k']

log(f"reading {LAT}")
a = ad.read_h5ad(LAT)
log(f"loaded {a.shape}")
rng = np.random.default_rng(0)
idx = np.sort(rng.choice(a.n_obs, min(N_SUB, a.n_obs), replace=False))
a = a[idx].copy()
a.obsm['X_scvi'] = a.X
log(f"subsampled to {a.shape}")

sc.pp.neighbors(a, n_neighbors=K, use_rep='X_scvi')
log("neighbors built")

# kNN-based mixing/purity from the connectivity graph
import scipy.sparse as sp
conn = a.obsp['distances']  # kNN (excludes self)
conn = (conn > 0).astype(int)
def same_frac(labels):
    codes = a.obs[labels].astype('category').cat.codes.to_numpy()
    same = 0; tot = 0
    indptr, indices = conn.indptr, conn.indices
    for i in range(conn.shape[0]):
        nb = indices[indptr[i]:indptr[i+1]]
        if len(nb)==0: continue
        same += int((codes[nb]==codes[i]).sum()); tot += len(nb)
    return same/tot

ds_same = same_frac('dataset_slug')
ot_same = same_frac('organoid_type')
# baselines: expected same-fraction if neighbors were random (sum p_i^2)
def baseline(labels):
    p = a.obs[labels].value_counts(normalize=True).to_numpy()
    return float((p**2).sum())
log("=== integration metrics (200k subsample, k=30) ===")
log(f"  batch (dataset_slug): same-neighbor frac = {ds_same:.3f} | random baseline = {baseline('dataset_slug'):.3f}")
log(f"     -> closer to baseline = better mixing; >>baseline = batch still structures the space")
log(f"  biology (organoid_type): same-neighbor frac = {ot_same:.3f} | random baseline = {baseline('organoid_type'):.3f}")
log(f"     -> well above baseline = biology preserved")

log("computing UMAP")
sc.tl.umap(a)
try:
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(16, 7))
    sc.pl.umap(a, color='dataset_slug', ax=ax[0], show=False, legend_loc=None,
               title=f'by dataset (113) — want MIXED  [same-nbr={ds_same:.2f}]', size=3)
    sc.pl.umap(a, color='organoid_type', ax=ax[1], show=False,
               title=f'by organoid_type (33) — want SEPARATED  [same-nbr={ot_same:.2f}]', size=3)
    fig.tight_layout()
    out_png = resolve_path(cfg, cfg['outputs']['figures_dir']) / (
        f"scvi_umap_eval_{cfg.get('out_tag', 'v5')}.png"
    )
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=110)
    log(f"wrote {out_png}")
    stamp_provenance(
        out_png.with_suffix('.provenance.json'),
        cfg,
        {'input': str(LAT), 'figure': str(out_png), 'k': K, 'subsample': N_SUB},
    )
except Exception as e:
    log(f"(plot skipped: {e})")
log("DONE")
