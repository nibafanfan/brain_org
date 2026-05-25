#!/usr/bin/env python3
"""scIB with a BASELINE embedding (review #2 refinement): compare unintegrated
lognorm-PCA vs the scVI latent in the same Benchmarker run, so the scaled
aggregate scores are meaningful (single-embedding aggregates aren't).

RUN WITH THE ISOLATED VENV:
  /Users/eg/.venvs/scib/bin/python scripts/scib_metrics_baseline.py
"""
import time
import numpy as np, anndata as ad

ROOT = '/Users/eg/brain_organoid'
N_SUB = 300_000
SEED = 0
np.random.seed(SEED)
t0 = time.time()
def log(m): print(f"[{time.time()-t0:7.1f}s] {m}", flush=True)

lat = ad.read_h5ad(f'{ROOT}/data/scvi_latent_v5_full.h5ad')          # X = scVI latent, ordered like preprocessed
pre = ad.read_h5ad(f'{ROOT}/data/processed/atlas_v5_preprocessed.h5ad', backed='r')
assert list(lat.obs_names[:3]) == list(pre.obs_names[:3]) and lat.n_obs == pre.n_obs, \
    "scvi latent and preprocessed must share row order (positional join for PCA)"
cal = ad.read_h5ad(f'{ROOT}/data/braun_transfer_full_calibrated.h5ad', backed='r')

rng = np.random.default_rng(SEED)
idx = np.sort(rng.choice(lat.n_obs, N_SUB, replace=False))
sub = lat[idx].copy()
sub.obs['CellClass_cal'] = cal.obs['CellClass_cal'].reindex(sub.obs_names).values
keep = sub.obs['CellClass_cal'].notna() & (sub.obs['CellClass_cal'] != 'Unknown')
sub = sub[keep.values].copy()
kept = idx[keep.values]

# embeddings: scVI latent + unintegrated lognorm-PCA(50) on the same cells
sub.obsm['X_scvi'] = np.asarray(sub.X)
import scipy.sparse as sp
from sklearn.decomposition import PCA
Xp = pre[kept].to_memory().X
Xp = Xp.toarray() if sp.issparse(Xp) else np.asarray(Xp)
sub.obsm['X_pca'] = PCA(n_components=50, random_state=SEED).fit_transform(Xp)
sub.obs['dataset_slug'] = sub.obs['dataset_slug'].astype('category')
sub.obs['CellClass_cal'] = sub.obs['CellClass_cal'].astype('category')
log(f"benchmarking {sub.n_obs:,} cells; embeddings: X_pca (unintegrated), X_scvi")

from scib_metrics.benchmark import Benchmarker
bm = Benchmarker(sub, batch_key='dataset_slug', label_key='CellClass_cal',
                 embedding_obsm_keys=['X_pca', 'X_scvi'], n_jobs=-1)
bm.benchmark()
res = bm.get_results(min_max_scale=False)
log("=== scIB: unintegrated PCA vs scVI latent (raw) ===")
log("\n" + res.T.to_string())
res.to_csv(f'{ROOT}/data/scib_metrics_baseline.tsv', sep='\t')
log("saved -> data/scib_metrics_baseline.tsv")
log("DONE")
