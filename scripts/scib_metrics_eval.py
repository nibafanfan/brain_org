#!/usr/bin/env python3
"""Review item #2: proper scIB metrics on the scVI latent (replaces the crude
same-neighbor ratio). Uses the CALIBRATED labels (CellClass_cal).

RUN WITH THE ISOLATED VENV (keeps scib-metrics/jax out of the scvi base env):
  /Users/eg/.venvs/scib/bin/python scripts/scib_metrics_eval.py

Batch correction (batch=dataset_slug): iLISI, kBET, graph connectivity, silhouette_batch, PCR.
Bio conservation  (label=CellClass_cal): cLISI, silhouette_label, isolated labels, NMI/ARI.
Outputs the full Benchmarker table + a stratified iLISI-per-CellClass supplement.
"""
import time
import numpy as np, pandas as pd, anndata as ad

ROOT = '/Users/eg/brain_organoid'
N_SUB = 300_000
t0 = time.time()
def log(m): print(f"[{time.time()-t0:7.1f}s] {m}", flush=True)

lat = ad.read_h5ad(f'{ROOT}/data/scvi_latent_v5_full.h5ad')
cal = ad.read_h5ad(f'{ROOT}/data/braun_transfer_full_calibrated.h5ad', backed='r')
lat.obs['CellClass_cal'] = cal.obs['CellClass_cal'].reindex(lat.obs_names).values
# drop abstained/unlabeled for bio metrics; need a clean categorical label
lat = lat[lat.obs['CellClass_cal'].notna() & (lat.obs['CellClass_cal'] != 'Unknown')].copy()
log(f"labeled latent: {lat.n_obs:,}")

rng = np.random.default_rng(0)
if lat.n_obs > N_SUB:
    lat = lat[np.sort(rng.choice(lat.n_obs, N_SUB, replace=False))].copy()
lat.obsm['X_scvi'] = np.asarray(lat.X)
lat.obs['dataset_slug'] = lat.obs['dataset_slug'].astype('category')
lat.obs['CellClass_cal'] = lat.obs['CellClass_cal'].astype('category')
log(f"subsample {lat.shape}; {lat.obs['dataset_slug'].nunique()} datasets, "
    f"{lat.obs['CellClass_cal'].nunique()} classes")

from scib_metrics.benchmark import Benchmarker
bm = Benchmarker(lat, batch_key='dataset_slug', label_key='CellClass_cal',
                 embedding_obsm_keys=['X_scvi'], n_jobs=-1)
bm.benchmark()
res = bm.get_results(min_max_scale=False)
log("=== scIB metrics (scVI latent, calibrated labels) ===")
log("\n" + res.T.to_string())
res.to_csv(f'{ROOT}/data/scib_metrics.tsv', sep='\t')
log("saved -> data/scib_metrics.tsv")
log("DONE")
