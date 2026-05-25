#!/usr/bin/env python3
"""scIB comparator: unintegrated PCA vs scVI vs scANVI (label-aware), same panel.
Decides the GO/NO-GO embedding-level thresholds (ΔiLISI/ΔkBET vs scVI, cLISI drop).

RUN WITH THE ISOLATED VENV:
  /Users/eg/.venvs/scib/bin/python scripts/scib_metrics_comparator.py
"""
import time
import numpy as np, anndata as ad
import scipy.sparse as sp
from sklearn.decomposition import PCA

ROOT = '/Users/eg/brain_organoid'
N_SUB = 300_000
SEED = 0
np.random.seed(SEED)
t0 = time.time()
def log(m): print(f"[{time.time()-t0:7.1f}s] {m}", flush=True)

scvi_lat = ad.read_h5ad(f'{ROOT}/data/scvi_latent_v5_full.h5ad')          # X = scVI latent
scanvi_lat = ad.read_h5ad(f'{ROOT}/data/scanvi_comparator_latent.h5ad')   # X = scANVI latent
pre = ad.read_h5ad(f'{ROOT}/data/processed/atlas_v5_preprocessed.h5ad', backed='r')
cal = ad.read_h5ad(f'{ROOT}/data/braun_transfer_full_calibrated.h5ad', backed='r')
assert scvi_lat.n_obs == scanvi_lat.n_obs == pre.n_obs and \
    list(scvi_lat.obs_names[:3]) == list(scanvi_lat.obs_names[:3]) == list(pre.obs_names[:3]), \
    "latents + preprocessed must share row order"

rng = np.random.default_rng(SEED)
idx = np.sort(rng.choice(scvi_lat.n_obs, N_SUB, replace=False))
sub = scvi_lat[idx].copy()
sub.obs['CellClass_cal'] = cal.obs['CellClass_cal'].reindex(sub.obs_names).values
keep = sub.obs['CellClass_cal'].notna() & (sub.obs['CellClass_cal'] != 'Unknown')
kept = idx[keep.values]
sub = sub[keep.values].copy()

sub.obsm['X_scvi'] = np.asarray(sub.X)
sub.obsm['X_scanvi'] = np.asarray(scanvi_lat[kept].X)
Xp = pre[kept].to_memory().X
Xp = Xp.toarray() if sp.issparse(Xp) else np.asarray(Xp)
sub.obsm['X_pca'] = PCA(n_components=50, random_state=SEED).fit_transform(Xp)
sub.obs['dataset_slug'] = sub.obs['dataset_slug'].astype('category')
sub.obs['CellClass_cal'] = sub.obs['CellClass_cal'].astype('category')
log(f"benchmarking {sub.n_obs:,} cells; embeddings: X_pca, X_scvi, X_scanvi")

from scib_metrics.benchmark import Benchmarker
bm = Benchmarker(sub, batch_key='dataset_slug', label_key='CellClass_cal',
                 embedding_obsm_keys=['X_pca', 'X_scvi', 'X_scanvi'], n_jobs=-1)
bm.benchmark()
res = bm.get_results(min_max_scale=False)
log("=== scIB: PCA vs scVI vs scANVI (raw) ===")
log("\n" + res.T.to_string())
res.to_csv(f'{ROOT}/data/scib_metrics_comparator.tsv', sep='\t')

# GO/NO-GO embedding-level check vs scVI
def g(emb, metric): return float(res.loc[emb, metric])
d_ilisi = g('X_scanvi', 'iLISI') - g('X_scvi', 'iLISI')
d_kbet = g('X_scanvi', 'KBET') - g('X_scvi', 'KBET')
d_clisi = g('X_scvi', 'cLISI') - g('X_scanvi', 'cLISI')   # positive = scANVI lost biology
log(f"\nGO/NO-GO (vs scVI): ΔiLISI={d_ilisi:+.4f} (need >=+0.05) | "
    f"ΔkBET={d_kbet:+.4f} (need >=+0.05) | cLISI_drop={d_clisi:+.4f} (need <=0.02)")
go = (d_ilisi >= 0.05) and (d_kbet >= 0.05) and (d_clisi <= 0.02)
log(f"EMBEDDING-LEVEL VERDICT: {'GO (proceed to OOD/confidence guardrails)' if go else 'NO-GO -> freeze scVI'}")
log("saved -> data/scib_metrics_comparator.tsv")
log("DONE")
