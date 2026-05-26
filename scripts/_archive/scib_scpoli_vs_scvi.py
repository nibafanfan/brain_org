#!/usr/bin/env python3
"""Pilot go/no-go: scIB metrics for scPoli vs scVI on the SAME pilot cells.
Does scPoli mix datasets better (iLISI/kBET up) without collapsing biology (cLISI)?

RUN WITH THE scib VENV:
  /Users/eg/.venvs/scib/bin/python scripts/scib_scpoli_vs_scvi.py
"""
import sys, time
from pathlib import Path
import numpy as np, anndata as ad
sys.path.insert(0, str(Path(__file__).resolve().parent))
from atlas_common import load_config

t0 = time.time()
def log(m): print(f"[{time.time()-t0:7.1f}s] {m}", flush=True)
cfg = load_config()

sp_lat = ad.read_h5ad(cfg.root / 'data/scpoli_pilot_latent.h5ad')   # X = scPoli latent (pilot cells)
names = sp_lat.obs_names
scvi = ad.read_h5ad(cfg.latent, backed='r')
pos = {n: i for i, n in enumerate(scvi.obs_names)}
rows = np.array([pos[n] for n in names])
scvi_lat = np.asarray(scvi[rows].to_memory().X)

# real X (preprocessed lognorm HVG, same cells) so the Benchmarker's internal PCA works
pre = ad.read_h5ad(cfg.preprocessed, backed='r')
a = ad.AnnData(pre[rows].to_memory().X)
a.obs_names = names
a.obsm['X_scvi'] = scvi_lat
a.obsm['X_scpoli'] = np.asarray(sp_lat.X)
a.obs['dataset_slug'] = sp_lat.obs['dataset_slug'].astype('category')
lab = sp_lat.obs['CellClass_cal'].astype(str)
a.obs['CellClass_cal'] = lab.values
a = a[(lab != 'Unknown').values].copy()
a.obs['CellClass_cal'] = a.obs['CellClass_cal'].astype('category')
log(f"benchmarking {a.n_obs:,} cells; {a.obs['dataset_slug'].nunique()} datasets")

from scib_metrics.benchmark import Benchmarker
bm = Benchmarker(a, batch_key='dataset_slug', label_key='CellClass_cal',
                 embedding_obsm_keys=['X_scvi', 'X_scpoli'], n_jobs=-1)
bm.benchmark()
res = bm.get_results(min_max_scale=False)
log("=== scIB: scVI vs scPoli (pilot) ===")
log("\n" + res.T.to_string())
res.to_csv(cfg.root / 'data/scib_scpoli_vs_scvi.tsv', sep='\t')

def g(emb, m): return float(res.loc[emb, m])
d_il, d_kb = g('X_scpoli', 'iLISI') - g('X_scvi', 'iLISI'), g('X_scpoli', 'KBET') - g('X_scvi', 'KBET')
d_cl = g('X_scvi', 'cLISI') - g('X_scpoli', 'cLISI')
log(f"\nvs scVI: ΔiLISI={d_il:+.4f} ΔkBET={d_kb:+.4f} cLISI_drop={d_cl:+.4f}")
log("GO if scPoli materially improves iLISI/kBET without cLISI collapse (drop<=0.02)")
log("DONE")
