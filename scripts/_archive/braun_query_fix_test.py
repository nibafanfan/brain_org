#!/usr/bin/env python3
"""Isolate the query-side 100%-RG collapse using the EXISTING good reference
model (data/braun_scanvi_full, dia­gnosed healthy: 88% ref accuracy).

On a 200k organoid subsample, for each surgery-epoch setting, compare:
  (A) scANVI classifier head  q_model.predict()
  (B) kNN on the joint latent  (Braun ref latent -> query latent)
Both should yield a sane spread (NOT ~100% Radial glia). This tells us whether
the fix is "fewer surgery epochs", "use kNN", or both.

  /opt/homebrew/Caskroom/miniforge/base/bin/python3.13 scripts/braun_query_fix_test.py
"""
import os
os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")
import sys, time
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad, scvi
from lightning.pytorch.callbacks import Callback
from sklearn.neighbors import KNeighborsClassifier
sys.path.insert(0, str(Path(__file__).resolve().parent))
from atlas_common import load_config, sym2ens as _sym2ens, model_genes, reindex_braun

cfg = load_config()
ROOT = cfg.root
MDIR = cfg.braun_scanvi_model
ACC = 'mps'
QUERY_N = 200_000
SURGERY_EPOCHS = [40, 10, 0]     # 40 reproduces the collapse; 10/0 are candidate fixes
t0 = time.time()
def log(m): print(f"[{time.time()-t0:7.1f}s] {m}", flush=True)

vn = model_genes(MDIR)
sym2ens = _sym2ens(cfg.canonical)

# ---- reference (reindex Braun to model genes) + load model
braun = reindex_braun(cfg.braun, vn)
braun.obs['CellClass'] = braun.obs['CellClass'].astype(str)
braun.obs['donor_id'] = braun.obs['donor_id'].astype(str)
scanvi = scvi.model.SCANVI.load(str(MDIR), adata=braun)
log(f"reference ready: {braun.shape}; model loaded")

# reference latent + labels for kNN (subsample for speed)
rng = np.random.default_rng(0)
ridx = np.sort(rng.choice(braun.n_obs, 300_000, replace=False))
ref_lat = scanvi.get_latent_representation(braun[ridx].copy())
ref_cls = braun.obs['CellClass'].to_numpy()[ridx]
log(f"ref latent for kNN: {ref_lat.shape}")

# ---- organoid query subsample, mapped to model genes
atlas_b = ad.read_h5ad(ROOT/'data/atlas_v5_full.h5ad', backed='r')
atlas_ens = np.array([sym2ens.get(s, '') for s in atlas_b.var_names])
ens_to_col = {}
for c, e in enumerate(atlas_ens):
    if e and e not in ens_to_col:
        ens_to_col[e] = c
qidx = np.sort(rng.choice(atlas_b.n_obs, QUERY_N, replace=False))
sub = atlas_b[qidx].to_memory()
cols = [ens_to_col[g] for g in vn if g in ens_to_col]
present = [g for g in vn if g in ens_to_col]
import scipy.sparse as sp
Xq = sub.X[:, cols].astype('float32')
query0 = ad.AnnData(X=Xq, obs=sub.obs.copy())
query0.var_names = present
log(f"query built: {query0.shape} ({len(present)}/{len(vn)} model genes present)")

def make_query():
    q = query0.copy()
    q.layers['counts'] = q.X.copy()
    q.obs['CellClass'] = 'Unknown'
    q.obs['donor_id'] = q.obs['tech_sample'].astype(str)
    scvi.model.SCANVI.prepare_query_anndata(q, scanvi)
    return q

knn = KNeighborsClassifier(n_neighbors=30, n_jobs=-1).fit(ref_lat, ref_cls)

def dist(arr, top=6):
    vc = pd.Series(arr).value_counts(normalize=True)
    return " ".join(f"{k}={v*100:.0f}%" for k, v in vc.head(top).items())

for ep in SURGERY_EPOCHS:
    q = make_query()
    qm = scvi.model.SCANVI.load_query_data(q, scanvi)
    if ep > 0:
        qm.train(max_epochs=ep, plan_kwargs={'weight_decay': 0.0},
                 accelerator=ACC, batch_size=512, enable_progress_bar=False)
    head = np.asarray(qm.predict())
    lat = qm.get_latent_representation()
    knn_pred = knn.predict(lat)
    log(f"=== surgery_epochs={ep} ===")
    log(f"   (A) scANVI head : {dist(head)}")
    log(f"   (B) kNN-latent  : {dist(knn_pred)}")
log("DONE")
