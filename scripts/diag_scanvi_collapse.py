#!/usr/bin/env python3
"""Localize the 100%-Radial-glia collapse: does the saved scANVI classifier
collapse on the BRAUN REFERENCE itself (known labels)?

Reindexes Braun directly to the model's stored 2006 genes (read from model.pt),
loads data/braun_scanvi_full, predicts on a reference subsample, crosstabs vs truth.

  ~100% RG on Braun too  -> classifier collapsed in TRAINING (fix scANVI config)
  Braun predicted well   -> collapse is QUERY/surgery-specific (use kNN-on-latent)
"""
import os
os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")
import sys, time
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad, scvi
sys.path.insert(0, str(Path(__file__).resolve().parent))
from atlas_common import load_config, model_genes, reindex_braun

cfg = load_config()
ROOT = cfg.root
MDIR = cfg.braun_scanvi_model
t0 = time.time()
def log(m): print(f"[{time.time()-t0:7.1f}s] {m}", flush=True)

vn = model_genes(MDIR)
log(f"model expects {len(vn)} genes")

braun = reindex_braun(cfg.braun, vn)             # reindex to model genes + counts layer
log(f"Braun reindexed: {braun.shape}")
braun.obs['CellClass'] = braun.obs['CellClass'].astype(str)
braun.obs['donor_id'] = braun.obs['donor_id'].astype(str)
log(f"reindexed reference: {braun.shape}")

log("=== TRUE Braun CellClass distribution ===")
for k, v in braun.obs['CellClass'].value_counts(normalize=True).items():
    log(f"   {k:24} {v*100:5.1f}%")

model = scvi.model.SCANVI.load(str(MDIR), adata=braun)
log("scANVI model loaded")

rng = np.random.default_rng(0)
idx = np.sort(rng.choice(braun.n_obs, 50000, replace=False))
sub = braun[idx].copy()
pred = np.asarray(model.predict(sub))
log("=== PREDICTED CellClass on Braun reference subsample (50k) ===")
for k, v in pd.Series(pred).value_counts(normalize=True).items():
    log(f"   {k:24} {v*100:5.1f}%")
true = sub.obs['CellClass'].to_numpy()
log(f"reference prediction accuracy vs true labels: {(pred==true).mean():.3f}")
# is the latent itself class-separable? (validates kNN-on-latent fix)
lat = model.get_latent_representation(sub)
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import cross_val_score
sc_acc = cross_val_score(KNeighborsClassifier(30), lat, true, cv=3).mean()
log(f"latent kNN 3-fold CV accuracy on Braun CellClass: {sc_acc:.3f}  "
    f"(high => kNN-on-latent is a sound fix)")
log("DONE")
