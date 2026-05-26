#!/usr/bin/env python3
"""Finalize Braun label transfer (Stage 2 only) reusing the diagnosed-healthy
reference model data/braun_scanvi_full. FIX: the scANVI classifier head
collapses to one class on the organoid query, so transfer BOTH CellClass and
Region by kNN on the joint scArches latent (the latent generalizes; the head
does not). Surgery ~15 epochs is enough (latent stable 10 vs 40).

Memory-safe: streams the 4M-cell query, keeping only the model's genes.
  /opt/homebrew/Caskroom/miniforge/base/bin/python3.13 scripts/braun_transfer_finalize.py
"""
import os
os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")
import sys, time
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad, scvi
import scipy.sparse as sp
from sklearn.neighbors import KNeighborsClassifier
sys.path.insert(0, str(Path(__file__).resolve().parent))
from atlas_common import load_config, sym2ens as _sym2ens, model_genes, reindex_braun, read_atlas_genes

cfg = load_config()
ROOT = cfg.root
MDIR = cfg.braun_scanvi_model
ACC = 'mps'
SURGERY_EPOCHS = 15
KNN_REF_N = 400_000          # Braun cells used to fit the kNN label transfer
OUT = ROOT / 'data/braun_transfer_full_knn.h5ad'
t0 = time.time()
def log(m): print(f"[{time.time()-t0:8.1f}s] {m}", flush=True)

vn = model_genes(MDIR)
s2e = _sym2ens(cfg.canonical)
log(f"model genes: {len(vn)}")

# ---- reference + model
braun = reindex_braun(cfg.braun, vn)
braun.obs['CellClass'] = braun.obs['CellClass'].astype(str)
braun.obs['Region'] = braun.obs['Region'].astype(str)
braun.obs['donor_id'] = braun.obs['donor_id'].astype(str)
scanvi = scvi.model.SCANVI.load(str(MDIR), adata=braun)
log(f"reference + model ready: {braun.shape}")

# ---- reference latent (subsample) to fit kNN transfer
rng = np.random.default_rng(0)
ridx = np.sort(rng.choice(braun.n_obs, min(KNN_REF_N, braun.n_obs), replace=False))
ref_lat = scanvi.get_latent_representation(braun[ridx].copy())
knn_cls = KNeighborsClassifier(30, n_jobs=-1).fit(ref_lat, braun.obs['CellClass'].to_numpy()[ridx])
knn_reg = KNeighborsClassifier(30, n_jobs=-1).fit(ref_lat, braun.obs['Region'].to_numpy()[ridx])
log(f"kNN fit on {len(ridx):,} ref cells")

# ---- full query, streamed, model genes only (memory-safe)
atlas_b = ad.read_h5ad(ROOT/'data/atlas_v5_full.h5ad', backed='r')
Xq, present = read_atlas_genes(atlas_b, vn, s2e)        # memory-safe chunked read (all rows)
query = ad.AnnData(X=Xq, obs=atlas_b.obs.copy())
query.var_names = present
query.layers['counts'] = query.X.copy()
query.obs['CellClass'] = 'Unknown'
query.obs['Region'] = 'Unknown'
query.obs['donor_id'] = query.obs['tech_sample'].astype(str)
log(f"query built: {query.shape} ({len(present)}/{len(vn)} model genes present)")

# ---- scArches surgery -> joint latent
scvi.model.SCANVI.prepare_query_anndata(query, scanvi)
qm = scvi.model.SCANVI.load_query_data(query, scanvi)
log(f"surgery training query ({SURGERY_EPOCHS} ep)…")
qm.train(max_epochs=SURGERY_EPOCHS, plan_kwargs={'weight_decay': 0.0},
         accelerator=ACC, batch_size=512, enable_progress_bar=False)
q_lat = qm.get_latent_representation()
log(f"query latent: {q_lat.shape}")

# ---- kNN transfer (CellClass + Region) with confidence
def transfer(knn, X):
    proba = knn.predict_proba(X)
    return knn.classes_[proba.argmax(1)], proba.max(1)
cc_pred, cc_conf = transfer(knn_cls, q_lat)
rg_pred, rg_conf = transfer(knn_reg, q_lat)
log("=== CellClass (kNN-on-latent) ===")
for k, v in pd.Series(cc_pred).value_counts(normalize=True).items():
    log(f"   {k:24} {v*100:5.1f}%")
log(f"   mean conf {cc_conf.mean():.3f} | frac>0.6 {(cc_conf>0.6).mean():.2f}")
log("=== Region (kNN-on-latent) ===")
for k, v in pd.Series(rg_pred).value_counts(normalize=True).items():
    log(f"   {k:24} {v*100:5.1f}%")
log(f"   mean conf {rg_conf.mean():.3f}")

# crosstab vs organoid cell_type_original (meaningful labels)
NULLISH = {'unknown', 'nan', 'none', '', 'na'}
col = 'cell_type_original'
if col in query.obs:
    m = ~query.obs[col].astype(str).str.lower().isin(NULLISH)
    if m.sum():
        log(f"=== CellClass vs organoid {col} ({m.sum():,} labeled) ===")
        log("\n" + pd.crosstab(query.obs.loc[m.values, col],
                               pd.Series(cc_pred, index=query.obs_names)[m.values]).to_string())

out = ad.AnnData(X=q_lat, obs=query.obs.copy())
out.obs['CellClass_pred'] = cc_pred; out.obs['CellClass_conf'] = cc_conf
out.obs['Region_pred'] = rg_pred;    out.obs['Region_conf'] = rg_conf
out.write_h5ad(OUT)
log(f"saved -> {OUT}")
log("DONE")
