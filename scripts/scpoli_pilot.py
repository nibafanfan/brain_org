#!/usr/bin/env python3
"""scPoli integration PILOT (HNOCA's method) on an organoid-atlas subsample, to test
whether it mixes datasets better than scVI before committing to a full run.

condition = bio_sample, cell-type labels = CellClass_cal (semi-supervised; 'Unknown'
= abstained). Saves the scPoli latent for a side-by-side scIB comparison vs scVI.

RUN WITH THE scPoli VENV:
  /Users/eg/.venvs/scpoli/bin/python scripts/scpoli_pilot.py
"""
import anndata as ad
ad.read = ad.read_h5ad                      # shim: scarches 0.6.1 imports the removed alias
import sys, time
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
sys.path.insert(0, str(Path(__file__).resolve().parent))
from atlas_common import load_config
from scarches.models.scpoli import scPoli

N_SUB = 120_000
N_EPOCHS, PRETRAIN = 12, 8
t0 = time.time()
def log(m): print(f"[{time.time()-t0:7.1f}s] {m}", flush=True)

cfg = load_config()
rng = np.random.default_rng(cfg.defaults['seed'])
pre = ad.read_h5ad(cfg.preprocessed)
cal = ad.read_h5ad(cfg.calibrated, backed='r')
idx = np.sort(rng.choice(pre.n_obs, min(N_SUB, pre.n_obs), replace=False))
a = pre[idx].copy()
a.X = a.layers['counts'].copy()             # scPoli (nb recon) needs raw counts in X
a.obs['bio_sample'] = a.obs['bio_sample'].astype(str)
a.obs['CellClass_cal'] = cal.obs['CellClass_cal'].reindex(a.obs_names).fillna('Unknown').astype(str).values
log(f"subsample {a.shape}; {a.obs['bio_sample'].nunique()} conditions; "
    f"{a.obs['CellClass_cal'].nunique()} cell types")

model = scPoli(adata=a, condition_keys=['bio_sample'], cell_type_keys=['CellClass_cal'],
               embedding_dims=5, recon_loss='nb', unknown_ct_names=['Unknown'])
log("scPoli built; training…")
model.train(n_epochs=N_EPOCHS, pretraining_epochs=PRETRAIN, eta=5)
model.save(str(cfg.root / 'data/scpoli_pilot_model'), overwrite=True)   # save BEFORE latent (retry-safe)
log("trained + model saved")

# manual latent extraction (scarches get_latent indexes X with a torch tensor -> numpy
# crash on this version; replicate it with numpy batch indices).
import torch
device = next(model.model.parameters()).device
X = a.X
clabels = []
for cond in model.condition_keys_:
    vals = a.obs[cond].values
    lab = np.zeros(len(vals))
    for condition, label in model.model.condition_encoders[cond].items():
        lab[vals == condition] = label
    clabels.append(lab)
cmat = torch.tensor(np.asarray(clabels), device=device).T
lat_parts = []
for b in np.array_split(np.arange(X.shape[0]), max(1, X.shape[0] // 512)):
    xb = X[b, :]
    xb = xb.toarray() if sp.issparse(xb) else np.asarray(xb)
    xb = torch.tensor(xb, device=device).float()
    z = model.model.get_latent(xb, cmat[b, :], True)
    lat_parts.append(z.cpu().detach())
lat = torch.cat(lat_parts).numpy()
log(f"scPoli latent {lat.shape}")

out = ad.AnnData(X=np.asarray(lat).astype('float32'), obs=a.obs.copy())
out.uns['scpoli'] = {'n_sub': int(len(idx)), 'n_epochs': N_EPOCHS, 'pretrain': PRETRAIN,
                     'condition': 'bio_sample', 'cell_type': 'CellClass_cal'}
out.write_h5ad(cfg.root / 'data/scpoli_pilot_latent.h5ad')
log(f"saved -> data/scpoli_pilot_latent.h5ad ({out.shape})")
log("DONE")
