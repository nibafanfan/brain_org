#!/usr/bin/env python3
"""Review item #1: calibrated kNN label transfer (replaces argmax/unweighted).

Reuses the EXISTING joint latents (no surgery re-run): query latent from
data/braun_transfer_full_knn.h5ad (X), Braun reference latent from the
braun_scanvi_full model. Adds:
  - distance-weighted votes
  - class-BALANCED reference sampling (prior correction -> rare-class recall)
  - abstention: max posterior < TAU -> 'Unknown'
  - OOD diagnostic: query kNN-distance vs the reference's own in-distribution dist
Reports rare-class (Immune/Vascular) recovery vs the old argmax run, abstention
rate, and a held-out Braun accuracy/recall check (balanced+weighted vs naive).
"""
import sys, time
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad, scvi, torch
from sklearn.neighbors import NearestNeighbors
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _provenance import stamp

ROOT = '/Users/eg/brain_organoid'
MDIR = f'{ROOT}/data/braun_scanvi_full'
K = 30
CAP = 6000          # max ref cells per CellClass (balanced sampling)
TAU = 0.40          # abstention threshold on max posterior
EPS = 1e-6
t0 = time.time()
def log(m): print(f"[{time.time()-t0:7.1f}s] {m}", flush=True)
rng = np.random.default_rng(0)

# ---- Braun reference latent (balanced subset only -> cheaper)
vn = list(torch.load(f'{MDIR}/model.pt', map_location='cpu', weights_only=False)['var_names'])
braun = ad.read_h5ad(f'{ROOT}/data/raw/braun_2023/braun_all.h5ad')[:, vn].copy()
braun.layers['counts'] = braun.X.copy()
braun.obs['CellClass'] = braun.obs['CellClass'].astype(str)
braun.obs['donor_id'] = braun.obs['donor_id'].astype(str)
scanvi = scvi.model.SCANVI.load(MDIR, adata=braun)

cls_all = braun.obs['CellClass'].to_numpy()
bal_idx = np.concatenate([
    rng.choice(np.where(cls_all == c)[0], min(CAP, (cls_all == c).sum()), replace=False)
    for c in np.unique(cls_all)])
rng.shuffle(bal_idx)
ref_lat = scanvi.get_latent_representation(braun[bal_idx].copy())
ref_lab = cls_all[bal_idx]
classes = np.unique(ref_lab)
log(f"balanced ref: {len(bal_idx):,} cells, {len(classes)} classes "
    f"(per-class {pd.Series(ref_lab).value_counts().min()}-{pd.Series(ref_lab).value_counts().max()})")

# ---- weighted-kNN posterior helper (manual: one neighbor search)
nn = NearestNeighbors(n_neighbors=K, n_jobs=-1).fit(ref_lat)
ref_codes = pd.Categorical(ref_lab, categories=classes).codes

def weighted_posterior(Q):
    dist, idx = nn.kneighbors(Q)
    w = 1.0 / (dist + EPS)
    codes = ref_codes[idx]                                  # (n, K)
    proba = np.zeros((Q.shape[0], len(classes)), np.float32)
    for c in range(len(classes)):
        proba[:, c] = np.where(codes == c, w, 0.0).sum(1)
    proba /= proba.sum(1, keepdims=True)
    return proba, dist

# ---- held-out Braun check: balanced+weighted vs naive(unbalanced+unweighted)
te = rng.permutation(len(bal_idx))
cut = int(0.8 * len(te)); tr, ts = te[:cut], te[cut:]
nn_tr = NearestNeighbors(n_neighbors=K, n_jobs=-1).fit(ref_lat[tr])
d, ii = nn_tr.kneighbors(ref_lat[ts])
codes_tr = ref_codes[tr][ii]
# weighted
wp = np.zeros((len(ts), len(classes)), np.float32)
for c in range(len(classes)):
    wp[:, c] = np.where(codes_tr == c, 1/(d+EPS), 0).sum(1)
pred_w = classes[wp.argmax(1)]
# unweighted
up = np.zeros((len(ts), len(classes)), np.float32)
for c in range(len(classes)):
    up[:, c] = (codes_tr == c).sum(1)
pred_u = classes[up.argmax(1)]
true = ref_lab[ts]
def recall(pred, lab):
    return {c: round((pred[true == c] == c).mean(), 2) for c in ['Immune', 'Vascular', 'Oligo', 'Neuron']}
log(f"held-out Braun acc: weighted={np.mean(pred_w==true):.3f} unweighted={np.mean(pred_u==true):.3f}")
log(f"  rare-class recall weighted:   {recall(pred_w, true)}")
log(f"  rare-class recall unweighted: {recall(pred_u, true)}")

# ---- apply to organoid query (reuse saved latent)
tr_h5 = ad.read_h5ad(f'{ROOT}/data/braun_transfer_full_knn.h5ad')
Q = np.asarray(tr_h5.X)
old_pred = tr_h5.obs['CellClass_pred'].to_numpy()
log(f"query latent {Q.shape}")
proba, dist = weighted_posterior(Q)
conf = proba.max(1)
new_pred = classes[proba.argmax(1)].astype(object)
abstain = conf < TAU
new_pred[abstain] = 'Unknown'

# OOD: query mean-kNN-distance vs reference in-distribution (ref self kNN dist)
d_ref, _ = nn.kneighbors(ref_lat[rng.choice(len(ref_lat), 20000, replace=False)])
ref_p95 = np.percentile(d_ref.mean(1), 95)
ood = dist.mean(1) > ref_p95
log(f"abstain (conf<{TAU}): {abstain.mean()*100:.1f}% | OOD (dist>ref p95): {ood.mean()*100:.1f}%")

log("=== CellClass distribution: OLD argmax vs NEW calibrated ===")
old_d = pd.Series(old_pred).value_counts(normalize=True)
new_d = pd.Series(new_pred).value_counts(normalize=True)
cmp = pd.DataFrame({'old_%': (old_d*100).round(2), 'new_%': (new_d*100).round(2)}).fillna(0)
log("\n" + cmp.sort_values('new_%', ascending=False).to_string())

out = ad.AnnData(X=Q, obs=tr_h5.obs.copy())
out.obs['CellClass_cal'] = new_pred
out.obs['CellClass_cal_conf'] = conf
out.obs['abstain'] = abstain
out.obs['ood'] = ood
stamp(out, __file__, {'K': K, 'CAP': CAP, 'TAU': TAU, 'EPS': EPS,
                      'ref_p95_dist': float(ref_p95), 'n_ref_balanced': int(len(bal_idx))})
out.write_h5ad(f'{ROOT}/data/braun_transfer_full_calibrated.h5ad')
log("saved -> data/braun_transfer_full_calibrated.h5ad")
log("DONE")
