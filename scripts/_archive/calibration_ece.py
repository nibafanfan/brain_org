#!/usr/bin/env python3
"""Calibration reliability + ECE for the transfer confidence (review refinement).

Organoid cells lack ground truth, so we calibrate on Braun (known labels): fit the
SAME weighted + class-balanced kNN used by braun_transfer_calibrated on a train split,
predict a held-out Braun test split, and compare the max-posterior CONFIDENCE to the
empirical ACCURACY per confidence bin. Reports a reliability table + ECE.
  ECE = sum_bins (n_bin/N) * |acc_bin - conf_bin|   (lower = better calibrated)
"""
import argparse, sys, time
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad, scvi
from sklearn.neighbors import NearestNeighbors
sys.path.insert(0, str(Path(__file__).resolve().parent))
from atlas_common import load_config, model_genes, reindex_braun, write_sidecar

ap = argparse.ArgumentParser(); ap.add_argument('--root', default=None); args = ap.parse_args()
cfg = load_config(root=args.root)
K, CAP, EPS, NBINS = 30, 6000, 1e-6, 10
t0 = time.time()
def log(m): print(f"[{time.time()-t0:7.1f}s] {m}", flush=True)
rng = np.random.default_rng(0)

vn = model_genes(cfg.braun_scanvi_model)
braun = reindex_braun(cfg.braun, vn)
braun.obs['CellClass'] = braun.obs['CellClass'].astype(str)
braun.obs['donor_id'] = braun.obs['donor_id'].astype(str)
scanvi = scvi.model.SCANVI.load(str(cfg.braun_scanvi_model), adata=braun)

cls_all = braun.obs['CellClass'].to_numpy()
bal_idx = np.concatenate([
    rng.choice(np.where(cls_all == c)[0], min(CAP, (cls_all == c).sum()), replace=False)
    for c in np.unique(cls_all)])
rng.shuffle(bal_idx)
lat = scanvi.get_latent_representation(braun[bal_idx].copy())
lab = cls_all[bal_idx]
classes = np.unique(lab)
codes = pd.Categorical(lab, categories=classes).codes
log(f"balanced Braun: {len(bal_idx):,} cells, {len(classes)} classes")

# train/test split, weighted-kNN posterior on test (mirrors calibrated transfer)
perm = rng.permutation(len(bal_idx)); cut = int(0.8 * len(perm))
tr, te = perm[:cut], perm[cut:]
nn = NearestNeighbors(n_neighbors=K, n_jobs=-1).fit(lat[tr])
dist, idx = nn.kneighbors(lat[te])
w = 1.0 / (dist + EPS)
ctr = codes[tr][idx]
proba = np.zeros((len(te), len(classes)), np.float32)
for c in range(len(classes)):
    proba[:, c] = np.where(ctr == c, w, 0.0).sum(1)
proba /= proba.sum(1, keepdims=True)
conf = proba.max(1)
pred = proba.argmax(1)
correct = (pred == codes[te]).astype(float)
log(f"overall test accuracy {correct.mean():.3f}; mean confidence {conf.mean():.3f}")

# reliability bins + ECE
edges = np.linspace(0, 1, NBINS + 1)
rows, ece = [], 0.0
for i in range(NBINS):
    m = (conf >= edges[i]) & (conf < edges[i + 1] if i < NBINS - 1 else conf <= edges[i + 1])
    n = int(m.sum())
    if n == 0:
        rows.append((f"[{edges[i]:.1f},{edges[i+1]:.1f})", 0, np.nan, np.nan)); continue
    acc, cf = correct[m].mean(), conf[m].mean()
    ece += n / len(te) * abs(acc - cf)
    rows.append((f"[{edges[i]:.1f},{edges[i+1]:.1f})", n, round(cf, 3), round(acc, 3)))
df = pd.DataFrame(rows, columns=['conf_bin', 'n', 'mean_conf', 'accuracy'])
log("=== reliability (Braun held-out test; weighted balanced kNN) ===")
log("\n" + df.to_string(index=False))
log(f"\nECE = {ece:.4f}  (lower=better; <0.05 well-calibrated)")
out_tsv = f'{cfg.root}/data/calibration_ece.tsv'
df.to_csv(out_tsv, sep='\t', index=False)
write_sidecar(out_tsv, __file__, {'K': K, 'CAP': CAP, 'NBINS': NBINS,
                                  'test_accuracy': float(correct.mean()),
                                  'mean_conf': float(conf.mean()), 'ECE': float(ece)})
log(f"saved -> {out_tsv}")
log("DONE")
