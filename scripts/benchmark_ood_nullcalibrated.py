#!/usr/bin/env python3
"""Null-calibrated OOD (review #1 refinement). The earlier OOD used the reference's
own self-distance p95 as the threshold; here we build a proper in-distribution null:
split Braun ref into train/test, compute test->train kNN-distance (cells we KNOW are
in-distribution), and set the OOD threshold at p95/p99 of that null. Then query
cells beyond it are OOD. Also a CellClass-stratified per-class null.
"""
import sys, time
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad, scvi, torch
from sklearn.neighbors import NearestNeighbors
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _provenance import write_sidecar

ROOT = '/Users/eg/brain_organoid'
N_TRAIN, N_TEST, N_QUERY = 300_000, 50_000, 500_000
K = 30
t0 = time.time()
def log(m): print(f"[{time.time()-t0:7.1f}s] {m}", flush=True)
rng = np.random.default_rng(0)

vn = list(torch.load(f'{ROOT}/data/braun_scanvi_full/model.pt',
                     map_location='cpu', weights_only=False)['var_names'])
braun = ad.read_h5ad(f'{ROOT}/data/raw/braun_2023/braun_all.h5ad')[:, vn].copy()
braun.layers['counts'] = braun.X.copy()
braun.obs['CellClass'] = braun.obs['CellClass'].astype(str)
braun.obs['donor_id'] = braun.obs['donor_id'].astype(str)
scanvi = scvi.model.SCANVI.load(f'{ROOT}/data/braun_scanvi_full', adata=braun)

sel = rng.choice(braun.n_obs, N_TRAIN + N_TEST, replace=False)
tr_idx, te_idx = sel[:N_TRAIN], sel[N_TRAIN:]
ref_lat = scanvi.get_latent_representation(braun[sel].copy())
tr_lat, te_lat = ref_lat[:N_TRAIN], ref_lat[N_TRAIN:]
tr_cls = braun.obs['CellClass'].to_numpy()[tr_idx]
te_cls = braun.obs['CellClass'].to_numpy()[te_idx]
log(f"ref latent: train {tr_lat.shape} test {te_lat.shape}")

nn = NearestNeighbors(n_neighbors=K, n_jobs=-1).fit(tr_lat)
te_d = nn.kneighbors(te_lat)[0].mean(1)            # in-distribution null
p95, p99 = np.percentile(te_d, [95, 99])
log(f"in-distribution null (test->train mean kNN-dist): p95={p95:.3f} p99={p99:.3f}")

cal = ad.read_h5ad(f'{ROOT}/data/braun_transfer_full_calibrated.h5ad', backed='r')
qi = np.sort(rng.choice(cal.n_obs, N_QUERY, replace=False))
q_lat = np.asarray(cal[qi].to_memory().X)
q_cls = cal.obs['CellClass_cal'].to_numpy()[qi]
q_d = nn.kneighbors(q_lat)[0].mean(1)
log("=== query OOD fraction (null-calibrated) ===")
log(f"  > p95: {(q_d > p95).mean()*100:.1f}%   > p99: {(q_d > p99).mean()*100:.1f}%")
log(f"  (prior ref-self-p95 method reported 78.3%)")

# CellClass-stratified per-class null thresholds
log("=== CellClass-stratified OOD (per-class test->train p95) ===")
rows = []
for c in np.unique(te_cls):
    m_te = te_cls == c
    if m_te.sum() < 100:
        continue
    nn_c = NearestNeighbors(n_neighbors=K, n_jobs=-1).fit(tr_lat[tr_cls == c]) \
        if (tr_cls == c).sum() > K else None
    if nn_c is None:
        continue
    null_c = nn_c.kneighbors(te_lat[m_te])[0].mean(1)
    thr95, thr99 = np.percentile(null_c, [95, 99])
    m_q = q_cls == c
    if m_q.sum() == 0:
        ood95 = ood99 = np.nan
    else:
        qd = nn_c.kneighbors(q_lat[m_q])[0].mean(1)        # query->train dist, computed once
        ood95 = (qd > thr95).mean() * 100
        ood99 = (qd > thr99).mean() * 100
    rows.append((c, int(m_q.sum()), round(thr95, 3), round(thr99, 3),
                 round(ood95, 1) if m_q.sum() else None,
                 round(ood99, 1) if m_q.sum() else None))
strat = pd.DataFrame(rows, columns=['CellClass', 'n_query', 'class_p95_thr', 'class_p99_thr',
                                    'ood_pct_p95', 'ood_pct_p99'])
log("\n" + strat.to_string(index=False))
out_tsv = f'{ROOT}/data/ood_nullcalibrated.tsv'
strat.to_csv(out_tsv, sep='\t', index=False)
write_sidecar(out_tsv, __file__, {'N_TRAIN': N_TRAIN, 'N_TEST': N_TEST, 'N_QUERY': N_QUERY,
                                  'K': K, 'global_p95': float(p95), 'global_p99': float(p99),
                                  'ood_p95_pct': float((q_d > p95).mean()*100),
                                  'ood_p99_pct': float((q_d > p99).mean()*100)})
log(f"saved -> {out_tsv}")
log("DONE")
