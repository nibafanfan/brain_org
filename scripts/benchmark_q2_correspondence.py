#!/usr/bin/env python3
"""Q2 Part 1: do organoid cell types transcriptomically correspond to primary?

Pseudobulk per CellClass (mean raw counts -> CP10K -> log1p) on the 2006 model
genes, for organoid (300k subsample) and Braun. Correlate organoid-class vs
Braun-class. Diagonal dominance (organoid-C best matches Braun-C) = correspondence.
"""
import argparse, sys, time
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad
import scipy.sparse as sp
sys.path.insert(0, str(Path(__file__).resolve().parent))
from atlas_common import load_config, sym2ens, model_genes, read_atlas_genes, reindex_braun, cp10k_log

ap = argparse.ArgumentParser(); ap.add_argument('--root', default=None); args = ap.parse_args()
cfg = load_config(root=args.root)
ROOT = cfg.root
N_SUB = cfg.defaults['n_sub']
t0 = time.time()
def log(m): print(f"[{time.time()-t0:7.1f}s] {m}", flush=True)

vn = model_genes(cfg.braun_scanvi_model)
s2e = sym2ens(cfg.canonical)

def _dense(X):
    return X.toarray() if sp.issparse(X) else np.asarray(X)

# ---- Braun pseudobulk
braun = reindex_braun(cfg.braun, vn)
braun.obs['CellClass'] = braun.obs['CellClass'].astype(str)
bdf = pd.DataFrame(_dense(braun.X), columns=vn, index=braun.obs_names)
bdf['__c'] = braun.obs['CellClass'].values
b_pb = cp10k_log(bdf.groupby('__c').mean())
log(f"Braun pseudobulk: {b_pb.shape}")

# ---- organoid pseudobulk (subsample), via shared memory-safe reader
tr = ad.read_h5ad(ROOT / 'data/braun_transfer_full_knn.h5ad', backed='r')
cc = tr.obs['CellClass_pred']
atlas = ad.read_h5ad(cfg.atlas_full, backed='r')
rng = np.random.default_rng(0)
idx = np.sort(rng.choice(atlas.n_obs, N_SUB, replace=False))
oX, present = read_atlas_genes(atlas, vn, s2e, row_idx=idx)
obs_names = atlas.obs_names[idx]
odf = pd.DataFrame(_dense(oX), columns=present, index=obs_names)
odf['__c'] = cc.reindex(obs_names).values
o_pb = cp10k_log(odf.groupby('__c').mean())
log(f"organoid pseudobulk: {o_pb.shape}")

# ---- correlation matrix (organoid rows x Braun cols)
shared_cls = [c for c in o_pb.index if c in b_pb.index]
corr = pd.DataFrame(index=o_pb.index, columns=b_pb.index, dtype=float)
for oc in o_pb.index:
    for bc in b_pb.index:
        corr.loc[oc, bc] = np.corrcoef(o_pb.loc[oc], b_pb.loc[bc])[0, 1]
corr = corr.astype(float)
corr.to_csv(ROOT / 'data/q2_correspondence_corr.tsv', sep='\t')

log("=== organoid CellClass -> best-matching Braun CellClass ===")
for oc in o_pb.index:
    row = corr.loc[oc]
    best = row.idxmax()
    diag = row[oc] if oc in row.index else np.nan
    flag = "OK" if best == oc else f"!! best={best}"
    log(f"  {oc:14} self-corr={diag:.3f}  argmax={best:14} (r={row[best]:.3f})  {flag}")
log(f"\ndiagonal dominance: {sum(corr.loc[c].idxmax()==c for c in shared_cls)}/{len(shared_cls)} "
    f"organoid classes best-match their own Braun class")
log("full matrix -> data/q2_correspondence_corr.tsv")
log("DONE")
