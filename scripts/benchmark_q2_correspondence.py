#!/usr/bin/env python3
"""Q2 Part 1: do organoid cell types transcriptomically correspond to primary?

Pseudobulk per CellClass (mean raw counts -> CP10K -> log1p) on the 2006 model
genes, for organoid (300k subsample) and Braun. Correlate organoid-class vs
Braun-class. Diagonal dominance (organoid-C best matches Braun-C) = correspondence.
"""
import time
import numpy as np, pandas as pd, anndata as ad, torch

ROOT = '/Users/eg/brain_organoid'
N_SUB = 300_000
t0 = time.time()
def log(m): print(f"[{time.time()-t0:7.1f}s] {m}", flush=True)

vn = list(torch.load(f'{ROOT}/data/braun_scanvi_full/model.pt',
                     map_location='cpu', weights_only=False)['var_names'])
can = pd.read_csv(f'{ROOT}/data/reference/hnoca_var_canonical.tsv', sep='\t')
sym2ens = {s: e for s, e in zip(can['hgnc_symbol'].astype(str), can['ensembl'].astype(str))
           if isinstance(e, str) and e.startswith('ENSG')}

def pseudobulk(counts_df_by_class):
    """rows=class, cols=gene mean raw counts -> CP10K + log1p."""
    cp = counts_df_by_class.div(counts_df_by_class.sum(1), axis=0) * 1e4
    return np.log1p(cp)

# ---- Braun pseudobulk
braun = ad.read_h5ad(f'{ROOT}/data/raw/braun_2023/braun_all.h5ad')[:, vn]
braun.obs['CellClass'] = braun.obs['CellClass'].astype(str)
import scipy.sparse as sp
bX = braun.X.tocsr() if sp.issparse(braun.X) else braun.X
bdf = pd.DataFrame(bX.toarray() if sp.issparse(bX) else np.asarray(bX),
                   columns=vn, index=braun.obs_names)
bdf['__c'] = braun.obs['CellClass'].values
b_pb = pseudobulk(bdf.groupby('__c').mean())
log(f"Braun pseudobulk: {b_pb.shape}")

# ---- organoid pseudobulk (subsample)
tr = ad.read_h5ad(f'{ROOT}/data/braun_transfer_full_knn.h5ad', backed='r')
cc = tr.obs['CellClass_pred']
atlas = ad.read_h5ad(f'{ROOT}/data/atlas_v5_full.h5ad', backed='r')
ens_to_col = {}
for c, e in enumerate(sym2ens.get(s, '') for s in atlas.var_names):
    if e and e not in ens_to_col:
        ens_to_col[e] = c
cols = [ens_to_col[g] for g in vn]               # all 2006 present (verified earlier)
rng = np.random.default_rng(0)
idx = np.sort(rng.choice(atlas.n_obs, N_SUB, replace=False))
sub = atlas[idx].to_memory()
oX = sub.X[:, cols]
odf = pd.DataFrame(oX.toarray() if sp.issparse(oX) else np.asarray(oX),
                   columns=vn, index=sub.obs_names)
odf['__c'] = cc.reindex(sub.obs_names).values
o_pb = pseudobulk(odf.groupby('__c').mean())
log(f"organoid pseudobulk: {o_pb.shape}")

# ---- correlation matrix (organoid rows x Braun cols)
shared_cls = [c for c in o_pb.index if c in b_pb.index]
corr = pd.DataFrame(index=o_pb.index, columns=b_pb.index, dtype=float)
for oc in o_pb.index:
    for bc in b_pb.index:
        corr.loc[oc, bc] = np.corrcoef(o_pb.loc[oc], b_pb.loc[bc])[0, 1]
corr = corr.astype(float)
corr.to_csv(f'{ROOT}/data/q2_correspondence_corr.tsv', sep='\t')

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
