#!/usr/bin/env python3
"""Single benchmark rerun on CALIBRATED labels (CellClass_cal) with deltas vs the
old argmax labels (CellClass_pred). Recomputes Q1 (coverage) and Q2 (pseudobulk
correspondence) for BOTH label sets in one pass. Q3 aggregates by Braun labels on
the organoid side, so it is unaffected by organoid-label calibration (noted only).
Integration quality is now reported by scIB (data/scib_metrics.tsv), not here.
"""
import argparse, sys, time
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad
import scipy.sparse as sp
from scipy.stats import mannwhitneyu
sys.path.insert(0, str(Path(__file__).resolve().parent))
from atlas_common import load_config, sym2ens, model_genes, read_atlas_genes, cp10k_log

ap = argparse.ArgumentParser(); ap.add_argument('--root', default=None); args = ap.parse_args()
cfg = load_config(root=args.root)
ROOT = cfg.root
N_SUB = cfg.defaults['n_sub']
MIN_FRAC, MIN_CELLS = 0.01, 20
NULLISH = {'Unknown', 'nan', 'none', ''}
t0 = time.time()
def log(m): print(f"[{time.time()-t0:7.1f}s] {m}", flush=True)

cal = ad.read_h5ad(cfg.calibrated)
o = cal.obs
o['ml'] = o['multi_lineage'].astype(str).map(
    {'1': 'multi', 'True': 'multi', '0': 'single', 'False': 'single', 'No': 'single'})
dep_ml = o.groupby('dataset_slug', observed=True)['ml'].agg(lambda s: s.mode().iloc[0])

# ---------- Q1: per-deposit presence, multi vs single, old vs new label ----------
def q1_diff(label):
    sub = o[~o[label].astype(str).isin(NULLISH)]
    ct = pd.crosstab(sub['dataset_slug'], sub[label])
    present = (ct.div(ct.sum(1), axis=0) >= MIN_FRAC) & (ct >= MIN_CELLS)
    present = present.reindex(dep_ml.index).fillna(False)
    rate_m = present[dep_ml == 'multi'].mean() * 100
    rate_s = present[dep_ml == 'single'].mean() * 100
    return (rate_m - rate_s).round(0)

q1 = pd.DataFrame({'old (pred)': q1_diff('CellClass_pred'),
                   'new (cal)': q1_diff('CellClass_cal')})
q1['delta'] = (q1['new (cal)'] - q1['old (pred)']).round(0)
log("=== Q1: multi-minus-single deposit-presence (pp), old vs calibrated ===")
log("\n" + q1.sort_values('new (cal)', ascending=False).to_string())

# ---------- Q2: pseudobulk correspondence vs Braun, old vs new label ----------
vn = model_genes(cfg.braun_scanvi_model)
s2e = sym2ens(cfg.canonical)
pb = cp10k_log

braun = ad.read_h5ad(cfg.braun)[:, vn]
braun.obs['CellClass'] = braun.obs['CellClass'].astype(str)
bX = braun.X.toarray() if sp.issparse(braun.X) else np.asarray(braun.X)
bdf = pd.DataFrame(bX, columns=vn); bdf['__c'] = braun.obs['CellClass'].values
b_pb = pb(bdf.groupby('__c').mean())
log(f"Braun pseudobulk {b_pb.shape}")

atlas = ad.read_h5ad(cfg.atlas_full, backed='r')
rng = np.random.default_rng(0)
idx = np.sort(rng.choice(atlas.n_obs, N_SUB, replace=False))
oX, present = read_atlas_genes(atlas, vn, s2e, row_idx=idx)
obs_names = atlas.obs_names[idx]
odf = pd.DataFrame(oX.toarray() if sp.issparse(oX) else np.asarray(oX),
                   columns=present, index=obs_names)

def self_corr(label):
    lab = cal.obs[label].reindex(obs_names)
    m = ~lab.astype(str).isin(NULLISH)
    d = odf[m.values].copy(); d['__c'] = lab[m.values].values
    o_pb = pb(d.groupby('__c').mean())
    res = {}
    diag_ok = 0; shared = [c for c in o_pb.index if c in b_pb.index]
    for oc in shared:
        cors = {bc: np.corrcoef(o_pb.loc[oc], b_pb.loc[bc])[0, 1] for bc in b_pb.index}
        best = max(cors, key=cors.get)
        res[oc] = (round(cors[oc], 3), best == oc)
        diag_ok += best == oc
    return res, diag_ok, len(shared)

old_r, old_ok, old_n = self_corr('CellClass_pred')
new_r, new_ok, new_n = self_corr('CellClass_cal')
rows = []
for c in sorted(set(old_r) | set(new_r)):
    rows.append((c, old_r.get(c, ('-',))[0], new_r.get(c, ('-',))[0]))
q2 = pd.DataFrame(rows, columns=['CellClass', 'old self-corr', 'new self-corr'])
log("=== Q2: organoid-vs-Braun self-correspondence, old vs calibrated ===")
log("\n" + q2.to_string(index=False))
log(f"diagonal dominance: old {old_ok}/{old_n}, new {new_ok}/{new_n}")
log("\nNOTE: Q3 (reference-coverage) aggregates by BRAUN labels on the organoid "
    "side; organoid-label calibration does not change it. Integration quality: "
    "see scIB metrics (data/scib_metrics.tsv) — iLISI 0.015, cLISI 0.962.")
log("DONE")
