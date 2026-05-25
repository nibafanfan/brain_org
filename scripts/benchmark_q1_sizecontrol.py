#!/usr/bin/env python3
"""Q1 robustness: re-check the multi>single rare-lineage enrichment with EQUAL
detection power. Subsample every deposit to FLOOR cells (drop smaller ones), so
presence rates can't be driven by multi-lineage deposits being larger.
"""
import numpy as np, pandas as pd, anndata as ad

ROOT = '/Users/eg/brain_organoid'
FLOOR = 8000          # cells per deposit after downsampling
MIN_FRAC, MIN_CELLS, CONF = 0.01, 20, 0.5
rng = np.random.default_rng(0)
def log(m): print(m, flush=True)

o = ad.read_h5ad(f'{ROOT}/data/braun_transfer_full_knn.h5ad', backed='r').obs.copy()
o['ml'] = o['multi_lineage'].astype(str).map(
    {'1': 'multi', 'True': 'multi', '0': 'single', 'False': 'single', 'No': 'single'})
o = o[o['CellClass_conf'].to_numpy() >= CONF]
dep_ml = o.groupby('dataset_slug', observed=True)['ml'].agg(lambda s: s.mode().iloc[0])

# keep deposits with >= FLOOR confident cells, downsample each to exactly FLOOR
sizes = o.groupby('dataset_slug', observed=True).size()
keep = sizes[sizes >= FLOOR].index
idx = []
for d in keep:
    rows = np.where((o['dataset_slug'] == d).to_numpy())[0]
    idx.append(rng.choice(rows, FLOOR, replace=False))
od = o.iloc[np.concatenate(idx)]
dml = dep_ml.reindex(keep)
log(f"deposits kept (>= {FLOOR} conf cells): {dml.value_counts().to_dict()} "
    f"(of {dep_ml.value_counts().to_dict()})")

ct = pd.crosstab(od['dataset_slug'], od['CellClass_pred'])
frac = ct.div(ct.sum(1), axis=0)
present = ((frac >= MIN_FRAC) & (ct >= MIN_CELLS)).reindex(dml.index)
n = present.sum(1)
m, s = n[dml == 'multi'], n[dml == 'single']
from scipy.stats import mannwhitneyu
log(f"\ndistinct CellClasses/deposit (equal power): multi {m.mean():.2f} | single {s.mean():.2f} "
    f"| MWU p={mannwhitneyu(m, s, alternative='greater').pvalue:.3f}")

rate_m = present[dml == 'multi'].mean() * 100
rate_s = present[dml == 'single'].mean() * 100
tab = pd.DataFrame({'multi_%': rate_m.round(0), 'single_%': rate_s.round(0)})
tab['diff'] = (tab['multi_%'] - tab['single_%']).round(0)
log("per-type deposit-presence (equal power):")
log("\n" + tab.sort_values('diff', ascending=False).to_string())
log("DONE")
