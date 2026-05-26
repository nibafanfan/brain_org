#!/usr/bin/env python3
"""Benchmark Q1: do multi-lineage protocols show broader cell-type / region
coverage than single-lineage?

Per deposit (dataset_slug), count how many distinct Braun CellClasses and
Regions are 'present' (>= MIN_FRAC of the deposit AND >= MIN_CELLS, among
confident calls). Compare multi- vs single-lineage deposits, and report which
specific lineages multi-lineage protocols add. Repertoire-based -> robust to the
residual batch structure in the embedding.
"""
import numpy as np, pandas as pd, anndata as ad
from scipy.stats import mannwhitneyu

ROOT = '/Users/eg/brain_organoid'
MIN_FRAC = 0.01
MIN_CELLS = 20
CONF = 0.5
def log(m): print(m, flush=True)

o = ad.read_h5ad(f'{ROOT}/data/braun_transfer_full_knn.h5ad', backed='r').obs.copy()
o['ml'] = o['multi_lineage'].astype(str).map(
    {'1': 'multi', 'True': 'multi', '0': 'single', 'False': 'single', 'No': 'single'})
dep_ml = o.groupby('dataset_slug', observed=True)['ml'].agg(lambda s: s.mode().iloc[0])
log(f"deposits: {dep_ml.value_counts().to_dict()} (total {dep_ml.size})")

conf_mask = o['CellClass_conf'].to_numpy() >= CONF
oc = o[conf_mask]
log(f"confident cells (conf>={CONF}): {len(oc):,}/{len(o):,}")

def coverage(label_col):
    # present types per deposit: >=MIN_FRAC of deposit AND >=MIN_CELLS
    ct = pd.crosstab(oc['dataset_slug'], oc[label_col])
    frac = ct.div(ct.sum(1), axis=0)
    present = (frac >= MIN_FRAC) & (ct >= MIN_CELLS)
    n_per_dep = present.sum(1)
    return present, n_per_dep

for label in ['CellClass_pred', 'Region_pred']:
    present, n = coverage(label)
    n = n.reindex(dep_ml.index)
    m, s = n[dep_ml == 'multi'], n[dep_ml == 'single']
    u, p = mannwhitneyu(m, s, alternative='greater')
    log(f"\n=== {label}: distinct types covered per deposit ===")
    log(f"  multi  (n={len(m)}): mean {m.mean():.2f}  median {m.median():.0f}")
    log(f"  single (n={len(s)}): mean {s.mean():.2f}  median {s.median():.0f}")
    log(f"  Mann-Whitney U (multi>single): p={p:.4f}")
    # which specific types are enriched in multi-lineage deposits?
    pm = present.reindex(dep_ml.index)
    rate_m = pm[dep_ml == 'multi'].mean()
    rate_s = pm[dep_ml == 'single'].mean()
    tab = pd.DataFrame({'multi_%': (rate_m*100).round(0),
                        'single_%': (rate_s*100).round(0)})
    tab['diff'] = (tab['multi_%'] - tab['single_%']).round(0)
    log("  per-type deposit-presence rate (% of deposits containing it):")
    log("\n" + tab.sort_values('diff', ascending=False).to_string())

# confound check: deposit size
sz = oc.groupby('dataset_slug', observed=True).size().reindex(dep_ml.index)
log(f"\nconfound check - median cells/deposit: multi {sz[dep_ml=='multi'].median():,.0f} "
    f"| single {sz[dep_ml=='single'].median():,.0f}")
log("(if multi deposits are much larger, some coverage gain is detection power, not biology)")
log("DONE")
