#!/usr/bin/env python3
"""Review item #4: cluster-level marker gating (replaces fragile per-cell score_genes).

Leiden-cluster the query on the scVI latent, score rare-lineage marker panels PER
CLUSTER (cluster means are ambient-robust), gate whole clusters, and report
prevalence + transferred-label agreement (CellClass_cal = Braun cross-ref) + mean
confidence + OOD fraction. OOD-stratified throughout (78% of cells are OOD vs Braun).
"""
import sys, time
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad, scanpy as sc
import scipy.sparse as sp
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _provenance import write_sidecar

ROOT = '/Users/eg/brain_organoid'
N_SUB = 400_000
RES = 2.0
CLUSTER_GATE = 0.30          # cluster-mean score_genes to flag a lineage cluster
t0 = time.time()
def log(m): print(f"[{time.time()-t0:7.1f}s] {m}", flush=True)
rng = np.random.default_rng(0)

PANELS = {
    'microglia':   ['PTPRC', 'P2RY12', 'CX3CR1', 'AIF1', 'C1QA', 'C1QB', 'C1QC',
                    'CSF1R', 'TMEM119', 'TYROBP', 'FCER1G', 'CD68'],
    'endothelium': ['CLDN5', 'PECAM1', 'CDH5', 'VWF', 'EGFL7'],
    'oligo':       ['MBP', 'PLP1', 'MOG'],
}

cal = ad.read_h5ad(f'{ROOT}/data/braun_transfer_full_calibrated.h5ad')   # X=latent
atlas = ad.read_h5ad(f'{ROOT}/data/atlas_v5_full.h5ad', backed='r')
assert list(cal.obs_names[:3]) == list(atlas.obs_names[:3]) and cal.n_obs == atlas.n_obs, \
    "calibrated and atlas row order must match (positional join)"
idx = np.sort(rng.choice(atlas.n_obs, N_SUB, replace=False))

asub = atlas[idx].to_memory()                      # expression (symbols) for scoring
a = ad.AnnData(X=asub.X.astype('float32'), obs=cal.obs.iloc[idx].copy())
a.var_names = atlas.var_names
a.obsm['X_scvi'] = np.asarray(cal.X)[idx]
sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
for lin, genes in PANELS.items():
    g = [x for x in genes if x in a.var_names]
    sc.tl.score_genes(a, g, score_name=f's_{lin}')
log(f"scored panels on {a.n_obs:,} cells")

sc.pp.neighbors(a, use_rep='X_scvi', n_neighbors=15)
sc.tl.leiden(a, resolution=RES, flavor='igraph', n_iterations=2, directed=False)
nC = a.obs['leiden'].nunique()
log(f"Leiden: {nC} clusters (res={RES})")

a.obs['ood'] = a.obs['ood'].astype(bool)
g = a.obs.groupby('leiden', observed=True)
tab = pd.DataFrame({
    'n': g.size(),
    'pct': (g.size() / a.n_obs * 100).round(2),
    's_microglia': g['s_microglia'].mean().round(3),
    's_endothelium': g['s_endothelium'].mean().round(3),
    's_oligo': g['s_oligo'].mean().round(3),
    'dom_CellClass_cal': g['CellClass_cal'].agg(lambda s: s.mode().iloc[0]),
    'dom_frac': g['CellClass_cal'].agg(lambda s: s.value_counts(normalize=True).iloc[0]).round(2),
    'mean_conf': g['CellClass_cal_conf'].mean().round(3),
    'pct_ood': (g['ood'].mean() * 100).round(0),
    'pct_multi': g['multi_lineage'].apply(
        lambda s: s.astype(str).isin(['1', 'True']).mean() * 100).round(0),
})

log("=== rare-lineage CLUSTERS (cluster-mean marker score > %.2f) ===" % CLUSTER_GATE)
for lin in PANELS:
    hit = tab[tab[f's_{lin}'] > CLUSTER_GATE].sort_values(f's_{lin}', ascending=False)
    if len(hit) == 0:
        log(f"  {lin}: no cluster above gate")
        continue
    tot = int(hit['n'].sum())
    log(f"  {lin}: {len(hit)} cluster(s), {tot} cells ({tot/a.n_obs*100:.2f}%)")
    log("\n" + hit[[f's_{lin}', 'n', 'pct', 'dom_CellClass_cal', 'dom_frac',
                    'mean_conf', 'pct_ood', 'pct_multi']].to_string())

# OOD-stratified rare-lineage prevalence (cluster-gated)
log("\n=== OOD-stratified rare-lineage prevalence (cluster-gated cells) ===")
gated = {lin: tab.index[tab[f's_{lin}'] > CLUSTER_GATE] for lin in PANELS}
for lin, cl in gated.items():
    mask = a.obs['leiden'].isin(cl).to_numpy()
    if mask.sum() == 0:
        continue
    ood = a.obs['ood'].to_numpy()
    log(f"  {lin}: in-dist {int((mask & ~ood).sum())} | OOD {int((mask & ood).sum())} "
        f"(OOD frac {(mask & ood).sum()/mask.sum()*100:.0f}%)")

out_tsv = f'{ROOT}/data/q2_clustergate.tsv'
tab.sort_values('pct', ascending=False).to_csv(out_tsv, sep='\t')
write_sidecar(out_tsv, __file__, {'N_SUB': N_SUB, 'RES': RES, 'CLUSTER_GATE': CLUSTER_GATE,
                                  'n_clusters': int(nC)})
log(f"saved -> {out_tsv} (+ provenance)")
log("DONE")
