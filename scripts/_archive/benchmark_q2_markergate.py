#!/usr/bin/env python3
"""Q2 Part 2: find rare added lineages by canonical markers (independent of the
Braun transfer, which gives them ~0%), then check (a) how many cells, (b) which
protocols (multi vs single), (c) how many the transfer MISSED, (d) whether the
gated cells transcriptomically match the Braun primary counterpart.
"""
import argparse, sys, time
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad, scanpy as sc
import scipy.sparse as sp
sys.path.insert(0, str(Path(__file__).resolve().parent))
from atlas_common import load_config, sym2ens, model_genes, read_atlas_genes

ap = argparse.ArgumentParser(); ap.add_argument('--root', default=None); args = ap.parse_args()
cfg = load_config(root=args.root)
ROOT = cfg.root
N_SUB = 500_000
t0 = time.time()
def log(m): print(f"[{time.time()-t0:7.1f}s] {m}", flush=True)

# Strict, lineage-SPECIFIC markers. Dropped contaminants: A2M/FLT1/KDR (astrocyte/
# progenitor-expressed) from endothelium; SOX10/OLIG1/OLIG2 (ventral-progenitor /
# OPC, not mature oligo) from oligo.
MARKERS = {
    'microglia':   ['PTPRC', 'P2RY12', 'CX3CR1', 'AIF1', 'C1QA', 'C1QB', 'C1QC',
                    'CSF1R', 'TMEM119', 'TYROBP', 'FCER1G', 'CD68'],
    'endothelium': ['CLDN5', 'PECAM1', 'CDH5', 'VWF', 'EGFL7'],
    'oligo':       ['MBP', 'PLP1', 'MOG'],
}
LINEAGE_TO_BRAUN = {'microglia': 'Immune', 'endothelium': 'Vascular', 'oligo': 'Oligo'}

vn = model_genes(cfg.braun_scanvi_model)
s2e = sym2ens(cfg.canonical)
mk_ens = {lin: [s2e[s] for s in syms if s in s2e and s2e[s] in set(vn)]
          for lin, syms in MARKERS.items()}
for lin, g in mk_ens.items():
    log(f"{lin}: {len(g)}/{len(MARKERS[lin])} markers in 2006-gene panel")

# ---- organoid subsample, 2006 genes, normalized
tr = ad.read_h5ad(ROOT / 'data/braun_transfer_full_knn.h5ad', backed='r')
atlas = ad.read_h5ad(cfg.atlas_full, backed='r')
rng = np.random.default_rng(0)
idx = np.sort(rng.choice(atlas.n_obs, N_SUB, replace=False))
oX, present = read_atlas_genes(atlas, vn, s2e, row_idx=idx)
a = ad.AnnData(X=oX, obs=atlas.obs.iloc[idx].copy())
a.var_names = present
a.obs['CellClass_pred'] = tr.obs['CellClass_pred'].reindex(a.obs_names).values
a.obs['ml'] = a.obs['multi_lineage'].astype(str).map(
    {'1': 'multi', 'True': 'multi', '0': 'single', 'False': 'single', 'No': 'single'})
sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
log(f"organoid subsample ready: {a.shape}")

for lin, genes in mk_ens.items():
    sc.tl.score_genes(a, genes, score_name=f'score_{lin}')

# ---- Braun pseudobulk for the matched classes
braun = ad.read_h5ad(f'{ROOT}/data/raw/braun_2023/braun_all.h5ad')[:, vn]
braun.obs['CellClass'] = braun.obs['CellClass'].astype(str)
sc.pp.normalize_total(braun, target_sum=1e4); sc.pp.log1p(braun)
bX = braun.X.toarray() if sp.issparse(braun.X) else np.asarray(braun.X)
bdf = pd.DataFrame(bX, columns=vn); bdf['__c'] = braun.obs['CellClass'].values
b_pb = bdf.groupby('__c').mean()
log("Braun pseudobulk ready")

aX = a.X.toarray() if sp.issparse(a.X) else np.asarray(a.X)
adf = pd.DataFrame(aX, columns=vn, index=a.obs_names)

log("\n=== marker-gated rare lineages (organoid, 500k subsample) ===")
for lin in MARKERS:
    s = a.obs[f'score_{lin}'].to_numpy()
    for thr in (0.5, 1.0):
        gate = s > thr
        n = int(gate.sum())
        if n == 0:
            log(f"  {lin:12} thr>{thr}: 0 cells"); continue
        ml = a.obs['ml'][gate].value_counts()
        # of gated cells, how many did the transfer label as the matched Braun class?
        braun_cls = LINEAGE_TO_BRAUN[lin]
        hit = (a.obs['CellClass_pred'][gate] == braun_cls).mean()
        # transcriptomic correspondence: gated-cell pseudobulk vs Braun matched class
        # (both are means of log1p-CP10K over the 2006 genes -> comparable)
        r = np.corrcoef(adf[gate].mean(0), b_pb.loc[braun_cls])[0, 1]
        log(f"  {lin:12} thr>{thr}: n={n:>5} ({n/N_SUB*100:.2f}%) "
            f"| multi={ml.get('multi',0)} single={ml.get('single',0)} "
            f"| transfer-labeled-{braun_cls}={hit*100:.0f}% "
            f"| corr-to-Braun-{braun_cls}={r:.3f}")
log("DONE")
