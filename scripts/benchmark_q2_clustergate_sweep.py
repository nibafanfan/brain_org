#!/usr/bin/env python3
"""Review #4 robustness sweep: resolutions {1.0,1.5,2.0} x gate thresholds
{0.25,0.30,0.35} x 2 subsample seeds, with per-cluster Braun centroid-correlation
(top class + margin). Confirms the cluster-gating conclusions are stable:
microglia surfaces (centroid->Immune), endothelium doesn't, "oligo"-marker cluster
maps to Neural crest. Long-format TSV + compact summary.
"""
import argparse, sys, time
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad, scanpy as sc
import scipy.sparse as sp
sys.path.insert(0, str(Path(__file__).resolve().parent))
from atlas_common import load_config, sym2ens, model_genes, ens_to_col, cp10k_log, write_sidecar

ap = argparse.ArgumentParser(); ap.add_argument('--root', default=None); args = ap.parse_args()
cfg = load_config(root=args.root)
ROOT = cfg.root
N_SUB = 400_000
RES = [1.0, 1.5, 2.0]
THR = [0.25, 0.30, 0.35]
SEEDS = [0, 1]
t0 = time.time()
def log(m): print(f"[{time.time()-t0:7.1f}s] {m}", flush=True)

PANELS = {
    'microglia':   ['PTPRC', 'P2RY12', 'CX3CR1', 'AIF1', 'C1QA', 'C1QB', 'C1QC',
                    'CSF1R', 'TMEM119', 'TYROBP', 'FCER1G', 'CD68'],
    'endothelium': ['CLDN5', 'PECAM1', 'CDH5', 'VWF', 'EGFL7'],
    'oligo':       ['MBP', 'PLP1', 'MOG'],
}

vn = model_genes(cfg.braun_scanvi_model)
s2e = sym2ens(cfg.canonical)

# Braun pseudobulk (2006 genes) for centroid-correlation
braun = ad.read_h5ad(cfg.braun)[:, vn]
bX = braun.X.toarray() if sp.issparse(braun.X) else np.asarray(braun.X)
b_pb = cp10k_log(pd.DataFrame(bX, columns=vn).assign(__c=braun.obs['CellClass'].astype(str).values)
                 .groupby('__c').mean())
log(f"Braun pseudobulk {b_pb.shape}")

cal = ad.read_h5ad(cfg.calibrated)
atlas = ad.read_h5ad(cfg.atlas_full, backed='r')
e2c = ens_to_col(atlas.var_names, s2e)
cols2006 = [e2c[g] for g in vn]

def centroid_corr(pb_vec):
    cors = {bc: np.corrcoef(pb_vec, b_pb.loc[bc])[0, 1] for bc in b_pb.index}
    s = pd.Series(cors).sort_values(ascending=False)
    return s.index[0], round(s.iloc[0] - s.iloc[1], 3)      # top class, margin

rows = []
for seed in SEEDS:
    rng = np.random.default_rng(seed)
    idx = np.sort(rng.choice(atlas.n_obs, N_SUB, replace=False))
    asub = atlas[idx].to_memory()
    a = ad.AnnData(X=asub.X.astype('float32'), obs=cal.obs.iloc[idx].copy())
    a.var_names = atlas.var_names
    a.obsm['X_scvi'] = np.asarray(cal.X)[idx]
    a.layers['counts'] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    for lin, genes in PANELS.items():
        sc.tl.score_genes(a, [x for x in genes if x in a.var_names], score_name=f's_{lin}')
    sc.pp.neighbors(a, use_rep='X_scvi', n_neighbors=15, random_state=seed)
    counts2006 = a.layers['counts'][:, cols2006]
    counts2006 = counts2006.toarray() if sp.issparse(counts2006) else np.asarray(counts2006)
    for res in RES:
        sc.tl.leiden(a, resolution=res, flavor='igraph', n_iterations=2,
                     directed=False, random_state=seed, key_added='leiden')
        cl = a.obs['leiden'].to_numpy()
        for lin in PANELS:
            sc_lin = a.obs[f's_{lin}'].groupby(a.obs['leiden'], observed=True).mean()
            for thr in THR:
                hits = sc_lin.index[sc_lin > thr]
                mask = np.isin(cl, hits)
                n = int(mask.sum())
                # centroid-corr of the union of gated clusters
                if n > 0:
                    pbv = np.log1p(counts2006[mask].mean(0) / counts2006[mask].mean(0).sum() * 1e4)
                    top, margin = centroid_corr(pbv)
                else:
                    top, margin = '-', np.nan
                rows.append((seed, res, lin, thr, len(hits), n, round(n / a.n_obs * 100, 3),
                             top, margin))
        log(f"  seed={seed} res={res}: done")
    del a, asub, counts2006

df = pd.DataFrame(rows, columns=['seed', 'res', 'lineage', 'thr', 'n_clusters', 'n_cells',
                                 'pct', 'centroid_top_braun', 'centroid_margin'])
out_tsv = f'{ROOT}/data/q2_clustergate_sweep.tsv'
df.to_csv(out_tsv, sep='\t', index=False)
write_sidecar(out_tsv, __file__, {'N_SUB': N_SUB, 'RES': RES, 'THR': THR, 'SEEDS': SEEDS})
log("=== stability summary (thr=0.30): n_cells (pct) + centroid top->margin, per seed x res ===")
for lin in PANELS:
    log(f"\n{lin}:")
    sub = df[(df.lineage == lin) & (df.thr == 0.30)]
    for _, r in sub.iterrows():
        log(f"  seed{r.seed} res{r.res}: {r.n_cells:>5} ({r.pct:.2f}%) "
            f"clusters={r.n_clusters} centroid->{r.centroid_top_braun} (margin {r.centroid_margin})")
log(f"\nsaved -> {out_tsv}")
log("DONE")
