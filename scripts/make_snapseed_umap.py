#!/usr/bin/env python3
"""HNOCA-true snapseed annotation UMAP: cluster the scVI latent, score each
snapseed cell-class marker panel PER CLUSTER (z-scored panel mean), assign each
cluster its winning class, and color the UMAP by that marker-derived label.
Cross-checks the snapseed labels against the Braun-transfer CellClass_cal.

Layout is unchanged (scVI latent); this only adds a marker-driven coloring +
independent annotation. No retraining.
"""
import sys, time
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad, scanpy as sc, yaml
import scipy.sparse as sp
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _provenance import write_sidecar

ROOT = Path('/Users/eg/brain_organoid')
# committed marker panel (docs/) preferred so the repo is self-contained; the
# identical data/reference copy is gitignored under data/*.
YAML = ROOT / 'docs/snapseed_markers_hnocadownload.yaml'
if not YAML.exists():
    YAML = ROOT / 'data/reference/hnoca_snapseed_markers.yaml'
N_SUB, RES = 300_000, 2.0
t0 = time.time()
def log(m): print(f"[{time.time()-t0:7.1f}s] {m}", flush=True)

# --- parse snapseed panels: {cell_class: [genes]}
def leaves(node):
    g = set()
    if isinstance(node, dict):
        for k, v in node.items():
            if k in ('marker_genes', 'markers', 'positive', 'genes') and isinstance(v, list):
                g |= {x for x in v if isinstance(x, str)}
            else:
                g |= leaves(v)
    elif isinstance(node, list):
        for x in node:
            g |= ({x} if isinstance(x, str) else leaves(x))
    return g
panels = {cls: sorted(leaves(node)) for cls, node in yaml.safe_load(open(YAML)).items()}
allg = sorted(set().union(*panels.values()))
log(f"snapseed: {len(panels)} classes, {len(allg)} marker genes")

# --- latent subsample + Braun labels
lat = ad.read_h5ad(ROOT / 'data/scvi_latent_v5_full.h5ad')
cal = ad.read_h5ad(ROOT / 'data/braun_transfer_full_calibrated.h5ad', backed='r')
rng = np.random.default_rng(0)
idx = np.sort(rng.choice(lat.n_obs, N_SUB, replace=False))
a = lat[idx].copy()
a.obsm['X_scvi'] = np.asarray(a.X)
a.obs['CellClass_cal'] = cal.obs['CellClass_cal'].reindex(a.obs_names).values

# --- marker expression (full atlas), normalized by library size, z-scored
atlas = ad.read_h5ad(ROOT / 'data/atlas_v5_full.h5ad', backed='r')
pos = pd.Series(np.arange(atlas.n_obs), index=atlas.obs_names)
mk = [g for g in allg if g in set(atlas.var_names)]
ge = atlas[pos.loc[a.obs_names].to_numpy()].to_memory()[:, mk]
M = ge.X.toarray() if sp.issparse(ge.X) else np.asarray(ge.X)
ncounts = a.obs['n_counts'].to_numpy()[:, None] if 'n_counts' in a.obs else M.sum(1, keepdims=True)
M = np.log1p(M / np.clip(ncounts, 1, None) * 1e4)
Z = (M - M.mean(0)) / (M.std(0) + 1e-9)
zdf = pd.DataFrame(Z, columns=mk)
log(f"markers z-scored: {len(mk)}/{len(allg)} present")

# --- cluster + per-cluster panel scoring -> snapseed label
sc.pp.neighbors(a, use_rep='X_scvi', n_neighbors=15, random_state=0)
sc.tl.leiden(a, resolution=RES, flavor='igraph', n_iterations=2, directed=False, random_state=0)
sc.tl.umap(a, random_state=0)
log(f"Leiden {a.obs['leiden'].nunique()} clusters + UMAP done")

cell_panel = pd.DataFrame(
    {cls: zdf[[g for g in genes if g in mk]].mean(1) for cls, genes in panels.items()})
cell_panel['leiden'] = a.obs['leiden'].to_numpy()
clust_score = cell_panel.groupby('leiden', observed=True).mean()
clust_label = clust_score.idxmax(1)
a.obs['snapseed'] = a.obs['leiden'].map(clust_label).astype('category')
# Braun label assigned by MAJORITY VOTE per Leiden cluster -> tiles the map like snapseed
braun_majority = a.obs.groupby('leiden', observed=True)['CellClass_cal'].agg(
    lambda s: s.mode().iloc[0] if len(s.mode()) else 'Unknown')
a.obs['Braun_clustervote'] = a.obs['leiden'].map(braun_majority).astype('category')
log("snapseed + Braun-majority per-cluster labels assigned")

# --- figure: snapseed (cluster) | Braun per-cell (soft) | Braun majority (cluster)
fig, ax = plt.subplots(1, 3, figsize=(30, 8))
sc.pl.umap(a, color='snapseed', ax=ax[0], show=False, size=3,
           title='snapseed (marker-derived, per-cluster)', legend_loc='right margin')
sc.pl.umap(a, color='CellClass_cal', ax=ax[1], show=False, size=3,
           title='Braun CellClass (per-cell kNN, soft)', legend_loc='right margin')
sc.pl.umap(a, color='Braun_clustervote', ax=ax[2], show=False, size=3,
           title='Braun CellClass (majority vote per cluster)', legend_loc='right margin')
fig.suptitle('scVI-latent UMAP — same layout; cluster-tiled vs per-cell labels', fontsize=14)
fig.tight_layout(); fig.savefig(ROOT / 'data/scvi_umap_snapseed.png', dpi=110)
log("wrote data/scvi_umap_snapseed.png")

# cluster-level cross-check: snapseed vs Braun-majority (both one-label-per-cluster)
cl = pd.DataFrame({'snapseed': clust_label, 'braun_majority': braun_majority})
log("=== per-cluster snapseed vs Braun-majority (one row per Leiden cluster) ===")
log("\n" + cl.sort_values('braun_majority').to_string())

# --- cross-check crosstab (row-normalized %)
ct = pd.crosstab(a.obs['snapseed'], a.obs['CellClass_cal'], normalize='index') * 100
out_tsv = ROOT / 'data/snapseed_vs_braun_crosstab.tsv'
ct.round(1).to_csv(out_tsv, sep='\t')
write_sidecar(out_tsv, __file__, {'N_SUB': N_SUB, 'RES': RES, 'n_markers': len(mk),
                                  'n_clusters': int(a.obs['leiden'].nunique())})
log("=== snapseed label -> top Braun CellClass (row %, cross-check) ===")
for s in ct.index:
    top = ct.loc[s].sort_values(ascending=False)
    log(f"  {s:28} -> {top.index[0]} {top.iloc[0]:.0f}%  (next {top.index[1]} {top.iloc[1]:.0f}%)")
log(f"saved crosstab -> {out_tsv}")
log("DONE")
