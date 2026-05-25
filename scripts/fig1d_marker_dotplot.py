#!/usr/bin/env python3
"""HNOCA Fig.1d-style marker dot plot: the snapseed curated marker panel (the SAME
genes HNOCA used, grouped by cell class) x our transferred CellClass — dot size =
fraction of cells expressing, color = mean expression. Validates the annotation and
matches HNOCA's marker selection (vs the ad-hoc 7-gene UMAP overlay set).

  python scripts/fig1d_marker_dotplot.py [--n-sub 400000]
"""
import argparse, sys, time
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad, scanpy as sc, yaml
import scipy.sparse as sp
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parent))
from atlas_common import load_config

# snapseed cell-class -> Braun CellClass (for ordering the dotplot rows sensibly)
CLASS_ORDER = ['Radial glia', 'Neuronal IPC', 'Neuroblast', 'Neuron', 'Glioblast', 'Oligo',
               'Immune', 'Vascular', 'Fibroblast', 'Neural crest', 'Placodes', 'Erythrocyte']
t0 = time.time()
def log(m): print(f"[{time.time()-t0:7.1f}s] {m}", flush=True)

ap = argparse.ArgumentParser()
ap.add_argument('--root', default=None)
ap.add_argument('--n-sub', type=int, default=400_000)
ap.add_argument('--out-tag', default='v5')
args = ap.parse_args()
cfg = load_config(root=args.root)
OUT = cfg.root / 'data/poster'; OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(cfg.defaults['seed'])

# snapseed panels: {cell_class: [genes]}
yml = cfg.snapseed_yaml if cfg.snapseed_yaml.exists() else cfg.root / 'data/reference/hnoca_snapseed_markers.yaml'
def leaves(node):
    g = set()
    if isinstance(node, dict):
        for k, v in node.items():
            g |= ({x for x in v if isinstance(x, str)} if k in ('marker_genes', 'markers', 'positive', 'genes')
                  and isinstance(v, list) else leaves(v))
    elif isinstance(node, list):
        for x in node:
            g |= ({x} if isinstance(x, str) else leaves(x))
    return g
panels = {cls: sorted(leaves(node)) for cls, node in yaml.safe_load(open(yml)).items()}

# read snapseed marker expression (symbols) for a subsample + join CellClass_cal
lat = ad.read_h5ad(cfg.latent)
cal = ad.read_h5ad(cfg.calibrated, backed='r')
atlas = ad.read_h5ad(cfg.atlas_full, backed='r')
allg = sorted(set().union(*panels.values()))
mk = [g for g in allg if g in set(atlas.var_names)]
panels = {c: [g for g in gs if g in mk] for c, gs in panels.items()}
panels = {c: gs for c, gs in panels.items() if gs}
log(f"snapseed: {len(panels)} classes, {len(mk)} marker genes present")

idx = np.sort(rng.choice(lat.n_obs, min(args.n_sub, lat.n_obs), replace=False))
names = lat.obs_names[idx]
Xm = atlas[idx].to_memory()[:, mk].X
Xm = Xm.toarray() if sp.issparse(Xm) else np.asarray(Xm)
a = ad.AnnData(X=Xm.astype('float32')); a.var_names = mk
a.obs_names = names
a.obs['CellClass'] = cal.obs['CellClass_cal'].reindex(names).astype(str).values
a = a[a.obs['CellClass'] != 'Unknown'].copy()
sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
order = [c for c in CLASS_ORDER if c in set(a.obs['CellClass'])]
a.obs['CellClass'] = pd.Categorical(a.obs['CellClass'], categories=order, ordered=True)
log(f"dotplot input {a.shape}; classes {order}")

dp = sc.pl.dotplot(a, panels, groupby='CellClass', standard_scale='var',
                   cmap='viridis', show=False, return_fig=True,
                   title='snapseed marker panel x CellClass (HNOCA Fig1d-style)')
dp.savefig(OUT / f'fig1d_marker_dotplot_{args.out_tag}.png', dpi=150, bbox_inches='tight')
# matrixplot variant (mean expression heatmap) as a companion
mp = sc.pl.matrixplot(a, panels, groupby='CellClass', standard_scale='var',
                      cmap='viridis', show=False, return_fig=True)
mp.savefig(OUT / f'fig1d_marker_matrix_{args.out_tag}.png', dpi=150, bbox_inches='tight')
log(f"wrote fig1d_marker_dotplot_{args.out_tag}.png + matrix -> {OUT}")
log("DONE")
