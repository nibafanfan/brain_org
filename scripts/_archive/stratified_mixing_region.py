#!/usr/bin/env python3
"""Finer stratification: CellClass x Region. Separates true batch effect from
subtype/regional biology left over in the coarse 12-class CellClass.

If e.g. (Neuron, Telencephalon) cells from different datasets STILL don't mix
(ratio stays ~10-15x), the residual is genuine batch effect. If the ratio drops
a lot vs Neuron-overall (15.8x), most of that was regional/subtype heterogeneity.
"""
import argparse, sys, time
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad
sys.path.insert(0, str(Path(__file__).resolve().parent))
from atlas_common import load_config, same_frac, baseline

ap = argparse.ArgumentParser()
ap.add_argument('--root', default=None)
args = ap.parse_args()
cfg = load_config(root=args.root)
ROOT = cfg.root
N_SUB = cfg.defaults['n_sub']
CAP = 40_000
MIN_CELLS = 1500
MIN_DATASETS = 10        # mixing only assessable if the group spans many datasets
K = cfg.defaults['knn_k']
BATCH = 'dataset_slug'
t0 = time.time()
def log(m): print(f"[{time.time()-t0:7.1f}s] {m}", flush=True)

lat = ad.read_h5ad(cfg.latent)
tr = ad.read_h5ad(ROOT / 'data/braun_transfer_full_knn.h5ad', backed='r')
trobs = tr.obs[['CellClass_pred', 'Region_pred']].reindex(lat.obs_names)
lat.obs['CellClass'] = trobs['CellClass_pred'].values
lat.obs['Region'] = trobs['Region_pred'].values
lat = lat[lat.obs['CellClass'].notna()].copy()
log(f"joined CellClass+Region; {lat.n_obs:,} labeled")

rng = np.random.default_rng(0)
if lat.n_obs > N_SUB:
    lat = lat[np.sort(rng.choice(lat.n_obs, N_SUB, replace=False))].copy()
X = np.asarray(lat.X)
ds = lat.obs[BATCH].astype('category')
dscodes = ds.cat.codes.to_numpy()
cc = lat.obs['CellClass'].to_numpy()
rg = lat.obs['Region'].to_numpy()
log(f"subsample {X.shape}")

# CellClass-only ratios (reference, from prior run)
cc_only = {'Radial glia': 15.5, 'Neuron': 15.8, 'Neuroblast': 12.9,
           'Glioblast': 9.6, 'Neuronal IPC': 8.7, 'Fibroblast': 7.0}

rows = []
pairs = pd.Series(list(zip(cc, rg))).value_counts()
for (c, r), n in pairs.items():
    if n < MIN_CELLS:
        continue
    mask = (cc == c) & (rg == r)
    codes = dscodes[mask]
    nds = np.unique(codes).size
    if nds < MIN_DATASETS:
        continue
    Xc = X[mask]
    if n > CAP:
        sel = rng.choice(n, CAP, replace=False)
        Xc, codes = Xc[sel], codes[sel]
    sf, bl = same_frac(Xc, codes), baseline(codes)
    rows.append((c, r, int(n), nds, round(sf, 3), round(bl, 3),
                 round(sf / bl, 1), cc_only.get(c, np.nan)))
    log(f"  {c:14}x {r:14} n={int(n):>6} nds={nds:>3} ratio={sf/bl:5.1f}x "
        f"(CellClass-only {cc_only.get(c,'?')}x)")

df = pd.DataFrame(rows, columns=['CellClass', 'Region', 'n', 'n_datasets',
                                 'same_frac', 'baseline', 'ratio', 'cellclass_only_ratio'])
df = df.sort_values('n', ascending=False)
log("\n=== CellClass x Region mixing (ratio vs CellClass-only) ===")
log("\n" + df.to_string(index=False))
df.to_csv(ROOT / 'data/stratified_mixing_region.tsv', sep='\t', index=False)
log("saved -> data/stratified_mixing_region.tsv")
log("DONE")
