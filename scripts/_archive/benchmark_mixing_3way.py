#!/usr/bin/env python3
"""Decisive mixing test: within (CellClass x Region x age-bin), do DATASETS mix?

Controls for the three biological axes that could explain non-mixing (cell type,
region, developmental age). If same-dataset kNN fraction stays >> the within-triplet
baseline, the residual structure is genuine TECHNICAL batch (integration incomplete).
If it drops toward baseline, the global non-mixing was biology (protocol composition /
maturation) and the integration is fine where it matters.
"""
import sys, time
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad
sys.path.insert(0, str(Path(__file__).resolve().parent))
from atlas_common import load_config, same_frac, baseline

N_SUB, CAP, MIN_CELLS, MIN_DATASETS = 300_000, 30_000, 800, 5
AGE_BINS = [0, 30, 60, 90, 120, 180, np.inf]
t0 = time.time()
def log(m): print(f"[{time.time()-t0:7.1f}s] {m}", flush=True)
cfg = load_config(); rng = np.random.default_rng(cfg.defaults['seed'])

lat = ad.read_h5ad(cfg.latent)
cal = ad.read_h5ad(cfg.calibrated, backed='r')
lat.obs['CellClass'] = cal.obs['CellClass_cal'].reindex(lat.obs_names).values
lat.obs['Region'] = cal.obs['Region_pred'].reindex(lat.obs_names).values
lat.obs['age'] = pd.to_numeric(lat.obs['age_days'], errors='coerce')
lat = lat[lat.obs['CellClass'].notna() & (lat.obs['CellClass'] != 'Unknown') & lat.obs['age'].notna()].copy()
if lat.n_obs > N_SUB:
    lat = lat[np.sort(rng.choice(lat.n_obs, N_SUB, replace=False))].copy()
X = np.asarray(lat.X)
cc = lat.obs['CellClass'].to_numpy(); rg = lat.obs['Region'].to_numpy()
ab = pd.cut(lat.obs['age'], AGE_BINS).astype(str).to_numpy()
ds = lat.obs['dataset_slug'].astype('category'); dscodes = ds.cat.codes.to_numpy()
log(f"subsample {X.shape}; {ds.cat.categories.size} datasets")

rows = []
trip = pd.Series(list(zip(cc, rg, ab))).value_counts()
for (c, r, a), n in trip.items():
    if n < MIN_CELLS:
        continue
    m = (cc == c) & (rg == r) & (ab == a)
    codes = dscodes[m]
    nds = np.unique(codes).size
    if nds < MIN_DATASETS:
        continue
    Xc = X[m]
    if n > CAP:
        sel = rng.choice(n, CAP, replace=False); Xc, codes = Xc[sel], codes[sel]
    sf, bl = same_frac(Xc, codes), baseline(codes)
    rows.append((c, r, a, int(n), nds, round(sf, 3), round(bl, 3), round(sf / bl, 1)))

df = pd.DataFrame(rows, columns=['CellClass', 'Region', 'age_bin', 'n', 'n_datasets',
                                 'same_frac', 'baseline', 'ratio']).sort_values('n', ascending=False)
df.to_csv(cfg.root / 'data/mixing_celltype_age_region.tsv', sep='\t', index=False)
log(f"=== datasets mixing WITHIN (CellClass x Region x age-bin) : {len(df)} assessable triplets ===")
log("\n" + df.head(20).to_string(index=False))
# cell-weighted summary
w = df['n'].to_numpy()
log(f"\ncell-weighted median ratio = {np.median(np.repeat(df['ratio'].values, (w/ w.min()).astype(int))):.1f}x")
log(f"ratio: min {df['ratio'].min()}x  median {df['ratio'].median()}x  max {df['ratio'].max()}x  "
    f"(1x=datasets mix; >>1x=residual technical batch)")
log("vs prior: CellClass-only ~15x, CellClass×Region ~8-13x")
log("DONE")
