#!/usr/bin/env python3
"""Cell-type-stratified batch mixing on the scVI v5 atlas latent.

The global kNN metric conflates batch with biology (a dataset can look isolated
just because it uniquely contains some cell type). This controls for that: WITHIN
each transferred Braun CellClass, build a kNN graph among only that type's cells
and ask what fraction of neighbors come from the same dataset. Near the within-
type baseline (sum p_i^2) = well mixed across labs; >> baseline = residual batch.

  same-frac ~ baseline  -> e.g. radial glia from Lab A mix with radial glia Lab B
  same-frac >> baseline -> same cell type still splits by dataset (batch effect)
"""
import time
import numpy as np, pandas as pd, anndata as ad
from sklearn.neighbors import NearestNeighbors

ROOT = '/Users/eg/brain_organoid'
N_SUB = 300_000
CAP = 40_000          # max cells per class for the kNN (speed)
MIN_CELLS = 1500      # skip classes too rare to assess
K = 30
BATCH = 'dataset_slug'
t0 = time.time()
def log(m): print(f"[{time.time()-t0:7.1f}s] {m}", flush=True)

lat = ad.read_h5ad(f'{ROOT}/data/scvi_latent_v5_full.h5ad')
tr = ad.read_h5ad(f'{ROOT}/data/braun_transfer_full_knn.h5ad', backed='r')
log(f"latent {lat.shape}")
# join CellClass_pred onto latent by barcode
cc = tr.obs['CellClass_pred'].reindex(lat.obs_names)
lat.obs['CellClass'] = cc.values
lat = lat[lat.obs['CellClass'].notna()].copy()
log(f"joined CellClass; {lat.n_obs:,} cells labeled")

rng = np.random.default_rng(0)
if lat.n_obs > N_SUB:
    idx = np.sort(rng.choice(lat.n_obs, N_SUB, replace=False))
    lat = lat[idx].copy()
X = np.asarray(lat.X)
cls = lat.obs['CellClass'].to_numpy()
ds = lat.obs[BATCH].astype('category')
log(f"subsample {X.shape}; {ds.cat.categories.size} datasets total")

def same_frac(Xc, codes, k=K):
    k = min(k, len(codes) - 1)
    nn = NearestNeighbors(n_neighbors=k + 1).fit(Xc)
    _, idx = nn.kneighbors(Xc)
    idx = idx[:, 1:]                       # drop self
    return float((codes[idx] == codes[:, None]).mean())

def baseline(codes):
    p = pd.Series(codes).value_counts(normalize=True).to_numpy()
    return float((p ** 2).sum())

rows = []
for c in pd.Series(cls).value_counts().index:
    mask = cls == c
    n = int(mask.sum())
    if n < MIN_CELLS:
        continue
    Xc = X[mask]
    codes = ds.cat.codes.to_numpy()[mask]
    if n > CAP:                            # cap for speed, keep dataset proportions
        sel = rng.choice(n, CAP, replace=False)
        Xc, codes = Xc[sel], codes[sel]
    nds = np.unique(codes).size
    sf, bl = same_frac(Xc, codes), baseline(codes)
    rows.append((c, n, nds, sf, bl, sf / bl if bl else np.nan))
    log(f"  {c:18} n={n:>7} ndatasets={nds:>3} same={sf:.3f} base={bl:.3f} ratio={sf/bl:5.1f}x")

df = pd.DataFrame(rows, columns=['CellClass', 'n', 'n_datasets', 'same_frac', 'baseline', 'ratio'])
# overall (ignoring class) for comparison with prior global numbers
allcodes = ds.cat.codes.to_numpy()
sel = rng.choice(len(allcodes), min(CAP, len(allcodes)), replace=False)
g_sf = same_frac(X[sel], allcodes[sel]); g_bl = baseline(allcodes[sel])
log("")
log(f"GLOBAL (not stratified): same={g_sf:.3f} base={g_bl:.3f} ratio={g_sf/g_bl:.1f}x")
log("=== stratified summary (lower ratio = better cross-dataset mixing within type) ===")
log("\n" + df.round(3).to_string(index=False))
df.to_csv(f'{ROOT}/data/stratified_mixing.tsv', sep='\t', index=False)
log(f"saved -> data/stratified_mixing.tsv")
log("DONE")
