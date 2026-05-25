#!/usr/bin/env python3
"""Q3: which primary (Braun) populations do organoids fail to recapitulate -
even multi-lineage ones?

Reference-coverage in the joint scArches latent: pool a balanced Braun + organoid
sample; for each Braun cell, fraction of its kNN that are organoid. ~0.5 = well
covered (organoids reach that latent region); ~0 = GAP. Aggregate by Braun
CellClass / Region / Age, separately for multi- vs single-lineage organoid cells.
Gaps persisting under multi-lineage = the answer. Expect (from Q2) older ages and
vascular/oligo to be poorly covered.
"""
import time
import numpy as np, pandas as pd, anndata as ad, scvi, torch
from sklearn.neighbors import NearestNeighbors

ROOT = '/Users/eg/brain_organoid'
MDIR = f'{ROOT}/data/braun_scanvi_full'
N = 80_000          # per source in each balanced pool
K = 30
t0 = time.time()
def log(m): print(f"[{time.time()-t0:7.1f}s] {m}", flush=True)
rng = np.random.default_rng(0)

# ---- Braun latent + labels
vn = list(torch.load(f'{MDIR}/model.pt', map_location='cpu', weights_only=False)['var_names'])
braun = ad.read_h5ad(f'{ROOT}/data/raw/braun_2023/braun_all.h5ad')[:, vn].copy()
braun.layers['counts'] = braun.X.copy()
for c in ['CellClass', 'Region', 'donor_id']:
    braun.obs[c] = braun.obs[c].astype(str)
scanvi = scvi.model.SCANVI.load(MDIR, adata=braun)
bidx = np.sort(rng.choice(braun.n_obs, N, replace=False))
bsub = braun[bidx].copy()
b_lat = scanvi.get_latent_representation(bsub)
b_obs = bsub.obs.copy()
b_obs['Age'] = pd.to_numeric(b_obs['Age'], errors='coerce')
b_obs['AgeBin'] = pd.cut(b_obs['Age'], [0, 8, 10, 12, 99],
                         labels=['<8pcw', '8-10', '10-12', '>12pcw'])
log(f"Braun latent {b_lat.shape}")

# ---- organoid latent + lineage
tr = ad.read_h5ad(f'{ROOT}/data/braun_transfer_full_knn.h5ad')
ml = tr.obs['multi_lineage'].astype(str).map(
    {'1': 'multi', 'True': 'multi', '0': 'single', 'False': 'single', 'No': 'single'})
oX = np.asarray(tr.X)
log(f"organoid latent {oX.shape}; multi={int((ml=='multi').sum()):,} single={int((ml=='single').sum()):,}")

def coverage(org_mask):
    oi = np.where(org_mask.to_numpy())[0]
    oi = rng.choice(oi, N, replace=False)
    pool = np.vstack([b_lat, oX[oi]])
    is_org = np.r_[np.zeros(N, bool), np.ones(N, bool)]
    nn = NearestNeighbors(n_neighbors=K + 1).fit(pool)
    _, idx = nn.kneighbors(pool[:N])           # query = Braun cells only
    return is_org[idx[:, 1:]].mean(1)          # frac organoid neighbors per Braun cell

cov_multi = coverage(ml == 'multi')
cov_single = coverage(ml == 'single')
log("coverage computed (multi + single)")

def summarize(by):
    g = pd.DataFrame({'by': b_obs[by].values, 'multi': cov_multi, 'single': cov_single})
    t = g.groupby('by', observed=True)[['multi', 'single']].mean()
    t['n'] = g.groupby('by', observed=True).size()
    t['multi_adds'] = (t['multi'] - t['single']).round(3)
    return t.round(3).sort_values('multi')

for by in ['CellClass', 'Region', 'AgeBin']:
    log(f"\n=== organoid coverage of Braun {by} (0.5=balanced/covered, ~0=GAP) ===")
    log("\n" + summarize(by).to_string())
out = pd.DataFrame({'CellClass': b_obs['CellClass'].values, 'Region': b_obs['Region'].values,
                    'AgeBin': b_obs['AgeBin'].values, 'cov_multi': cov_multi, 'cov_single': cov_single})
out.to_csv(f'{ROOT}/data/q3_coverage.tsv', sep='\t', index=False)
log("saved -> data/q3_coverage.tsv")
log("DONE")
