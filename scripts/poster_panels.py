#!/usr/bin/env python3
"""Poster Fig.1-style panels, built directly from the split atlas files (no
assembly step) — same machinery as make_snapseed_umap.

UMAP scatter + marker panels use a subsample (can't plot 4M); the per-sample
composition bars + metadata summary use ALL cells (accurate proportions). Writes
panel-data TSVs, PNG/PDF figures, and a manifest. Deterministic (seed).

  python scripts/poster_panels.py [--n-sub 300000] [--out-tag v5]
"""
import argparse, json, subprocess, sys, time
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad, scanpy as sc
import scipy.sparse as sp
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parent))
from atlas_common import load_config, git_sha

MARKERS = ['SOX2', 'EOMES', 'DCX', 'NEUROD2', 'AQP4', 'P2RY12', 'CLDN5']
LEGEND_TOP_N = 15
t0 = time.time()
def log(m): print(f"[{time.time()-t0:7.1f}s] {m}", flush=True)

ap = argparse.ArgumentParser()
ap.add_argument('--root', default=None)
ap.add_argument('--n-sub', type=int, default=300_000)
ap.add_argument('--out-tag', default='v5')
args = ap.parse_args()
cfg = load_config(root=args.root)
OUT = cfg.root / 'data/poster'; OUT.mkdir(parents=True, exist_ok=True)
TAG = args.out_tag
rng = np.random.default_rng(cfg.defaults['seed'])
outputs = []
def topn(s, n=LEGEND_TOP_N):
    keep = s.value_counts().index[:n]
    return s.where(s.isin(keep), 'Other').astype('category')

# ---------- full-atlas obs for tables (all cells) ----------
cal = ad.read_h5ad(cfg.calibrated, backed='r')
o = cal.obs.copy()
o['age_days'] = pd.to_numeric(o['age_days'], errors='coerce')
o['ml'] = o['multi_lineage'].astype(str).map(
    {'1': 'multi', 'True': 'multi', '0': 'single', 'False': 'single', 'No': 'single'})
log(f"full atlas obs: {len(o):,} cells")

# metadata summary (datasets, samples, ages, protocols)
meta = {
    'n_cells': len(o), 'n_datasets': o['dataset_slug'].nunique(),
    'n_bio_samples': o['bio_sample'].nunique(),
    'n_organoid_types': o['organoid_type'].nunique(),
    'n_protocols': o['protocol'].nunique(),
    'multi_cells_pct': round((o['ml'] == 'multi').mean() * 100, 1),
    'age_median_days': round(o['age_days'].median(), 1),
    'age_min_days': round(o['age_days'].min(), 1), 'age_max_days': round(o['age_days'].max(), 1),
}
pd.Series(meta).to_csv(OUT / f'metadata_summary_{TAG}.tsv', sep='\t', header=['value'])
outputs.append(f'metadata_summary_{TAG}.tsv')
log("metadata: " + ", ".join(f"{k}={v}" for k, v in meta.items()))

# per-bio_sample cell-class composition (all cells), grouped by dataset, ordered by age
samp = (o.groupby('bio_sample', observed=True)
          .agg(dataset=('dataset_slug', 'first'), age=('age_days', 'median'),
               n=('CellClass_cal', 'size')))
comp = pd.crosstab(o['bio_sample'], o['CellClass_cal'], normalize='index')
comp.columns = [str(c) for c in comp.columns]      # drop CategoricalIndex (breaks .join)
comp = comp.join(samp).sort_values(['dataset', 'age'])
comp.to_csv(OUT / f'sample_composition_{TAG}.tsv', sep='\t')
outputs.append(f'sample_composition_{TAG}.tsv')
log(f"sample composition: {comp.shape[0]} bio_samples x {comp.shape[1]-3} classes")

# ---------- subsample for UMAP scatter + markers ----------
lat = ad.read_h5ad(cfg.latent)
idx = np.sort(rng.choice(lat.n_obs, min(args.n_sub, lat.n_obs), replace=False))
a = lat[idx].copy(); a.obsm['X_scvi'] = np.asarray(a.X)
names = a.obs_names
for c in ['CellClass_cal', 'Region_pred']:
    a.obs[c] = cal.obs[c].reindex(names).values
a.obs['age_days'] = pd.to_numeric(a.obs['age_days'], errors='coerce')
a.obs['organoid_type_topN'] = topn(a.obs['organoid_type'].astype(str))
# markers from full-gene atlas (symbols), lib-size normalized
atlas = ad.read_h5ad(cfg.atlas_full, backed='r')
mk = [m for m in MARKERS if m in set(atlas.var_names)]
Xm = atlas[idx].to_memory()[:, mk].X
Xm = (Xm.toarray() if sp.issparse(Xm) else np.asarray(Xm))
nc = a.obs['n_counts'].to_numpy()[:, None] if 'n_counts' in a.obs else Xm.sum(1, keepdims=True)
Xm = np.log1p(Xm / np.clip(nc, 1, None) * 1e4)
for i, g in enumerate(mk):
    a.obs[g] = Xm[:, i]
sc.pp.neighbors(a, use_rep='X_scvi', n_neighbors=15, random_state=cfg.defaults['seed'])
sc.tl.umap(a, random_state=cfg.defaults['seed'])
# Leiden cluster -> cluster-tiled categorical labels (majority vote per cluster) so the
# CellClass/Region panels read as crisp tiles (like scvi_umap_snapseed), not per-cell speckle.
sc.tl.leiden(a, resolution=2.0, flavor='igraph', n_iterations=2, directed=False,
             random_state=cfg.defaults['seed'])
for c in ['CellClass_cal', 'Region_pred']:
    maj = a.obs.groupby('leiden', observed=True)[c].agg(
        lambda s: s.mode().iloc[0] if len(s.mode()) else (s.iloc[0] if len(s) else 'NA'))
    a.obs[c + '_tiled'] = a.obs['leiden'].map(maj).astype('category')
log(f"UMAP + Leiden ({a.obs['leiden'].nunique()} clusters); markers {mk}")

# umap point table (per-cell labels + tiled labels)
pts = pd.DataFrame(a.obsm['X_umap'][:, :2], columns=['UMAP1', 'UMAP2'], index=names)
for c in (['CellClass_cal', 'CellClass_cal_tiled', 'Region_pred', 'Region_pred_tiled',
           'leiden', 'age_days', 'organoid_type', 'protocol'] + mk):
    if c in a.obs:
        pts[c] = a.obs[c].values
pts.to_csv(OUT / f'umap_points_{TAG}.tsv.gz', sep='\t', compression='gzip')
outputs.append(f'umap_points_{TAG}.tsv.gz')

# shuffle plotting order so dense classes don't bury sparse ones
a = a[rng.permutation(a.n_obs)].copy()
PT = 6  # larger points read cleaner on big panels

# Panel set 1: integrated UMAP by label/age/protocol — large panels (one per cell)
fig, ax = plt.subplots(2, 2, figsize=(28, 24))
sc.pl.umap(a, color='CellClass_cal_tiled', ax=ax[0, 0], show=False, size=PT, legend_loc='right margin', frameon=False, title='CellClass (cluster-tiled)')
sc.pl.umap(a, color='Region_pred_tiled', ax=ax[0, 1], show=False, size=PT, legend_loc='right margin', frameon=False, title='Region (cluster-tiled)')
sc.pl.umap(a, color='age_days', ax=ax[1, 0], show=False, size=PT, cmap='viridis', frameon=False, title='organoid age (days)')
sc.pl.umap(a, color='organoid_type_topN', ax=ax[1, 1], show=False, size=PT, legend_loc='right margin', frameon=False, title=f'organoid_type (top {LEGEND_TOP_N})')
fig.suptitle('Integrated organoid atlas (scVI latent UMAP)', fontsize=18, y=1.0)
fig.tight_layout(rect=[0, 0, 1, 0.97])
for ext in ('png', 'pdf'):
    fig.savefig(OUT / f'fig1_integrated_umap_{TAG}.{ext}', dpi=150)
outputs += [f'fig1_integrated_umap_{TAG}.png', f'fig1_integrated_umap_{TAG}.pdf']
plt.close(fig)

# Panel set 2: marker grid — bigger panels + per-gene p99 contrast (true CP10K stays correct)
ncol = 3; nrow = int(np.ceil(len(mk) / ncol))
fig, ax = plt.subplots(nrow, ncol, figsize=(6.5 * ncol, 5.5 * nrow)); ax = np.atleast_1d(ax).ravel()
for i, g in enumerate(mk):
    sc.pl.umap(a, color=g, ax=ax[i], show=False, size=PT, cmap='viridis', frameon=False,
               vmin=0, vmax='p99', title=g)            # p99 vmax -> use the dynamic range
for j in range(len(mk), len(ax)):
    ax[j].axis('off')
fig.suptitle('Marker-gene expression (log1p CP10K, vmax=p99)', fontsize=15); fig.tight_layout()
fig.savefig(OUT / f'fig1_markers_{TAG}.png', dpi=150); outputs.append(f'fig1_markers_{TAG}.png')
plt.close(fig)

# Panel set 3: per-sample composition stacked bars (all cells)
classes = [c for c in comp.columns if c not in ('dataset', 'age', 'n')]
fig, axb = plt.subplots(figsize=(20, 6))
bottom = np.zeros(len(comp))
cmap = plt.get_cmap('tab20')
for k, cls in enumerate(classes):
    axb.bar(range(len(comp)), comp[cls].values, bottom=bottom, width=1.0,
            color=cmap(k % 20), label=cls)
    bottom += comp[cls].values
axb.set(xlim=(0, len(comp)), ylim=(0, 1), xlabel='biological sample (grouped by dataset, ordered by age)',
        ylabel='cell-class fraction', title='Per-sample cell-class composition')
axb.legend(ncol=6, fontsize=7, loc='upper center', bbox_to_anchor=(0.5, -0.12))
fig.tight_layout()
fig.savefig(OUT / f'fig1_sample_composition_{TAG}.png', dpi=120, bbox_inches='tight')
outputs.append(f'fig1_sample_composition_{TAG}.png')
plt.close(fig)

manifest = {'git_sha': git_sha(), 'script': __file__, 'tag': TAG, 'n_sub': int(len(idx)),
            'seed': cfg.defaults['seed'], 'markers': mk,
            'inputs': {'latent': str(cfg.latent), 'calibrated': str(cfg.calibrated),
                       'atlas_full': str(cfg.atlas_full)},
            'outputs': sorted(outputs), 'timestamp': time.strftime('%Y-%m-%dT%H:%M:%S')}
(OUT / f'figure1_manifest_{TAG}.json').write_text(json.dumps(manifest, indent=2))
log(f"wrote {len(outputs)} outputs + manifest -> {OUT}")
log("DONE")
