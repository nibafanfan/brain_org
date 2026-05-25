#!/usr/bin/env python3
"""Braun 2023 -> organoid atlas label transfer via scvi-native scArches.

Stage 1  train scVI -> scANVI on the Braun fetal-brain reference (label_key=CellClass).
Stage 2  scArches surgery-map the organoid query; predict CellClass + soft-prob confidence.
Stage 3  kNN-transfer Region from Braun onto the query in the joint scANVI latent.
Stage 4  sanity report + save query obs/latent.

Gene bridge: organoid var_names are HGNC symbols; Braun is Ensembl. Map via
data/reference/hnoca_var_canonical.tsv (hgnc_symbol -> ensembl), intersect with Braun.

Run with the miniforge BASE python (scvi 1.4.1, MPS):
  /opt/homebrew/Caskroom/miniforge/base/bin/python3.13 scripts/braun_label_transfer.py --pilot

--pilot subsamples (ref-n/query-n) and uses few epochs to validate end-to-end fast.
Drop --pilot for the full run (all cells; bump epochs).
"""
import os
os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")

import argparse, time
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scvi
from lightning.pytorch.callbacks import Callback
from sklearn.neighbors import KNeighborsClassifier

ROOT = Path('/Users/eg/brain_organoid')
BRAUN = ROOT / 'data/raw/braun_2023/braun_all.h5ad'
ATLAS = ROOT / 'data/atlas_v5_full.h5ad'
CANON = ROOT / 'data/reference/hnoca_var_canonical.tsv'

t0 = time.time()
def log(m): print(f"[{time.time()-t0:8.1f}s] {m}", flush=True)


class EpochLogger(Callback):
    KEYS = ("elbo_train", "elbo_validation", "reconstruction_loss_train")
    def on_train_epoch_end(self, trainer, pl_module):
        m = trainer.callback_metrics
        parts = [f"{k}={float(m[k]):.1f}" for k in self.KEYS if k in m]
        print(f"  [epoch {trainer.current_epoch+1:>3}/{trainer.max_epochs}] "
              + " ".join(parts), flush=True)


ap = argparse.ArgumentParser()
ap.add_argument('--pilot', action='store_true')
ap.add_argument('--ref-n', type=int, default=100_000)
ap.add_argument('--query-n', type=int, default=200_000)
ap.add_argument('--ref-epochs', type=int, default=40)
ap.add_argument('--scanvi-epochs', type=int, default=20)
ap.add_argument('--surgery-epochs', type=int, default=40)
ap.add_argument('--n-hvg', type=int, default=2000)
ap.add_argument('--accelerator', default='mps')
ap.add_argument('--batch-size', type=int, default=512)
ap.add_argument('--out-tag', default='pilot')
ap.add_argument('--ref-batch-key', default='donor_id')
ap.add_argument('--query-batch-key', default='tech_sample')
ap.add_argument('--label-key', default='CellClass')
ap.add_argument('--region-key', default='Region')
args = ap.parse_args()

UNLAB = 'Unknown'
train_kw = dict(accelerator=args.accelerator, batch_size=args.batch_size,
                enable_progress_bar=False, callbacks=[EpochLogger()])

# Curated rare-lineage markers force-added to the Braun HVG panel. These cell
# types exist in the atlas (~1-3%: microglia P2RY12/TMEM119, endothelial
# CLDN5/PECAM1) but their specific markers are too low-variance to be selected as
# HVGs, so scANVI can't resolve Immune/Vascular without them (the pilot's 0%).
RARE_PANEL = [
    # microglia / immune / myeloid
    'PTPRC', 'AIF1', 'CX3CR1', 'P2RY12', 'TMEM119', 'C1QA', 'C1QB', 'C1QC',
    'CSF1R', 'TYROBP', 'FCER1G', 'CD68', 'ITGAM', 'CD74', 'LAPTM5',
    # endothelial / vascular
    'PECAM1', 'CLDN5', 'CDH5', 'FLT1', 'KDR', 'VWF', 'CD34', 'A2M', 'ESAM',
    'EGFL7', 'EMCN',
    # mural / pericyte
    'PDGFRB', 'RGS5', 'ACTA2', 'NOTCH3', 'KCNJ8',
    # erythrocyte
    'HBB', 'HBA1', 'HBA2', 'ALAS2', 'GYPA',
    # fibroblast
    'COL1A1', 'COL1A2', 'COL3A1', 'DCN', 'LUM',
    # oligodendrocyte
    'SOX10', 'OLIG1', 'OLIG2', 'MBP', 'PLP1', 'MOG',
]

log(f"scvi {scvi.__version__} | pilot={args.pilot} | out-tag={args.out_tag}")

# ---------------------------------------------------------------- gene bridge
can = pd.read_csv(CANON, sep='\t')
sym2ens = {s: e for s, e in zip(can['hgnc_symbol'].astype(str), can['ensembl'].astype(str))
           if isinstance(e, str) and e.startswith('ENSG')}
log(f"bridge: {len(sym2ens)} symbol->ensembl pairs")

# ------------------------------------------------------------ Stage 1: Braun
braun_b = ad.read_h5ad(BRAUN, backed='r')
rng = np.random.default_rng(0)
n_ref = min(args.ref_n, braun_b.n_obs) if args.pilot else braun_b.n_obs
ridx = np.sort(rng.choice(braun_b.n_obs, n_ref, replace=False)) if n_ref < braun_b.n_obs \
       else np.arange(braun_b.n_obs)
braun = braun_b[ridx].to_memory()
log(f"Braun loaded: {braun.shape}")

# shared gene panel in Ensembl space
atlas_syms = ad.read_h5ad(ATLAS, backed='r').var_names
atlas_ens = {sym2ens[s] for s in atlas_syms if s in sym2ens}
shared = [g for g in braun.var_names if g in atlas_ens]
log(f"shared genes (Braun ∩ atlas-mapped): {len(shared)}")
braun = braun[:, shared].copy()

# HVG via scanpy 'seurat' flavor (seurat_v3 banned: scikit-misc ABI conflict)
braun.layers['counts'] = braun.X.copy()
sc.pp.normalize_total(braun, target_sum=1e4)
sc.pp.log1p(braun)
sc.pp.highly_variable_genes(braun, n_top_genes=args.n_hvg, flavor='seurat')
# force-include curated rare-lineage markers that are present in the shared gene
# space but too low-variance to be HVGs — needed to resolve Immune/Vascular/etc.
panel_ens = {sym2ens[s] for s in RARE_PANEL if s in sym2ens}
forced = braun.var_names.isin(panel_ens) & ~braun.var['highly_variable'].values
braun.var.loc[forced, 'highly_variable'] = True
log(f"forced {int(forced.sum())} rare-lineage markers into HVG panel "
    f"({len(panel_ens)} of {len(RARE_PANEL)} bridged to ensembl & in shared space)")
braun = braun[:, braun.var.highly_variable].copy()
braun.X = braun.layers['counts'].copy()
log(f"HVG-subset Braun -> {braun.shape}")

# drop tiny label/batch categories that break stratified validation split
braun.obs[args.label_key] = braun.obs[args.label_key].astype(str)
braun.obs[args.ref_batch_key] = braun.obs[args.ref_batch_key].astype(str)

scvi.model.SCVI.setup_anndata(braun, layer='counts', batch_key=args.ref_batch_key)
ref = scvi.model.SCVI(braun, n_layers=2, n_latent=30, gene_likelihood='nb')
log(f"training reference scVI ({args.ref_epochs} ep)…")
ref.train(max_epochs=args.ref_epochs, **train_kw)

log(f"training reference scANVI on {args.label_key} ({args.scanvi_epochs} ep)…")
scanvi = scvi.model.SCANVI.from_scvi_model(ref, unlabeled_category=UNLAB,
                                           labels_key=args.label_key)
scanvi.train(max_epochs=args.scanvi_epochs, **train_kw)
mdir = ROOT / f'data/braun_scanvi_{args.out_tag}'
scanvi.save(str(mdir), overwrite=True)
log(f"saved reference scANVI -> {mdir}")

# ----------------------------------------------------------- Stage 2: query
# Memory-safe load: keep ONLY the ~2000 reference genes. The full atlas X is
# ~92 GB in memory (4M x 36842 @ 7.7% density) and copying it OOMs even on
# 192 GB. So map the reference gene panel back to atlas columns and read just
# those, streaming row-chunks for the full run (atlas X is CSR / row-major).
import scipy.sparse as sp
ref_genes = list(braun.var_names)                       # ensembl: the scANVI panel
atlas_b = ad.read_h5ad(ATLAS, backed='r')
atlas_ens = np.array([sym2ens.get(s, '') for s in atlas_b.var_names])
ref_set = set(ref_genes)
keep_cols, keep_ens, seen = [], [], set()
for c, e in enumerate(atlas_ens):
    if e in ref_set and e not in seen:
        seen.add(e); keep_cols.append(c); keep_ens.append(e)
keep_cols = np.array(keep_cols)
log(f"query: {len(keep_cols)}/{len(ref_genes)} reference genes present in atlas "
    f"(missing get zero-padded by prepare_query_anndata)")

n_q = min(args.query_n, atlas_b.n_obs) if args.pilot else atlas_b.n_obs
if args.pilot and n_q < atlas_b.n_obs:
    qidx = np.sort(rng.choice(atlas_b.n_obs, n_q, replace=False))
    Xq = atlas_b[qidx].to_memory().X[:, keep_cols]
    qobs = atlas_b.obs.iloc[qidx].copy()
else:
    qobs = atlas_b.obs.copy()
    parts, seen_n = [], 0
    for i, item in enumerate(atlas_b.chunked_X(250_000)):
        chunk = item[0] if isinstance(item, tuple) else item  # yields (X, start, end)
        parts.append(chunk[:, keep_cols]); seen_n += chunk.shape[0]
        if i % 4 == 0:
            log(f"  query chunk {i}: {seen_n:,}/{atlas_b.n_obs:,} cells read")
    Xq = sp.vstack(parts).tocsr(); del parts
Xq = Xq.astype('float32')
query = ad.AnnData(X=Xq, obs=qobs)
query.var_names = keep_ens
query.layers['counts'] = query.X.copy()                 # small now (~7 GB)
query.obs[args.label_key] = UNLAB
query.obs[args.ref_batch_key] = query.obs[args.query_batch_key].astype(str)
log(f"query built: {query.shape}")

# align to reference gene panel (zero-pads missing), then surgery
scvi.model.SCANVI.prepare_query_anndata(query, scanvi)
q_model = scvi.model.SCANVI.load_query_data(query, scanvi)
log(f"surgery training query ({args.surgery_epochs} ep)…")
q_model.train(max_epochs=args.surgery_epochs, plan_kwargs={'weight_decay': 0.0}, **train_kw)

pred = q_model.predict()
soft = q_model.predict(soft=True)
query.obs['CellClass_pred'] = np.asarray(pred)
query.obs['CellClass_conf'] = soft.to_numpy().max(1) if hasattr(soft, 'to_numpy') else np.asarray(soft).max(1)
query.obsm['X_scanvi'] = q_model.get_latent_representation()
log("CellClass predicted (+ confidence)")

# ------------------------------------------------- Stage 3: Region via kNN
ref_lat = scanvi.get_latent_representation()
knn = KNeighborsClassifier(n_neighbors=30, n_jobs=-1)
knn.fit(ref_lat, braun.obs[args.region_key].astype(str).to_numpy())
rprob = knn.predict_proba(query.obsm['X_scanvi'])
query.obs['Region_pred'] = knn.classes_[rprob.argmax(1)]
query.obs['Region_conf'] = rprob.max(1)
log("Region kNN-transferred (+ confidence)")

# ------------------------------------------------------- Stage 4: report
log("=== TRANSFERRED CellClass (organoid query) ===")
vc = query.obs['CellClass_pred'].value_counts(normalize=True)
for k, v in vc.items():
    log(f"   {k:24} {v*100:5.1f}%")
log(f"   mean CellClass confidence: {query.obs['CellClass_conf'].mean():.3f} "
    f"| frac > 0.8: {(query.obs['CellClass_conf']>0.8).mean():.2f}")
log("=== TRANSFERRED Region ===")
for k, v in query.obs['Region_pred'].value_counts(normalize=True).items():
    log(f"   {k:24} {v*100:5.1f}%")
log(f"   mean Region confidence: {query.obs['Region_conf'].mean():.3f}")

# crosstab vs organoid's own annotation — pick the col with the most MEANINGFUL
# labels (cell_type is 100% 'unknown'; cell_type_original has ~47k real labels)
NULLISH = {'unknown', 'nan', 'none', '', 'na'}
def meaningful_n(col):
    if col not in query.obs: return -1
    return (~query.obs[col].astype(str).str.lower().isin(NULLISH)).sum()
# NB: exclude 'cell_type_origin' — it's the pluripotent line (esc/ipsc), NOT a
# cell type, so it can't validate CellClass even though it's fully populated.
cands = sorted(('annotation', 'cell_type', 'cell_type_original'),
               key=meaningful_n, reverse=True)
best = cands[0]
if meaningful_n(best) > 0:
    m = ~query.obs[best].astype(str).str.lower().isin(NULLISH)
    log(f"=== CellClass_pred vs organoid {best} ({m.sum():,} labeled cells) ===")
    ct = pd.crosstab(query.obs.loc[m, best], query.obs.loc[m, 'CellClass_pred'])
    log("\n" + ct.to_string())

# --- Region_pred vs finalized organoid_type (100% coverage) — headline cross-check.
# organoid_type is the strongest sheet ground truth: a Midbrain organoid should
# map to Midbrain region, a Cortical/Cerebral one to Telencephalon/Forebrain.
if 'organoid_type' in query.obs:
    ot = query.obs['organoid_type'].astype(str)
    m = ~ot.str.lower().isin(NULLISH)
    if m.sum() > 0:
        log(f"=== Region_pred vs organoid_type ({m.sum():,} cells, row-normalized %) ===")
        ct = (pd.crosstab(query.obs.loc[m, 'organoid_type'],
                          query.obs.loc[m, 'Region_pred'], normalize='index') * 100)
        top_ot = ot[m].value_counts().head(12).index
        log("\n" + ct.loc[ct.index.intersection(top_ot)].round(1).to_string())

# --- transfer confidence stratified by annotation provenance (gsm vs deposit),
# so benchmark claims aren't inflated/deflated by coarse deposit-level fallback.
if 'annotation_level' in query.obs:
    log("=== confidence by annotation_level (gsm = authoritative, deposit = coarse) ===")
    for lvl, sub in query.obs.groupby('annotation_level', observed=True):
        log(f"   {lvl:8} n={len(sub):>9,} | CellClass conf {sub['CellClass_conf'].mean():.3f}"
            f" | Region conf {sub['Region_conf'].mean():.3f}")

# save transfer + the columns needed to stratify / cross-check downstream
keep_obs = ['CellClass_pred', 'CellClass_conf', 'Region_pred', 'Region_conf',
            args.query_batch_key]
for extra in ('organoid_type', 'annotation_level', 'cell_type_original', 'gsm'):
    if extra in query.obs and extra not in keep_obs:
        keep_obs.append(extra)
out = ad.AnnData(X=query.obsm['X_scanvi'], obs=query.obs[keep_obs].copy())
op = ROOT / f'data/braun_transfer_{args.out_tag}.h5ad'
out.write_h5ad(op)
log(f"saved -> {op}")
log("DONE")
