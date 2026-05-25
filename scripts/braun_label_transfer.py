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
atlas_b = ad.read_h5ad(ATLAS, backed='r')
n_q = min(args.query_n, atlas_b.n_obs) if args.pilot else atlas_b.n_obs
qidx = np.sort(rng.choice(atlas_b.n_obs, n_q, replace=False)) if n_q < atlas_b.n_obs \
       else np.arange(atlas_b.n_obs)
query = atlas_b[qidx].to_memory()
log(f"query loaded: {query.shape}")

# rename query genes symbol->ensembl, dedup, set counts layer
ens = np.array([sym2ens.get(s, '') for s in query.var_names])
keep = ens != ''
query = query[:, keep].copy()
query.var_names = ens[keep]
query = query[:, ~query.var_names.duplicated()].copy()
query.layers['counts'] = query.X.copy()
query.obs[args.label_key] = UNLAB
query.obs[args.ref_batch_key] = query.obs[args.query_batch_key].astype(str)
log(f"query renamed to ensembl -> {query.shape}")

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

# crosstab vs organoid's own annotation if present
for col in ('annotation', 'cell_type', 'cell_type_original'):
    if col in query.obs and query.obs[col].notna().any():
        log(f"=== CellClass_pred vs organoid {col} (top rows) ===")
        ct = pd.crosstab(query.obs[col], query.obs['CellClass_pred'])
        log("\n" + ct.head(12).to_string())
        break

out = ad.AnnData(X=query.obsm['X_scanvi'],
                 obs=query.obs[['CellClass_pred', 'CellClass_conf',
                                'Region_pred', 'Region_conf',
                                args.query_batch_key]].copy())
op = ROOT / f'data/braun_transfer_{args.out_tag}.h5ad'
out.write_h5ad(op)
log(f"saved -> {op}")
log("DONE")
