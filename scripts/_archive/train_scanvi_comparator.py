#!/usr/bin/env python3
"""Scoped label-aware integration comparator (Codex GO): scANVI on the organoid
atlas, semi-supervised by the calibrated Braun CellClass labels, to test whether
it integrates datasets better than plain scVI (the 'secret weapon').

Reuses the trained scVI model (scvi_model_v5_full) via SCANVI.from_scvi_model,
labels = CellClass_cal (abstained/Unknown -> unlabeled_category), batch_key as in
the scVI setup (tech_sample). Saves a comparator latent for scIB eval.

  /opt/homebrew/Caskroom/miniforge/base/bin/python3.13 scripts/train_scanvi_comparator.py
"""
import os
os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")
import sys, time
from pathlib import Path
import numpy as np, anndata as ad, scvi
from lightning.pytorch.callbacks import Callback
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _provenance import stamp

ROOT = Path('/Users/eg/brain_organoid')
PRE = ROOT / 'data/processed/atlas_v5_preprocessed.h5ad'
SCVI_DIR = ROOT / 'data/scvi_model_v5_full'
CAL = ROOT / 'data/braun_transfer_full_calibrated.h5ad'
EPOCHS = 20
ACC = 'mps'
UNLAB = 'Unknown'
t0 = time.time()
def log(m): print(f"[{time.time()-t0:8.1f}s] {m}", flush=True)


class EpochLogger(Callback):
    KEYS = ("elbo_train", "elbo_validation")
    def on_train_epoch_end(self, tr, pl):
        m = tr.callback_metrics
        log("  [ep %d/%d] %s" % (tr.current_epoch + 1, tr.max_epochs,
            " ".join(f"{k}={float(m[k]):.1f}" for k in self.KEYS if k in m)))


adata = ad.read_h5ad(PRE)
log(f"preprocessed atlas {adata.shape}")
cal = ad.read_h5ad(CAL, backed='r')
adata.obs['CellClass_cal'] = cal.obs['CellClass_cal'].reindex(adata.obs_names).fillna(UNLAB).values
n_lab = int((adata.obs['CellClass_cal'] != UNLAB).sum())
log(f"labels joined: {n_lab:,}/{adata.n_obs:,} labeled ({adata.obs['CellClass_cal'].nunique()} classes incl '{UNLAB}')")

model = scvi.model.SCVI.load(str(SCVI_DIR), adata=adata)
log("loaded scVI model")
scanvi = scvi.model.SCANVI.from_scvi_model(model, unlabeled_category=UNLAB,
                                           labels_key='CellClass_cal')
log(f"training scANVI ({EPOCHS} ep)…")
scanvi.train(max_epochs=EPOCHS, accelerator=ACC, batch_size=512,
             enable_progress_bar=False, callbacks=[EpochLogger()])
scanvi.save(str(ROOT / 'data/scanvi_comparator_model'), overwrite=True)

lat = scanvi.get_latent_representation()
out = ad.AnnData(X=lat.astype('float32'), obs=adata.obs.copy())
stamp(out, __file__, {'EPOCHS': EPOCHS, 'batch_size': 512, 'n_labeled': n_lab,
                      'from_model': str(SCVI_DIR)})
out.write_h5ad(ROOT / 'data/scanvi_comparator_latent.h5ad')
log(f"saved comparator latent {lat.shape} -> data/scanvi_comparator_latent.h5ad")
log("DONE")
