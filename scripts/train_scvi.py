#!/usr/bin/env python3
"""Step 4: train scVI on the preprocessed organoid atlas.

Run with the miniforge BASE env python (3.13: scvi 1.4.1, torch 2.10, MPS):
  /opt/homebrew/Caskroom/miniforge/base/bin/python scripts/train_scvi.py [opts]

  --subsample N     train on N random cells (pilot/validation)
  --max-epochs E    epochs (default: scvi heuristic)
  --batch-key K     obs batch column (default tech_sample)
  --accelerator A   mps | cpu | auto  (default mps)
"""
import os
# MPS fallback: route any op Metal doesn't implement to the CPU instead of
# crashing the run. MUST be set before torch is imported (directly or via scvi).
os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")

import argparse, time, sys
from pathlib import Path
from brain_organoid.config import add_common_overrides, load_config, merge_cli_overrides, resolve_path, require_file, stamp_provenance
import numpy as np
import anndata as ad
import scvi
import torch
from lightning.pytorch.callbacks import Callback


class EpochLogger(Callback):
    """Print one clean line of metrics per epoch (ELBO + losses).

    Replaces the live tqdm bar, which garbles in a redirected/tmux terminal on
    long runs. Reads whatever scvi has logged into trainer.callback_metrics.
    """

    KEYS = ("elbo_train", "elbo_validation", "reconstruction_loss_train",
            "kl_local_train", "train_loss_epoch", "validation_loss")

    def on_train_epoch_end(self, trainer, pl_module):
        m = trainer.callback_metrics
        parts = []
        for k in self.KEYS:
            if k in m:
                try:
                    parts.append(f"{k}={float(m[k]):.2f}")
                except (TypeError, ValueError):
                    pass
        ep = trainer.current_epoch + 1
        total = trainer.max_epochs
        print(f"[epoch {ep:>3}/{total}] " + " ".join(parts), flush=True)


t0 = time.time()
def log(m): print(f"[{time.time()-t0:8.1f}s] {m}", flush=True)

ap = add_common_overrides(argparse.ArgumentParser())
ap.add_argument('--subsample', type=int, default=0)
ap.add_argument('--max-epochs', type=int, default=0)
ap.add_argument('--batch-key', default='tech_sample')
ap.add_argument('--accelerator', default='mps')
ap.add_argument('--out-tag', default='full')
ap.add_argument('--early-stopping', action='store_true')
ap.add_argument('--batch-size', type=int, default=2048,
                help='minibatch size; 128 starves the GPU at this scale')
ap.add_argument('--num-workers', type=int, default=0,
                help='DataLoader workers. WARNING on macOS (spawn): >0 pickles '
                     'the in-memory dataset into each worker -> N+1x RAM. '
                     'Leave 0 unless you have RAM headroom; data is already in '
                     'memory so this rarely helps on MPS anyway.')
args = ap.parse_args()
cfg = merge_cli_overrides(load_config(args.config), args)
out_tag = args.out_tag_override or args.out_tag
root = resolve_path(cfg, cfg.get('root','.'))
in_rel = args.in_override or cfg['paths'].get('input') or f"{cfg['outputs']['processed_dir']}/atlas_{cfg.get('out_tag','v5')}_preprocessed.h5ad"
IN = require_file(resolve_path(cfg, in_rel), 'preprocessed atlas input')

log(f"torch {torch.__version__} | mps={torch.backends.mps.is_available()} | scvi {scvi.__version__}")
log(f"reading {IN}")
adata = ad.read_h5ad(IN)
log(f"loaded {adata.shape}")

if args.subsample and args.subsample < adata.n_obs:
    rng = np.random.default_rng(0)
    idx = np.sort(rng.choice(adata.n_obs, args.subsample, replace=False))
    adata = adata[idx].copy()
    log(f"subsampled to {adata.shape}")

scvi.model.SCVI.setup_anndata(adata, layer='counts', batch_key=args.batch_key)
log(f"setup_anndata: layer=counts, batch_key={args.batch_key} "
    f"({adata.obs[args.batch_key].nunique()} batches)")

model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood='nb')
log(f"model: {model}")

kw = dict(
    accelerator=args.accelerator,
    devices=1,
    # (1) Large minibatch: feed the 60-core GPU and steady the gradient over
    #     the 505 batches. (2) trainer_kwargs (**kw beyond named args) go to the
    #     Lightning Trainer -> drop the live tqdm bar, print epoch-level metrics.
    batch_size=args.batch_size,
    enable_progress_bar=False,
    callbacks=[EpochLogger()],
)
# (3) num_workers is plumbed through datasplitter_kwargs -> the DataLoader.
#     scvi has no `dataloader_kwargs` arg. persistent_workers avoids respawning
#     (and re-pickling) every epoch under macOS spawn.
if args.num_workers > 0:
    kw['datasplitter_kwargs'] = dict(
        num_workers=args.num_workers,
        persistent_workers=True,
    )
if args.max_epochs:
    kw['max_epochs'] = args.max_epochs
if args.early_stopping:
    kw['early_stopping'] = True
    kw['early_stopping_patience'] = 10
log(f"training start: {kw}")
te = time.time()
model.train(**kw)
log(f"training done in {time.time()-te:.0f}s")

mdir = resolve_path(cfg, cfg['outputs']['models_dir']) / f'scvi_model_{out_tag}'
model.save(str(mdir), overwrite=True)
log(f"saved model -> {mdir}")

adata.obsm['X_scvi'] = model.get_latent_representation()
lat = resolve_path(cfg, cfg['outputs']['models_dir']) / f'scvi_latent_{out_tag}.h5ad'
out = ad.AnnData(X=adata.obsm['X_scvi'], obs=adata.obs.copy())
out.write_h5ad(lat)
log(f"saved latent ({adata.obsm['X_scvi'].shape}) -> {lat}")
stamp_provenance(lat.with_suffix('.provenance.json'), cfg, {'input': str(IN), 'latent': str(lat), 'model_dir': str(mdir), 'out_tag': out_tag})
log("DONE")
