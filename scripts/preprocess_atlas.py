#!/usr/bin/env python3
"""Step 3: preprocess atlas_v4_full into a scVI-ready file.

RECOMPUTES our own HVGs on the atlas (does NOT reuse HNOCA's list), per directive.

Feature selection: analytic Pearson residuals (sctransform-style), the
sc-best-practices recommendation. Operates on RAW counts, so this is single-pass:
  read full raw -> HVG on raw counts -> keep raw counts of HVGs (counts layer)
  -> normalize_total(1e4 on full-gene totals) + log1p on the HVG subset (X).

Settings: pearson_residuals, n_top_genes=3000, batch_key='bio_sample'.
"""
import os
# Use the server's cores for the underlying BLAS / scipy / numba ops.
# Must be set BEFORE numpy is imported.
N_THREADS = 16
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMBA_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, str(N_THREADS))
import time, gc, threading
import numpy as np
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
sc.settings.n_jobs = N_THREADS

FULL      = 'data/atlas_v5_full.h5ad'
OUT       = 'data/processed/atlas_v5_preprocessed.h5ad'
HVG_TSV   = 'data/atlas_v5_hvg.tsv'
N_HVG     = 3000
FLAVOR    = 'cell_ranger'        # dispersion-binning, no loess (can't throw singular-matrix)
BATCH_KEY = 'bio_sample'
MIN_BATCH = 50                   # drop bio_samples with fewer cells (sampling noise)
TARGET    = 1e4

t0 = time.time()
def log(m): print(f"[{time.time()-t0:8.1f}s] {m}", flush=True)

log(f"reading full atlas {FULL} (~77 GB)")
A = ad.read_h5ad(FULL)
log(f"loaded {A.shape}, X={type(A.X).__name__} {A.X.dtype}")
assert sp.issparse(A.X)
A.X = A.X.tocsr()

log("computing per-cell totals over ALL genes (normalize denominator)")
totals = np.asarray(A.X.sum(axis=1)).ravel().astype(np.float64)

# --- FATAL-DATA FIX: drop zero-count cells before HVG ---------------------
# ~63,714 cells have total count == 0 in the canonical 36,842-gene space
# (genes lost during projection). All-zero cells contribute no information and
# distort per-gene mean/dispersion, so remove them before HVG + normalization.
n_zero = int((totals == 0).sum())
keep = totals > 0
log(f"FILTER: dropping {n_zero} zero-count cells of {A.n_obs} "
    f"({100*n_zero/A.n_obs:.2f}%)")
A = A[keep].copy()
gc.collect()
log(f"after zero-cell filter: {A.shape}")

# --- prune micro-batches: bio_samples with <MIN_BATCH cells distort per-bin
#     dispersion averages via technical sampling noise -----------------------
vc = A.obs[BATCH_KEY].value_counts()
small = list(vc.index[vc < MIN_BATCH])
if small:
    keepb = ~A.obs[BATCH_KEY].isin(small)
    log(f"PRUNE: removing {len(small)} bio_samples <{MIN_BATCH} cells "
        f"({int((~keepb).sum())} cells)")
    A = A[keepb.values].copy(); gc.collect()
log(f"after micro-batch prune: {A.shape} | {A.obs[BATCH_KEY].nunique()} batches")

# --- preserve raw counts, then log-normalize X for cell_ranger HVG --------
# cell_ranger uses a binned dispersion lookup (no loess), so it cannot throw
# the singular-matrix error seurat_v3 hit; it requires log-normalized input.
log("layers['counts']=raw; normalize_total(1e4)+log1p -> X")
A.layers['counts'] = A.X.copy()
sc.pp.normalize_total(A, target_sum=TARGET)
sc.pp.log1p(A)

# --- progress heartbeat (HVG batch loop has no per-iteration hook) --------
def _heartbeat(stop_evt, label):
    s = time.time()
    while not stop_evt.wait(15):
        print(f"   ...{label} running {time.time()-s:6.0f}s", flush=True)

log(f"HVG: highly_variable_genes(flavor={FLAVOR}, n_top_genes={N_HVG}, "
    f"batch_key={BATCH_KEY}) on lognorm | n_jobs={sc.settings.n_jobs}, "
    f"{A.obs[BATCH_KEY].nunique()} batches")
_stop = threading.Event()
_hb = threading.Thread(target=_heartbeat, args=(_stop, 'HVG'), daemon=True)
_hb.start()
try:
    sc.pp.highly_variable_genes(
        A, flavor=FLAVOR, n_top_genes=N_HVG, batch_key=BATCH_KEY)
finally:
    _stop.set(); _hb.join()
mask = np.asarray(A.var['highly_variable'].values)
log(f"selected {int(mask.sum())} HVGs")

log("writing HVG table + subsetting to HVGs (X=lognorm, layers['counts']=raw)")
A.var[mask].to_csv(HVG_TSV, sep='\t')
out = A[:, mask].copy()
del A
gc.collect()
out.uns['log1p'] = {'base': None}
out.uns['hvg_source'] = f'recomputed on atlas (NOT HNOCA): flavor={FLAVOR}, n_top_genes={N_HVG}, batch_key={BATCH_KEY}'
out.uns['normalization'] = 'normalize_total target_sum=1e4 on full 36842-gene totals, then log1p'
out.uns['qc_filter'] = f'dropped zero-count cells + bio_samples <{MIN_BATCH} cells'
log(f"output {out.shape}, layers={list(out.layers.keys())}, counts nnz={out.layers['counts'].nnz:,}")

log(f"writing {OUT} (gzip)")
out.write_h5ad(OUT, compression='gzip')
log("write complete; verifying")

B = ad.read_h5ad(OUT, backed='r')
xs = B.X[:2000]; cs = B.layers['counts'][:2000]
log(f"verify shape={B.shape} layers={list(B.layers.keys())}")
log(f"  X (lognorm) max={xs.max():.3f} | counts max={cs.max():.0f}")
log(f"  is_control all True: {bool(B.obs['is_control'].all())} | datasets: {B.obs['dataset_slug'].nunique()}")
log("VERIFY OK — done")
