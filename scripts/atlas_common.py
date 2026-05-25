"""Shared helpers for the benchmark pipeline (review item #5 de-dup).

Consolidates the logic that was copy-pasted across the benchmark/transfer scripts:
central config, the HGNC-symbol<->Ensembl gene bridge, model-gene lookup, the
memory-safe chunked atlas reader, Braun reindexing, and the pseudobulk / kNN-mixing
metrics. Import these instead of re-implementing inline.

  from atlas_common import load_config, sym2ens, model_genes, ens_to_col, \
                          read_atlas_genes, reindex_braun, cp10k_log, same_frac, baseline
Provenance helpers are re-exported from _provenance.
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp

from _provenance import stamp, write_sidecar, git_sha  # re-export

CHUNK = 250_000


# ---------------------------------------------------------------- config
class Config:
    """Resolved config: attribute access to absolute Paths + defaults."""
    def __init__(self, d, root):
        self.root = Path(root)
        self._paths = {k: self.root / v for k, v in d.get('paths', {}).items()}
        self.defaults = d.get('defaults', {})

    def __getattr__(self, name):
        if name in self.__dict__.get('_paths', {}):
            return self._paths[name]
        raise AttributeError(name)

    def path(self, name):
        return self._paths[name]


def load_config(config_path=None, root=None):
    """Load config.yaml; `root` (arg > ATLAS_ROOT env > yaml) overrides path base."""
    import yaml
    here = Path(__file__).resolve().parent
    config_path = Path(config_path) if config_path else here.parent / 'config.yaml'
    d = yaml.safe_load(open(config_path))
    root = root or os.environ.get('ATLAS_ROOT') or d.get('root') or here.parent
    return Config(d, root)


# ---------------------------------------------------------------- gene bridge
def sym2ens(canonical_path):
    """HGNC symbol -> Ensembl id map (atlas var_names are symbols; Braun is Ensembl)."""
    can = pd.read_csv(canonical_path, sep='\t')
    return {s: e for s, e in zip(can['hgnc_symbol'].astype(str), can['ensembl'].astype(str))
            if isinstance(e, str) and e.startswith('ENSG')}


def model_genes(model_dir):
    """The Ensembl gene list a saved scvi/scanvi model expects (from model.pt)."""
    import torch
    pt = Path(model_dir) / 'model.pt' if Path(model_dir).is_dir() else Path(model_dir)
    return list(torch.load(pt, map_location='cpu', weights_only=False)['var_names'])


def ens_to_col(atlas_var_names, s2e):
    """Map Ensembl id -> column index in the symbol-keyed atlas (first occurrence)."""
    m = {}
    for c, s in enumerate(atlas_var_names):
        e = s2e.get(s, '')
        if e and e not in m:
            m[e] = c
    return m


# ---------------------------------------------------------------- io
def read_atlas_genes(atlas_backed, ens_genes, s2e, row_idx=None, chunk=CHUNK):
    """Memory-safe read of the full-gene atlas restricted to `ens_genes` (Ensembl).

    Returns (X_csr float32, present_ens_list). Streams row-chunks (atlas X is CSR)
    so the ~92 GB full matrix is never materialized. If row_idx given, reads just
    those rows in memory; else streams all rows.
    """
    e2c = ens_to_col(atlas_backed.var_names, s2e)
    present = [g for g in ens_genes if g in e2c]
    cols = np.array([e2c[g] for g in present])
    if row_idx is not None:
        X = atlas_backed[row_idx].to_memory().X[:, cols]
    else:
        parts = []
        for item in atlas_backed.chunked_X(chunk):
            chunkX = item[0] if isinstance(item, tuple) else item
            parts.append(chunkX[:, cols])
        X = sp.vstack(parts).tocsr()
    return X.astype('float32'), present


def reindex_braun(braun_path, ens_genes):
    """Load Braun reindexed to `ens_genes` with raw counts in X and layers['counts']."""
    import anndata as ad
    b = ad.read_h5ad(braun_path)[:, list(ens_genes)].copy()
    b.layers['counts'] = b.X.copy()
    return b


# ---------------------------------------------------------------- metrics
def cp10k_log(counts_by_group):
    """Pseudobulk normalize: rows=group, cols=gene mean counts -> CP10K -> log1p."""
    return np.log1p(counts_by_group.div(counts_by_group.sum(1), axis=0) * 1e4)


def baseline(codes):
    """Expected same-neighbor fraction under random mixing: sum p_i^2."""
    p = pd.Series(codes).value_counts(normalize=True).to_numpy()
    return float((p ** 2).sum())


def same_frac(X, codes, k=30):
    """Mean fraction of each cell's k nearest neighbors sharing its `codes` label."""
    from sklearn.neighbors import NearestNeighbors
    codes = np.asarray(codes)
    k = min(k, len(codes) - 1)
    nn = NearestNeighbors(n_neighbors=k + 1).fit(X)
    idx = nn.kneighbors(X)[1][:, 1:]            # drop self
    return float((codes[idx] == codes[:, None]).mean())
