"""Lightweight run-provenance stamping (review item #5-lite).

stamp(adata, script, params) -> writes adata.uns['provenance']; or
write_sidecar(path, script, params) -> writes <path>.provenance.json next to a TSV.
Records git SHA, script name+args, timestamp, python + key library versions.
"""
import json
import subprocess
import sys
import time
from pathlib import Path

REPO = '/Users/eg/brain_organoid'
_LIBS = ['numpy', 'scipy', 'pandas', 'anndata', 'scanpy', 'sklearn', 'scvi', 'torch',
         'scib_metrics', 'jax']


def git_sha(repo=REPO):
    try:
        return subprocess.check_output(['git', '-C', repo, 'rev-parse', 'HEAD'],
                                       text=True, stderr=subprocess.DEVNULL).strip()
    except Exception:
        return 'unknown'


def _lib_versions():
    import importlib
    out = {}
    for m in _LIBS:
        try:
            # str(): some libs (e.g. torch) return a non-str version object that
            # anndata cannot serialize into .uns -> coerce to plain str.
            out[m] = str(getattr(importlib.import_module(m), '__version__', '?'))
        except Exception:
            pass
    return out


def record(script, params):
    return {
        'script': str(script),
        'git_sha': git_sha(),
        'params': {k: (str(v) if not isinstance(v, (int, float, str, bool, type(None))) else v)
                   for k, v in (params or {}).items()},
        'timestamp': time.strftime('%Y-%m-%dT%H:%M:%S'),
        'python': sys.version.split()[0],
        'libs': _lib_versions(),
    }


def stamp(adata, script, params=None):
    """Attach provenance to an AnnData's .uns (serializable for h5ad)."""
    adata.uns['provenance'] = record(script, params)
    return adata


def write_sidecar(out_path, script, params=None):
    """Write <out_path>.provenance.json next to a non-AnnData output (e.g. TSV)."""
    p = Path(str(out_path) + '.provenance.json')
    p.write_text(json.dumps(record(script, params), indent=2))
    return p
