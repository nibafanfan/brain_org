#!/usr/bin/env python3
"""Smoke test for atlas_common (review #5). Verifies each shared helper against
real data slices and tiny synthetic cases — no full-matrix reads. Run:
  /opt/homebrew/Caskroom/miniforge/base/bin/python3.13 scripts/test_atlas_common.py
"""
import sys, time
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad
sys.path.insert(0, str(Path(__file__).resolve().parent))
import atlas_common as ac

t0 = time.time()
def ok(m): print(f"[{time.time()-t0:6.1f}s] PASS {m}", flush=True)

cfg = ac.load_config()
assert cfg.canonical.exists() and cfg.atlas_full.exists() and cfg.braun_scanvi_model.exists()
ok(f"config resolves paths under {cfg.root} (n_sub default={cfg.defaults['n_sub']})")

# gene bridge
s2e = ac.sym2ens(cfg.canonical)
assert len(s2e) > 30000 and s2e.get('TSPAN6', '').startswith('ENSG')
ok(f"sym2ens: {len(s2e):,} symbol->ensembl pairs")

# model genes
vn = ac.model_genes(cfg.braun_scanvi_model)
assert len(vn) == 2006 and all(g.startswith('ENSG') for g in vn[:5])
ok(f"model_genes: {len(vn)} genes (braun_scanvi_full)")

# ens->col map + chunked read (2 chunks only, cheap)
atlas = ad.read_h5ad(cfg.atlas_full, backed='r')
e2c = ac.ens_to_col(atlas.var_names, s2e)
assert len(set(vn) & set(e2c)) == len(vn), "all model genes should map to atlas cols"
ok(f"ens_to_col: {len(vn)}/{len(vn)} model genes present in atlas")

ridx = np.arange(5000)                            # tiny row slice -> in-memory path
X, present = ac.read_atlas_genes(atlas, vn, s2e, row_idx=ridx)
assert X.shape == (5000, len(vn)) and str(X.dtype) == 'float32'
ok(f"read_atlas_genes(row_idx): {X.shape} dtype={X.dtype}")

# metrics: cp10k_log on a tiny known matrix
df = pd.DataFrame({'gA': [10.0, 0.0], 'gB': [0.0, 5.0]}, index=['c1', 'c2'])
out = ac.cp10k_log(df)
assert np.isclose(out.loc['c1', 'gA'], np.log1p(1e4)) and np.isclose(out.loc['c1', 'gB'], 0.0)
ok("cp10k_log: row-normalized log1p CP10K correct")

# same_frac / baseline on synthetic: two tight clusters w/ pure labels -> ~1.0
rng = np.random.default_rng(0)
Xc = np.vstack([rng.normal(0, 0.01, (200, 5)), rng.normal(10, 0.01, (200, 5))])
codes = np.r_[np.zeros(200, int), np.ones(200, int)]
sf, bl = ac.same_frac(Xc, codes, k=10), ac.baseline(codes)
assert sf > 0.95 and np.isclose(bl, 0.5), (sf, bl)
ok(f"same_frac={sf:.3f} (pure clusters ~1), baseline={bl:.2f} (balanced=0.5)")

# provenance re-exports present
assert callable(ac.stamp) and callable(ac.write_sidecar) and isinstance(ac.git_sha(), str)
ok("provenance re-exports (stamp/write_sidecar/git_sha)")
print("ALL PASS")
