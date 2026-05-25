#!/usr/bin/env python3
"""Regenerate data/manifest.tsv from the actual data/processed/*.h5ad files.

rebuild_atlas.py --slugs rewrites the whole manifest from only the processed
slugs (footgun), which clobbered it to 5 rows. The h5ads on disk are the real
truth, so rebuild the manifest by scanning them. Excludes non-deposit files and
DROP_SLUGS (e.g. gse168323 = HNOCA overlap).
"""
import anndata as ad
import pandas as pd
from pathlib import Path

PROC = Path('/Users/eg/brain_organoid/data/processed')
OUT = Path('/Users/eg/brain_organoid/data/manifest.tsv')
EXCLUDE = {'atlas_v4_preprocessed'}      # the concatenated atlas, not a deposit
DROP_SLUGS = {'gse168323'}               # HNOCA dataset -> dropped from this atlas

rows = []
for p in sorted(PROC.glob('*.h5ad')):
    slug = p.stem
    if slug in EXCLUDE or slug in DROP_SLUGS:
        print(f"  exclude {slug}", flush=True); continue
    a = ad.read_h5ad(p, backed='r')
    o = a.obs
    acc = str(o['accession'].iloc[0]) if 'accession' in o and a.n_obs else slug.upper()
    n_ctrl = int(o['is_control'].sum()) if 'is_control' in o else a.n_obs
    ot = o['organoid_type'].mode().iloc[0] if 'organoid_type' in o and a.n_obs else ''
    nsamp = int(o['bio_sample'].nunique()) if 'bio_sample' in o else (int(o['gsm'].nunique()) if 'gsm' in o else 1)
    rows.append({'slug': slug, 'accession': acc, 'path': str(p),
                 'n_cells': a.n_obs, 'n_genes': a.n_vars, 'n_samples': nsamp,
                 'n_control_samples': nsamp, 'n_control_cells': n_ctrl,
                 'organoid_type': ot, 'filter': 'disk_scan_2026-05-24', 'status': 'ok'})
    a.file.close()

df = pd.DataFrame(rows)
df.to_csv(OUT, sep='\t', index=False)
print(f"\nWROTE {OUT}: {len(df)} rows | total cells {df.n_cells.sum():,}", flush=True)
