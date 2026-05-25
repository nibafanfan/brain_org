#!/usr/bin/env python3
"""Rebuild gse296775_strada control-only from per-sub-sample 10x trios.

WHY: this deposit (STRADA / mTORopathy forebrain organoids) packs MULTIPLE
sub-samples under each GSM, mixing isogenic control (_C1/_C2) with STRADA-mutant
(_M1/_M2, ± rapamycin). The GSM-level filter in rebuild_atlas.py can't separate
them, and load_10x_mtx_per_gsm groups all files by GSM prefix and reads a single
(mismatched) trio. So we select control sub-samples explicitly by filename.

Control = sub-sample whose name contains _C<digit>. Mutants (_M*) are dropped.
Region (Dorsal/Ventral) and age (D35/D77) are parsed from the sub-sample name.

Run with: /opt/homebrew/Caskroom/miniforge/base/bin/python
"""
import re
import sys
from pathlib import Path

import anndata as ad
import pandas as pd

ROOT = Path('/Users/eg/brain_organoid')
sys.path.insert(0, str(ROOT / 'scripts'))
from rebuild_atlas import _read_trio, apply_qc  # noqa: E402
from migrate_to_hnoca_schema import (  # noqa: E402
    load_canonical, map_to_hgnc, attach_canonical_var, migrate_obs, sanitize_for_h5,
)

SLUG = 'gse296775_strada'
ACCESSION = 'GSE296775'
RAW = ROOT / 'data' / 'raw' / 'gse296775'
OUT = ROOT / 'data' / 'processed' / f'{SLUG}.h5ad'

CONTROL_RE = re.compile(r'_C\d')          # _C1, _C2 = isogenic control
MUTANT_RE = re.compile(r'_M\d')
GSM_RE = re.compile(r'(GSM\d+)')


def main():
    canonical, e2h = load_canonical()
    mtx_files = sorted(f for f in RAW.iterdir() if f.name.endswith('_matrix.mtx.gz'))
    parts = []
    for mtx in mtx_files:
        stem = mtx.name[:-len('_matrix.mtx.gz')]          # e.g. GSM8976366_Dorsal_C1_D35_rep1-mRNA
        if not CONTROL_RE.search(stem):
            continue                                       # skip mutant / non-control
        bc = RAW / f'{stem}_barcodes.tsv.gz'
        ft = RAW / f'{stem}_features.tsv.gz'
        if not (bc.exists() and ft.exists()):
            print(f"  !! missing trio sibling for {stem}", flush=True); continue
        gsm = GSM_RE.match(stem).group(1)
        subsample = stem[len(gsm) + 1:].replace('-mRNA', '')   # Dorsal_C1_D35_rep1
        a = _read_trio(bc, ft, mtx, sample_id=subsample.lower(), gsm=gsm)
        region = 'Dorsal Forebrain' if 'Dorsal' in subsample else ('Ventral Forebrain' if 'Ventral' in subsample else 'Forebrain')
        age = 35 if 'D35' in subsample else (77 if 'D77' in subsample else -1)
        a.obs['bio_sample'] = subsample
        a.obs['organoid_type'] = region
        a.obs['organoid_age_days'] = age
        parts.append(a)
        print(f"  + {gsm} {subsample}: {a.n_obs} cells | {region} d{age}", flush=True)

    a = ad.concat(parts, axis=0, join='outer', fill_value=0, merge='same')
    a.obs_names_make_unique()
    print(f"control sub-samples: {len(parts)} | total {a.n_obs} cells", flush=True)

    # --- standard pipeline --- (apply_qc subsets rows, map_to_hgnc subsets genes;
    # both preserve obs of kept cells, so per-cell fields ride along.)
    a, qc = apply_qc(a, canonical, SLUG)
    a, mp = map_to_hgnc(a, canonical, e2h)
    a = attach_canonical_var(a, canonical)
    print(f"QC {qc['pre_qc']}->{qc['post_qc']} | genes {mp['in']}->{mp['kept']} ({mp['ns']})", flush=True)

    # stash per-cell fields migrate_obs would clobber (bio_sample<-gsm, age<-NaN)
    a.obs['_bio'] = a.obs['bio_sample'].astype(str)
    a.obs['_age'] = a.obs['organoid_age_days']
    a.obs['accession'] = ACCESSION
    a.obs['dataset_slug'] = SLUG
    a.obs['multi_lineage'] = 'False'        # each organoid is single-region (dorsal OR ventral)
    a.obs['is_control'] = True              # only _C sub-samples loaded
    a.obs['dataset_filter'] = 'subsample_control_C_2026-05-24'
    a = migrate_obs(a, SLUG)                 # organoid_type is NOT clobbered; bio_sample/age are
    a.obs['bio_sample'] = a.obs['_bio']
    a.obs['organoid_age_days'] = a.obs['_age']
    a.obs['tech_sample'] = (SLUG + '_' + a.obs['_bio'].astype(str)).astype(str)
    del a.obs['_bio'], a.obs['_age']
    a.obs['protocol'] = 'Conventional'
    a.obs['cell_type'] = 'unknown'
    a.obs['annot_source'] = 'subsample_parse_2026-05-24'

    a = sanitize_for_h5(a)
    a.write_h5ad(OUT, compression='gzip')
    print(f"WROTE {OUT} shape={a.shape}", flush=True)
    print(f"  organoid_type: {a.obs['organoid_type'].value_counts().to_dict()}", flush=True)
    print(f"  bio_sample (sub-samples): {a.obs['bio_sample'].nunique()} | is_control all: {bool(a.obs['is_control'].all())}", flush=True)


if __name__ == '__main__':
    main()
