#!/usr/bin/env python3
"""Reuse archived deposits that are safe (HGNC namespace + valid pct_mito).

For each SAFE_REUSE row in data/archive_reuse_audit.tsv:
  1. Load archived h5ad from _archive/processed_v2_pre_rebuild_2026-05-17 or
     fall back to _archive/processed_pre_rebuild_2026-05-17.
  2. For INCLUDE_EXPLICIT decisions, filter cells to gsm whitelist.
  3. Ensure obs schema matches new rebuild (is_control=True, dataset_filter
     tag, plus HNOCA-style columns via migrate_obs).
  4. Re-attach canonical var (in case archive has different var schema).
  5. Write to data/processed/<slug>.h5ad and append to manifest.

For reuse_redo_qc rows: same but recompute pct_mito after re-mapping to HGNC.
"""
import argparse
import json
import sys
import time
import traceback
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path('/Users/eg/brain_organoid')
sys.path.insert(0, str(ROOT / 'scripts'))
from migrate_to_hnoca_schema import (  # noqa
    load_canonical, map_to_hgnc, attach_canonical_var, migrate_obs,
    sanitize_for_h5, add_lengthnorm_layer, SMART_SEQ2,
    normalize_ensembl, detect_namespace,
)
from rebuild_atlas import build_mt_mask, apply_qc  # noqa

AUDIT = ROOT / 'data' / 'archive_reuse_audit.tsv'
CONFIG = ROOT / 'data' / 'rebuild_config.tsv'
ARCHIVE_V2 = ROOT / 'data' / '_archive' / 'processed_v2_pre_rebuild_2026-05-17'
ARCHIVE_V1 = ROOT / 'data' / '_archive' / 'processed_pre_rebuild_2026-05-17'
OUT_DIR = ROOT / 'data' / 'processed'
MANIFEST = ROOT / 'data' / 'manifest.tsv'
LOG_PATH = ROOT / 'data' / 'reuse_log.json'


def process_one(slug, audit_row, cfg_row, canonical, ensembl_to_hgnc):
    """Reuse one archived deposit."""
    verdict = audit_row['verdict']
    decision = audit_row.get('decision', '')

    # Find archive path
    arch = ARCHIVE_V2 if audit_row['archive'] == 'v2' else ARCHIVE_V1
    src = arch / f'{slug}.h5ad'
    if not src.exists():
        return {'slug': slug, 'status': 'no_archive_file'}

    t0 = time.time()
    a = ad.read_h5ad(src)
    n_loaded = a.n_obs

    # GSM filter for INCLUDE_EXPLICIT
    if decision == 'INCLUDE_EXPLICIT':
        gsms = cfg_row.get('gsms', '')
        if gsms:
            gsm_set = set(g.strip().upper() for g in gsms.split(',') if g.strip())
            if 'gsm' in a.obs.columns:
                keep_mask = a.obs['gsm'].astype(str).str.upper().isin(gsm_set).values
                a = a[keep_mask].copy()
                if a.n_obs == 0:
                    return {'slug': slug, 'status': 'zero_post_gsm_filter',
                            'n_loaded': n_loaded, 'gsm_set_size': len(gsm_set)}
            else:
                return {'slug': slug, 'status': 'no_gsm_col_in_archive',
                        'n_loaded': n_loaded}

    # If reuse_redo_qc: re-map to HGNC + recompute QC (ensures MT mask is canonical)
    if verdict == 'reuse_redo_qc':
        a, map_stats = map_to_hgnc(a, canonical, ensembl_to_hgnc)
        if a is None:
            return {'slug': slug, 'status': 'fail_map_post_archive', **map_stats}
        a, qc_stats = apply_qc(a, canonical, slug)
        if a.n_obs == 0:
            return {'slug': slug, 'status': 'zero_post_qc_redo', **qc_stats}

    # Ensure obs fields match rebuild output
    a.obs['accession'] = cfg_row.get('accession', slug.upper())
    a.obs['dataset_slug'] = slug
    a.obs['organoid_type'] = str(cfg_row.get('organoid_type', '') or '')
    a.obs['multi_lineage'] = str(cfg_row.get('multi_lineage', '') or '')
    a.obs['is_control'] = True
    a.obs['dataset_filter'] = f'reuse_{decision}_2026-05-17'

    # Apply HNOCA obs migration (idempotent — adds _original cols if missing)
    a = migrate_obs(a, slug)
    a = attach_canonical_var(a, canonical)
    if cfg_row.get('is_smartseq2', False) or slug in SMART_SEQ2:
        a = add_lengthnorm_layer(a)
    a = sanitize_for_h5(a)

    out_path = OUT_DIR / f'{slug}.h5ad'
    a.write_h5ad(out_path, compression='gzip')

    return {
        'slug': slug,
        'status': 'ok',
        'verdict': verdict,
        'decision': decision,
        'cells_loaded': n_loaded,
        'final_shape': list(a.shape),
        'path': str(out_path),
        'seconds': round(time.time() - t0, 1),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--slugs', nargs='*')
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    canonical, ensembl_to_hgnc = load_canonical()

    audit = pd.read_csv(AUDIT, sep='\t', dtype=str).fillna('')
    audit = audit[audit['verdict'].isin(['SAFE_REUSE', 'reuse_redo_qc'])].copy()
    # Dedup duplicate slug rows
    audit = audit.drop_duplicates(subset=['slug'])
    if args.slugs:
        audit = audit[audit['slug'].isin(args.slugs)]
    print(f'reuse candidates: {len(audit)}', flush=True)

    config = pd.read_csv(CONFIG, sep='\t', dtype=str).fillna('')
    slug_to_cfg = {r['slug']: r for _, r in config.iterrows() if r.get('slug')}

    results = []
    for i, row in audit.iterrows():
        slug = row['slug']
        cfg = slug_to_cfg.get(slug)
        if cfg is None:
            print(f'  {slug}: no config row', flush=True)
            results.append({'slug': slug, 'status': 'no_config_row'})
            continue
        print(f'\n[{len(results)+1}/{len(audit)}] {slug} ({row["verdict"]}, {cfg.get("decision","")})', flush=True)
        try:
            r = process_one(slug, row, cfg, canonical, ensembl_to_hgnc)
        except Exception as e:
            r = {'slug': slug, 'status': 'err', 'error': str(e)}
            traceback.print_exc()
        print(f'  -> {r.get("status")} {r.get("final_shape","")} {r.get("seconds","")}s', flush=True)
        results.append(r)

    # Update manifest with reused deposits
    if not args.dry_run:
        existing = {}
        if MANIFEST.exists():
            df = pd.read_csv(MANIFEST, sep='\t', dtype=str).fillna('')
            for _, r in df.iterrows():
                existing[r['slug']] = r.to_dict()
        for r in results:
            if r.get('status') != 'ok':
                continue
            slug = r['slug']
            cfg = slug_to_cfg.get(slug, {})
            existing[slug] = {
                'slug': slug,
                'accession': cfg.get('accession', slug.upper()),
                'path': r['path'],
                'n_cells': r['final_shape'][0],
                'n_genes': r['final_shape'][1],
                'n_samples': '',
                'n_control_samples': '',
                'n_control_cells': r['final_shape'][0],
                'organoid_type': cfg.get('organoid_type', ''),
                'filter': f'reuse_archive_{r["verdict"]}',
                'status': 'ok',
            }
        df_out = pd.DataFrame(list(existing.values()))
        df_out.to_csv(MANIFEST, sep='\t', index=False)
        print(f'\nmanifest updated: {len(df_out)} rows', flush=True)

    with open(LOG_PATH, 'w') as f:
        json.dump(results, f, indent=2, default=str)

    n_ok = sum(1 for r in results if r.get('status') == 'ok')
    print(f'\n===== REUSE SUMMARY =====')
    print(f'attempted: {len(results)}  ok: {n_ok}  failed: {len(results)-n_ok}')
    cells = sum(r['final_shape'][0] for r in results if r.get('status') == 'ok')
    print(f'cells reused: {cells:,}')


if __name__ == '__main__':
    main()
