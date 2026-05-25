#!/usr/bin/env python3
"""B2: migrate per-deposit AnnDatas to HNOCA-compatible schema.

For each `data/processed/<slug>.h5ad`, produce `data/processed_v2/<slug>.h5ad`
with:
  - var.index = HGNC symbols matching HNOCA's 36,842-gene set
  - var['ensembl'], var['gene_symbol'], var['gene_length'], var['mt'],
    var['highly_variable'], var['highly_variable_nbatches'] from canonical
  - obs with HNOCA-style harmonization columns (cell_type/cell_type_original
    pattern + bio_sample, tech_sample, batch, individual)
  - X = raw integer counts (unchanged)
  - layers['counts_lengthnorm'] for Smart-seq2 deposits only

Skips deposits whose path doesn't exist or whose var can't be mapped.
"""
import argparse
import json
import re
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

_ENSG_RE = re.compile(r'ENSG\d{11}')


def normalize_ensembl(v):
    """Extract a clean ENSG\\d{11} ID from any decorated string.

    Handles: ENSG00000123456 (clean), ENSG00000123456.5 (versioned),
    hg19_ENSG00000123456 (build prefix), GRCh38_2020A_v5_ENSG00000123456
    (CellRanger ref prefix). Returns None if no ENSG ID present.
    """
    m = _ENSG_RE.search(str(v))
    return m.group(0) if m else None

ROOT = Path('/Users/eg/brain_organoid')
SRC_DIR = ROOT / 'data' / 'processed'
DST_DIR = ROOT / 'data' / 'processed_v2'
MANIFEST = ROOT / 'data' / 'manifest.tsv'
MANIFEST_V2 = ROOT / 'data' / 'manifest_v2.tsv'
CANONICAL = ROOT / 'data' / 'reference' / 'hnoca_var_canonical.tsv'

# Smart-seq2 deposits get length-normalized layer; all others are 10x
SMART_SEQ2 = {'gse195692', 'gse185052', 'gse75140', 'gse124299', 'gse115011_pollen'}


def load_canonical():
    df = pd.read_csv(CANONICAL, sep='\t', index_col=0)
    df.index.name = 'hgnc_symbol'
    # Build ensembl → hgnc map
    ensembl_to_hgnc = {}
    for hgnc, row in df.iterrows():
        eid = row['ensembl']
        if isinstance(eid, str) and eid.startswith('ENSG'):
            ensembl_to_hgnc[eid] = hgnc
    return df, ensembl_to_hgnc


def detect_namespace(var_index):
    """Decide if var_names look like Ensembl IDs (decorated or not) vs HGNC symbols.

    Inclusive: matches any string containing ENSG\\d{11} so prefixed/versioned
    deposits get classified as Ensembl.
    """
    sample = [str(v) for v in var_index[:50]]
    n_ensg = sum(1 for v in sample if _ENSG_RE.search(v))
    return 'ensembl' if n_ensg > len(sample) / 2 else 'symbol'


def map_to_hgnc(adata, canonical, ensembl_to_hgnc):
    """Reindex adata.var to HNOCA's HGNC namespace. Drops unmappable genes."""
    ns = detect_namespace(adata.var_names)
    var_names = adata.var_names.astype(str)
    if ns == 'ensembl':
        # Normalize: strip any prefix/version → bare ENSG\d{11}, then map ENSG → HGNC
        hgnc_aligned = pd.Series(
            [ensembl_to_hgnc.get(normalize_ensembl(v)) for v in var_names],
            index=var_names,
        )
    else:
        # Already symbols — pass through; only keep symbols present in HNOCA
        hgnc_aligned = pd.Series([v if v in canonical.index else None for v in var_names], index=var_names)
    keep_mask = hgnc_aligned.notna().to_numpy()
    n_in = adata.shape[1]
    n_kept = int(keep_mask.sum())
    if n_kept == 0:
        return None, {'in': n_in, 'kept': 0, 'ns': ns, 'reason': 'no genes mapped to HNOCA'}
    # Slice to mappable genes; HGNC names may be duplicated (multiple Ensembl → same HGNC).
    # Drop duplicates by HGNC BEFORE make_unique so we keep only canonical-resolvable names.
    new_names = hgnc_aligned[keep_mask].values
    is_first = ~pd.Series(new_names).duplicated(keep='first').to_numpy()
    a2 = adata[:, keep_mask][:, is_first].copy()
    a2.var_names = pd.Index(new_names[is_first])
    return a2, {'in': n_in, 'kept': a2.shape[1], 'ns': ns}


def attach_canonical_var(adata, canonical):
    """Add HNOCA-style var columns (ensembl, gene_symbol, gene_length, mt, hv flags)."""
    keep_cols = ['ensembl', 'gene_symbol', 'gene_length', 'mt',
                 'highly_variable', 'highly_variable_rank', 'highly_variable_nbatches']
    canonical_sub = canonical[[c for c in keep_cols if c in canonical.columns]]
    aligned = canonical_sub.reindex(adata.var_names)
    for c in aligned.columns:
        adata.var[c] = aligned[c].values
    return adata


def migrate_obs(adata, slug):
    """Add HNOCA-style harmonization columns with _original preservation."""
    obs = adata.obs.copy()

    # Add _original versions of fields HNOCA harmonizes
    pairs = [
        ('cell_type', 'unknown'),
        ('cell_line', obs.get('cell_line', pd.Series(['unknown'] * len(obs), index=obs.index)).astype(str).iloc[0] if 'cell_line' in obs.columns else 'unknown'),
        ('disease', 'healthy'),
        ('assay_sc', 'Smart-seq2' if slug in SMART_SEQ2 else "10x 3' v3"),
        ('organ', 'brain'),
        ('organism', 'Homo sapiens'),
        ('sex', 'unknown'),
        ('ethnicity', 'unknown'),
        ('suspension_type', 'cell'),
    ]
    for col, default in pairs:
        orig_col = f'{col}_original'
        if col in obs.columns:
            obs[orig_col] = obs[col].astype(str)
        else:
            obs[orig_col] = str(default)
            obs[col] = str(default)

    # Sample hierarchy
    obs['bio_sample'] = obs.get('gsm', pd.Series(['pool'] * len(obs), index=obs.index)).astype(str)
    obs['tech_sample'] = (str(slug) + '_' + obs['bio_sample']).astype(str)
    obs['batch'] = obs.get('batch', pd.Series([str(slug)] * len(obs), index=obs.index)).astype(str)
    obs['individual'] = obs.get('cell_line', pd.Series(['unknown'] * len(obs), index=obs.index)).astype(str)

    # Other HNOCA columns
    if 'organoid_type' in obs.columns:
        obs['assay_differentiation'] = obs['organoid_type'].astype(str)
    obs['organoid_age_days'] = obs.get('protocol_age_days', pd.Series([np.nan] * len(obs), index=obs.index))
    obs['publication'] = obs.get('accession', pd.Series([slug.upper()] * len(obs), index=obs.index)).astype(str)
    obs['sample_source'] = 'organoid'
    obs['state_exact'] = 'unknown'

    # Disease: refine from existing genotype/diagnosis columns
    if 'genotype' in obs.columns:
        gt = obs['genotype'].astype(str).str.lower()
        is_wt = gt.str.contains('wt|wild|isogenic|ctrl|corrected|parental', regex=True, na=False)
        obs['disease'] = np.where(is_wt, 'healthy', obs['genotype'].astype(str))
    elif 'diagnosis' in obs.columns:
        obs['disease'] = obs['diagnosis'].astype(str)

    adata.obs = obs
    return adata


def sanitize_for_h5(adata):
    """Coerce object/category obs and var columns to str so h5py vlen-string can write them.

    HDF5 fails on mixed-type object columns (NaN + strings) with
    'Can\\'t implicitly convert non-string objects to strings'.
    Applies to both adata.obs and adata.var.
    """
    for df_name in ('obs', 'var'):
        df = getattr(adata, df_name)
        for c in list(df.columns):
            ser = df[c]
            if ser.dtype == object or str(ser.dtype) == 'category':
                df[c] = ser.astype(str).fillna('').replace({'nan': '', 'None': ''})
    return adata


# Backwards-compatible alias
sanitize_obs_for_h5 = sanitize_for_h5


def add_lengthnorm_layer(adata):
    """For Smart-seq2 deposits: counts_lengthnorm = X / gene_length * mean(gene_length)."""
    if 'gene_length' not in adata.var.columns:
        return adata
    lengths = adata.var['gene_length'].astype(float).values
    if np.isnan(lengths).any():
        # use median for NaN
        median = np.nanmedian(lengths)
        lengths = np.where(np.isnan(lengths), median, lengths)
    mean_len = np.mean(lengths)
    # X is cells × genes; divide each column by gene length, multiply by mean
    X = adata.X
    if sp.issparse(X):
        # Build scaling diag matrix
        scale = mean_len / lengths
        scale_diag = sp.diags(scale, format='csr')
        adata.layers['counts_lengthnorm'] = X @ scale_diag
    else:
        scale = mean_len / lengths
        adata.layers['counts_lengthnorm'] = X * scale[np.newaxis, :]
    return adata


def migrate_one(slug, src_path, canonical, ensembl_to_hgnc, dry_run=False):
    a = ad.read_h5ad(src_path)
    a2, stats = map_to_hgnc(a, canonical, ensembl_to_hgnc)
    if a2 is None:
        return None, {'status': 'fail', **stats}
    a2 = attach_canonical_var(a2, canonical)
    a2 = migrate_obs(a2, slug)
    a2 = sanitize_obs_for_h5(a2)
    if slug in SMART_SEQ2:
        a2 = add_lengthnorm_layer(a2)
    return a2, {'status': 'ok',
                'shape_in': a.shape, 'shape_out': a2.shape,
                'genes_kept': stats['kept'], 'genes_in': stats['in'],
                'namespace': stats['ns'],
                'has_lengthnorm': 'counts_lengthnorm' in a2.layers}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dry-run', action='store_true', help="Don't write output files")
    ap.add_argument('--slugs', nargs='*', help='Limit to specific slugs (for testing)')
    args = ap.parse_args()

    DST_DIR.mkdir(parents=True, exist_ok=True)
    canonical, ensembl_to_hgnc = load_canonical()
    print(f'canonical: {canonical.shape[0]} genes; {len(ensembl_to_hgnc)} ensembl→hgnc mappings', flush=True)

    manifest = pd.read_csv(MANIFEST, sep='\t')
    if args.slugs:
        manifest = manifest[manifest['slug'].isin(args.slugs)]
    print(f'manifest rows: {len(manifest)}', flush=True)

    results = []
    t0 = time.time()
    for i, row in manifest.iterrows():
        slug = row['slug']
        src = Path(row['path'])
        if not src.exists():
            print(f'  [{i+1}/{len(manifest)}] {slug}: SKIP (missing src {src})', flush=True)
            results.append({'slug': slug, 'status': 'missing_src'})
            continue
        try:
            a2, stats = migrate_one(slug, src, canonical, ensembl_to_hgnc, dry_run=args.dry_run)
            if a2 is None:
                print(f'  [{i+1}/{len(manifest)}] {slug}: FAIL {stats}', flush=True)
                results.append({'slug': slug, **stats})
                continue
            if not args.dry_run:
                dst = DST_DIR / f'{slug}.h5ad'
                a2.write_h5ad(dst, compression='gzip')
                stats['out'] = str(dst)
            print(f'  [{i+1}/{len(manifest)}] {slug}: {stats["shape_in"]} → {stats["shape_out"]} ({stats["genes_kept"]:,}/{stats["genes_in"]:,} genes, {stats["namespace"]}, lengthnorm={stats["has_lengthnorm"]})', flush=True)
            results.append({'slug': slug, **stats})
        except Exception as e:
            import traceback
            print(f'  [{i+1}/{len(manifest)}] {slug}: ERR {e}', flush=True)
            results.append({'slug': slug, 'status': 'err', 'error': str(e)})

    # Build manifest_v2.tsv (merge with existing rows so --slugs runs don't clobber)
    if not args.dry_run:
        # Load existing manifest_v2 if present (rows from prior runs)
        existing_rows = {}
        if MANIFEST_V2.exists():
            with open(MANIFEST_V2) as f:
                header_line = f.readline()
                for line in f:
                    cols = line.rstrip('\n').split('\t')
                    if cols and cols[0]:
                        existing_rows[cols[0]] = line.rstrip('\n')

        # Compose updated rows for slugs migrated in this run
        for r in results:
            if r.get('status') != 'ok':
                continue
            slug = r['slug']
            try:
                v2 = ad.read_h5ad(DST_DIR / f'{slug}.h5ad', backed='r')
                acc = v2.obs['accession'].iloc[0] if 'accession' in v2.obs.columns else slug.upper()
                n_cells = v2.shape[0]; n_genes = v2.shape[1]
                n_samples = int(v2.obs['sample_id'].nunique()) if 'sample_id' in v2.obs.columns else 1
                n_ctrl = int(v2.obs['is_control'].sum()) if 'is_control' in v2.obs.columns else 0
                n_ctrl_s = int(v2.obs[v2.obs['is_control']]['sample_id'].nunique()) if 'is_control' in v2.obs.columns and 'sample_id' in v2.obs.columns else 0
                ot = v2.obs['organoid_type'].iloc[0] if 'organoid_type' in v2.obs.columns else 'unknown'
                v2.file.close()
                existing_rows[slug] = '\t'.join(map(str, [slug, acc, str(DST_DIR / f'{slug}.h5ad'),
                                                          n_cells, n_genes, n_samples, n_ctrl_s, n_ctrl,
                                                          ot, 'hnoca-schema-v2', 'ok']))
            except Exception as e:
                print(f'  manifest_v2 row failed for {slug}: {e}', flush=True)

        # Sort by slug and write
        with open(MANIFEST_V2, 'w') as f:
            f.write('slug\taccession\tpath\tn_cells\tn_genes\tn_samples\tn_control_samples\tn_control_cells\torganoid_type\tfilter\tstatus\n')
            for slug in sorted(existing_rows):
                f.write(existing_rows[slug] + '\n')
        print(f'\nWROTE {MANIFEST_V2}: {len(existing_rows)} rows total', flush=True)

    # Summary
    elapsed = time.time() - t0
    n_ok = sum(1 for r in results if r.get('status') == 'ok')
    n_fail = sum(1 for r in results if r.get('status') not in ('ok', None))
    print(f'\ndone in {elapsed:.1f}s: {n_ok} ok, {n_fail} failed/skipped', flush=True)

    # Stats: dropped-gene report
    if results:
        gene_keep_pct = [100 * r['genes_kept']/max(1, r['genes_in']) for r in results if r.get('status') == 'ok']
        if gene_keep_pct:
            print(f'gene retention: mean {np.mean(gene_keep_pct):.1f}%, min {min(gene_keep_pct):.1f}%, max {max(gene_keep_pct):.1f}%', flush=True)

    # Write log
    log_path = Path('/tmp/migrate_log.json')
    with open(log_path, 'w') as f:
        json.dump([{k: (str(v) if isinstance(v, tuple) else v) for k, v in r.items()} for r in results], f, indent=2)
    print(f'log: {log_path}', flush=True)


if __name__ == '__main__':
    main()
