#!/usr/bin/env python3
"""Inject the 8 finalized annotation fields onto every concat temp, keyed by GSM.

Source of truth = brain_organoid_GSMannotations.xlsx 'GSM Annotations' sheet
(100% complete). Each cell gets every field, so the columns are universal across
all temps and survive concat_on_disk's column intersection.

Per-cell gsm is read from obs['gsm']; gse251684_striato_nigral has no gsm column
but its obs['sample_id'] (STR / STR_SN) maps to GSM groups with identical
annotations, so we assign a representative GSM per group.

Usage:
  scripts/inject_finalized_annotations.py --dry-run   # report match rates only
  scripts/inject_finalized_annotations.py             # write columns into temps
"""
import argparse, glob, os, re, sys
import numpy as np, pandas as pd, h5py
from anndata._io.specs import read_elem, write_elem

ROOT = '/Users/eg/brain_organoid'
SHEET = f'{ROOT}/data/brain_organoid_GSMannotations.xlsx'
TMPDIR = f'{ROOT}/data/_concat_tmp'

# sheet column -> obs column
FIELD_MAP = {
    'Cell Type Origin': 'cell_type_origin',
    'Age (Days in Vitro)': 'age_days',
    'Organoid Type': 'organoid_type',
    'Protocol': 'protocol',
    'Unguided?': 'unguided',
    'Multi-Lineage?': 'multi_lineage',
    'Vascularized?': 'vascularized',
    'Slice?': 'slice',
}
OBS_FIELDS = list(FIELD_MAP.values())

# gse251684: sample_id group -> representative GSM (annotations identical within group)
SAMPLEID_GSM = {
    'gse251684_striato_nigral': {'STR': 'GSM7985991', 'STR_SN': 'GSM7985993'},
}


def load_gsm_annotations():
    df = pd.read_excel(SHEET, sheet_name='GSM Annotations')
    ann = {}
    for _, r in df.iterrows():
        g = str(r['GSM']).strip().upper()
        if not g.startswith('GSM'):
            continue
        ann[g] = {obs: ('' if pd.isna(r[sheet]) else str(r[sheet]).strip())
                  for sheet, obs in FIELD_MAP.items()}
    return ann


def read_obs_col(f, col):
    if col not in f['obs']:
        return None
    n = f[f'obs/{col}']
    if isinstance(n, h5py.Group):
        cats = [x.decode() if isinstance(x, bytes) else str(x) for x in n['categories'][:]]
        codes = n['codes'][:]
        return np.array([cats[i] if 0 <= i < len(cats) else '' for i in codes], dtype=object)
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in n[:]], dtype=object)


def per_cell_gsm(f, slug, n):
    g = read_obs_col(f, 'gsm')
    if g is not None:
        return np.array([str(x).strip().upper() for x in g], dtype=object)
    # fallback: sample_id -> representative GSM
    if slug in SAMPLEID_GSM:
        sid = read_obs_col(f, 'sample_id')
        mp = SAMPLEID_GSM[slug]
        return np.array([mp.get(str(s), '') for s in sid], dtype=object)
    return np.full(n, '', dtype=object)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dry-run', action='store_true')
    args = ap.parse_args()

    ann = load_gsm_annotations()
    print(f'finalized GSM annotations loaded: {len(ann)} GSMs', flush=True)

    temps = sorted(glob.glob(f'{TMPDIR}/*.h5ad'))
    rows = []
    for p in temps:
        slug = os.path.basename(p).replace('.h5ad', '')
        with h5py.File(p, 'r' if args.dry_run else 'a') as f:
            n = f['X'].attrs.get('shape', None)
            n = int(n[0]) if n is not None else read_obs_col(f, 'n_counts').shape[0]
            gsm = per_cell_gsm(f, slug, n)
            matched = np.array([g in ann for g in gsm])
            rate = matched.mean() if n else 0.0
            rows.append((slug, n, int(matched.sum()), round(100 * rate, 1)))
            if args.dry_run:
                continue
            # build each field column with fill-don't-clobber semantics:
            #   - cell's gsm matches a finalized GSM -> authoritative sheet value
            #   - else -> keep the deposit's existing obs value (series-level
            #     deposits + gse297594 carry valid deposit-level annotation);
            #     'unknown' only if no existing value exists.
            matched_mask = np.array([g in ann for g in gsm])
            cols = {}
            for obs in OBS_FIELDS:
                existing = read_obs_col(f, obs)
                out = np.empty(n, dtype=object)
                for i, g in enumerate(gsm):
                    if g in ann and ann[g][obs] != '':
                        out[i] = ann[g][obs]
                    elif existing is not None and str(existing[i]) not in ('', 'nan', 'None'):
                        out[i] = str(existing[i])
                    else:
                        out[i] = 'unknown'
                cols[obs] = out
            # provenance flag: 'gsm' = authoritative per-GSM sheet value,
            # 'deposit' = fell back to existing deposit-level annotation.
            cols['annotation_level'] = np.where(matched_mask, 'gsm', 'deposit').astype(object)
            # also (re)write a clean per-cell gsm + the provenance flag
            write_names = ['gsm'] + OBS_FIELDS + ['annotation_level']
            arrs = {'gsm': np.array([g if g else 'unknown' for g in gsm], object),
                    'annotation_level': cols['annotation_level'],
                    **{obs: cols[obs] for obs in OBS_FIELDS}}
            for name in write_names:
                if f'obs/{name}' in f:
                    del f[f'obs/{name}']
                write_elem(f, f'obs/{name}', pd.Categorical(pd.Series(arrs[name]).astype(str)))
            # ensure column-order attr includes the new names
            co = f['obs'].attrs.get('column-order')
            if co is not None:
                names = [c.decode() if isinstance(c, bytes) else c for c in co]
                for name in write_names:
                    if name not in names:
                        names.append(name)
                f['obs'].attrs['column-order'] = np.array(names, dtype=h5py.special_dtype(vlen=str))

    rep = pd.DataFrame(rows, columns=['slug', 'cells', 'matched', 'pct'])
    bad = rep[rep['pct'] < 99.0].sort_values('pct')
    print(f"\n{'deposit':34}{'cells':>10}{'matched':>10}{'pct':>7}")
    for _, r in bad.iterrows():
        print(f"{r['slug']:34}{r['cells']:>10}{r['matched']:>10}{r['pct']:>7}")
    print(f"\n{len(rep)} deposits | {len(bad)} with <99% GSM match")
    print(f"total cells: {rep['cells'].sum():,} | matched: {rep['matched'].sum():,} "
          f"({100*rep['matched'].sum()/rep['cells'].sum():.1f}%)")
    if not args.dry_run:
        print('\nWROTE annotation columns into all temps.')


if __name__ == '__main__':
    sys.exit(main())
