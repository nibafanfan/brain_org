#!/usr/bin/env python3
"""Patch finalized per-GSM annotations into existing processed h5ads (obs only).

The finalized labels (brain_organoid_GSMannotations.xlsx, sheet 'GSM Annotations')
are per-GSM *metadata* — they don't change counts/QC/genes, so deposits whose
control-cell MEMBERSHIP is unchanged only need an obs refresh, not a step-1
reload (see scripts/membership_diff.py / data/membership_diff_2026-05-24.tsv).

This patches: MATCH deposits (per-cell by gsm) + the 8 INCLUDE_ALL_SAMPLES
series-level deposits (deposit-level, since they have no per-cell gsm).
SKIPS: the 10 MEMBERSHIP_CHANGED (full rebuild), gse168323 (pending decision),
gse297594_mecp2 (already fixed), atlas_v4_preprocessed (not a deposit).

Default is --dry-run (reports per-deposit changes, writes nothing).
Use --apply to rewrite obs in place. counts/var/layers are untouched.
Run with: /opt/homebrew/Caskroom/miniforge/base/bin/python
"""
import argparse
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

ROOT = Path('/Users/eg/brain_organoid')
XLSX = ROOT / 'data' / 'brain_organoid_GSMannotations.xlsx'
DIFF = ROOT / 'data' / 'membership_diff_2026-05-24.tsv'
PROC = ROOT / 'data' / 'processed'
GSM_RE = re.compile(r'GSM\d+')
ANNOT_TAG = 'GSMannotations_2026-05-24'

SKIP_SLUGS = {'gse168323', 'gse297594_mecp2', 'atlas_v4_preprocessed'}

# finalized sheet column -> obs column we write
FIELD_MAP = {
    'Organoid Type': 'organoid_type',
    'Protocol': 'protocol',
    'Cell Type Origin': 'cell_origin',
}
BOOL_FIELDS = {'Multi-Lineage?': 'multi_lineage', 'Vascularized?': 'vascularized', 'Slice?': 'slice'}


def norm_bool(v):
    s = str(v).strip().lower()
    if s in ('1', '1.0', 'yes', 'true', 'y'): return True
    if s in ('0', '0.0', 'no', 'false', 'n', '', 'nan'): return False
    return False


def norm_age(v):
    m = re.findall(r'\d+', str(v))
    return int(m[0]) if m else np.nan


def build_lookups():
    ann = pd.read_excel(XLSX, sheet_name='GSM Annotations')
    ann['GSM'] = ann['GSM'].astype(str).str.strip()
    ann['GEO'] = ann['GEO'].astype(str).str.strip().str.upper()
    ann = ann[ann['GSM'].str.fullmatch(GSM_RE)].copy()

    def row_to_fields(r):
        d = {}
        for src, dst in FIELD_MAP.items():
            d[dst] = str(r[src]).strip() if pd.notna(r[src]) else ''
        for src, dst in BOOL_FIELDS.items():
            d[dst] = norm_bool(r[src])
        d['organoid_age_days'] = norm_age(r['Age (Days in Vitro)'])
        d['guided'] = 'unguided' if norm_bool(r['Unguided?']) else 'guided'
        return d

    per_gsm = {r['GSM']: row_to_fields(r) for _, r in ann.iterrows()}

    # deposit-level aggregate: one value per accession (single value if uniform, else 'mixed:a|b')
    per_acc = {}
    fields = list(FIELD_MAP.values()) + list(BOOL_FIELDS.values()) + ['organoid_age_days', 'guided']
    for acc, grp in ann.groupby('GEO'):
        agg = {}
        recs = [row_to_fields(r) for _, r in grp.iterrows()]
        for f in fields:
            vals = [rec[f] for rec in recs]
            uniq = list(dict.fromkeys([v for v in vals if str(v) not in ('', 'nan', 'None')]))
            if len(uniq) == 1:
                agg[f] = uniq[0]
            elif f == 'organoid_age_days':
                agg[f] = np.nan  # ambiguous numeric -> leave blank
            elif len(uniq) == 0:
                agg[f] = ''
            else:
                agg[f] = 'mixed:' + '|'.join(map(str, uniq))
        per_acc[acc] = agg
    return per_gsm, per_acc, fields


def patch_one(path, per_gsm, per_acc, fields, apply):
    a = ad.read_h5ad(path)
    acc = str(a.obs['accession'].iloc[0]).upper() if 'accession' in a.obs and a.n_obs else ''
    gsms = a.obs['gsm'].astype(str) if 'gsm' in a.obs else pd.Series([''] * a.n_obs, index=a.obs_names)
    is_gsm = gsms.str.fullmatch(GSM_RE)
    mode = 'per_gsm' if is_gsm.mean() > 0.5 and any(g in per_gsm for g in gsms[is_gsm].unique()) else 'per_deposit'

    new_cols = {f: [] for f in fields}
    if mode == 'per_gsm':
        for g in gsms:
            rec = per_gsm.get(g, per_acc.get(acc, {}))
            for f in fields:
                new_cols[f].append(rec.get(f, '' if f != 'organoid_age_days' else np.nan))
    else:
        rec = per_acc.get(acc, {})
        for f in fields:
            new_cols[f] = [rec.get(f, '' if f != 'organoid_age_days' else np.nan)] * a.n_obs

    # summarize change
    summary = {}
    for f in fields:
        col = pd.Series(new_cols[f])
        if f == 'organoid_age_days':
            summary[f] = f"{np.nanmin(col.values):.0f}-{np.nanmax(col.values):.0f}" if col.notna().any() else 'NA'
        else:
            u = [x for x in col.unique() if str(x) not in ('', 'nan')]
            summary[f] = '|'.join(map(str, u[:3])) + ('…' if len(u) > 3 else '')

    if apply:
        for f in fields:
            a.obs[f] = new_cols[f]
        a.obs['annot_source'] = ANNOT_TAG
        # keep obs h5-safe (mirror sanitize: object cols -> str)
        for c in a.obs.columns:
            if a.obs[c].dtype == object or str(a.obs[c].dtype) == 'category':
                a.obs[c] = a.obs[c].astype(str).replace({'nan': '', 'None': ''})
        a.write_h5ad(path, compression='gzip')
    return mode, summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--apply', action='store_true', help='rewrite obs in place (default: dry-run)')
    ap.add_argument('--slugs', nargs='*', help='limit to these slugs')
    args = ap.parse_args()

    per_gsm, per_acc, fields = build_lookups()
    print(f"lookups: {len(per_gsm)} GSMs, {len(per_acc)} accessions", flush=True)

    diff = pd.read_csv(DIFF, sep='\t')
    targets = set(diff[diff.status.isin(['MATCH', 'NON_GSM_KEYED'])]['slug']) - SKIP_SLUGS
    if args.slugs:
        targets = set(args.slugs)
    print(f"{'APPLY' if args.apply else 'DRY-RUN'}: {len(targets)} deposits to patch\n", flush=True)

    for slug in sorted(targets):
        p = PROC / f'{slug}.h5ad'
        if not p.exists():
            print(f"  {slug}: MISSING file, skip", flush=True); continue
        mode, summ = patch_one(p, per_gsm, per_acc, fields, args.apply)
        print(f"  {slug:30} [{mode}] type={summ['organoid_type']!s:24} "
              f"age={summ['organoid_age_days']!s:8} multi={summ['multi_lineage']} "
              f"protocol={summ['protocol']!s:.30}", flush=True)
    print(f"\n{'WROTE obs in place.' if args.apply else 'dry-run only — rerun with --apply to write.'}", flush=True)


if __name__ == '__main__':
    main()
