#!/usr/bin/env python3
"""Recover per-cell finalized annotations for deposits whose combined matrix
encodes the sample in the cell barcode (preserved in the atlas obs index).

Handles 3 deposits the GSM-join missed because their obs['gsm'] was series-level:
  gse281535_brainstem : obs index has timepoint tag D20-N..D60-N -> sheet GSM Label
  gse227640           : obs index has Term/Line/Condition/Week -> sheet GSM Label
  gse195690           : single sample -> apply the one sheet row uniformly

For matched cells, writes the 8 finalized fields + gsm from the sheet and sets
annotation_level='gsm'. Run with --dry-run to validate mapping coverage first.
"""
import argparse, re, sys
import pandas as pd, numpy as np, h5py
from anndata._io.specs import read_elem, write_elem

SHEET = '/Users/eg/brain_organoid/data/brain_organoid_GSMannotations.xlsx'
TMP = '/Users/eg/brain_organoid/data/_concat_tmp'
FIELD_MAP = {'Cell Type Origin': 'cell_type_origin', 'Age (Days in Vitro)': 'age_days',
             'Organoid Type': 'organoid_type', 'Protocol': 'protocol', 'Unguided?': 'unguided',
             'Multi-Lineage?': 'multi_lineage', 'Vascularized?': 'vascularized', 'Slice?': 'slice'}
OBS_FIELDS = list(FIELD_MAP.values())


def row_to_ann(r):
    d = {obs: ('unknown' if pd.isna(r[s]) else str(r[s]).strip()) for s, obs in FIELD_MAP.items()}
    d['gsm'] = str(r['GSM']).strip().upper()
    return d


def acc_rows(df, accs):
    pat = '|'.join(accs)
    return df[df['GEO'].astype(str).str.contains(pat, case=False, na=False)]


def norm_factors(s):
    """(term, line, cond, week) from either obs tag or sheet GSM Label of GSE227640."""
    low = s.lower()
    term = 'long' if 'long' in low else 'short'
    line = 'H1' if re.search(r'h1', low) else ('L13234' if '13234' in low else '?')
    cond = 'LIF' if 'lif' in low else 'control'
    m = re.search(r'wk?(?:eek)?(\d+)', low)
    week = m.group(1) if m else '?'
    return (term, line, cond, week)


def get_index(f):
    i = f['obs'].attrs.get('_index', '_index'); i = i.decode() if isinstance(i, bytes) else i
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in f[f'obs/{i}'][:]])


def build_keys(slug, idx, df):
    """Return per-cell annotation dict-or-None array for the given deposit."""
    n = len(idx)
    if slug == 'gse281535_brainstem':
        rows = acc_rows(df, ['GSE281535'])
        lut = {str(r['GSM Label']).strip(): row_to_ann(r) for _, r in rows.iterrows()}
        out = []
        for x in idx:
            m = re.search(r'__(D\d+-N)_', x)
            out.append(lut.get(m.group(1)) if m else None)
        return out
    if slug == 'gse227640':
        rows = acc_rows(df, ['GSE227640'])
        lut = {norm_factors(str(r['GSM Label'])): row_to_ann(r) for _, r in rows.iterrows()}
        out = []
        for x in idx:
            m = re.match(r'gse227640_gse227640_(.+?)_[ACGT]{8,}$', x)
            out.append(lut.get(norm_factors(m.group(1))) if m else None)
        return out
    if slug == 'gse195690':
        rows = acc_rows(df, ['GSE195690'])
        ann = row_to_ann(rows.iloc[0])
        return [ann] * n
    raise ValueError(slug)


def read_col(f, c):
    if c not in f['obs']: return None
    nd = f[f'obs/{c}']
    if isinstance(nd, h5py.Group):
        cats = [x.decode() if isinstance(x, bytes) else str(x) for x in nd['categories'][:]]
        codes = nd['codes'][:]
        return np.array([cats[i] if 0 <= i < len(cats) else 'unknown' for i in codes], dtype=object)
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in nd[:]], dtype=object)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--dry-run', action='store_true')
    args = ap.parse_args()
    df = pd.read_excel(SHEET, sheet_name='GSM Annotations')
    for slug in ['gse281535_brainstem', 'gse227640', 'gse195690']:
        p = f'{TMP}/{slug}.h5ad'
        with h5py.File(p, 'r' if args.dry_run else 'a') as f:
            idx = get_index(f)
            anns = build_keys(slug, idx, df)
            matched = sum(a is not None for a in anns)
            print(f"{slug:22} {len(idx):>7} cells | mapped {matched} ({100*matched/len(idx):.1f}%)")
            if args.dry_run:
                continue
            cur = {c: read_col(f, c) for c in OBS_FIELDS + ['gsm', 'annotation_level']}
            for c in OBS_FIELDS + ['gsm', 'annotation_level']:
                arr = cur[c].copy() if cur[c] is not None else np.array(['unknown'] * len(idx), dtype=object)
                for i, a in enumerate(anns):
                    if a is None:
                        continue
                    arr[i] = a.get(c) if c in a else ('gsm' if c == 'annotation_level' else arr[i])
                if c == 'annotation_level':
                    arr = np.array(['gsm' if anns[i] is not None else arr[i] for i in range(len(idx))], dtype=object)
                if f'obs/{c}' in f:
                    del f[f'obs/{c}']
                write_elem(f, f'obs/{c}', pd.Categorical(pd.Series(arr).astype(str)))
    if not args.dry_run:
        print("wrote per-cell annotations into the 3 temps.")


if __name__ == '__main__':
    sys.exit(main())
