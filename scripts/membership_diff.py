#!/usr/bin/env python3
"""Diff the finalized GSM annotations vs what's actually built in the atlas.

Decides, per deposit, whether the NEXT rebuild needs a step-1 reload (cell
membership changed) or just an obs-patch (only labels changed). Membership truth
= set of obs['gsm'] in each data/processed/<slug>.h5ad; finalized truth = the
'GSM Annotations' sheet of brain_organoid_GSMannotations.xlsx grouped by GEO.

Writes data/membership_diff_2026-05-24.tsv.
"""
import re
from pathlib import Path

import anndata as ad
import pandas as pd

ROOT = Path('/Users/eg/brain_organoid')
XLSX = ROOT / 'data' / 'brain_organoid_GSMannotations.xlsx'
PROC = ROOT / 'data' / 'processed'
OUT = ROOT / 'data' / 'membership_diff_2026-05-24.tsv'

GSM_RE = re.compile(r'GSM\d+')

# ---- finalized control GSMs per accession ----
ann = pd.read_excel(XLSX, sheet_name='GSM Annotations')
ann['GSM'] = ann['GSM'].astype(str).str.strip()
ann['GEO'] = ann['GEO'].astype(str).str.strip().str.upper()
ann = ann[ann['GSM'].str.match(GSM_RE)]
final_by_acc = ann.groupby('GEO')['GSM'].apply(lambda s: set(s)).to_dict()
print(f"finalized: {len(ann)} GSM rows across {len(final_by_acc)} accessions", flush=True)

# ---- as-built GSMs per processed deposit ----
rows = []
built_accs = set()
for p in sorted(PROC.glob('*.h5ad')):
    slug = p.stem
    a = ad.read_h5ad(p, backed='r')
    obs = a.obs
    acc = str(obs['accession'].iloc[0]).upper() if 'accession' in obs.columns and a.n_obs else ''
    built = set(str(g) for g in obs['gsm'].unique()) if 'gsm' in obs.columns else set()
    built_gsm = {g for g in built if GSM_RE.fullmatch(g)}
    a.file.close()
    built_accs.add(acc)
    fin = final_by_acc.get(acc, set())
    non_gsm_keyed = len(built_gsm) == 0  # series-level loaders key by GSE, can't diff by GSM
    removed = built_gsm - fin     # in atlas but no longer finalized -> drop these cells
    to_add = fin - built_gsm      # finalized but not in atlas -> add these cells
    if non_gsm_keyed:
        status = 'NON_GSM_KEYED'
    elif not fin:
        status = 'NOT_IN_FINALIZED'   # built but accession absent from finalized sheet
    elif not removed and not to_add:
        status = 'MATCH'             # obs-patch only
    else:
        status = 'MEMBERSHIP_CHANGED'  # step-1 rebuild needed
    rows.append({
        'slug': slug, 'accession': acc, 'status': status,
        'n_built': len(built_gsm), 'n_finalized': len(fin),
        'n_removed': len(removed), 'n_to_add': len(to_add),
        'removed': ','.join(sorted(removed))[:200],
        'to_add': ','.join(sorted(to_add))[:200],
    })

df = pd.DataFrame(rows).sort_values(['status', 'slug'])
df.to_csv(OUT, sep='\t', index=False)

# finalized accessions with NO built deposit at all = brand-new builds
new_accs = sorted(set(final_by_acc) - built_accs)

print("\n=== status counts (built deposits) ===", flush=True)
print(df['status'].value_counts().to_string(), flush=True)
print(f"\nfinalized accessions NOT built at all (new step-1 builds): {len(new_accs)}", flush=True)
print('  ' + ', '.join(new_accs), flush=True)
print(f"\nMEMBERSHIP_CHANGED deposits (need step-1 rebuild):", flush=True)
print(df[df.status=='MEMBERSHIP_CHANGED'][['slug','n_built','n_finalized','n_removed','n_to_add']].to_string(index=False), flush=True)
print(f"\nNOT_IN_FINALIZED (built but absent from finalized sheet — verify/drop?):", flush=True)
print(df[df.status=='NOT_IN_FINALIZED'][['slug','accession','n_built']].to_string(index=False), flush=True)
print(f"\nwrote {OUT}", flush=True)
