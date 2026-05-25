#!/usr/bin/env python3
"""Review item #3: de-circularize Q2 correspondence.

The transfer labels (CellClass_cal) were assigned by latent proximity to Braun,
and the latent was built from the 2006 scANVI feature genes -> correlating
organoid vs Braun pseudobulk on those SAME genes inflates the diagonal. Here we
recompute correspondence on genes HELD OUT of the transfer feature set (shared
atlas∩Braun genes minus the 2006), with bootstrap CIs over cells. If the diagonal
still dominates on held-out genes, the correspondence is real, not circular.
"""
import argparse, sys, time
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad
import scipy.sparse as sp
sys.path.insert(0, str(Path(__file__).resolve().parent))
from atlas_common import (load_config, sym2ens, model_genes, ens_to_col,
                          read_atlas_genes, cp10k_log, write_sidecar)

ap = argparse.ArgumentParser(); ap.add_argument('--root', default=None); args = ap.parse_args()
cfg = load_config(root=args.root)
ROOT = cfg.root
N_SUB = 200_000
N_GENES = 4000          # random held-out genes for the pseudobulk correlation
B = 200                 # bootstrap resamples
CAP = 15_000            # cells/class cap for bootstrap speed
NULLISH = {'Unknown', 'nan', 'none', ''}
# on-feature self-corr (calibrated rerun, 2006 transfer genes) for side-by-side
ON_FEATURE = {'Neuron': 0.932, 'Neuroblast': 0.930, 'Glioblast': 0.900, 'Radial glia': 0.878,
              'Neuronal IPC': 0.910, 'Oligo': 0.883, 'Fibroblast': 0.841, 'Neural crest': 0.810,
              'Immune': 0.764, 'Vascular': 0.760, 'Placodes': 0.681, 'Erythrocyte': 0.659}
t0 = time.time()
def log(m): print(f"[{time.time()-t0:7.1f}s] {m}", flush=True)
rng = np.random.default_rng(0)

vn = set(model_genes(cfg.braun_scanvi_model))
s2e = sym2ens(cfg.canonical)

braun_all = ad.read_h5ad(cfg.braun)
atlas = ad.read_h5ad(cfg.atlas_full, backed='r')
e2c = ens_to_col(atlas.var_names, s2e)
shared = [g for g in braun_all.var_names if g in e2c]
held_out = sorted(set(shared) - vn)                 # genes NOT used by the transfer
sample = list(rng.choice(held_out, min(N_GENES, len(held_out)), replace=False))
log(f"shared={len(shared)} held_out={len(held_out)} -> sampled {len(sample)} held-out genes")

# Braun pseudobulk on held-out genes
braun = braun_all[:, sample]
bX = braun.X.toarray() if sp.issparse(braun.X) else np.asarray(braun.X)
bdf = pd.DataFrame(bX, columns=sample); bdf['__c'] = braun_all.obs['CellClass'].astype(str).values
b_pb = cp10k_log(bdf.groupby('__c').mean())
log(f"Braun pseudobulk (held-out): {b_pb.shape}")

# organoid pseudobulk on held-out genes, by calibrated label
cal = ad.read_h5ad(cfg.calibrated, backed='r')
idx = np.sort(rng.choice(atlas.n_obs, N_SUB, replace=False))
oX, _present = read_atlas_genes(atlas, sample, s2e, row_idx=idx)
obs_names = atlas.obs_names[idx]
oX = oX.toarray() if sp.issparse(oX) else np.asarray(oX)
lab = cal.obs['CellClass_cal'].reindex(obs_names).astype(str).to_numpy()
keep = ~pd.Series(lab).isin(NULLISH).to_numpy()
oX, lab = oX[keep], lab[keep]
log(f"organoid cells (labeled): {oX.shape[0]:,}")

classes = [c for c in b_pb.index if c in set(lab)]
rows = []
for c in classes:
    m = lab == c
    Xc = oX[m]
    if Xc.shape[0] > CAP:
        Xc = Xc[rng.choice(Xc.shape[0], CAP, replace=False)]
    pb_c = np.log1p(Xc.mean(0) / Xc.mean(0).sum() * 1e4)
    cors = {bc: np.corrcoef(pb_c, b_pb.loc[bc])[0, 1] for bc in b_pb.index}
    best = max(cors, key=cors.get)
    # bootstrap CI of self-corr
    boots = []
    n = Xc.shape[0]
    for _ in range(B):
        bs = Xc[rng.integers(0, n, n)]
        v = np.log1p(bs.mean(0) / bs.mean(0).sum() * 1e4)
        boots.append(np.corrcoef(v, b_pb.loc[c])[0, 1])
    lo, hi = np.percentile(boots, [2.5, 97.5])
    rows.append((c, n, round(cors[c], 3), round(lo, 3), round(hi, 3),
                 best, round(cors[best], 3), best == c, ON_FEATURE.get(c, np.nan)))

df = pd.DataFrame(rows, columns=['CellClass', 'n', 'heldout_selfcorr', 'ci_lo', 'ci_hi',
                                 'argmax', 'argmax_corr', 'diag_ok', 'on_feature_selfcorr'])
df['circular_inflation'] = (df['on_feature_selfcorr'] - df['heldout_selfcorr']).round(3)
log("=== Q2 correspondence on HELD-OUT genes (de-circularized) ===")
log("\n" + df.to_string(index=False))
log(f"diagonal dominance (held-out): {int(df['diag_ok'].sum())}/{len(df)}")
out_tsv = f'{ROOT}/data/q2_heldout_correspondence.tsv'
df.to_csv(out_tsv, sep='\t', index=False)
write_sidecar(out_tsv, __file__, {'N_SUB': N_SUB, 'N_GENES': N_GENES, 'B': B, 'CAP': CAP,
                                  'n_held_out_genes': len(held_out)})
log(f"saved -> {out_tsv} (+ provenance sidecar)")
log("DONE")
