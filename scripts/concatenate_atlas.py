#!/usr/bin/env python3
"""Concatenate the 113-deposit atlas into a single AnnData (atlas step 2).

Strategy (memory-safe):
  1. Project each deposit's X onto the HNOCA canonical 36,842-gene set via a
     sparse projection matrix (X @ P). Genes not in canonical are dropped;
     canonical genes absent in a deposit are zero-filled. Duplicate source
     symbols mapping to one canonical gene are summed.
  2. Write each reindexed deposit to a temp h5ad with *identical* var (the full
     canonical var table) and slug-prefixed, unique obs_names.
  3. concat_on_disk the temp files into one atlas h5ad without a giant
     in-memory vstack.

Resumable: existing temp files are reused. Use --limit N to test on N deposits.
"""
import argparse
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path("/Users/eg/brain_organoid")
CANON_TSV = ROOT / "data/reference/hnoca_var_canonical.tsv"
MANIFEST = ROOT / "data/manifest.tsv"
TMPDIR = ROOT / "data/_concat_tmp"


def load_canonical():
    canon = pd.read_csv(CANON_TSV, sep="\t")
    canon["hgnc_symbol"] = canon["hgnc_symbol"].astype(str)
    canon = canon.set_index("hgnc_symbol")
    canon.index.name = None
    return canon


def project_deposit(a, canon_idx, canon_var):
    """Return AnnData reindexed onto the canonical gene set."""
    src = a.var.index.astype(str)
    rows, cols = [], []
    for i, g in enumerate(src):
        j = canon_idx.get(g)
        if j is not None:
            rows.append(i)
            cols.append(j)
    NG = len(canon_var)
    P = sparse.csr_matrix(
        (np.ones(len(rows), dtype=np.float32), (rows, cols)),
        shape=(len(src), NG),
    )
    X = a.X
    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)
    newX = (X.tocsr().astype(np.float32) @ P).tocsr()
    obs = a.obs.copy()
    # multi_lineage is per-cell and arrives as bool in some deposits, categorical
    # in others; concat_on_disk can't merge mixed dtypes. Coerce to string-
    # categorical (like organoid_type) so all temps share one dtype.
    if "multi_lineage" in obs.columns:
        obs["multi_lineage"] = pd.Categorical(obs["multi_lineage"].astype(str))
    nv = ad.AnnData(X=newX, obs=obs, var=canon_var.copy())
    return nv, len(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None,
                    help="only process first N deposits (for testing)")
    ap.add_argument("--out", default=str(ROOT / "data/atlas_v4_full.h5ad"))
    ap.add_argument("--reindex-only", action="store_true",
                    help="build temp files but skip final concat")
    args = ap.parse_args()

    TMPDIR.mkdir(parents=True, exist_ok=True)
    canon_var = load_canonical()
    canon_idx = {g: i for i, g in enumerate(canon_var.index)}
    manifest = pd.read_csv(MANIFEST, sep="\t")
    if args.limit:
        manifest = manifest.head(args.limit)

    tmp_files = []
    n = len(manifest)
    for k, (_, r) in enumerate(manifest.iterrows(), 1):
        slug = r["slug"]
        out = TMPDIR / f"{slug}.h5ad"
        if out.exists():
            tmp_files.append(str(out))
            print(f"[{k}/{n}] {slug}: temp exists, skip", flush=True)
            continue
        a = ad.read_h5ad(r["path"])
        nv, matched = project_deposit(a, canon_idx, canon_var)
        nv.obs_names = slug + "_" + nv.obs_names.astype(str)
        nv.write_h5ad(out)
        tmp_files.append(str(out))
        print(f"[{k}/{n}] {slug}: {a.shape[0]} cells, "
              f"{matched}/{a.shape[1]} genes matched canonical -> {out.name}",
              flush=True)
        del a, nv

    if args.reindex_only:
        print(f"reindex-only done: {len(tmp_files)} temp files in {TMPDIR}")
        return

    from anndata.experimental import concat_on_disk
    print(f"concat_on_disk: {len(tmp_files)} files -> {args.out}", flush=True)
    # merge=None: var index preserved (identical across files); merging the var
    # annotation columns triggers an anndata 0.10.9 writer bug, so we re-attach
    # the canonical var metadata directly afterward.
    concat_on_disk(tmp_files, args.out, axis=0, join="inner", merge=None)

    # Guard: the concat must yield the full canonical gene axis. A 0-gene (or
    # short) result means the join collapsed (mismatched temp var) -> abort BEFORE
    # the destructive re-attach below, so we never overwrite var on a bad concat.
    # (An interrupted re-attach is what produced the earlier 0-gene atlas_v5.)
    import h5py
    from anndata._io.specs import read_elem, write_elem
    NG = len(canon_var)
    with h5py.File(args.out, "r") as f:
        var_index = pd.Index(read_elem(f["var"]).index.astype(str))
    if len(var_index) != NG:
        raise SystemExit(
            f"ABORT: concat produced {len(var_index)} genes, expected {NG}. "
            f"Temp files' var is not aligned; leaving {args.out} untouched.")

    # Re-attach canonical var metadata atomically: write under a temp key, verify,
    # then swap. If the process dies mid-write, the original var stays intact.
    print("re-attaching canonical var metadata...", flush=True)
    new_var = canon_var.loc[var_index]
    assert len(new_var) == NG, f"canon_var.loc lost rows: {len(new_var)} != {NG}"
    with h5py.File(args.out, "a") as f:
        if "var_new" in f:
            del f["var_new"]
        write_elem(f, "var_new", new_var)
        assert len(pd.Index(read_elem(f["var_new"]).index)) == NG, "var_new readback bad"
        del f["var"]
        f.move("var_new", "var")

    # Final verification of the written file.
    a = ad.read_h5ad(args.out, backed="r")
    n_cells, n_genes = a.shape
    a.file.close()
    assert n_genes == NG, f"final atlas has {n_genes} genes, expected {NG}"
    print(f"VERIFIED: {n_cells} cells x {n_genes} genes", flush=True)
    print("done.", flush=True)


if __name__ == "__main__":
    sys.exit(main())
