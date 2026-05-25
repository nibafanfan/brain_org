#!/usr/bin/env python3
"""Build Fig.1-style integrated atlas panels from standardized inputs.

Outputs deterministic TSV/PNG/PDF artifacts plus a JSON manifest.
"""

from __future__ import annotations

import argparse
import json
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import anndata as ad
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


DEFAULT_SEED = 0


@dataclass
class RunConfig:
    h5ad: str
    outdir: str
    prefix: str
    seed: int
    max_cells_plot: int
    age_col: str
    protocol_col: str
    dataset_col: str
    sample_col: str
    class_col: str
    region_col: str
    umap_key: str
    latent_key: str
    add_advanced: bool
    marker_genes: list[str]


def _first_existing(frame: pd.DataFrame, candidates: Iterable[str], required: bool = True) -> str | None:
    for c in candidates:
        if c in frame.columns:
            return c
    if required:
        raise KeyError(f"None of columns found: {list(candidates)}")
    return None


def _ensure_umap(a: ad.AnnData, seed: int, latent_key: str, umap_key: str) -> None:
    if umap_key in a.obsm:
        return
    import scanpy as sc

    if latent_key in a.obsm:
        use_rep = latent_key
    elif "X_scvi" in a.obsm:
        use_rep = "X_scvi"
    else:
        raise KeyError(f"Missing both obsm['{umap_key}'] and latent key '{latent_key}'")
    sc.pp.neighbors(a, use_rep=use_rep, random_state=seed)
    sc.tl.umap(a, random_state=seed)


def _write(df: pd.DataFrame, p: Path) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(p, sep="\t", index=False)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--outdir", default="data/figure1")
    ap.add_argument("--prefix", default="figure1")
    ap.add_argument("--seed", type=int, default=DEFAULT_SEED)
    ap.add_argument("--max-cells-plot", type=int, default=250000)
    ap.add_argument("--age-col", default="organoid_age_days")
    ap.add_argument("--protocol-col", default="protocol")
    ap.add_argument("--dataset-col", default="dataset_slug")
    ap.add_argument("--sample-col", default="bio_sample")
    ap.add_argument("--class-col", default="CellClass_cal")
    ap.add_argument("--region-col", default="organoid_type")
    ap.add_argument("--umap-key", default="X_umap")
    ap.add_argument("--latent-key", default="X_scvi")
    ap.add_argument("--add-advanced", action="store_true")
    ap.add_argument("--marker-genes", nargs="*", default=["SOX2", "EOMES", "DCX", "NEUROD2"])
    args = ap.parse_args()

    cfg = RunConfig(
        h5ad=args.h5ad,
        outdir=args.outdir,
        prefix=args.prefix,
        seed=args.seed,
        max_cells_plot=args.max_cells_plot,
        age_col=args.age_col,
        protocol_col=args.protocol_col,
        dataset_col=args.dataset_col,
        sample_col=args.sample_col,
        class_col=args.class_col,
        region_col=args.region_col,
        umap_key=args.umap_key,
        latent_key=args.latent_key,
        add_advanced=args.add_advanced,
        marker_genes=args.marker_genes,
    )

    outdir = Path(cfg.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(cfg.seed)

    a = ad.read_h5ad(cfg.h5ad)
    obs = a.obs.copy()

    dataset_col = _first_existing(obs, [cfg.dataset_col, "dataset", "study", "geo_accession"])
    sample_col = _first_existing(obs, [cfg.sample_col, "sample_id", "gsm"])
    class_col = _first_existing(obs, [cfg.class_col, "cell_class", "cell_type", "cell_type_origin"])
    region_col = _first_existing(obs, [cfg.region_col, "region", "brain_region"], required=False)
    protocol_col = _first_existing(obs, [cfg.protocol_col, "protocol_type", "technology"], required=False)
    age_col = _first_existing(obs, [cfg.age_col, "age_days", "protocol_age_days"], required=False)

    # Panel 1: pipeline overview schematic summary
    steps = pd.DataFrame(
        {
            "step_order": [1, 2, 3, 4, 5],
            "step": [
                "Input annotated atlas",
                "Reuse latent/labels",
                "Compute/reuse embedding",
                "Aggregate sample-level stats",
                "Export panel tables/figures",
            ],
            "n_cells": [a.n_obs, a.n_obs, a.n_obs, a.n_obs, a.n_obs],
            "n_features": [a.n_vars, a.n_vars, a.n_vars, a.n_vars, a.n_vars],
        }
    )
    tsv_pipeline = outdir / f"{cfg.prefix}_pipeline_overview.tsv"
    _write(steps, tsv_pipeline)

    # Panel 2: metadata summary
    meta = pd.DataFrame(
        {
            "metric": ["n_cells", "n_genes", "n_datasets", "n_samples", "n_cell_classes"],
            "value": [
                int(a.n_obs),
                int(a.n_vars),
                int(obs[dataset_col].nunique(dropna=True)),
                int(obs[sample_col].nunique(dropna=True)),
                int(obs[class_col].nunique(dropna=True)),
            ],
        }
    )
    if age_col:
        age_vals = pd.to_numeric(obs[age_col], errors="coerce")
        meta = pd.concat([meta, pd.DataFrame({"metric": ["age_min", "age_median", "age_max"], "value": [age_vals.min(), age_vals.median(), age_vals.max()]})], ignore_index=True)
    tsv_meta = outdir / f"{cfg.prefix}_metadata_summary.tsv"
    _write(meta, tsv_meta)

    # subset for plotting deterministically
    if a.n_obs > cfg.max_cells_plot:
        pick = np.sort(rng.choice(a.n_obs, size=cfg.max_cells_plot, replace=False))
        ap = a[pick].copy()
    else:
        ap = a.copy()

    _ensure_umap(ap, cfg.seed, cfg.latent_key, cfg.umap_key)
    um = pd.DataFrame(ap.obsm[cfg.umap_key], columns=["umap1", "umap2"], index=ap.obs_names)
    um["cell_id"] = um.index
    for c in [dataset_col, sample_col, class_col]:
        um[c] = ap.obs[c].astype(str).values
    if region_col:
        um[region_col] = ap.obs[region_col].astype(str).values
    if protocol_col:
        um[protocol_col] = ap.obs[protocol_col].astype(str).values
    if age_col:
        um[age_col] = pd.to_numeric(ap.obs[age_col], errors="coerce").values

    marker_present = [g for g in cfg.marker_genes if g in ap.var_names]
    if marker_present:
        X = ap[:, marker_present].X
        if hasattr(X, "toarray"):
            X = X.toarray()
        for j, g in enumerate(marker_present):
            um[f"expr_{g}"] = X[:, j]

    tsv_umap = outdir / f"{cfg.prefix}_umap_points.tsv"
    _write(um.reset_index(drop=True), tsv_umap)

    # Panel 4: sample-level stacked composition
    comp = (
        obs[[dataset_col, sample_col, class_col] + ([age_col] if age_col else [])]
        .assign(n=1)
        .groupby([dataset_col, sample_col, class_col], dropna=False, observed=False)["n"]
        .sum()
        .reset_index()
    )
    total = comp.groupby([dataset_col, sample_col], dropna=False)["n"].sum().rename("total").reset_index()
    comp = comp.merge(total, on=[dataset_col, sample_col], how="left")
    comp["fraction"] = comp["n"] / comp["total"]

    ages = None
    if age_col and age_col in obs.columns:
        ages = pd.to_numeric(obs[[sample_col, age_col]].drop_duplicates()[age_col], errors="coerce")
        age_map = obs[[sample_col, age_col]].drop_duplicates().set_index(sample_col)[age_col]
        comp[age_col] = pd.to_numeric(comp[sample_col].map(age_map), errors="coerce")

    tsv_comp = outdir / f"{cfg.prefix}_sample_composition.tsv"
    _write(comp, tsv_comp)

    # Plot multi-panel summary
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.ravel()
    axes[0].axis("off")
    axes[0].set_title("Pipeline summary")
    txt = f"cells={a.n_obs:,}\ngenes={a.n_vars:,}\ndatasets={obs[dataset_col].nunique()}\nsamples={obs[sample_col].nunique()}"
    axes[0].text(0.1, 0.6, txt, fontsize=12)

    if age_col:
        av = pd.to_numeric(obs[age_col], errors="coerce")
        axes[1].hist(av.dropna(), bins=30)
        axes[1].set_title("Sample age distribution")
    else:
        axes[1].text(0.2, 0.5, "No age column")
        axes[1].axis("off")

    if protocol_col:
        obs[protocol_col].astype(str).value_counts().head(15).plot(kind="bar", ax=axes[2])
        axes[2].set_title("Protocol counts")
    else:
        axes[2].text(0.2, 0.5, "No protocol column")
        axes[2].axis("off")

    axes[3].scatter(um["umap1"], um["umap2"], c=pd.Categorical(um[class_col]).codes, s=2, alpha=0.8)
    axes[3].set_title("Integrated UMAP by cell class")

    ccol = age_col if age_col and age_col in um.columns else class_col
    color = pd.to_numeric(um[ccol], errors="coerce") if ccol == age_col else pd.Categorical(um[ccol]).codes
    axes[4].scatter(um["umap1"], um["umap2"], c=color, s=2, alpha=0.8)
    axes[4].set_title(f"UMAP by {ccol}")

    top_samples = (comp[[dataset_col, sample_col, "total"]].drop_duplicates().sort_values([dataset_col, age_col if age_col else sample_col], na_position="last").head(25))
    comp_top = comp.merge(top_samples[[dataset_col, sample_col]], on=[dataset_col, sample_col], how="inner")
    pivot = comp_top.pivot_table(index=sample_col, columns=class_col, values="fraction", fill_value=0)
    pivot.plot(kind="bar", stacked=True, ax=axes[5], legend=False, width=0.9)
    axes[5].set_title("Sample composition (top 25 samples)")
    axes[5].set_ylabel("fraction")

    fig.tight_layout()
    png = outdir / f"{cfg.prefix}_panels.png"
    pdf = outdir / f"{cfg.prefix}_panels.pdf"
    fig.savefig(png, dpi=150)
    fig.savefig(pdf)
    plt.close(fig)

    extra_outputs: list[str] = []
    if cfg.add_advanced:
        import scanpy as sc

        sc.pp.neighbors(ap, use_rep=cfg.latent_key if cfg.latent_key in ap.obsm else None, random_state=cfg.seed)
        sc.tl.diffmap(ap)
        d1 = pd.DataFrame(ap.obsm["X_diffmap"][:, :3], columns=["dc1", "dc2", "dc3"], index=ap.obs_names)
        d1 = d1.reset_index().rename(columns={"index": "cell_id"})
        tsv_dc = outdir / f"{cfg.prefix}_diffusion_components.tsv"
        _write(d1, tsv_dc)
        extra_outputs.append(str(tsv_dc))

        fig2, ax2 = plt.subplots(1, 2, figsize=(10, 4))
        ax2[0].scatter(um["umap1"], um["umap2"], c=d1["dc1"].values, s=2)
        ax2[0].set_title("Diffusion component 1")

        if "transition_matrix" in ap.uns:
            T = ap.uns["transition_matrix"]
            if hasattr(T, "toarray"):
                T = T.toarray()
            step = max(1, len(um) // 1000)
            idx = np.arange(0, len(um), step)
            target = T[idx].argmax(axis=1)
            dx = um["umap1"].to_numpy()[target] - um["umap1"].to_numpy()[idx]
            dy = um["umap2"].to_numpy()[target] - um["umap2"].to_numpy()[idx]
            ax2[1].quiver(um["umap1"].to_numpy()[idx], um["umap2"].to_numpy()[idx], dx, dy, angles="xy", scale_units="xy", scale=1)
            ax2[1].set_title("Transition arrows")
        else:
            ax2[1].text(0.2, 0.5, "No transition_matrix in uns")
            ax2[1].axis("off")
        fig2.tight_layout()
        p2 = outdir / f"{cfg.prefix}_advanced_panels.png"
        fig2.savefig(p2, dpi=150)
        plt.close(fig2)
        extra_outputs.append(str(p2))

    commit_sha = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    manifest = {
        "sources": [str(Path(cfg.h5ad))],
        "commit_sha": commit_sha,
        "config": asdict(cfg),
        "outputs": [str(tsv_pipeline), str(tsv_meta), str(tsv_umap), str(tsv_comp), str(png), str(pdf)] + extra_outputs,
    }
    manifest_path = outdir / "figure1_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
