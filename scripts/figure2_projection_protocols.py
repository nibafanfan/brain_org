#!/usr/bin/env python3
"""Reproducible Figure 2 pipeline: projection panels + protocol benchmarking.

Inputs
------
A calibrated transfer AnnData file containing at least:
- predictions/confidence columns for mapped class and region
- metadata columns: protocol, unguided, vascularized, multi_lineage,
  age_days, annotation_level

Outputs
-------
- Projection panels (UMAP maps, age-stage heatmap, protocol-region bars)
- Coverage/presence diagnostics (max-presence distributions, under-represented states)
- Protocol comparison tables (age-stratified composition, rare-lineage prevalence,
  OOD/confidence/abstention summaries)
- Covariate-adjusted protocol effect-size tables
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.linear_model import LogisticRegression


DEFAULT_INPUT = "data/transfer/braun_transfer_full_calibrated.h5ad"


def _first_present(options: list[str], columns: pd.Index, *, required: bool = True) -> str | None:
    for col in options:
        if col in columns:
            return col
    if required:
        raise ValueError(f"None of columns found: {options}")
    return None


def _ensure_obs_columns(adata: ad.AnnData, required_cols: list[str]) -> None:
    missing = [c for c in required_cols if c not in adata.obs.columns]
    if missing:
        raise ValueError(f"Missing required .obs columns: {missing}")


def _safe_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    return series.astype(str).str.lower().isin({"1", "true", "t", "yes", "y"})


def _abstention_summary(df: pd.DataFrame, confidence_col: str, threshold: float) -> pd.DataFrame:
    out = (
        df.assign(abstained=df[confidence_col] < threshold)
        .groupby("protocol", observed=True)
        .agg(
            n_cells=(confidence_col, "size"),
            mean_confidence=(confidence_col, "mean"),
            median_confidence=(confidence_col, "median"),
            abstention_rate=("abstained", "mean"),
        )
        .reset_index()
    )
    return out


def _model_protocol_effects(
    df: pd.DataFrame,
    protocol_col: str,
    target_col: str,
    age_col: str,
    unguided_col: str,
    vascularized_col: str,
    annotation_level_col: str,
    out_csv: Path,
) -> pd.DataFrame:
    model_df = df[
        [protocol_col, target_col, age_col, unguided_col, vascularized_col, annotation_level_col]
    ].dropna()
    y = (model_df[target_col].astype(str) == "1").astype(int)

    x = pd.DataFrame(
        {
            "age_days": pd.to_numeric(model_df[age_col], errors="coerce"),
            "unguided": _safe_bool(model_df[unguided_col]).astype(int),
            "vascularized": _safe_bool(model_df[vascularized_col]).astype(int),
            "annotation_level": model_df[annotation_level_col].astype(str),
            "protocol": model_df[protocol_col].astype(str),
        }
    ).dropna()
    y = y.loc[x.index]

    x_design = pd.get_dummies(x, columns=["annotation_level", "protocol"], drop_first=True)
    if y.nunique() < 2:
        effect = pd.DataFrame(columns=["term", "coef", "odds_ratio", "target"])
        effect.to_csv(out_csv, index=False)
        return effect

    lr = LogisticRegression(max_iter=1000)
    lr.fit(x_design, y)

    effect = pd.DataFrame({"term": x_design.columns, "coef": lr.coef_[0]})
    effect["odds_ratio"] = np.exp(effect["coef"])
    effect["target"] = target_col
    effect = effect.sort_values("odds_ratio", ascending=False)
    effect.to_csv(out_csv, index=False)
    return effect


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=DEFAULT_INPUT)
    ap.add_argument("--outdir", default="outputs/figure2_projection_protocols")
    ap.add_argument("--abstain-threshold", type=float, default=0.5)
    ap.add_argument("--rare-threshold", type=float, default=0.01)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    figdir = outdir / "figures"
    tabdir = outdir / "tables"
    figdir.mkdir(parents=True, exist_ok=True)
    tabdir.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.input)
    required_metadata = [
        "protocol",
        "unguided",
        "vascularized",
        "multi_lineage",
        "age_days",
        "annotation_level",
    ]
    _ensure_obs_columns(adata, required_metadata)

    cls_col = _first_present(["CellClass_pred", "cell_class_pred", "mapped_class"], adata.obs.columns)
    reg_col = _first_present(["Region_pred", "region_pred", "mapped_region"], adata.obs.columns)
    conf_col = _first_present(
        ["CellClass_conf", "cell_class_conf", "confidence", "calibrated_confidence"], adata.obs.columns
    )
    stage_col = _first_present(["primary_stage", "PrimaryStage", "stage_pred"], adata.obs.columns, required=False)

    umap_key = "X_umap" if "X_umap" in adata.obsm else None

    df = adata.obs.copy()
    df["protocol"] = df["protocol"].astype(str)
    df["age_days"] = pd.to_numeric(df["age_days"], errors="coerce")

    # 1) mapped cell classes and regions on organoid UMAP
    if umap_key is not None:
        umap = pd.DataFrame(adata.obsm[umap_key][:, :2], columns=["UMAP1", "UMAP2"], index=adata.obs_names)
        plot_df = umap.join(df[[cls_col, reg_col]])

        for col, fname in [(cls_col, "umap_mapped_cell_class.png"), (reg_col, "umap_mapped_region.png")]:
            plt.figure(figsize=(8, 7))
            sns.scatterplot(
                data=plot_df.sample(min(200_000, len(plot_df)), random_state=0),
                x="UMAP1",
                y="UMAP2",
                hue=col,
                s=3,
                linewidth=0,
                alpha=0.7,
                legend=False,
            )
            plt.title(f"Organoid UMAP colored by {col}")
            plt.tight_layout()
            plt.savefig(figdir / fname, dpi=220)
            plt.close()

    # 2) age-to-primary-stage matching heatmap
    if stage_col is not None:
        age_bins = pd.cut(df["age_days"], bins=[0, 30, 60, 90, 120, 180, 365, np.inf], include_lowest=True)
        hm = pd.crosstab(age_bins, df[stage_col], normalize="index")
        hm.to_csv(tabdir / "age_to_primary_stage_heatmap_values.csv")

        plt.figure(figsize=(10, 4))
        sns.heatmap(hm, cmap="viridis")
        plt.title("Age-to-primary-stage matching")
        plt.tight_layout()
        plt.savefig(figdir / "age_to_primary_stage_matching_heatmap.png", dpi=220)
        plt.close()

    # 3) regional proportion bars by dataset/protocol
    if "dataset" not in df.columns:
        df["dataset"] = "organoid"
    region_prop = (
        df.groupby(["dataset", "protocol", reg_col], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    region_prop["prop"] = region_prop["n_cells"] / region_prop.groupby(["dataset", "protocol"], observed=True)[
        "n_cells"
    ].transform("sum")
    region_prop.to_csv(tabdir / "regional_proportion_by_dataset_protocol.csv", index=False)

    plt.figure(figsize=(11, 5))
    sns.barplot(data=region_prop, x="protocol", y="prop", hue=reg_col)
    plt.xticks(rotation=45, ha="right")
    plt.title("Regional proportion by protocol")
    plt.tight_layout()
    plt.savefig(figdir / "regional_proportions_by_protocol.png", dpi=220)
    plt.close()

    # 4) max presence score distributions by class/region
    max_presence = (
        df.groupby([cls_col, reg_col], observed=True)[conf_col].max().rename("max_presence").reset_index()
    )
    max_presence.to_csv(tabdir / "max_presence_by_class_region.csv", index=False)

    plt.figure(figsize=(9, 5))
    sns.histplot(data=max_presence, x="max_presence", hue=cls_col, bins=30, element="step", stat="density")
    plt.title("Max presence score distribution by class/region")
    plt.tight_layout()
    plt.savefig(figdir / "max_presence_distributions.png", dpi=220)
    plt.close()

    # 5) under-represented primary states visualization
    cls_prop = df[cls_col].value_counts(normalize=True).rename("prop").reset_index(names=["cell_class"])
    under = cls_prop[cls_prop["prop"] < args.rare_threshold].copy()
    under.to_csv(tabdir / "underrepresented_primary_states.csv", index=False)
    if not under.empty:
        plt.figure(figsize=(8, max(3, 0.3 * len(under))))
        sns.barplot(data=under.sort_values("prop"), x="prop", y="cell_class", color="#bc5090")
        plt.title("Under-represented primary states")
        plt.tight_layout()
        plt.savefig(figdir / "underrepresented_primary_states.png", dpi=220)
        plt.close()

    # Protocol comparison tables
    age_bin = pd.cut(df["age_days"], bins=[0, 30, 60, 90, 120, 180, 365, np.inf], include_lowest=True)
    age_comp = (
        df.assign(age_bin=age_bin)
        .groupby(["protocol", "age_bin", cls_col], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    age_comp["prop"] = age_comp["n_cells"] / age_comp.groupby(["protocol", "age_bin"], observed=True)[
        "n_cells"
    ].transform("sum")
    age_comp.to_csv(tabdir / "age_stratified_cell_class_composition.csv", index=False)

    rare_lineage = (
        df.assign(is_rare=df[cls_col].map(df[cls_col].value_counts(normalize=True)) < args.rare_threshold)
        .groupby("protocol", observed=True)["is_rare"]
        .mean()
        .rename("rare_lineage_prevalence")
        .reset_index()
    )
    rare_lineage.to_csv(tabdir / "rare_lineage_prevalence.csv", index=False)

    ood = _abstention_summary(df, conf_col, args.abstain_threshold)
    ood["ood_rate"] = (df[conf_col] < 0.2).groupby(df["protocol"]).mean().values
    ood.to_csv(tabdir / "ood_confidence_abstention_summary.csv", index=False)

    # Covariate-adjusted protocol effects
    df["multi_lineage_bin"] = _safe_bool(df["multi_lineage"]).astype(int).astype(str)
    _model_protocol_effects(
        df,
        protocol_col="protocol",
        target_col="multi_lineage_bin",
        age_col="age_days",
        unguided_col="unguided",
        vascularized_col="vascularized",
        annotation_level_col="annotation_level",
        out_csv=tabdir / "covariate_adjusted_protocol_effects_multi_lineage.csv",
    )

    print(f"Done. Wrote outputs to: {outdir}")


if __name__ == "__main__":
    main()
