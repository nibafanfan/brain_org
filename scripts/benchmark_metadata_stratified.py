#!/usr/bin/env python3
"""Metadata-aware protocol benchmarking from calibrated transfer output.

Uses per-cell finalized GSM annotations (cell_type_origin, age_days, organoid_type,
protocol, unguided, multi_lineage, vascularized, slice, annotation_level) to produce
stratified summaries for protocol comparison and projection-quality auditing.

Inputs:
  - data/braun_transfer_full_calibrated.h5ad
Outputs:
  - data/metadata_protocol_scorecard.tsv
  - data/metadata_strata_summary.tsv
"""
import argparse
import re
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad


def norm_boolish(v):
    s = str(v).strip().lower()
    if s in {"1", "true", "yes", "y"}:
        return "yes"
    if s in {"0", "false", "no", "n"}:
        return "no"
    if "sort of" in s:
        return "sort_of"
    if s in {"", "nan", "none", "unknown"}:
        return "unknown"
    return "unknown"


def norm_multi(v):
    s = str(v).strip().lower()
    if s in {"1", "true", "yes", "y"}:
        return "multi"
    if s in {"0", "false", "no", "n"}:
        return "single"
    if "sort of" in s:
        return "sort_of"
    return "unknown"


def parse_age_days(v):
    s = str(v).strip().lower()
    if s in {"", "nan", "none", "unknown"}:
        return np.nan
    nums = re.findall(r"\d+", s)
    if not nums:
        return np.nan
    vals = [float(x) for x in nums]
    return float(np.median(vals))


def age_bin(x):
    if pd.isna(x):
        return "unknown"
    x = float(x)
    if x <= 30:
        return "0-30"
    if x <= 60:
        return "31-60"
    if x <= 120:
        return "61-120"
    return "120+"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-h5ad", default="data/braun_transfer_full_calibrated.h5ad")
    ap.add_argument("--out-scorecard", default="data/metadata_protocol_scorecard.tsv")
    ap.add_argument("--out-strata", default="data/metadata_strata_summary.tsv")
    ap.add_argument("--cellclass-col", default="CellClass_cal")
    ap.add_argument("--conf-col", default="CellClass_cal_conf")
    ap.add_argument("--ood-col", default="ood")
    args = ap.parse_args()

    a = ad.read_h5ad(args.in_h5ad, backed="r")
    obs = a.obs.copy()

    required = [
        "protocol", "multi_lineage", "unguided", "vascularized", "slice",
        "annotation_level", "age_days", args.cellclass_col, args.conf_col, args.ood_col,
        "dataset_slug"
    ]
    missing = [c for c in required if c not in obs.columns]
    if missing:
        raise ValueError(f"missing required obs columns: {missing}")

    d = pd.DataFrame(index=obs.index)
    d["protocol"] = obs["protocol"].astype(str).str.strip().replace({"": "unknown", "nan": "unknown"})
    d["multi_lineage_norm"] = obs["multi_lineage"].map(norm_multi)
    d["unguided_norm"] = obs["unguided"].map(norm_boolish)
    d["vascularized_norm"] = obs["vascularized"].map(norm_boolish)
    d["slice_norm"] = obs["slice"].map(norm_boolish)
    d["annotation_level"] = obs["annotation_level"].astype(str).str.strip().replace({"": "unknown", "nan": "unknown"})
    d["age_days_num"] = obs["age_days"].map(parse_age_days)
    d["age_bin"] = d["age_days_num"].map(age_bin)
    d["CellClass_cal"] = obs[args.cellclass_col].astype(str)
    d["conf"] = pd.to_numeric(obs[args.conf_col], errors="coerce")
    d["ood"] = obs[args.ood_col].astype(bool)
    d["dataset_slug"] = obs["dataset_slug"].astype(str)

    # protocol-level scorecard
    grp = d.groupby("protocol", observed=True)
    prot = grp.agg(
        n_cells=("CellClass_cal", "size"),
        n_deposits=("dataset_slug", "nunique"),
        age_median=("age_days_num", "median"),
        age_iqr=("age_days_num", lambda x: np.nanpercentile(x, 75) - np.nanpercentile(x, 25)),
        ood_rate=("ood", "mean"),
        mean_conf=("conf", "mean"),
        gsm_level_rate=("annotation_level", lambda x: np.mean(x == "gsm")),
        multi_rate=("multi_lineage_norm", lambda x: np.mean(x == "multi")),
        vascularized_rate=("vascularized_norm", lambda x: np.mean(x == "yes")),
    ).reset_index()

    # repertoire breadth at protocol level: classes with >=0.5% prevalence and >=100 cells
    breadth = []
    for p, sub in grp:
        vc = sub["CellClass_cal"].value_counts()
        frac = vc / vc.sum()
        keep = vc[(vc >= 100) & (frac >= 0.005)]
        breadth.append((p, int(len(keep)), ";".join(keep.index.tolist()[:12])))
    breadth_df = pd.DataFrame(breadth, columns=["protocol", "repertoire_breadth", "top_classes"])
    prot = prot.merge(breadth_df, on="protocol", how="left")
    prot = prot.sort_values("n_cells", ascending=False)

    # strata summary for downstream benchmarking
    strata = d.groupby([
        "multi_lineage_norm", "unguided_norm", "vascularized_norm", "age_bin", "annotation_level"
    ], observed=True).agg(
        n_cells=("CellClass_cal", "size"),
        n_protocols=("protocol", "nunique"),
        n_deposits=("dataset_slug", "nunique"),
        ood_rate=("ood", "mean"),
        mean_conf=("conf", "mean"),
        immune_rate=("CellClass_cal", lambda x: np.mean(x == "Immune")),
        vascular_rate=("CellClass_cal", lambda x: np.mean(x == "Vascular")),
        oligo_rate=("CellClass_cal", lambda x: np.mean(x == "Oligo")),
    ).reset_index().sort_values("n_cells", ascending=False)

    Path(args.out_scorecard).parent.mkdir(parents=True, exist_ok=True)
    prot.to_csv(args.out_scorecard, sep='\t', index=False)
    strata.to_csv(args.out_strata, sep='\t', index=False)

    print(f"wrote {args.out_scorecard} ({len(prot)} protocols)")
    print(f"wrote {args.out_strata} ({len(strata)} strata rows)")


if __name__ == "__main__":
    main()
