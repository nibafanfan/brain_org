#!/usr/bin/env python3
"""Metadata-adjusted protocol effects from calibrated transfer output.

Companion to benchmark_metadata_stratified.py. Fits covariate-adjusted models so
headline effects (e.g., multi-lineage enrichment) are not purely confounded by
age/protocol/guidance/vascularization/provenance mix.

Outputs:
  - data/metadata_adjusted_logit.tsv
  - data/metadata_adjusted_linear.tsv
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
    nums = re.findall(r"\d+", s)
    if not nums:
        return np.nan
    vals = [float(x) for x in nums]
    return float(np.median(vals))


def prep_df(obs, cellclass_col, conf_col, ood_col):
    d = pd.DataFrame(index=obs.index)
    d["protocol"] = obs["protocol"].astype(str).str.strip().replace({"": "unknown", "nan": "unknown"})
    d["multi_lineage_norm"] = obs["multi_lineage"].map(norm_multi)
    d["unguided_norm"] = obs["unguided"].map(norm_boolish)
    d["vascularized_norm"] = obs["vascularized"].map(norm_boolish)
    d["annotation_level"] = obs["annotation_level"].astype(str).str.strip().replace({"": "unknown", "nan": "unknown"})
    d["age_days_num"] = obs["age_days"].map(parse_age_days)
    d["CellClass_cal"] = obs[cellclass_col].astype(str)
    d["conf"] = pd.to_numeric(obs[conf_col], errors="coerce")
    d["ood"] = obs[ood_col].astype(bool).astype(int)
    d = d.replace([np.inf, -np.inf], np.nan)
    return d


def fit_logit_models(d):
    import statsmodels.formula.api as smf
    outcomes = {
        "is_immune": "CellClass_cal == 'Immune'",
        "is_vascular": "CellClass_cal == 'Vascular'",
        "is_oligo": "CellClass_cal == 'Oligo'",
    }
    rows = []
    # keep only clear binary strata for multi-lineage
    dd = d[d["multi_lineage_norm"].isin(["multi", "single"])].copy()
    dd["is_multi"] = (dd["multi_lineage_norm"] == "multi").astype(int)

    for name, expr in outcomes.items():
        dd[name] = dd.eval(expr).astype(int)
        sub = dd.dropna(subset=[name, "age_days_num"]).copy()
        if sub[name].sum() < 50:
            continue
        # sample for speed if huge
        if len(sub) > 1_000_000:
            sub = sub.sample(1_000_000, random_state=0)
        f = (
            f"{name} ~ is_multi + age_days_num + C(protocol) + C(unguided_norm) + "
            "C(vascularized_norm) + C(annotation_level)"
        )
        m = smf.logit(formula=f, data=sub).fit(disp=0, maxiter=100)
        if "is_multi" in m.params.index:
            beta = float(m.params["is_multi"])
            se = float(m.bse["is_multi"])
            p = float(m.pvalues["is_multi"])
            rows.append({
                "outcome": name,
                "coef_is_multi": beta,
                "or_is_multi": float(np.exp(beta)),
                "se": se,
                "z": beta / se if se else np.nan,
                "pvalue": p,
                "n": int(len(sub)),
                "events": int(sub[name].sum()),
            })
    return pd.DataFrame(rows)


def fit_linear_models(d):
    import statsmodels.formula.api as smf
    rows = []
    dd = d[d["multi_lineage_norm"].isin(["multi", "single"])].copy()
    dd["is_multi"] = (dd["multi_lineage_norm"] == "multi").astype(int)
    for outcome in ["conf", "ood"]:
        sub = dd.dropna(subset=[outcome, "age_days_num"]).copy()
        if len(sub) > 1_000_000:
            sub = sub.sample(1_000_000, random_state=0)
        f = (
            f"{outcome} ~ is_multi + age_days_num + C(protocol) + C(unguided_norm) + "
            "C(vascularized_norm) + C(annotation_level)"
        )
        m = smf.ols(formula=f, data=sub).fit()
        if "is_multi" in m.params.index:
            beta = float(m.params["is_multi"])
            se = float(m.bse["is_multi"])
            p = float(m.pvalues["is_multi"])
            rows.append({
                "outcome": outcome,
                "coef_is_multi": beta,
                "se": se,
                "t": beta / se if se else np.nan,
                "pvalue": p,
                "n": int(len(sub)),
            })
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-h5ad", default="data/braun_transfer_full_calibrated.h5ad")
    ap.add_argument("--out-logit", default="data/metadata_adjusted_logit.tsv")
    ap.add_argument("--out-linear", default="data/metadata_adjusted_linear.tsv")
    ap.add_argument("--cellclass-col", default="CellClass_cal")
    ap.add_argument("--conf-col", default="CellClass_cal_conf")
    ap.add_argument("--ood-col", default="ood")
    args = ap.parse_args()

    a = ad.read_h5ad(args.in_h5ad, backed="r")
    obs = a.obs.copy()
    required = ["protocol", "multi_lineage", "unguided", "vascularized", "annotation_level", "age_days",
                args.cellclass_col, args.conf_col, args.ood_col]
    missing = [c for c in required if c not in obs.columns]
    if missing:
        raise ValueError(f"missing required obs columns: {missing}")

    d = prep_df(obs, args.cellclass_col, args.conf_col, args.ood_col)

    try:
        logit = fit_logit_models(d)
        linear = fit_linear_models(d)
    except ModuleNotFoundError as e:
        raise SystemExit("statsmodels is required for this script. Install it in your env.") from e

    Path(args.out_logit).parent.mkdir(parents=True, exist_ok=True)
    logit.to_csv(args.out_logit, sep='\t', index=False)
    linear.to_csv(args.out_linear, sep='\t', index=False)
    print(f"wrote {args.out_logit} ({len(logit)} rows)")
    print(f"wrote {args.out_linear} ({len(linear)} rows)")


if __name__ == "__main__":
    main()
