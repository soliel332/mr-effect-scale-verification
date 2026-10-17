"""
POST HOC comparator analysis. Written before any benchmark output was inspected.

Frame: every release is presented to both methods as if its coefficients were log odds ratios —
the use case of CheckSumStats predict_lnor_sh and the default assumption of an MR analysis.
Status per release and method, over ALL releases (unconditional denominators):
  CheckSumStats  : not_assessable (< 3 variants after package filters) | alarm (slope outside
                   0.8-1.20, Haycock et al. rule) | no_alarm
  3-sig primary  : not_assessable (sufficiency abstention) | alarm (assigned scale != log odds)
                   | no_alarm
  3-sig secondary: as primary, plus alarm when distance to the assigned signature exceeds tau
Case-fraction-supplied records only (CheckSumStats uses case and control counts, so arm A4b has
no CheckSumStats analogue and is reported for the three-signature rule alone).
"""
import glob, os
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
TAU = float(open(os.path.join(HERE, "sim_tau.txt")).read().split("=")[-1])


def load():
    x = pd.concat([pd.read_parquet(f) for f in sorted(glob.glob(os.path.join(HERE, "posthoc_out", "*.parquet")))], ignore_index=True)
    return x


def verify(x):
    fr = pd.read_parquet(os.path.join(HERE, "sim_results_records.parquet"))
    fr = fr[fr.set == "test"]
    m = x.merge(fr, on=["cell", "rep", "pipeline", "p_supplied_as"], suffixes=("_n", "_o"), how="left")
    return dict(regenerated_records=len(x), matched_to_frozen=int(m.assigned_o.notna().sum()),
                label_mismatch=int((m.assigned_n != m.assigned_o).sum()),
                n_pass_mismatch=int((m.n_pass_n != m.n_pass_o).sum()),
                max_abs_diff_sigma_hat=float(np.nanmax(np.abs(m.sigma_hat_n - m.sigma_hat_o))),
                finiteness_mismatch=int((np.isfinite(m.sigma_hat_n) != np.isfinite(m.sigma_hat_o)).sum()))


def status(x):
    x = x.copy()
    x["css"] = np.where(x.css_slope.isna(), "not_assessable",
                np.where((x.css_slope > 1.20) | (x.css_slope < 0.8), "alarm", "no_alarm"))
    ab = x.assigned.str.startswith("abstain")
    x["sig_primary"] = np.where(ab, "not_assessable", np.where(x.assigned != "log_odds", "alarm", "no_alarm"))
    x["sig_secondary"] = np.where(ab, "not_assessable",
                         np.where((x.assigned != "log_odds") | (x.d_min > TAU), "alarm", "no_alarm"))
    x["scale_deviation"] = x.true_scale != "log_odds"
    x["nonscale_deviation"] = x.arm.isin(["A5", "A6", "A7", "A8"])
    x["condition"] = np.select([x.scale_deviation & x.nonscale_deviation, x.scale_deviation, x.nonscale_deviation],
                               ["scale + non-scale", "scale only", "non-scale only"], "none")
    lam = 1 / (x.p_real * (1 - x.p_real))
    ls = np.log(np.abs(x.css_slope))
    cand = np.c_[np.log(lam), 0.5 * np.log(lam), np.zeros(len(x))]
    lab = np.array(["linear_probability", "standardised", "log_odds"])
    dist = np.abs(ls.to_numpy()[:, None] - cand)
    dist = np.where(np.isfinite(dist), dist, np.inf)
    x["css_implied_scale"] = np.where(x.css_slope.isna() | (x.css_slope <= 0), "none", lab[np.argmin(dist, axis=1)])
    return x


def tables(x):
    y = x[x.p_supplied_as == "case_fraction"]
    rows = []
    keys = ["arm", "level", "p_target", "K", "pipeline", "true_scale", "condition"]
    for k, g in y.groupby(keys):
        r = dict(zip(keys, k)); r["n"] = len(g)
        for mth in ("css", "sig_primary", "sig_secondary"):
            for s in ("alarm", "no_alarm", "not_assessable"):
                r[f"{mth}_{s}"] = (g[mth] == s).mean()
        r["css_slope_median"] = g.css_slope.median()
        r["css_slope_q10"] = g.css_slope.quantile(.10); r["css_slope_q90"] = g.css_slope.quantile(.90)
        r["css_n_used_median"] = g.css_n_used.median()
        cls = g[g.assigned.isin(["linear_probability", "standardised", "log_odds"])]
        r["sig_label_correct"] = (g.assigned == g.true_scale).mean()
        r["css_implied_label_correct"] = (g.css_implied_scale == g.true_scale).mean()
        rows.append(r)
    return pd.DataFrame(rows)


def a7(x):
    y = x[(x.arm == "A7") & (x.p_supplied_as == "case_fraction")]
    out = []
    for (lvl, p, pl), g in y.groupby(["level", "p_target", "pipeline"]):
        fr = np.concatenate([np.array(s.split(";"), float) for s in g.varN_fracs if s])
        r = dict(level=lvl, p_target=p, pipeline=pl, releases=len(g), selected_variants=len(fr),
                 frac_variants_full_N=(fr == 1).mean(), variantN_frac_median=np.median(fr),
                 variantN_frac_p10=np.quantile(fr, .1), variantN_frac_min=fr.min(),
                 css_slope_median=g.css_slope.median(), css_slope_q10=g.css_slope.quantile(.1),
                 css_slope_min=g.css_slope.min(), frac_css_slope_below_0_8=(g.css_slope < 0.8).mean(),
                 frac_css_slope_at_or_below_0_30=(g.css_slope <= 0.30).mean())
        if lvl == "nhet_matched":
            sc = np.round(fr * 15).astype(int)
            r["study_count_median"] = float(np.median(sc)); r["study_count_min"] = int(sc.min())
            r["frac_variants_15_studies"] = float((sc == 15).mean())
        out.append(r)
    return pd.DataFrame(out)


if __name__ == "__main__":
    x = load()
    v = verify(x); pd.Series(v).to_csv(os.path.join(HERE, "posthoc_regeneration_check.csv"), header=["value"])
    x = status(x)
    T = tables(x); T.to_csv(os.path.join(HERE, "posthoc_comparator_by_cell.csv"), index=False)
    A = a7(x); A.to_csv(os.path.join(HERE, "posthoc_A7_variantN_and_slope.csv"), index=False)
    x.drop(columns=["varN_fracs"]).to_parquet(os.path.join(HERE, "posthoc_records.parquet"), index=False)
    print(pd.Series(v).to_string())
