"""
POST HOC. Regenerates every comparator number reported in COMPARATOR_ADJUDICATION.md and in the
manuscript/Online Resource from the two post hoc record-level files:

    posthoc_records.parquet      log-odds frame  (posthoc_checksumstats.py -> posthoc_compare.py)
    posthoc_sd_records.parquet   declared-SD frame (posthoc_sd.py)

Writes:
    posthoc_logodds_anomaly_paired.csv       CheckSumStats slope alarm vs three-signature adequacy alarm
    posthoc_logodds_attribution_paired.csv   slope-implied scale vs three-signature primary assignment
    posthoc_declaredSD_comparison.csv        declared-SD frame, both methods
    posthoc_factor_quantiles.csv             sqrt(1 + z^2/N) and sd_est/sigma_hat quantiles
    posthoc_scalar_summaries.csv             label agreement, slope-sigma relationship, default-ncontrol effect

Run from the repository root:  python sim/posthoc_comparator_tables.py
"""
import os
import numpy as np, pandas as pd

H = os.path.dirname(os.path.abspath(__file__))
TAU = float(open(os.path.join(H, "sim_tau.txt")).read().split("=")[-1])
LAB = np.array(["linear_probability", "standardised", "log_odds"])


def grp(pl):
    return np.where(pl == "P1", "P1 linear", np.where(pl == "P2", "P2 standardised", "P3-P7 log odds"))


def setting(x):
    return np.where(x.arm == "core", "core K=" + x.K.astype(str), x.arm + " " + x.level)


def nearest(logs, p):
    L = 0.5 * np.log(p * (1 - p))
    d = np.abs(np.asarray(logs)[:, None] - np.c_[L, 0 * L, -L])
    d = np.where(np.isfinite(d), d, np.inf)
    return LAB[np.argmin(d, 1)], d.min(1)


def logodds_frame():
    x = pd.read_parquet(os.path.join(H, "posthoc_records.parquet"))
    x = x[x.p_supplied_as == "case_fraction"].copy()
    x["grp"] = grp(x.pipeline); x["setting"] = setting(x)
    x["css_alarm"] = x.css == "alarm"; x["css_na"] = x.css == "not_assessable"
    x["adq_alarm"] = x.sig_secondary == "alarm"; x["sig_na"] = x.sig_primary == "not_assessable"
    x["css_lab"] = x.css_implied_scale == x.true_scale; x["sig_lab"] = x.assigned == x.true_scale
    x["css_wronglab"] = (x.css_implied_scale != "none") & (x.css_implied_scale != x.true_scale)
    x["sig_wronglab"] = (~x.assigned.str.startswith("abstain")) & (x.assigned != x.true_scale)
    x["css_nolab"] = x.css_implied_scale == "none"; x["sig_nolab"] = x.assigned.str.startswith("abstain")
    A = x.groupby(["setting", "grp"])[["css_alarm", "adq_alarm", "css_na", "sig_na"]].mean().reset_index()
    B = x.groupby(["setting", "grp"])[["css_lab", "sig_lab", "css_wronglab", "sig_wronglab", "css_nolab", "sig_nolab"]].mean().reset_index()
    A.to_csv(os.path.join(H, "posthoc_logodds_anomaly_paired.csv"), index=False)
    B.to_csv(os.path.join(H, "posthoc_logodds_attribution_paired.csv"), index=False)
    return len(x), x.cell.nunique()


def sd_frame():
    x = pd.read_parquet(os.path.join(H, "posthoc_sd_records.parquet"))
    x["grp"] = grp(x.pipeline); x["setting"] = setting(x)
    ab = x.assigned.str.startswith("abstain")
    x["sig_prim_alarm"] = (~ab) & (x.assigned != "standardised")
    x["sig_adq_alarm"] = (~ab) & ((x.assigned != "standardised") | (x.d_min > TAU))
    na = x.css_slope_declaredSD.isna()
    x["css_sd_alarm"] = (~na) & ((x.css_slope_declaredSD > 1.2) | (x.css_slope_declaredSD < 0.8))
    x["css_selfnorm_alarm"] = (~x.css_slope_selfnormalised.isna()) & ((x.css_slope_selfnormalised > 1.2) | (x.css_slope_selfnormalised < 0.8))
    S = x.groupby(["setting", "grp"]).agg(css_declSD_alarm=("css_sd_alarm", "mean"), css_selfnorm_alarm=("css_selfnorm_alarm", "mean"),
                                          css_na=("css_slope_declaredSD", lambda s: s.isna().mean()),
                                          sig_primary_alarm=("sig_prim_alarm", "mean"), sig_adequacy_alarm=("sig_adq_alarm", "mean"),
                                          sig_na=("assigned", lambda s: s.str.startswith("abstain").mean()),
                                          css_slope_median=("css_slope_declaredSD", "median")).reset_index()
    S.to_csv(os.path.join(H, "posthoc_declaredSD_comparison.csv"), index=False)
    ok = x.sigma_hat.notna()
    x["r_total"] = x.css_sd_est_totalN / x.sigma_hat
    x["r_default"] = x.css_sd_est_default_ncontrol / x.sigma_hat
    q = [0, .5, .99, 1]
    F = pd.DataFrame({"per_release_median_factor": x.factor_median[ok].quantile(q).values,
                      "per_release_max_variant_factor": x.factor_max[ok].quantile(q).values,
                      "release_ratio_sd_est_totalN_over_sigma_hat": x.r_total[ok].quantile(q).values,
                      "release_ratio_sd_est_default_ncontrol_over_sigma_hat": x.r_default[ok].quantile(q).values},
                     index=["min", "median", "p99", "max"])
    F.to_csv(os.path.join(H, "posthoc_factor_quantiles.csv"))
    same = np.where(ab, x.css_sd_est_nearest.str.startswith("abstain"), x.css_sd_est_nearest == x.assigned)
    lab_def, _ = nearest(np.log(x.css_sd_est_default_ncontrol.to_numpy()), x.p_real.to_numpy())
    core = x[(x.arm == "core") & ok]
    acc_def = (pd.Series(lab_def, index=x.index)[core.index] == core.true_scale).groupby(core.p_target).mean()
    pos = ok & x.css_slope_declaredSD.notna() & (x.css_slope_declaredSD > 0)
    sc = {"records": len(x), "test_cells": x.cell.nunique(),
          "label_agreement_sd_est_totalN_vs_sigma_hat": float(np.mean(same)),
          "label_disagreements": int((~same).sum()),
          "corr_log_declSD_slope_vs_minus_log_sigma_hat": float(np.corrcoef(np.log(x.css_slope_declaredSD[pos]), -np.log(x.sigma_hat[pos]))[0, 1]),
          "median_declSD_slope_times_sigma_hat": float((x.css_slope_declaredSD * x.sigma_hat)[pos].median()),
          "selfnormalised_alarm_P1_core": float(core[core.pipeline == "P1"].pipe(lambda g: ((g.css_slope_selfnormalised > 1.2) | (g.css_slope_selfnormalised < 0.8)).mean()))}
    for p, v in acc_def.items():
        sc[f"core_label_accuracy_default_ncontrol_p{p}"] = float(v)
    pd.Series(sc).to_csv(os.path.join(H, "posthoc_scalar_summaries.csv"), header=["value"])
    return sc


if __name__ == "__main__":
    n, c = logodds_frame()
    sc = sd_frame()
    print("log-odds frame records (case fraction supplied):", n, "| cells:", c)
    for k, v in sc.items():
        print(f"{k}: {v}")
    print(pd.read_csv(os.path.join(H, "posthoc_factor_quantiles.csv"), index_col=0).round(5).to_string())
