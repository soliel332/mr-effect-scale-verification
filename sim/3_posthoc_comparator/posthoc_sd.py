"""
POST HOC declared-SD comparator. Same deterministic subset as posthoc_checksumstats.py (first 200
test replicates of each of the 86 test cells). Frozen modules imported, never modified.

CheckSumStats commit ba670767, R/make_plot_predlnor.R, predict_beta_sd():
    z = beta/se ; maf = folded eaf ; n = dat[, sample_size]   (DEFAULT sample_size = "ncontrol")
    beta_sd = z / sqrt(2 maf (1-maf) (n + z^2))
    sd_est  = median(beta / beta_sd) = median( se * sqrt(2 maf(1-maf)(n + z^2)) )
Manuscript eq. (1): sigma_hat = median( se * sqrt(2 N f (1-f)) ).
Per-variant ratio sd_est-term / sigma_hat-term = sqrt(1 + z^2/N) when n = N.

Declared-SD use of make_plot_pred_effect(): pred_beta = "beta_sd", beta = "beta", sd_est = NULL
(documented setting when "the reported effect size is already standardised"), no ncase column
(otherwise the log-OR filter format_data_predlnor_sh() is triggered). Slope = OLS with intercept
of beta_sd on reported beta. Alarm outside 0.8-1.20 (Haycock et al. rule, not in the code).
The default call (sd_est estimated from the same data) rescales the reported beta by sd_est and is
reported separately as the self-normalising setting.
"""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np, pandas as pd
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sim_scale as ss, posthoc_checksumstats as pc

OUTD = os.path.join(HERE, "posthoc_sd_out")


def ols_slope(y, x):
    if len(x) < 3:
        return np.nan
    return np.linalg.lstsq(np.c_[np.ones(len(x)), x], y, rcond=None)[0][1]


def one(cell, rep, seq):
    rng = np.random.default_rng(seq)
    K = cell["K"]; g = dict(cell); g["K"] = 3 * K
    d = ss.generate(rng, g)
    out, fl = ss.analyse(d)
    lp = ss.pipeline_logp(out, fl["logp_spa"])
    N = d["N_rep"]; ncase = d["p_real"] * N; nctrl = N - ncase
    rows = []
    for pl in ss.PIPELINES:
        idx, npass = ss.select(lp[pl], K)
        b = out[pl][0][idx]; se = out[pl][1][idx]; f = d["f_rep"][idx]
        r = ss.classify(se, f, N, d["p_real"])
        z = b / se; maf = np.where(f > 0.5, 1 - f, f)
        ok = np.isfinite(se) & (se > 0) & np.isfinite(f) & (f > 0) & (f < 1)
        fac = np.sqrt(1 + z[ok] ** 2 / N)
        sd_N = np.median(se[ok] * np.sqrt(2 * maf[ok] * (1 - maf[ok]) * (N + z[ok] ** 2))) if ok.sum() else np.nan
        sd_ctrl = np.median(se[ok] * np.sqrt(2 * maf[ok] * (1 - maf[ok]) * (nctrl + z[ok] ** 2))) if ok.sum() else np.nan
        bsd = z[ok] / np.sqrt(2 * maf[ok] * (1 - maf[ok]) * (N + z[ok] ** 2))
        r.update(cell=cell["cell"], arm=cell["arm"], level=cell["level"], p_target=cell["p"], K=K, rep=rep,
                 pipeline=pl, true_scale=ss.TRUE_SCALE[pl], p_real=d["p_real"], n_pass=npass,
                 z2overN_median=float(np.median(z[ok] ** 2 / N)) if ok.sum() else np.nan,
                 factor_median=float(np.median(fac)) if ok.sum() else np.nan,
                 factor_max=float(fac.max()) if ok.sum() else np.nan,
                 css_sd_est_totalN=sd_N, css_sd_est_default_ncontrol=sd_ctrl,
                 css_slope_declaredSD=ols_slope(bsd, b[ok]) if ok.sum() >= 3 else np.nan,
                 css_slope_selfnormalised=ols_slope(bsd, b[ok] / sd_N) if ok.sum() >= 3 else np.nan)
        # three-signature rule applied to CheckSumStats' own estimator
        if np.isfinite(sd_N) and ok.sum() >= 5:
            L = 0.5 * np.log(d["p_real"] * (1 - d["p_real"]))
            dist = np.abs(np.log(sd_N) - np.array([L, 0.0, -L]))
            r["css_sd_est_nearest"] = ss.SCALES[np.argmin(dist)]; r["css_sd_est_dmin"] = float(dist.min())
        else:
            r["css_sd_est_nearest"] = "abstain_sufficiency"; r["css_sd_est_dmin"] = np.nan
        rows.append(r)
    return rows


def run(task):
    cell, a, b, seqs = task
    fn = os.path.join(OUTD, f"{cell['cell']}__{a:04d}.parquet")
    if os.path.exists(fn):
        return fn
    rows = []
    for i, s in zip(range(a, b), seqs):
        rows += one(cell, i, s)
    pd.DataFrame(rows).to_parquet(fn + ".tmp", index=False); os.replace(fn + ".tmp", fn)
    return fn


if __name__ == "__main__":
    os.makedirs(OUTD, exist_ok=True)
    T = [(c, a, b, s) for (c, a, b, s) in pc.task_list()]
    T.sort(key=lambda t: -t[0]["K"] - (20 if t[0].get("design") == "casecontrol" else 0))
    print("test cells:", len({t[0]["cell"] for t in T}), "| tasks:", len(T), flush=True)
    done = 0
    with Pool(10) as p:
        for _ in p.imap_unordered(run, T):
            done += 1
            if done % 25 == 0 or done == len(T):
                print(f"{done}/{len(T)}", flush=True)
    print("complete")
