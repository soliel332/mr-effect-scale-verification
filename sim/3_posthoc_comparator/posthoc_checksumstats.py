"""
POST HOC benchmark — not part of the frozen protocol.

Compares the three-signature rule with the CheckSumStats `predict_lnor_sh` implementation on the
SAME releases, regenerated deterministically from the frozen seeds. The frozen modules
sim_scale.py and run_sim.py are imported, never modified; sim/out/ is read, never written.

CheckSumStats source: github.com/MRCIEU/CheckSumStats, commit
ba670767d42a9f6041b65e4a155b64a3ac096434 (2023-05-02), file R/make_plot_predlnor.R, functions
predict_lnor_sh(), format_data_predlnor_sh(), make_plot_pred_effect(). Reproduced here:
  * predicted log OR = argmin over a signed grid of |x / se(x) - z|, with se(x) the closed-form
    expression in that function (Harrison approach), inputs z, folded MAF, ncase, ncontrol;
    the package's 1e-6 grid is reproduced by a 1e-3 grid refined at 1e-6 around the minimum;
  * filters of format_data_predlnor_sh(): minor-allele count >= 50 in cases and in controls,
    |predicted log OR| <= 1.999;
  * slope = coefficient of stats::lm(predicted ~ reported), i.e. OLS WITH an intercept, unweighted.
Not reproduced: exclusion of variants in the package's 1000 Genomes MAF reference set (rsid-based,
not applicable to simulated variants).
The decision rule "slope > 1.20 or < 0.8" is taken from Haycock et al. 2023 (main text); the
package code reports the slope and draws the plot but applies no threshold.
"""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np, pandas as pd
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sim_scale as ss, run_sim as rs

REPS = 200        # first 200 test replicates of every test cell
OUTD = os.path.join(HERE, "posthoc_out")


def _se_sh(x, n0, n1, n2, N, odds):
    m = n1 + 2 * n2
    E0 = np.exp(x * m / N); E1 = np.exp(x * (m - N) / N); E2 = np.exp(x * (m - 2 * N) / N)
    t0 = n0 * odds * E0 / (odds + E0) ** 2; t1 = n1 * odds * E1 / (odds + E1) ** 2; t2 = n2 * odds * E2 / (odds + E2) ** 2
    A = t0 + t1 + t2; B = t1 + 4 * t2; C = (t1 + 2 * t2) ** 2
    return np.sqrt(A / (A * B - C))


def predict_lnor_sh(z, eaf, ncase, ncontrol):
    maf = np.where(eaf > 0.5, 1 - eaf, eaf)
    N = ncase + ncontrol; odds = ncase / (N - ncase)
    n0, n1, n2 = N * maf ** 2, N * 2 * maf * (1 - maf), N * (1 - maf) ** 2
    sgn = np.where(z >= 0, 1.0, -1.0)
    out = np.empty_like(z)
    coarse = np.arange(1, 6001) * 1e-3
    for j in range(len(z)):
        x = sgn[j] * coarse
        y = np.abs(x / _se_sh(x, n0[j], n1[j], n2[j], N, odds) - z[j])
        c = coarse[np.nanargmin(y)]
        fine = np.arange(max(c - 2e-3, 1e-6), c + 2e-3, 1e-6)
        x = sgn[j] * fine
        y = np.abs(x / _se_sh(x, n0[j], n1[j], n2[j], N, odds) - z[j])
        mn = np.nanmin(y)
        out[j] = x[y == mn].mean()
    return out, maf


def checksumstats(beta, se, eaf, ncase, ncontrol):
    z = beta / se
    pred, maf = predict_lnor_sh(z, eaf, ncase, ncontrol)
    keep = (maf * ncase * 2 >= 50) & (maf * ncontrol * 2 >= 50) & (np.abs(pred) <= 1.999)
    n = int(keep.sum())
    if n < 3:
        return dict(css_n_used=n, css_slope=np.nan, css_intercept=np.nan)
    X = np.c_[np.ones(n), beta[keep]]
    coef = np.linalg.lstsq(X, pred[keep], rcond=None)[0]
    return dict(css_n_used=n, css_slope=coef[1], css_intercept=coef[0])


def one(setname, cell, rep, seq):
    """Mirror of run_sim.one_replicate, returning the same records plus CheckSumStats fields."""
    rng = np.random.default_rng(seq)
    K = cell["K"]; gcell = dict(cell); gcell["K"] = 3 * K
    d = ss.generate(rng, gcell)
    out, fl = ss.analyse(d)
    lp = ss.pipeline_logp(out, fl["logp_spa"])
    sups = [("case_fraction", d["p_real"])]
    if cell.get("design") == "casecontrol":
        sups.append(("prevalence", cell["prev"]))
    nvar = d["m"].sum(0)
    ncase = d["p_real"] * d["N_rep"]; nctrl = (1 - d["p_real"]) * d["N_rep"]
    rows = []
    for pl in ss.PIPELINES:
        idx, npass = ss.select(lp[pl], K)
        se = out[pl][1][idx]; f = d["f_rep"][idx]; b = out[pl][0][idx]
        css = checksumstats(b, se, f, ncase, nctrl) if len(idx) >= 1 else dict(css_n_used=0, css_slope=np.nan, css_intercept=np.nan)
        vf = nvar[idx] / cell["N"]
        for supname, psup in sups:
            r = ss.classify(se, f, d["N_rep"], psup)
            r.update(set=setname, cell=cell["cell"], arm=cell["arm"], level=cell["level"], p_target=cell["p"],
                     K=K, rep=rep, pipeline=pl, true_scale=ss.TRUE_SCALE[pl], p_supplied_as=supname,
                     p_supplied=psup, p_real=d["p_real"], n_pass=npass, **css,
                     varN_frac_median=float(np.median(vf)) if len(vf) else np.nan,
                     varN_frac_min=float(vf.min()) if len(vf) else np.nan,
                     varN_frac_full=float((vf == 1).mean()) if len(vf) else np.nan,
                     varN_fracs=";".join(f"{v:.4f}" for v in vf))
            rows.append(r)
    return rows


def task_list():
    T = []
    for setname, cell, reps, seq in rs.seed_plan():
        if setname != "test":
            continue
        kids = seq.spawn(reps)
        for a in range(0, REPS, 50):
            T.append((cell, a, min(a + 50, REPS), kids[a:min(a + 50, REPS)]))
    return T


def run(task):
    cell, a, b, seqs = task
    fn = os.path.join(OUTD, f"{cell['cell']}__{a:04d}.parquet")
    if os.path.exists(fn):
        return fn
    rows = []
    for i, s in zip(range(a, b), seqs):
        rows += one("test", cell, i, s)
    pd.DataFrame(rows).to_parquet(fn + ".tmp", index=False); os.replace(fn + ".tmp", fn)
    return fn


if __name__ == "__main__":
    os.makedirs(OUTD, exist_ok=True)
    T = task_list(); T.sort(key=lambda t: -t[0]["K"] - (20 if t[0].get("design") == "casecontrol" else 0))
    done = 0
    with Pool(10) as p:
        for _ in p.imap_unordered(run, T):
            done += 1
            if done % 25 == 0 or done == len(T):
                print(f"{done}/{len(T)}", flush=True)
    print("complete")
