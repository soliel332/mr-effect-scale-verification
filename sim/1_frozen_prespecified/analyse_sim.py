"""
Analysis for SIMULATION_PROTOCOL.md, written before any result was inspected.
Implements the estimands of section 7, the adequacy threshold of section 3.6 (tuning set only),
and the pre-specified consequences R1-R5 of section 9.
"""
import glob, os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CRIT = 0.95


def load():
    fs = sorted(glob.glob(os.path.join(HERE, "out", "*.parquet")))
    return pd.concat([pd.read_parquet(f) for f in fs], ignore_index=True)


def tau_from_tuning(df):
    t = df[(df.set == "tuning") & df.pipeline.isin(["P1", "P2", "P3"]) & df.sigma_hat.notna()]
    col = {"linear_probability": "d_lin", "standardised": "d_std", "log_odds": "d_lo"}
    d_true = np.where(t.true_scale == "linear_probability", t.d_lin,
             np.where(t.true_scale == "standardised", t.d_std, t.d_lo))
    return float(np.quantile(d_true, 0.99))


def estimands(df, tau):
    x = df[df.set == "test"].copy()
    x["assigned_adequacy"] = np.where(x.assigned.str.startswith("abstain"), x.assigned,
                              np.where(x.d_min > tau, "abstain_adequacy", x.assigned))
    rows = []
    keys = ["arm", "level", "cell", "p_target", "K", "p_supplied_as", "pipeline", "true_scale"]
    for k, g in x.groupby(keys):
        r = dict(zip(keys, k)); n = len(g)
        for rule, col in (("primary", "assigned"), ("adequacy", "assigned_adequacy")):
            a = g[col]; nab = a.str.startswith("abstain").sum(); cls = a[~a.str.startswith("abstain")]
            acc = (cls == r["true_scale"]).mean() if len(cls) else np.nan
            r[f"{rule}_n"] = n
            r[f"{rule}_abstain"] = nab / n
            r[f"{rule}_accuracy"] = acc
            r[f"{rule}_mcse"] = np.sqrt(acc * (1 - acc) / max(len(cls), 1)) if len(cls) else np.nan
            for s in ("linear_probability", "standardised", "log_odds"):
                r[f"{rule}_to_{s}"] = (cls == s).mean() if len(cls) else np.nan
        r["abstain_sufficiency"] = (g.assigned == "abstain_sufficiency").mean()
        r["median_n_pass"] = g.n_pass.median()
        r["median_sigma_hat"] = g.sigma_hat.median()
        r["median_rel_sep"] = g.rel_sep.median()
        r["conv_logit"] = g.conv_logit.mean(); r["conv_firth"] = g.conv_firth.mean()
        r["spa_ok"] = g.spa_ok_frac.mean()
        rows.append(r)
    return pd.DataFrame(rows), x


def cell_pass(oc, rule="primary"):
    g = oc.groupby(["arm", "level", "cell", "p_target", "K", "p_supplied_as"])
    return g[f"{rule}_accuracy"].agg(lambda s: bool((s >= CRIT).all()) if s.notna().all() else False).rename("meets")


def claim_rules(oc):
    cp = cell_pass(oc).reset_index()
    core = cp[cp.arm == "core"]
    arm = lambda a, sup="case_fraction": cp[(cp.arm == a) & (cp.p_supplied_as == sup)]
    out = []
    r1 = core.meets.all() and all(arm(a).meets.all() for a in ("A1", "A2", "A3", "A4"))
    out.append(("R1 broad support", bool(r1),
                f"core cells meeting criterion {int(core.meets.sum())}/{len(core)}; "
                + "; ".join(f"{a} {int(arm(a).meets.sum())}/{len(arm(a))}" for a in ("A1", "A2", "A3", "A4"))))
    failing = core[~core.meets][["cell"]].cell.tolist() + sum(
        (arm(a)[~arm(a).meets].cell.tolist() for a in ("A1", "A2", "A3", "A4")), [])
    out.append(("R2 bounded support — failing cells", bool(len(failing) > 0), "; ".join(failing) or "none"))
    c56 = oc[(oc.arm == "core") & oc.pipeline.isin(["P5", "P6"])]
    f56 = c56[c56.primary_accuracy < CRIT]
    out.append(("R3 saddlepoint pipelines fail in a core cell", bool(len(f56) > 0),
                "; ".join(f"{r.pipeline} {r.cell} acc={r.primary_accuracy:.3f}" for r in f56.itertuples()) or "none"))
    mis = oc[oc.arm.isin(["A5", "A6", "A7", "A8"])].copy()
    mis["misattribution"] = 1 - mis.primary_accuracy
    bad = mis[mis.misattribution > 0.05]
    out.append(("R4 misattribution > 0.05", bool(len(bad) > 0),
                "; ".join(f"{r.level} {r.pipeline} p={r.p_target} mis={r.misattribution:.3f}" for r in bad.itertuples()) or "none"))
    a4a, a4b = arm("A4", "case_fraction"), arm("A4", "prevalence")
    out.append(("R5 prevalence fails where case fraction passes", bool(a4a.meets.all() and not a4b.meets.all()),
                f"case fraction {int(a4a.meets.sum())}/{len(a4a)}; prevalence {int(a4b.meets.sum())}/{len(a4b)}"))
    return pd.DataFrame(out, columns=["rule", "triggered", "detail"])


def separation_calibration(x):
    c = x[(x.arm == "core") & x.rel_sep.notna() & ~x.assigned.str.startswith("abstain")].copy()
    c["correct"] = c.assigned == c.true_scale
    c["decile"] = pd.qcut(c.rel_sep.rank(method="first"), 10, labels=False)
    return c.groupby("decile").agg(rel_sep_lo=("rel_sep", "min"), rel_sep_hi=("rel_sep", "max"),
                                   n=("correct", "size"), accuracy=("correct", "mean")).reset_index()


if __name__ == "__main__":
    df = load()
    tau = tau_from_tuning(df)
    oc, x = estimands(df, tau)
    oc.to_csv(os.path.join(HERE, "sim_operating_characteristics.csv"), index=False)
    cr = claim_rules(oc); cr.to_csv(os.path.join(HERE, "sim_claim_rules.csv"), index=False)
    sc = separation_calibration(x); sc.to_csv(os.path.join(HERE, "sim_separation_calibration.csv"), index=False)
    df.to_parquet(os.path.join(HERE, "sim_results_records.parquet"), index=False)
    with open(os.path.join(HERE, "sim_tau.txt"), "w") as f:
        f.write(f"adequacy threshold tau (99th percentile, tuning set, P1-P3 pooled over core grid) = {tau:.6f}\n")
    print("records:", len(df), "| tau:", round(tau, 4))
    print(cr.to_string(index=False))
