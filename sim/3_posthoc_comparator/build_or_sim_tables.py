"""
Builds Online Resource Tables S16-S20 from saved CSVs (no simulation is run, no frozen file is written).
S16, S17  - frozen, pre-specified simulation (from gate_unconditional_decomposition.csv)
S18-S20   - POST HOC comparator (from posthoc_*.csv written by posthoc_comparator_tables.py)
    python sim/build_or_sim_tables.py <out_dir>
"""
import os, sys
import numpy as np, pandas as pd
H = os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[1] if len(sys.argv) > 1 else "."
GROUP = {"P1": "P1 linear probability", "P2": "P2 standardised", "P3": "P3-P7 log odds", "P4": "P3-P7 log odds",
         "P5": "P3-P7 log odds", "P6": "P3-P7 log odds", "P7": "P3-P7 log odds"}
f3 = lambda v: "" if pd.isna(v) else f"{v:.3f}"

U = pd.read_csv(os.path.join(H, "gate_unconditional_decomposition.csv"))
U["wrong"] = U[["to_linear", "to_standardised", "to_log_odds"]].sum(1)
prim = U[U.rule.str.startswith("primary")]; adq = U[U.rule.str.startswith("secondary")]

# S16 core grid, primary rule, by K x pipeline, pooled over case fraction; adequacy abstention alongside
c = prim[prim.arm == "core"].groupby(["K", "pipeline"])[["correct", "wrong", "abstain_sufficiency"]].mean()
ca = adq[adq.arm == "core"].groupby(["K", "pipeline"])["abstain_adequacy"].mean()
S16 = c.join(ca).reset_index()
S16.columns = ["K", "pipeline", "correct", "wrong_any_scale", "abstain_sufficiency", "adequacy_alarm_on_correct_scale"]
for k in S16.columns[2:]: S16[k] = S16[k].map(f3)
S16.to_csv(os.path.join(OUT, "or1_tableS16_sim_core.csv"), index=False)

# S17 non-core arms, K = 10: by arm, level, supplied p, pipeline group; mean over case fraction and worst cell
x = prim[prim.arm != "core"].copy(); xa = adq[adq.arm != "core"].copy()
for d in (x, xa): d["group"] = d.pipeline.map(GROUP); d["p_supplied_as"] = d.p_supplied_as.fillna("case_fraction")
k = ["arm", "level", "p_supplied_as", "group"]
a = x.groupby(k).agg(primary_wrong_mean=("wrong", "mean"), primary_wrong_worst_cell=("wrong", "max"),
                     abstain_sufficiency_mean=("abstain_sufficiency", "mean"))
b = xa.groupby(k).agg(adequacy_wrong_worst_cell=("wrong", "max"), adequacy_alarm_mean=("abstain_adequacy", "mean"))
S17 = a.join(b).reset_index()
for kk in S17.columns[4:]: S17[kk] = S17[kk].map(f3)
S17.to_csv(os.path.join(OUT, "or1_tableS17_sim_arms.csv"), index=False)

# S18 post hoc, declared-log-odds frame: anomaly pairing and attribution pairing
A = pd.read_csv(os.path.join(H, "posthoc_logodds_anomaly_paired.csv"))
B = pd.read_csv(os.path.join(H, "posthoc_logodds_attribution_paired.csv"))
S18 = A.merge(B, on=["setting", "grp"])[["setting", "grp", "css_alarm", "css_na", "adq_alarm", "sig_na",
                                         "css_lab", "css_wronglab", "sig_lab", "sig_wronglab"]]
S18.columns = ["setting", "releases", "CSS_slope_alarm", "CSS_not_assessable", "3sig_adequacy_alarm", "3sig_not_assessable",
               "CSS_slope_implied_scale_correct", "CSS_slope_implied_scale_wrong", "3sig_primary_correct", "3sig_primary_wrong"]
for kk in S18.columns[2:]: S18[kk] = S18[kk].map(f3)
S18.to_csv(os.path.join(OUT, "or1_tableS18_posthoc_logodds_frame.csv"), index=False)

# S19 post hoc, declared-standardised frame
D = pd.read_csv(os.path.join(H, "posthoc_declaredSD_comparison.csv"))
S19 = D[["setting", "grp", "css_declSD_alarm", "css_na", "css_selfnorm_alarm", "sig_primary_alarm", "sig_adequacy_alarm", "sig_na", "css_slope_median"]].copy()
S19.columns = ["setting", "releases", "CSS_declaredSD_alarm", "CSS_not_assessable", "CSS_selfnormalised_alarm",
               "3sig_primary_alarm", "3sig_adequacy_alarm", "3sig_not_assessable", "CSS_slope_median"]
for kk in S19.columns[2:]: S19[kk] = S19[kk].map(f3)
S19.to_csv(os.path.join(OUT, "or1_tableS19_posthoc_declaredSD_frame.csv"), index=False)

# S20 post hoc, estimator agreement
F = pd.read_csv(os.path.join(H, "posthoc_factor_quantiles.csv"), index_col=0).T
F.index.name = "quantity"; F = F.reset_index()
for kk in ["min", "median", "p99", "max"]: F[kk] = F[kk].map(lambda v: f"{v:.5f}")
F.to_csv(os.path.join(OUT, "or1_tableS20_posthoc_estimator_agreement.csv"), index=False)
print({n: len(t) for n, t in [("S16", S16), ("S17", S17), ("S18", S18), ("S19", S19), ("S20", F)]})
