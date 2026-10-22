#!/usr/bin/env python3
"""Checks for the figure materials added in v1.2.

1. Main-text Figure 3 (figures/submitted/fig3_plotted_values.csv) plots exactly the four rows of
   case_study/bridge/bridge_estimates.csv, and its three non-bridge rows equal the v1.1 figure source
   data (data/repository/Tier2_scale_correction_source_data.csv).
2. Online Resource 1, Fig. S1: the IVW estimates of the original and rescaled Tier 2 analysis, recomputed
   from data/mr/, equal results/invariance_estimator_unit_test.csv; every leave-one-out estimate of the
   rescaled analysis is the original divided by lambda, with the same P value.
3. The leave-one-out estimates drawn by figures/code/make_mr_plots.R, recomputed here, equal the committed
   TwoSampleMR leave-one-out output in case_study/leave_one_out/ for all six analyses.
"""
import csv, math, os
import numpy as np
from scipy.stats import norm

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
F = lambda *p: os.path.join(ROOT, *p)
rows = lambda *p: list(csv.DictReader(open(F(*p), newline="", encoding="utf-8-sig")))
LAMBDA = 1 / (0.3450365567 * (1 - 0.3450365567))   # as in effect_scale_mr_analysis.R, section 4


def ivw(bx, by, sy):
    w = 1 / sy ** 2
    b = np.sum(w * bx * by) / np.sum(w * bx ** 2)
    sig = math.sqrt(np.sum(w * (by - b * bx) ** 2) / (len(bx) - 1))
    se = math.sqrt(1 / np.sum(w * bx ** 2)) * max(1.0, sig)
    return b, se, 2 * norm.sf(abs(b / se))


def harmonised(name):
    d = [r for r in rows("data", "mr", name + "_harmonised_kept.csv") if r["mr_keep"].upper() == "TRUE"]
    bx = np.array([float(r["beta.exposure"]) for r in d]); by = np.array([float(r["beta.outcome"]) for r in d])
    sy = np.array([float(r["se.outcome"]) for r in d]); s = np.sign(bx)
    return [r["SNP"] for r in d], bx * s, by * s, sy


# 1. Figure 3
fig3 = rows("figures", "submitted", "fig3_plotted_values.csv")
bridge = rows("case_study", "bridge", "bridge_estimates.csv")
assert len(fig3) == 4 == len(bridge)
for a, b in zip(fig3, bridge):
    for x, y in (("OR", "OR"), ("low", "CI_low"), ("high", "CI_high")):
        assert abs(float(a[x]) - float(b[y])) < 1e-12, (a["analysis"], x)
assert [r["post_hoc"] for r in fig3] == ["FALSE", "FALSE", "TRUE", "FALSE"]
old = rows("data", "repository", "Tier2_scale_correction_source_data.csv")
# The v1.1 source data used the multiplier 1.96 for the interval and the bridge output uses the normal
# quantile 1.959964; converting the latter to 1.96 reproduces the v1.1 values (identical to two decimals).
for a, b in zip([fig3[0], fig3[1], fig3[3]], old):
    lb = math.log(float(a["OR"])); se = (math.log(float(a["high"])) - lb) / norm.ppf(0.975)
    assert abs(float(a["OR"]) - float(b["OR"])) < 1e-9, b["analysis"]
    assert abs(math.exp(lb - 1.96 * se) - float(b["low"])) < 1e-6 and abs(math.exp(lb + 1.96 * se) - float(b["high"])) < 1e-6, b["analysis"]
    for x in ("OR", "low", "high"):
        assert f"{float(a[x]):.2f}" == f"{float(b[x]):.2f}", (b["analysis"], x)

# 2. Fig. S1
snp, bx, by, sy = harmonised("Tier_2_superseded_mixed_scale__UK_Biobank_asthma")
unit = {r["method"]: r for r in rows("results", "invariance_estimator_unit_test.csv")}
ivw_row = unit["Inverse variance weighted (multiplicative random effects)"]
b0, se0, p0 = ivw(bx, by, sy); b1, se1, p1 = ivw(bx * LAMBDA, by, sy)
assert abs(b0 - float(ivw_row["b_original"])) < 1e-10 and abs(se0 - float(ivw_row["se_original"])) < 1e-10
assert abs(b1 - float(ivw_row["b_rescaled"])) < 1e-10 and abs(se1 - float(ivw_row["se_rescaled"])) < 1e-10
for i in range(len(bx)):
    k = np.arange(len(bx)) != i
    lo = ivw(bx[k], by[k], sy[k]); lr = ivw(bx[k] * LAMBDA, by[k], sy[k])
    assert abs(lo[0] / lr[0] - LAMBDA) < 1e-9 and abs(lo[2] - lr[2]) < 1e-12, snp[i]
check = rows("figures", "online_resource_1", "mr_plots_check.csv")
assert all(r["ivw_matches_results"] == "TRUE" for r in check)

# 3. Leave-one-out against the committed TwoSampleMR output
n = 0
for name in ["Tier_1_affective_instability__FinnGen_asthma", "Tier_1_affective_instability__UK_Biobank_asthma",
             "Tier_2_help_seeking_corrected__FinnGen_asthma", "Tier_2_help_seeking_corrected__UK_Biobank_asthma",
             "Tier_3_major_depression_meta_analysis__FinnGen_asthma", "Tier_3_major_depression_meta_analysis__UK_Biobank_asthma"]:
    snp, bx, by, sy = harmonised(name)
    ref = {r["SNP"]: r for r in rows("case_study", "leave_one_out", name + "__leaveoneout.csv") if r["SNP"] != "All"}
    for i, s in enumerate(snp):
        k = np.arange(len(bx)) != i
        b, se, _ = ivw(bx[k], by[k], sy[k])
        assert abs(b - float(ref[s]["b"])) < 1e-10 and abs(se - float(ref[s]["se"])) < 1e-10, (name, s)
        n += 1
print(f"Figure validation passed: Figure 3 (4 rows), Fig. S1 (35 leave-one-out pairs), {n} leave-one-out estimates.")
