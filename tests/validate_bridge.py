#!/usr/bin/env python3
"""Recompute the post hoc bridge comparison in memory and check it against the committed outputs.

Does not write files, so the committed CSVs are not regenerated in CI. Floating-point values are
compared with a tolerance; counts are compared exactly.
"""
import importlib.util, os, sys
import pandas as pd
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
spec = importlib.util.spec_from_file_location("bridge", os.path.join(ROOT, "case_study", "bridge", "bridge_analysis.py"))
mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
out, meta = mod.compute()
ref = pd.read_csv(os.path.join(ROOT, "case_study", "bridge", "bridge_estimates.csv"))
assert list(out.analysis) == list(ref.analysis)
for col in ["OR", "CI_low", "CI_high", "P", "Q", "I2"]:
    d = (out[col] - ref[col]).abs().max()
    assert d < 1e-9 * max(1.0, ref[col].abs().max()), (col, d)
assert list(out.n_variants) == [35, 35, 35, 26]
m = pd.read_csv(os.path.join(ROOT, "case_study", "bridge", "bridge_variant_summary.csv")).iloc[0].to_dict()
assert {k: int(v) for k, v in m.items()} == meta, (m, meta)
assert meta == dict(variants_available=35, allele_aligned=35, palindromic=0,
                    genome_wide_significant_in_replacement=19, shared_with_replacement_instrument_set=7)
assert (out.Q / (out.n_variants - 1) > 1).all(), "dispersion floor would be active"
r = {round(x, 2) for x in out.OR}
assert r == {3.12, 1.29, 1.28, 1.07}, r                      # manuscript Sections 3.7-3.8; OR1 Table S22
assert round(out.CI_low[2], 2) == 1.12 and round(out.CI_high[2], 2) == 1.46

# The Python estimator must equal the pinned R (TwoSampleMR) estimator on the three rows R already produced,
# so the bridge differs from those analyses only in its input exposure coefficients.
import numpy as np
inv = pd.read_csv(os.path.join(ROOT, "results", "invariance_estimator_unit_test.csv"))
inv = inv[inv.method.str.contains("multiplicative random effects")].iloc[0]
fe = pd.read_csv(os.path.join(ROOT, "results", "MR_rerun_five_estimators.csv"))
rep = fe[fe.method.str.contains("multiplicative random effects") & (fe.model == "Help-seeking (corrected) | UK Biobank asthma")].iloc[0]
se = (np.log(out.CI_high) - np.log(out.OR)) / 1.959964
for k, (b_r, se_r) in {0: (inv.b_original, inv.se_original), 1: (inv.b_rescaled, inv.se_rescaled), 3: (rep.b, rep.se)}.items():
    assert abs(np.log(out.OR[k]) - b_r) < 1e-6 and abs(se[k] - se_r) < 1e-6, (k, np.log(out.OR[k]), b_r, se[k], se_r)
print("bridge comparison reproduced (estimator equal to the pinned R output on the original, rescaled and replacement analyses): 3.12 -> 1.29 (rescaled), 1.28 (same 35 variants, log-odds coefficients), 1.07 (26-variant set)")
