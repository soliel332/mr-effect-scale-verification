#!/usr/bin/env python3
"""Post hoc diagnostic bridge between the arithmetic rescaling and the replacement-release analysis.

usage: python case_study/bridge/bridge_analysis.py   (run from the repository root)

Question: does the move from the rescaled estimate (1.29) to the replacement-release estimate (1.07)
occur at the change of scale or at the change of instrument set?

The bridge keeps the 35 variants, outcome associations (UK Biobank asthma), outcome standard errors,
allele harmonisation and estimator of the original Tier 2 analysis, and substitutes only the exposure
coefficients of the replacement release (ebi-a-GCST90013959, log odds) for those of ukb-b-6991. The
replacement-release exposure coefficients for these 35 variants are taken from the calibration file
data/calibration/GP_visit_MRC_IEU_SNP_data.csv, in which the replacement release is the second
('outcome') side and alleles are already aligned to ukb-b-6991.

Estimator: first-order inverse-variance-weighted, multiplicative random effects, normal reference
distribution. This implementation floors the dispersion Q/(K-1) at 1; TwoSampleMR::mr_ivw_mre at the pinned
commit does not. The floor is inactive in all four analyses here (Q/(K-1) > 1; checked in
tests/validate_bridge.py), so the two should give identical estimates; this is checked in the pinned R
environment by case_study/bridge/verify_bridge_pinned.R. In this estimator the exposure standard errors
do not enter, so only the exposure coefficients are substituted.

This analysis was added post hoc, after the replacement-release result was known. It is a diagnostic
comparison of two operations, not a new primary MR estimate.
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
F = lambda *p: os.path.join(ROOT, *p)
P_CASE = 0.3450365567  # case fraction of ukb-b-6991, as in effect_scale_mr_analysis.R
LAMBDA = 1 / (P_CASE * (1 - P_CASE))  # 4.4250 (manuscript Section 3.7); same expression as the core script


def ivw(bx, by, sy):
    w = 1 / sy ** 2
    b = np.sum(w * bx * by) / np.sum(w * bx ** 2)
    res = (by - b * bx) / sy
    q = np.sum(res ** 2); k = len(by)
    phi = max(1.0, q / (k - 1))
    se = np.sqrt(phi / np.sum(w * bx ** 2))
    i2 = max(0.0, (q - (k - 1)) / q) if q > 0 else 0.0
    return b, se, 2 * norm.sf(abs(b / se)), q, i2


def row(label, variants, exposure, b, se, p, q, i2):
    return dict(analysis=label, n_variants=variants, exposure_coefficients=exposure,
                OR=np.exp(b), CI_low=np.exp(b - 1.959964 * se), CI_high=np.exp(b + 1.959964 * se),
                P=p, Q=q, I2=i2)


def compute():
    orig = pd.read_csv(F("data", "mr", "Tier_2_superseded_mixed_scale__UK_Biobank_asthma_harmonised_kept.csv"))
    repl = pd.read_csv(F("data", "mr", "Tier_2_help_seeking_corrected__UK_Biobank_asthma_harmonised_kept.csv"))
    cal = pd.read_csv(F("data", "calibration", "GP_visit_MRC_IEU_SNP_data.csv")).set_index("SNP")
    assert set(cal["id.outcome"]) == {"ebi-a-GCST90013959"}
    assert len(orig) == 35 and orig.SNP.is_unique
    missing = sorted(set(orig.SNP) - set(cal.index)); assert not missing, missing
    c = cal.loc[orig.SNP]
    same = (c["effect_allele.outcome"].values == orig["effect_allele.exposure"].values) & \
           (c["other_allele.outcome"].values == orig["other_allele.exposure"].values)
    assert same.all(), "allele orientation differs; harmonisation would be required"
    pal = [{a, b} in ({"A", "T"}, {"C", "G"}) for a, b in zip(orig["effect_allele.exposure"], orig["other_allele.exposure"])]
    bx_repl = c["beta.outcome"].values
    p_repl = c["pval.outcome"].values
    by, sy = orig["beta.outcome"].values, orig["se.outcome"].values
    rows = [
        row("Original analysis, linear-probability coefficients read as log odds", 35, "ukb-b-6991",
            *ivw(orig["beta.exposure"].values, by, sy)),
        row("Arithmetic re-expression (estimate divided by lambda)", 35, "ukb-b-6991 / lambda",
            *ivw(orig["beta.exposure"].values * LAMBDA, by, sy)),
        row("Bridge: same 35 variants, replacement-release exposure coefficients", 35, "ebi-a-GCST90013959",
            *ivw(bx_repl, by, sy)),
        row("Replacement release, own instrument set", len(repl), "ebi-a-GCST90013959",
            *ivw(repl["beta.exposure"].values, repl["beta.outcome"].values, repl["se.outcome"].values)),
    ]
    out = pd.DataFrame(rows)
    meta = dict(variants_available=int(len(orig) - len(missing)), allele_aligned=int(same.sum()),
                palindromic=int(sum(pal)), genome_wide_significant_in_replacement=int((p_repl < 5e-8).sum()),
                shared_with_replacement_instrument_set=int(len(set(orig.SNP) & set(repl.SNP))))
    return out, meta


def main():
    out, meta = compute()
    out.to_csv(F("case_study", "bridge", "bridge_estimates.csv"), index=False)
    pd.DataFrame([meta]).to_csv(F("case_study", "bridge", "bridge_variant_summary.csv"), index=False)
    print(out.round(4).to_string(index=False)); print(meta)


if __name__ == "__main__":
    main()
