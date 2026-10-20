#!/usr/bin/env python3
"""Recompute every claim in results/claim_evidence_manifest.csv from its evidence file.

Each claim has a recomputation below; a claim without one, or whose recomputed value differs from the
manifest, fails the test. The manuscript locations are listed in the manifest itself.
"""
import csv
import math
import os

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def f(rel):
    return os.path.join(ROOT, rel)


def check_file(rel, key, col="value"):
    return float(pd.read_csv(f(rel)).set_index(pd.read_csv(f(rel)).columns[0]).loc[key, col])


def core_grid(kind):
    u = pd.read_csv(f("sim/2_gate_analysis/gate_unconditional_decomposition.csv"))
    u = u[(u["arm"] == "core") & (u["rule"] == "primary (forced nearest)")]
    if kind == "n":
        return float(u["n"].sum())
    return float((u["n"] * (1 - u["correct"] - u["abstain_sufficiency"] - u["abstain_adequacy"])).sum().round(9))


def presso(col):
    p = pd.read_csv(f("case_study/pinned_mrpresso/phase2_mrpresso_pinned.csv")).set_index("model")
    return p.loc["Tier_2_help_seeking_corrected__UK_Biobank_asthma", col]


def source(col):
    r = pd.read_csv(f("scheme_c/source_version_check/Scheme_C_source_version_reverse_check.csv")).set_index("analysis_label")
    a = r.loc["ADHD"]
    return float(a["alleles_aligned_n"]) if col == "n" else float(round(a["beta_within_1e-4_fraction"] * a["alleles_aligned_n"]))


def summary(metric):
    s = pd.read_csv(f("results/repository_scale_summary.csv")).set_index("metric")
    return float(s.loc[metric, "value"])


def calib(i):
    return float(pd.read_csv(f("results/calibration_reanalysis.csv"))["observed_slope"].iloc[i])


def tier2(row):
    d = pd.read_csv(f("data/repository/Tier2_scale_correction_source_data.csv"))
    return float(d["OR"].iloc[row])


RECOMPUTE = {
    "Recovered stratified repository sample": lambda: summary("reconstructed_stratified_sample_size"),
    "Archived classifiable repository set": lambda: summary("classifiable_datasets_in_preserved_output"),
    "SD-labelled discordance count": lambda: summary("declared_SD_and_inferred_linear_probability"),
    "Frozen SD-label carrier count": lambda: float(len(pd.read_csv(f("data/repository/SD_label_carriers.csv")))),
    "Help-seeking calibration slope (MRC-IEU)": lambda: calib(0),
    "Help-seeking calibration slope (Neale)": lambda: calib(1),
    "Mood-swings calibration slope": lambda: calib(2),
    "Tier 2 mixed-scale UKB OR": lambda: tier2(0),
    "Tier 2 arithmetic-rescaling OR": lambda: tier2(1),
    "Tier 2 replacement-release UKB OR": lambda: tier2(2),
    "Core-grid simulated releases": lambda: core_grid("n"),
    "Core-grid wrong scale attributions": lambda: core_grid("wrong"),
    "Adequacy threshold tau": lambda: float(open(f("sim/1_frozen_prespecified/sim_tau.txt")).read().split("=")[-1]),
    "Comparator test cells": lambda: float(pd.read_csv(f("sim/3_posthoc_comparator/posthoc_sd_regeneration_check.csv"), index_col=0).loc["cells", "value"]),
    "Comparator regenerated records, declared-log-odds frame": lambda: check_file("sim/3_posthoc_comparator/posthoc_regeneration_check.csv", "regenerated_records"),
    "Comparator label mismatches, declared-log-odds frame": lambda: check_file("sim/3_posthoc_comparator/posthoc_regeneration_check.csv", "label_mismatch"),
    "Comparator regenerated records, declared-standardised frame": lambda: check_file("sim/3_posthoc_comparator/posthoc_sd_regeneration_check.csv", "records"),
    "Comparator label mismatches, declared-standardised frame": lambda: check_file("sim/3_posthoc_comparator/posthoc_sd_regeneration_check.csv", "label_mismatch"),
    "Tier 2 UKB MR-PRESSO outliers": lambda: float(presso("n_outliers")),
    "Tier 2 UKB MR-PRESSO outlier-corrected OR": lambda: float(np.exp(presso("corrected_estimate"))),
    "Tier 2 UKB MR-PRESSO outlier-corrected P": lambda: float(presso("corrected_p")),
    "Tier 2 UKB MR-PRESSO distortion P": lambda: float(presso("distortion_p")),
    "Tier 2 UKB minimum leave-one-out P": lambda: float(pd.read_csv(f("case_study/leave_one_out/Tier_2_help_seeking_corrected__UK_Biobank_asthma__leaveoneout.csv")).query("SNP != 'All'")["p"].min()),
    "Audit assignments changed under reported N x 0.8-1.2": lambda: float(pd.read_csv(f("audit/three_signature/or1_tableS14_N_perturbation_p_fixed.csv"))["changed"].sum()),
    "ADHD source-version effects within 1e-4": lambda: source("within"),
    "ADHD source-version overlapping variants": lambda: source("n"),
}

def prior(col):
    return float(pd.read_csv(f("prior_studies/reexpressed_estimates.csv")).set_index("exposure_accession").loc["ukb-b-6991", col])


RECOMPUTE.update({
    "Prior study (Yan 2024) case-fraction factor lambda": lambda: prior("lambda_full"),
    "Prior study (Yan 2024) re-expressed OR": lambda: prior("reexpressed_or_full"),
    "Prior study (Yan 2024) re-expressed CI lower": lambda: prior("reexpressed_ci_low_full"),
    "Prior study (Yan 2024) re-expressed CI upper": lambda: prior("reexpressed_ci_high_full"),
    "Prior study (Yan 2024) change in Wald z": lambda: prior("z_after") - prior("z_before"),
})


def bridge(col):
    t = pd.read_csv(f("case_study/bridge/bridge_estimates.csv"))
    return float(t.loc[t.analysis.str.startswith("Bridge"), col].iloc[0])


RECOMPUTE.update({
    "Bridge comparison OR (same 35 variants)": lambda: bridge("OR"),
    "Bridge comparison CI lower": lambda: bridge("CI_low"),
    "Bridge comparison CI upper": lambda: bridge("CI_high"),
    "Bridge variants allele-aligned": lambda: float(pd.read_csv(f("case_study/bridge/bridge_variant_summary.csv")).allele_aligned.iloc[0]),
})


rows = list(csv.DictReader(open(f("results/claim_evidence_manifest.csv"), newline="")))
missing = [r["claim"] for r in rows if r["claim"] not in RECOMPUTE]
assert not missing, f"claims without a recomputation: {missing}"
for r in rows:
    assert os.path.isfile(f(r["evidence_file"])), r["evidence_file"]
    got, want = RECOMPUTE[r["claim"]](), float(r["value"])
    assert math.isclose(got, want, rel_tol=1e-9, abs_tol=1e-12), (r["claim"], got, want)
print(f"claim manifest: all {len(rows)} claims recomputed from their evidence files")
