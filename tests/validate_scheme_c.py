#!/usr/bin/env python3
"""Claim-level validation for the frozen Scheme C submission outputs."""

from __future__ import annotations

import csv
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCHEME = ROOT / "scheme_c"


def rows(relative: str) -> list[dict[str, str]]:
    path = SCHEME / relative
    if not path.is_file():
        raise AssertionError(f"missing required file: {path}")
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def close(value: str, expected: float, tolerance: float) -> None:
    actual = float(value)
    if not math.isfinite(actual) or abs(actual - expected) > tolerance:
        raise AssertionError(f"expected {expected} +/- {tolerance}, got {actual}")


pairwise = rows("results/genomic/Scheme_C_GenomicSEM_pairwise_rg.csv")
assert len(pairwise) == 21, f"expected 21 pairwise correlations, got {len(pairwise)}"


def pair(a: str, b: str) -> dict[str, str]:
    for record in pairwise:
        if {record["trait_1"], record["trait_2"]} == {a, b}:
            return record
    raise AssertionError(f"missing pair {a}-{b}")


close(pair("T2", "T3")["rg"], 0.983, 0.002)
close(pair("T1", "NEURO")["rg"], 0.868, 0.002)
assert float(pair("T2", "T3")["ci_high"]) > 1.0

qc = rows("results/genomic/Scheme_C_univariate_ldsc_qc.csv")
assert len(qc) == 7, f"expected 7 univariate QC rows, got {len(qc)}"
z_field = next((name for name in qc[0] if name.lower().replace("_", "") in {"h2z", "z"}), None)
if z_field is None:
    raise AssertionError("could not identify the heritability-Z field")
assert min(float(record[z_field]) for record in qc) >= 10.0

contrasts = rows("results/genomic/Scheme_C_formal_tier_rg_contrasts.csv")
assert len(contrasts) == 12, f"expected 12 planned contrasts, got {len(contrasts)}"
adhd_t2_t3 = next(
    record for record in contrasts
    if record["comparator"] == "ADHD" and record["tier_a"] == "T2" and record["tier_b"] == "T3"
)
close(adhd_t2_t3["difference"], -0.065, 0.002)
assert float(adhd_t2_t3["p_bonferroni_12"]) < 0.002
asthma = [record for record in contrasts if record["comparator"] == "ASTHMA"]
assert len(asthma) == 3
assert all(float(record["p_bonferroni_12"]) >= 0.05 for record in asthma)

check = {r["analysis_label"]: r for r in rows("source_version_check/Scheme_C_source_version_reverse_check.csv")}
adhd = check["ADHD"]
assert int(adhd["hm3_snps_found"]) == 39 and int(adhd["alleles_aligned_n"]) == 39
assert round(float(adhd["beta_within_1e-4_fraction"]) * 39) == 1, "ADHD: 1 of 39 effects within 1e-4"
t3 = check["T3"]
assert int(t3["hm3_snps_found"]) == 55 and round(float(t3["beta_within_1e-4_fraction"]) * 55) == 53
retained = rows("source_version_check/retained_instrument_effects.csv")
assert set(retained[0]) == {"label", "SNP", "source_id", "effect_allele", "other_allele", "beta", "p"}
counts = {}
for r in retained:
    counts[r["label"]] = counts.get(r["label"], 0) + 1
for label in ("T1", "T2", "T3", "ADHD", "NEURO", "ASTHMA"):
    assert counts[label] == int(check[label]["retained_unique_snps"]), label

for relative in (
    "results/figures/Scheme_C_genetic_correlation_heatmap.png",
    "results/figures/Scheme_C_tier_profile_forest.png",
    "code/ldsc_input_manifest.tsv",
    "provenance/SOFTWARE_AND_INPUT_PROVENANCE.md",
):
    path = SCHEME / relative
    assert path.is_file() and path.stat().st_size > 0, f"missing or empty: {path}"

print("Scheme C validation passed: 7 traits, 21 correlations, 12 contrasts, and the source-version check.")
