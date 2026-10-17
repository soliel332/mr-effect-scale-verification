#!/usr/bin/env python3
"""Source-version check: compare the public genome-wide files used for the LDSC analysis with the
instrument-level effects retained from the original analyses.

The retained effects are read from scheme_c/source_version_check/retained_instrument_effects.csv
(SNP, alleles, beta and P value per source release), extracted once from the retained historical
analysis objects; no other historical object is needed. The public files are prepared by
download_scheme_c_inputs.sh and prepare_ldsc_inputs.py into ldsc_prepared_inputs/.
"""

from __future__ import annotations

import csv
import os
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(os.environ.get("SCHEME_C_ROOT", Path(__file__).resolve().parent))

ID_TO_LABEL = {
    "ebi-a-GCST90013923": "T1",
    "ebi-a-GCST90013959": "T2",
    "ieu-b-102": "T3",
    "ieu-a-1183": "ADHD",
    "ebi-a-GCST90029028": "NEURO",
}


def source_rows(label: str) -> pd.DataFrame:
    path = ROOT / "ldsc_prepared_inputs" / f"{label}.munging_input.tsv.gz"
    return pd.read_csv(path, sep="\t", dtype={"SNP": str})


def compare(label: str, retained: pd.DataFrame) -> dict[str, object]:
    current = source_rows(label)
    merged = retained.merge(current, on="SNP", how="left")
    found = merged[merged["A1"].notna()].copy()
    same = (found["effect_allele"] == found["A1"]) & (found["other_allele"] == found["A2"])
    swapped = (found["effect_allele"] == found["A2"]) & (found["other_allele"] == found["A1"])
    aligned = same | swapped
    expected_beta = np.where(same, found["retained_beta"], -found["retained_beta"])
    beta_abs_diff = np.abs(found["BETA"] - expected_beta)
    p_abs_diff = np.abs(found["P"] - found["retained_p"])
    beta_close = beta_abs_diff[aligned] <= 1e-4
    p_close = p_abs_diff[aligned] <= 1e-4
    beta_close_fraction = float(beta_close.mean()) if aligned.any() else np.nan
    p_close_fraction = float(p_close.mean()) if aligned.any() else np.nan
    if len(found) and aligned.all() and beta_close_fraction >= 0.99 and p_close_fraction >= 0.99:
        assessment = "exact_or_rounding_equivalent"
    elif len(found) and aligned.all() and beta_close_fraction >= 0.95:
        assessment = "mostly_equivalent_with_local_discrepancy"
    elif len(found) and aligned.all():
        assessment = "alleles_match_but_numeric_version_differs"
    else:
        assessment = "review"
    return {
        "analysis_label": label,
        "retained_unique_snps": int(retained["SNP"].nunique()),
        "hm3_snps_found": int(len(found)),
        "alleles_aligned_n": int(aligned.sum()),
        "alleles_aligned_fraction": float(aligned.mean()) if len(found) else np.nan,
        "beta_max_abs_difference": float(beta_abs_diff[aligned].max()) if aligned.any() else np.nan,
        "beta_median_abs_difference": float(beta_abs_diff[aligned].median()) if aligned.any() else np.nan,
        "beta_within_1e-4_fraction": beta_close_fraction,
        "p_max_abs_difference": float(p_abs_diff[aligned].max()) if aligned.any() else np.nan,
        "p_within_1e-4_fraction": p_close_fraction,
        "assessment": assessment,
    }


RETAINED = Path(__file__).resolve().parents[1] / "source_version_check" / "retained_instrument_effects.csv"


def _retained(label: str) -> pd.DataFrame:
    data = pd.read_csv(RETAINED, dtype={"SNP": str})
    subset = data[data["label"] == label].rename(columns={"beta": "retained_beta", "p": "retained_p"})
    return subset[["SNP", "effect_allele", "other_allele", "retained_beta", "retained_p"]].reset_index(drop=True)


def retained_exposures() -> dict[str, pd.DataFrame]:
    return {label: _retained(label) for label in ID_TO_LABEL.values()}


def retained_asthma() -> pd.DataFrame:
    return _retained("ASTHMA")


def main() -> None:
    rows = []
    for label, retained in retained_exposures().items():
        rows.append(compare(label, retained))
    rows.append(compare("ASTHMA", retained_asthma()))
    # ASD had no retained effect-size file among the historical analysis objects. Its public
    # PGC source is checked by checksum and format in the separate input audit.
    rows.append({
        "analysis_label": "ASD",
        "retained_unique_snps": 0,
        "hm3_snps_found": 0,
        "alleles_aligned_n": 0,
        "alleles_aligned_fraction": np.nan,
        "beta_max_abs_difference": np.nan,
        "beta_median_abs_difference": np.nan,
        "beta_within_1e-4_fraction": np.nan,
        "p_max_abs_difference": np.nan,
        "p_within_1e-4_fraction": np.nan,
        "assessment": "not_applicable_no_frozen_effect_file",
    })
    output = RETAINED.parent / "Scheme_C_source_version_reverse_check.csv"
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
