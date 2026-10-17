#!/usr/bin/env python3
"""Collect Scheme C LDSC logs into auditable machine-readable tables."""

from __future__ import annotations

import csv
import math
import os
import re
from pathlib import Path


ROOT = Path(os.environ.get("SCHEME_C_ROOT", Path(__file__).resolve().parent))
RESULTS = ROOT / "ldsc_results"
OUT = ROOT / "results"


def p_from_z(z: float) -> float:
    return math.erfc(abs(z) / math.sqrt(2.0))


def label(path_text: str) -> str:
    return Path(path_text).name.replace(".sumstats.gz", "")


def parse_rg() -> list[dict[str, object]]:
    rows: dict[tuple[str, str], dict[str, object]] = {}
    for log in sorted(RESULTS.glob("rg_*_full.log")):
        in_summary = False
        for raw in log.read_text(encoding="utf-8", errors="replace").splitlines():
            if raw.startswith("Summary of Genetic Correlation Results"):
                in_summary = True
                continue
            if not in_summary or ".sumstats.gz" not in raw:
                continue
            fields = raw.split()
            if len(fields) != 12 or not fields[0].endswith(".sumstats.gz"):
                continue
            trait1, trait2 = label(fields[0]), label(fields[1])
            rg, se, z = map(float, fields[2:5])
            key = tuple(sorted((trait1, trait2)))
            rows[key] = {
                "trait_1": trait1,
                "trait_2": trait2,
                "rg": rg,
                "se": se,
                "ci_low": rg - 1.959963984540054 * se,
                "ci_high": rg + 1.959963984540054 * se,
                "z": z,
                "p_two_sided": p_from_z(z),
                "h2_trait_2": float(fields[6]),
                "h2_trait_2_se": float(fields[7]),
                "h2_intercept_trait_2": float(fields[8]),
                "h2_intercept_trait_2_se": float(fields[9]),
                "gencov_intercept": float(fields[10]),
                "gencov_intercept_se": float(fields[11]),
                "source_log": log.name,
            }
    return [rows[key] for key in sorted(rows)]


def parse_h2() -> list[dict[str, object]]:
    patterns = {
        "h2": re.compile(r"Total Observed scale h2: (\S+) \(([^)]+)\)"),
        "mean_chi2": re.compile(r"Mean Chi\^2: (\S+)"),
        "intercept": re.compile(r"^Intercept: (\S+) \(([^)]+)\)$", re.MULTILINE),
        "ratio": re.compile(r"^Ratio: (\S+) \(([^)]+)\)$", re.MULTILINE),
        "n_snps": re.compile(r"After merging with regression SNP LD, ([0-9]+) SNPs remain"),
    }
    rows = []
    for log in sorted(RESULTS.glob("h2_*.log")):
        text = log.read_text(encoding="utf-8", errors="replace")
        h2 = patterns["h2"].search(text)
        mean_chi2 = patterns["mean_chi2"].search(text)
        intercept = patterns["intercept"].search(text)
        ratio = patterns["ratio"].search(text)
        n_snps = patterns["n_snps"].search(text)
        rows.append({
            "trait": log.stem.replace("h2_", ""),
            "h2_observed": float(h2.group(1)) if h2 else None,
            "h2_se": float(h2.group(2)) if h2 else None,
            "h2_z": (float(h2.group(1)) / float(h2.group(2))) if h2 else None,
            "mean_chi2": float(mean_chi2.group(1)) if mean_chi2 else None,
            "intercept": float(intercept.group(1)) if intercept else None,
            "intercept_se": float(intercept.group(2)) if intercept else None,
            "ratio": float(ratio.group(1)) if ratio else None,
            "ratio_se": float(ratio.group(2)) if ratio else None,
            "n_snps": int(n_snps.group(1)) if n_snps else None,
            "source_log": log.name,
        })
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise RuntimeError(f"No rows parsed for {path.name}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rg_rows = parse_rg()
    h2_rows = parse_h2()
    write_csv(OUT / "Scheme_C_pairwise_genetic_correlations.csv", rg_rows)
    write_csv(OUT / "Scheme_C_univariate_ldsc_qc.csv", h2_rows)
    print(f"Parsed {len(rg_rows)} unique genetic correlations and {len(h2_rows)} h2 results")


if __name__ == "__main__":
    main()
