#!/usr/bin/env python3
"""Create HapMap3-restricted, publication-auditable LDSC inputs.

Raw GWAS files stay in private_raw_inputs and are never intended for GitHub.
The output contains only the columns LDSC needs. Source-specific parsing and
filtering counts are written to a machine-readable audit table.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
from pathlib import Path


ANALYSIS_ROOT = Path(
    os.environ.get("SCHEME_C_ROOT", Path(__file__).resolve().parent)
)
RAW_DIR = ANALYSIS_ROOT / "private_raw_inputs"
PREP_DIR = ANALYSIS_ROOT / "ldsc_prepared_inputs"
HM3_FILE = Path(
    os.environ.get(
        "HM3_SNP_LIST",
        ANALYSIS_ROOT / "ldsc_reference_repo/GRCh37/data/hm3/w_hm3.snplist",
    )
)
MANIFEST_FILE = ANALYSIS_ROOT / "ldsc_input_manifest.tsv"


def open_text(file_path: Path):
    if file_path.suffix == ".gz":
        return gzip.open(file_path, "rt", encoding="utf-8", errors="replace")
    return file_path.open("rt", encoding="utf-8", errors="replace")


def load_hm3() -> set[str]:
    keep: set[str] = set()
    with HM3_FILE.open("rt", encoding="utf-8") as handle:
        next(handle)
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields and fields[0].startswith("rs"):
                keep.add(fields[0])
    return keep


def load_manifest() -> dict[str, dict[str, str]]:
    with MANIFEST_FILE.open("rt", encoding="utf-8") as handle:
        return {
            row["analysis_label"]: row
            for row in csv.DictReader(handle, delimiter="\t")
        }


def finite_float(value: str | None) -> float | None:
    if value is None:
        return None
    clean = value.strip()
    if clean in {"", "NA", "NaN", "nan", ".", "-"}:
        return None
    try:
        numeric = float(clean)
    except ValueError:
        return None
    return numeric if math.isfinite(numeric) else None


def standard_record(
    label: str, row: dict[str, str], reported_n: int
) -> tuple[str, str, str, float, float, int, str | None]:
    if label in {"T1", "T2", "NEURO", "ASTHMA"}:
        snp = row.get("hm_rsid", "")
        a1 = row.get("hm_effect_allele", "").upper()
        a2 = row.get("hm_other_allele", "").upper()
        beta = finite_float(row.get("hm_beta"))
        p_value = finite_float(row.get("p_value"))
        info = None
    elif label == "ADHD":
        snp = row.get("SNP", "")
        a1 = row.get("A1", "").upper()
        a2 = row.get("A2", "").upper()
        odds_ratio = finite_float(row.get("OR"))
        beta = math.log(odds_ratio) if odds_ratio and odds_ratio > 0 else None
        p_value = finite_float(row.get("P"))
        info = row.get("INFO")
    elif label == "ASD":
        snp = row.get("SNP", "")
        a1 = row.get("A1", "").upper()
        a2 = row.get("A2", "").upper()
        odds_ratio = finite_float(row.get("OR"))
        beta = math.log(odds_ratio) if odds_ratio and odds_ratio > 0 else None
        p_value = finite_float(row.get("P"))
        info = row.get("INFO")
    elif label == "T3":
        snp = row.get("MarkerName", "")
        a1 = row.get("A1", "").upper()
        a2 = row.get("A2", "").upper()
        beta = finite_float(row.get("LogOR"))
        p_value = finite_float(row.get("P"))
        info = None
    else:
        raise ValueError(f"Unsupported label: {label}")
    return snp, a1, a2, beta, p_value, reported_n, info


def prepare_one(label: str, metadata: dict[str, str], hm3: set[str]) -> dict[str, object]:
    raw_path = RAW_DIR / metadata["local_filename"]
    if not raw_path.exists():
        raise FileNotFoundError(raw_path)
    output_path = PREP_DIR / f"{label}.munging_input.tsv.gz"
    reported_n = int(metadata["reported_N"])
    counts = {
        "analysis_label": label,
        "accession": metadata["accession"],
        "raw_file": str(raw_path),
        "prepared_file": str(output_path),
        "rows_seen": 0,
        "rows_not_hm3": 0,
        "rows_invalid_alleles": 0,
        "rows_missing_statistic": 0,
        "rows_low_info": 0,
        "rows_duplicate_snp": 0,
        "rows_written": 0,
        "reported_N": reported_n,
    }
    seen: set[str] = set()
    with open_text(raw_path) as source, gzip.open(
        output_path, "wt", encoding="utf-8", newline=""
    ) as target:
        if label == "T3":
            fieldnames = source.readline().strip().split()
            rows = (dict(zip(fieldnames, line.strip().split())) for line in source)
        else:
            rows = csv.DictReader(source, delimiter="\t")
        writer = csv.writer(target, delimiter="\t", lineterminator="\n")
        writer.writerow(["SNP", "A1", "A2", "P", "BETA", "N"])
        for row in rows:
            counts["rows_seen"] += 1
            snp, a1, a2, beta, p_value, sample_n, info = standard_record(
                label, row, reported_n
            )
            if snp not in hm3:
                counts["rows_not_hm3"] += 1
                continue
            if a1 not in {"A", "C", "G", "T"} or a2 not in {"A", "C", "G", "T"}:
                counts["rows_invalid_alleles"] += 1
                continue
            if beta is None or p_value is None or p_value < 0 or p_value > 1:
                counts["rows_missing_statistic"] += 1
                continue
            info_numeric = finite_float(info)
            if info is not None and (info_numeric is None or info_numeric < 0.9):
                counts["rows_low_info"] += 1
                continue
            if snp in seen:
                counts["rows_duplicate_snp"] += 1
                continue
            seen.add(snp)
            writer.writerow([snp, a1, a2, max(p_value, 1e-300), beta, sample_n])
            counts["rows_written"] += 1
    return counts


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "labels",
        nargs="*",
        help="Subset of T1 T2 T3 ADHD ASD NEURO ASTHMA; default is all available",
    )
    args = parser.parse_args()
    PREP_DIR.mkdir(parents=True, exist_ok=True)
    hm3 = load_hm3()
    manifest = load_manifest()
    labels = args.labels or list(manifest)
    audits = []
    for label in labels:
        audit = prepare_one(label, manifest[label], hm3)
        audits.append(audit)
        print(json.dumps(audit, ensure_ascii=False, sort_keys=True), flush=True)
    audit_path = PREP_DIR / "preparation_audit.jsonl"
    existing = []
    if audit_path.exists():
        with audit_path.open("rt", encoding="utf-8") as handle:
            existing = [json.loads(line) for line in handle if line.strip()]
    merged = {row["analysis_label"]: row for row in existing}
    merged.update({row["analysis_label"]: row for row in audits})
    with audit_path.open("wt", encoding="utf-8") as handle:
        for label in manifest:
            if label in merged:
                handle.write(json.dumps(merged[label], ensure_ascii=False, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
