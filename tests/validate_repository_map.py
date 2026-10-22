#!/usr/bin/env python3
"""Validate the submission-to-repository crosswalk and public source manifest."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


mapping = read_csv(ROOT / "MANUSCRIPT_TO_REPOSITORY_MAP.csv")
assert len(mapping) == 20, f"expected 20 manuscript-artifact mappings, got {len(mapping)}"

missing: list[str] = []
for record in mapping:
    for value in record["repository_evidence"].split(";"):
        relative = value.strip()
        if relative and not (ROOT / relative).is_file():
            missing.append(f"{record['manuscript_location']}: {relative}")
assert not missing, "missing mapped evidence files:\n" + "\n".join(missing)

sources = read_csv(ROOT / "DATASET_SOURCE_MANIFEST.csv")
assert len(sources) == 12, f"expected 12 dataset-source records, got {len(sources)}"
assert len({row["dataset_id"] for row in sources}) == 12, "dataset identifiers must be unique"
assert all(row["access_url"].startswith("https://") for row in sources), "all source records need HTTPS access URLs"
assert all(row["raw_genomewide_redistributed"] == "no" for row in sources), "raw genome-wide files must not be marked as redistributed"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


# Files that are generated locally or restored from the archive assets are not in the manifest.
assets = {row["restore_to"] for row in read_csv(ROOT / "release_assets" / "RELEASE_ASSETS.csv")}


def generated(rel: str) -> bool:
    return rel in assets or rel.startswith("results/figures/") or rel.startswith("online_resource_2/Online_Resource_2_")


manifest = read_csv(ROOT / "FILE_MANIFEST.csv")
manifest_paths = {row["relative_path"] for row in manifest}
actual_paths = {
    path.relative_to(ROOT).as_posix()
    for path in ROOT.rglob("*")
    if path.is_file()
    and path.name not in {".DS_Store", "FILE_MANIFEST.csv"}
    and ".git" not in path.relative_to(ROOT).parts
    and "__pycache__" not in path.relative_to(ROOT).parts
    and not generated(path.relative_to(ROOT).as_posix())
}
assert manifest_paths == actual_paths, (
    f"manifest mismatch; missing={sorted(actual_paths - manifest_paths)}, "
    f"unexpected={sorted(manifest_paths - actual_paths)}"
)
for row in manifest:
    path = ROOT / row["relative_path"]
    assert path.stat().st_size == int(row["bytes"]), f"size mismatch: {row['relative_path']}"
    assert sha256(path) == row["sha256"], f"checksum mismatch: {row['relative_path']}"

print(
    f"Repository map validation passed: {len(mapping)} manuscript mappings, "
    f"12 source records, and {len(manifest)} manifest entries."
)
