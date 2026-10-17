#!/usr/bin/env python3
"""Restore the record-level simulation outputs, which are distributed as archive assets rather than
Git objects, and verify them against release_assets/RELEASE_ASSETS.csv.

usage:
  python release_assets/fetch_release_assets.py --from-dir DIR   # assets already downloaded to DIR
  python release_assets/fetch_release_assets.py --base-url URL   # download URL/<asset> for each asset
  python release_assets/fetch_release_assets.py --verify-only    # check assets already restored

Each asset is accepted only if its size and SHA-256 equal the recorded values. The archive location
(release or DOI) is stated in README.md under "Record-level simulation outputs".
"""
import argparse
import csv
import hashlib
import os
import shutil
import sys
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TABLE = os.path.join(ROOT, "release_assets", "RELEASE_ASSETS.csv")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--from-dir")
    g.add_argument("--base-url")
    g.add_argument("--verify-only", action="store_true")
    a = ap.parse_args()
    rows = list(csv.DictReader(open(TABLE, newline="")))
    bad = 0
    for r in rows:
        dst = os.path.join(ROOT, r["restore_to"])
        if a.from_dir:
            shutil.copyfile(os.path.join(a.from_dir, r["asset"]), dst)
        elif a.base_url:
            urllib.request.urlretrieve(a.base_url.rstrip("/") + "/" + r["asset"], dst)
        ok = os.path.isfile(dst) and os.path.getsize(dst) == int(r["bytes"]) and sha256(dst) == r["sha256"]
        print(("OK      " if ok else "FAILED  ") + r["restore_to"])
        bad += not ok
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
