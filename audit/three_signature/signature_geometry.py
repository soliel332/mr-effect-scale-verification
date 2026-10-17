#!/usr/bin/env python3
"""Three-signature geometry of the 41 classifiable audit records (Online Resource 1, Table S3b).

usage: python audit/three_signature/signature_geometry.py [OUT_CSV]

For each record, the observed per-allele scale sigma_obs is compared on the log scale with the
three signatures implied by its case fraction p: sqrt(p(1-p)) (linear probability), 1
(standardised) and 1/sqrt(p(1-p)) (log odds). The record is assigned to the nearest signature;
separation is the distance to the runner-up minus the distance to the assigned signature, and its
attainable ceiling is |0.5 log(p(1-p))|. Values are rounded to four decimals and rows are then ordered by relative separation and
identifier. Online Resource 1 Table S3b lists the same values; the 24 records with relative separation
1.0000 appear there in a different order. The assignments
must equal the frozen assignments in results/repository_scale_reanalysis.csv, which were made by
the earlier two-signature comparison; the standardised signature changes none of them.
"""
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "or1_tableS3b_signature_geometry.csv")
NAMES = np.array(["linear", "standardised", "logodds"])
FROZEN = {"linear_probability": "linear", "log_odds": "logodds"}


def main():
    rs = pd.read_csv(os.path.join(ROOT, "results", "repository_scale_reanalysis.csv"))
    p, sigma = rs["prev"].to_numpy(), rs["sigma_obs"].to_numpy()
    ceiling = np.abs(0.5 * np.log(p * (1 - p)))
    expected = np.stack([np.sqrt(p * (1 - p)), np.ones_like(p), 1 / np.sqrt(p * (1 - p))], 1)
    dist = np.abs(np.log(sigma)[:, None] - np.log(expected))
    order = np.argsort(dist, 1)
    rows = np.arange(len(p))
    assigned = NAMES[order[:, 0]]
    assert (assigned == rs["inferred_recomputed"].map(FROZEN).to_numpy()).all(), "assignment differs from the frozen audit"
    sep = dist[rows, order[:, 1]] - dist[rows, order[:, 0]]
    out = pd.DataFrame({
        "id": rs["id"], "case_fraction": p,
        "exp_linear": expected[:, 0], "exp_standardised": expected[:, 1], "exp_logodds": expected[:, 2],
        "d_linear": dist[:, 0], "d_standardised": dist[:, 1], "d_logodds": dist[:, 2],
        "runner_up": NAMES[order[:, 1]], "separation": sep, "max_separation": ceiling,
        "rel_separation": sep / ceiling, "boundary_low": np.exp(-ceiling / 2), "boundary_high": np.exp(ceiling / 2),
    }).round(4).sort_values(["rel_separation", "id"], kind="mergesort")
    out.to_csv(OUT, index=False)
    print(f"wrote {OUT}: {len(out)} records; assignments equal the frozen audit; runner-up "
          f"{', '.join(sorted(set(out['runner_up'])))}")


if __name__ == "__main__":
    main()
