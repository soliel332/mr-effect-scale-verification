#!/usr/bin/env python3
"""Independent recomputation of the leave-one-out IVW estimates from the retained harmonised files.

usage: python case_study/leave_one_out/check_leave_one_out.py

For each of the six corrected analyses, every variant is dropped in turn and the multiplicative
random-effects IVW estimate is recomputed (TwoSampleMR convention: residual standard error
floored at 1, normal reference distribution). The results must equal the distributed
*__leaveoneout.csv files, which were produced by phase2_expanded_mr.R with TwoSampleMR.
"""
import glob
import os

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
TOL = 1e-10


def ivw(bx, by, se):
    w = bx ** 2 / se ** 2
    b = np.sum(bx * by / se ** 2) / np.sum(w)
    k = len(bx)
    phi = np.sqrt(np.sum(((by - b * bx) / se) ** 2) / (k - 1)) if k > 1 else 1.0
    s = np.sqrt(1 / np.sum(w)) * max(1.0, phi)
    return b, s, 2 * stats.norm.sf(abs(b / s))


def main():
    files = sorted(glob.glob(os.path.join(HERE, "*__leaveoneout.csv")))
    assert len(files) == 6, len(files)
    worst = 0.0
    for f in files:
        model = os.path.basename(f).replace("__leaveoneout.csv", "")
        d = pd.read_csv(os.path.join(ROOT, "data", "mr", model + "_harmonised_kept.csv"))
        d = d[d["mr_keep"] == True]  # noqa: E712
        ref = pd.read_csv(f)
        ref = ref[ref["SNP"].isin(d["SNP"])]
        assert len(ref) == len(d), (model, len(ref), len(d))
        for row in ref.itertuples(index=False):
            k = d[d["SNP"] != row.SNP]
            b, s, p = ivw(k["beta.exposure"].values, k["beta.outcome"].values, k["se.outcome"].values)
            diff = max(abs(b - row.b), abs(s - row.se), abs(p - row.p))
            assert diff < TOL, (model, row.SNP, diff)
            worst = max(worst, diff)
        print(f"{model}: {len(ref)} leave-one-out estimates reproduced; minimum P = {ref['p'].min():.3f}")
    print(f"all six analyses reproduced; largest absolute difference {worst:.1e}")


if __name__ == "__main__":
    main()
