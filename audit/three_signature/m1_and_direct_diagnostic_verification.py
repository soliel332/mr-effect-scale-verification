"""
Numerical verification of the audit classification and the direct scale diagnostic (no new analysis).

  M1  - analysis-sample case fraction p = ncase/(ncase+ncontrol) versus the case-to-reported-N
        eligibility ratio ncase/N; three-signature reclassification of the 41 audit records with
        reported N perturbed by 0.8-1.2 while p is held fixed.
  A4  - direct scale diagnostic for ukb-b-6991 (44 variants) and ieu-b-102 (50 variants),
        recomputed from pre-existing variant-level input files and compared with the frozen
        output data/repository/case_study_scale_diagnostic.csv.

usage: python m1_and_direct_diagnostic_verification.py REPO_DIR SUPPORTING_DIR OUT_DIR
"""
import os, sys
import numpy as np, pandas as pd

REPO, SUP, OUT = sys.argv[1:4]
NAMES = np.array(["linear probability", "standardised", "log odds"])


def geometry(sigma, p):
    L = 0.5 * np.log(p * (1 - p))
    d = np.abs(np.log(np.asarray(sigma))[..., None] - np.stack([L, 0 * L, -L], -1))
    o = np.argsort(d, -1)
    best = np.take_along_axis(d, o[..., :1], -1)[..., 0]
    sep = np.take_along_axis(d, o[..., 1:2], -1)[..., 0] - best
    return NAMES[o[..., 0]], best, sep, d


# ---------------------------------------------------------------- M1
rs = pd.read_csv(os.path.join(REPO, "results/repository_scale_reanalysis.csv"))
s98 = pd.read_csv(os.path.join(REPO, "data/repository/repository_sample_98_reconstructed.csv"))
cls = s98[s98.classifiable_in_preserved_output == True]
assert len(cls) == 41 and set(cls.id) == set(rs.id)
assert (cls.N == cls.ncase + cls.ncontrol).all()
odd = s98[s98.N != s98.ncase + s98.ncontrol]
assert len(odd) == 3 and not odd.classifiable_in_preserved_output.any()
m = rs.merge(cls[["id", "ncase", "ncontrol"]], on="id")
assert np.allclose(m.prev, m.ncase / (m.ncase + m.ncontrol), rtol=0, atol=1e-12)
frozen = rs.inferred_recomputed.str.replace("_", " ")
rows = []
for f in (0.8, 0.9, 1.0, 1.1, 1.2):
    lab, best, sep, _ = geometry(rs.sigma_obs.values * np.sqrt(f), rs.prev.values)   # sigma-hat ~ sqrt(N); p fixed
    rows.append(dict(N_multiplier=f, n_records=len(rs), agree_with_frozen_assignment=int((lab == frozen).sum()),
                     changed=int((lab != frozen).sum()), min_separation=round(float(sep.min()), 4),
                     max_distance_to_assigned=round(float(best.max()), 4)))
pd.DataFrame(rows).to_csv(os.path.join(OUT, "or1_tableS14_N_perturbation_p_fixed.csv"), index=False)

# ---------------------------------------------------------------- A4
fz = pd.read_csv(os.path.join(REPO, "data/repository/case_study_scale_diagnostic.csv")).set_index("id")
cat = pd.read_csv(os.path.join(REPO, "data/repository/eligible_binary_european_catalogue.csv")).set_index("id")
out = []
for acc, fn in [("ukb-b-6991", "ukb-b-6991_direct_diagnostic_44_variants.csv"),
                ("ieu-b-102", "ieu-b-102_direct_diagnostic_50_variants.csv")]:
    d = pd.read_csv(os.path.join(SUP, fn))
    d = d[(d["id.exposure"] == acc) & d["se.exposure"].notna() & d["eaf.exposure"].between(0, 1, inclusive="neither")].drop_duplicates("SNP")
    nc, nco = int(cat.loc[acc, "ncase"]), int(cat.loc[acc, "ncontrol"])
    N = nc + nco
    assert N == int(cat.loc[acc, "N"]) == int(fz.loc[acc, "N"])
    p = nc / N
    sig = float(np.median(d["se.exposure"] * np.sqrt(2 * N * d["eaf.exposure"] * (1 - d["eaf.exposure"]))))
    lab, best, sep, dist = geometry(np.array([sig]), np.array([p]))
    assert len(d) == int(fz.loc[acc, "n_usable"]) and abs(sig - fz.loc[acc, "implied_sd"]) < 1e-9
    assert lab[0].replace(" ", "_") == fz.loc[acc, "inferred_scale"]
    L = abs(0.5 * np.log(p * (1 - p)))
    out.append(dict(release=acc, N=N, cases=nc, controls=nco, p=p, usable_variants=len(d), sigma_hat=sig,
                    sig_linear=np.sqrt(p * (1 - p)), sig_standardised=1.0, sig_log_odds=1 / np.sqrt(p * (1 - p)),
                    dist_linear=dist[0, 0], dist_standardised=dist[0, 1], dist_log_odds=dist[0, 2],
                    assigned=lab[0], separation=sep[0], attainable_ceiling=L,
                    adequacy_alarm_tau_0_1633=bool(best[0] > 0.163276)))
pd.DataFrame(out).to_csv(os.path.join(OUT, "direct_diagnostic_summary_full_precision.csv"), index=False)
print("all assertions passed")
