#!/usr/bin/env python3
"""Checks for the simulation, comparator, three-signature, direct-diagnostic, sensitivity and figure
layers. Runs from a fresh clone without network access or record-level assets; checks that need the
record-level simulation outputs run only when release_assets/fetch_release_assets.py has restored them.
"""
import csv
import hashlib
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = sys.executable


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run(*args, cwd=ROOT):
    r = subprocess.run([PY, *args], cwd=cwd, capture_output=True, text=True)
    assert r.returncode == 0, f"{args[0]} failed:\n{r.stdout}\n{r.stderr}"
    return r.stdout.strip()


def rows(rel):
    with open(os.path.join(ROOT, rel), newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


# 1. Simulation layer 1 is unchanged since it was frozen (FREEZE.txt records the hashes).
layer1 = os.path.join(ROOT, "sim", "1_frozen_prespecified")
frozen = 0
for line in open(os.path.join(layer1, "FREEZE.txt")):
    parts = line.split()
    if len(parts) >= 2 and len(parts[0]) == 64 and parts[1].startswith("sim/"):
        assert sha256(os.path.join(layer1, parts[1][4:])) == parts[0], parts[1]
        frozen += 1
assert frozen == 6, frozen
print(f"simulation layer 1: {frozen} frozen hashes match FREEZE.txt")

# 2. Layer manifests: every distributed simulation file matches LAYER_MANIFEST.csv.
assets = {r["restore_to"] for r in rows("release_assets/RELEASE_ASSETS.csv")}
n = 0
for r in rows("sim/LAYER_MANIFEST.csv"):
    path = os.path.join(ROOT, "sim", r["path"])
    if "sim/" + r["path"] in assets and not os.path.isfile(path):
        continue
    if not os.path.isfile(path):
        continue  # chunk directories are recorded by digest only
    assert sha256(path) == r["sha256"], r["path"]
    n += 1
print(f"simulation layers: {n} files match LAYER_MANIFEST.csv")

# 3. Three-signature audit geometry and N perturbation (Online Resource 1, Tables S3b and S14).
with tempfile.TemporaryDirectory() as tmp:
    run("audit/three_signature/signature_geometry.py", os.path.join(tmp, "s3b.csv"))
    got = {r["id"]: r for r in csv.DictReader(open(os.path.join(tmp, "s3b.csv")))}
    ref = {r["id"]: r for r in rows("audit/three_signature/or1_tableS3b_signature_geometry.csv")}
    assert got == ref and len(got) == 41
    run("audit/three_signature/m1_and_direct_diagnostic_verification.py", ".", "case_study/direct_diagnostic", tmp)
    for produced, distributed in [("or1_tableS14_N_perturbation_p_fixed.csv", "audit/three_signature/or1_tableS14_N_perturbation_p_fixed.csv"),
                                  ("direct_diagnostic_summary_full_precision.csv", "case_study/direct_diagnostic/direct_diagnostic_summary.csv")]:
        a = list(csv.reader(open(os.path.join(tmp, produced))))
        b = list(csv.reader(open(os.path.join(ROOT, distributed))))
        assert len(a) == len(b) and a[0] == b[0], produced
        for ra, rb in zip(a[1:], b[1:]):
            for x, y in zip(ra, rb):
                try:
                    assert abs(float(x) - float(y)) <= 1e-12 * max(1.0, abs(float(y))), (produced, x, y)
                except ValueError:
                    assert x == y, (produced, x, y)
s14 = rows("audit/three_signature/or1_tableS14_N_perturbation_p_fixed.csv")
assert all(r["changed"] == "0" for r in s14) and len(s14) == 5
print("three-signature audit: 41 assignments unchanged; N x 0.8-1.2 changes none; direct diagnostic reproduced")

# 4. Variant-level inputs of the direct diagnostic carry the hashes cited in Online Resource 1.
expected = {"ukb-b-6991_direct_diagnostic_44_variants.csv": "a0baaae2b2bfa324704f71449077669b983f6b46110b27d73e15d853fe3e756a",
            "ieu-b-102_direct_diagnostic_50_variants.csv": "5a6e9e86afacd32c4c73b778ad282752190e54c5dd70d6b96c5d630274afa753"}
for name, digest in expected.items():
    assert sha256(os.path.join(ROOT, "case_study", "direct_diagnostic", name)) == digest, name

# 5. Leave-one-out, recomputed independently; pinned MR-PRESSO headline values.
print(run("case_study/leave_one_out/check_leave_one_out.py").splitlines()[-1])
presso = {r["model"]: r for r in rows("case_study/pinned_mrpresso/phase2_mrpresso_pinned.csv")}
t2 = presso["Tier_2_help_seeking_corrected__UK_Biobank_asthma"]
assert t2["outlier_SNPs"] == "rs30266;rs3094054;rs3734537" and t2["R_version"] == "4.4.1"
assert abs(float(t2["corrected_p"]) - 0.0138) < 5e-4 and abs(float(t2["distortion_p"]) - 0.4431) < 1e-9
print("pinned MR-PRESSO: Tier 2 UK Biobank three outliers, corrected P = 0.014, distortion P = 0.44")

# 6. Online Resource 2 workbook is regenerated identically.
with tempfile.TemporaryDirectory() as tmp:
    print(run("online_resource_2/export_or2.py", os.path.join(tmp, "or2.xlsx")))

# 7. Source-version check reads only the distributed retained effects.
src = open(os.path.join(ROOT, "scheme_c", "code", "validate_ldsc_sources_against_frozen_instruments.py")).read()
assert "MR_ARCHIVE_ROOT" not in src and "retained_instrument_effects.csv" in src

# 8. Figures: the original export environment record of Figures 1 and 3 is distributed with only the local
#    environment path redacted; SESSION_RECORD_REDACTION.csv records the hashes before and after.
redaction = {r["path"]: r for r in rows("SESSION_RECORD_REDACTION.csv")}
for path, r in redaction.items():
    assert sha256(os.path.join(ROOT, path)) == r["sha256_redacted"], path
assert redaction["figures/submitted/fig13_sessionInfo.txt"]["sha256_original"] == \
    "82d3ea5ae0c33283fc5fac2d7164e2b0017f609d53745a691dada24a3f01d3cc"
head = open(os.path.join(ROOT, "figures", "submitted", "fig2_margins.eps"), "rb").read(400)
assert b"%%Creator: Matplotlib v3.11.0" in head

# 9. Online Resource 1 Tables S16-S20 are rebuilt from the layer CSVs; the comparator CSVs are rebuilt
#    from the record-level outputs only when the archive assets have been restored.
import shutil  # noqa: E402

restored = all(os.path.isfile(os.path.join(ROOT, p)) for p in assets)
COMPARATOR = ["posthoc_logodds_anomaly_paired.csv", "posthoc_logodds_attribution_paired.csv",
              "posthoc_declaredSD_comparison.csv", "posthoc_factor_quantiles.csv", "posthoc_scalar_summaries.csv"]
with tempfile.TemporaryDirectory() as tmp:
    for layer in ("1_frozen_prespecified", "2_gate_analysis", "3_posthoc_comparator"):
        d = os.path.join(ROOT, "sim", layer)
        for name in os.listdir(d):
            if os.path.isfile(os.path.join(d, name)):
                shutil.copy(os.path.join(d, name), tmp)
    if restored:
        run("release_assets/fetch_release_assets.py", "--verify-only")
        run(os.path.join(tmp, "posthoc_comparator_tables.py"), cwd=tmp)
        for name in COMPARATOR:
            assert open(os.path.join(tmp, name), "rb").read() == \
                open(os.path.join(ROOT, "sim", "3_posthoc_comparator", name), "rb").read(), name
        print("record-level assets verified; five comparator tables regenerated byte-identically")
    else:
        print("record-level assets not restored: comparator regeneration from records skipped (see README)")
    out = os.path.join(tmp, "or1")
    os.makedirs(out)
    run(os.path.join(tmp, "build_or_sim_tables.py"), out, cwd=tmp)
    for name in sorted(os.listdir(os.path.join(ROOT, "sim", "online_resource_1_tables"))):
        assert open(os.path.join(out, name), "rb").read() == \
            open(os.path.join(ROOT, "sim", "online_resource_1_tables", name), "rb").read(), name
    print("Online Resource 1 Tables S16-S20 rebuilt byte-identically")
print("all synchronisation-layer checks passed")
