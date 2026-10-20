#!/usr/bin/env python3
"""Check the re-expressed prior-study estimate reported in the Discussion."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "prior_studies"))
import reexpress_published_estimate as r

out = {o["exposure_accession"]: o for o in r.main()}
y = out["ukb-b-6991"]
assert y["lambda_"].startswith("4.4250"), y["lambda_"]
assert y["reexpressed_or"] == "1.42" and y["reexpressed_ci"] == "1.12-1.81", y
assert y["z_before"] == y["z_after"], (y["z_before"], y["z_after"])
assert abs(float(y["p_recovered_from_ci"]) - float(y["published_p"])) < 0.0005, y
print("prior-study re-expression check passed")
