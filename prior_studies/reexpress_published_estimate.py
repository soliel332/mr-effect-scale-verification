#!/usr/bin/env python3
"""Arithmetic re-expression of a published odds ratio whose exposure release was assessed in this study as
statistically consistent with linear-probability effects.

usage: python prior_studies/reexpress_published_estimate.py

Reads prior_studies/prior_study_provenance.csv and the frozen eligible-record catalogue, computes the
case-fraction factor lambda = 1 / (p (1 - p)) of the exposure release from its case and control counts,
and divides the published log odds ratio and its confidence limits by lambda. Because the coefficient and
its standard error are divided by the same factor, the Wald z statistic and P value are unchanged; the
script checks this by recovering the published P value from the published confidence interval.
The re-expression is conditional on that assessment. It is a re-expression of a published estimate on a
different scale, not a reanalysis of the underlying data.
"""
import csv, math, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
Z975 = 1.959963984540054


def norm_sf(x):
    return 0.5 * math.erfc(x / math.sqrt(2))


def main():
    cat = {r["id"]: r for r in csv.DictReader(open(os.path.join(ROOT, "data", "repository",
                                                                 "eligible_binary_european_catalogue.csv")))}
    rows = [r for r in csv.DictReader(open(os.path.join(ROOT, "prior_studies", "prior_study_provenance.csv")))
            if r["published_or"]]
    out = []
    for r in rows:
        c = cat[r["exposure_accession"]]
        p = float(c["ncase"]) / (float(c["ncase"]) + float(c["ncontrol"]))
        lam = 1.0 / (p * (1.0 - p))
        orr, lo, hi = float(r["published_or"]), float(r["published_ci_low"]), float(r["published_ci_high"])
        se = (math.log(hi) - math.log(lo)) / (2 * Z975)
        p_from_ci = 2 * norm_sf(abs(math.log(orr)) / se)
        new = [math.exp(math.log(v) / lam) for v in (orr, lo, hi)]
        out.append(dict(study=r["study"],
                        lambda_full=repr(lam), reexpressed_or_full=repr(new[0]), reexpressed_ci_low_full=repr(new[1]),
                        reexpressed_ci_high_full=repr(new[2]), exposure_accession=r["exposure_accession"], case_fraction=f"{p:.6f}",
                        lambda_=f"{lam:.6f}", published_or=r["published_or"], published_ci=f"{lo}-{hi}",
                        published_p=r["published_p"], p_recovered_from_ci=f"{p_from_ci:.4f}",
                        reexpressed_or=f"{new[0]:.2f}", reexpressed_ci=f"{new[1]:.2f}-{new[2]:.2f}",
                        z_before=f"{math.log(orr) / se:.4f}", z_after=f"{(math.log(orr) / lam) / (se / lam):.4f}"))
    path = os.path.join(ROOT, "prior_studies", "reexpressed_estimates.csv")
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
    for o in out:
        print(f"{o['study']}: OR {o['published_or']} ({o['published_ci']}) -> {o['reexpressed_or']} ({o['reexpressed_ci']}); "
              f"lambda {float(o['lambda_']):.4f}; z {o['z_before']} -> {o['z_after']}; P from CI {o['p_recovered_from_ci']}")
    return out


if __name__ == "__main__":
    main()
