#!/usr/bin/env python3
"""Write the instrument-level tables of Online Resource 2 from the retained harmonised files.

usage: python online_resource_2/export_or2.py [OUT_XLSX]

Reads data/mr/*_harmonised_kept.csv (variants with mr_keep == TRUE), writes one sheet per
analysis plus an index sheet, and checks that the inverse-variance-weighted odds ratios
recomputed from the written sheets equal results/MR_rerun_five_estimators.csv.

Values are written unmodified (full precision). Sample sizes that a release does not report per
variant are filled with the analysis sample size from Table 1 of the main text; allele
frequencies absent from a release are written as "NR"; P values of exactly zero are written as
released. The Word rendering of Online Resource 2 rounds these values for display only.
"""
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "online_resource_2", "Online_Resource_2_instrument_level_data.xlsx")

# Analysis-sample sizes from Table 1 of the main text, used only where a release reports none.
N_FINNGEN_ASTHMA = "156,078 (20,629 cases; 135,449 controls)"
N_IEU_B_102 = ("not stated in release metadata; 170,756 cases and 329,443 controls "
               "(Table 1 of the main text)")

# table, sheet, exposure label, file stem, outcome label
ANALYSES = [
    ("S2-2", "T1 AffInstab - UKB asthma", "Tier 1: affective instability (mood swings)",
     "Tier_1_affective_instability__UK_Biobank_asthma", "Asthma, UK Biobank"),
    ("S2-3", "T1 AffInstab - FinnGen asthma", "Tier 1: affective instability (mood swings)",
     "Tier_1_affective_instability__FinnGen_asthma", "Asthma, FinnGen"),
    ("S2-4", "T2 HelpSeek - UKB asthma", "Tier 2: help-seeking distress, replacement release",
     "Tier_2_help_seeking_corrected__UK_Biobank_asthma", "Asthma, UK Biobank"),
    ("S2-5", "T2 HelpSeek - FinnGen asthma", "Tier 2: help-seeking distress, replacement release",
     "Tier_2_help_seeking_corrected__FinnGen_asthma", "Asthma, FinnGen"),
    ("S2-6", "T3 Depression - UKB asthma", "Tier 3: depression meta-analysis",
     "Tier_3_major_depression_meta_analysis__UK_Biobank_asthma", "Asthma, UK Biobank"),
    ("S2-7", "T3 Depression - FinnGen asthma", "Tier 3: depression meta-analysis",
     "Tier_3_major_depression_meta_analysis__FinnGen_asthma", "Asthma, FinnGen"),
    ("S2-8", "T2 original mixed-scale - UKB",
     "Tier 2: help-seeking distress, original release, superseded analysis",
     "Tier_2_superseded_mixed_scale__UK_Biobank_asthma", "Asthma, UK Biobank"),
]
COLUMNS = ["SNP", "Chromosome", "Position", "Effect allele", "Other allele", "EAF (exposure)",
           "EAF (outcome)", "Beta (exposure)", "SE (exposure)", "P (exposure)", "N (exposure)",
           "Beta (outcome, log odds)", "SE (outcome)", "P (outcome)", "N (outcome)", "F statistic",
           "Palindromic", "Outcome proxy used", "Outcome proxy SNP"]
README_ROWS = ['Online Resource 2 — Instrument-level data for the depression–asthma case study', None, 'Contents',
    'One table per Mendelian randomization analysis, listing every variant retained after harmonisation, with its exposure and outcome associations. Table S2-1 is an index of the analyses. Instrument strength and the estimates themselves are reported in Online Resource 1 and the main text and are not repeated here.', None, 'Source',
    'Values are the harmonised, retained variants from the analysis repository (data/mr/*_harmonised_kept.csv), unmodified. Variants were selected and harmonised as described in Section 2.7 of the main text. Summary statistics were obtained from IEU OpenGWAS under the identifiers listed in Table 1 of the main text.', None, 'Columns',
    'Effect allele / Other allele: alleles after harmonisation; exposure and outcome betas refer to the same effect allele.',
    'Beta, SE, P (exposure): exposure association as released. Beta (outcome): log odds of asthma.',
    'F statistic: F = (β/SE)², from the exposure association.',
    'Palindromic: A/T or C/G variant, strand inferred from allele frequencies (retained only where inference was possible). Allele frequencies are absent from the retrieved summary statistics of the REGENIE releases (ebi-a-GCST90013923, ebi-a-GCST90013959 and ebi-a-GCST90013938) and are shown as NR.',
    'Outcome proxy used / Outcome proxy SNP: FinnGen lookups in which the instrument was replaced by an LD proxy in the outcome data.', None, 'Scale of the exposure betas',
    'The six main analyses use exposure releases reporting log-odds effects. Table S2-8 is the original, superseded analysis of the help-seeking release ukb-b-6991, whose coefficients are on the linear-probability scale although its metadata declare standardised units. It is included because the main text reports its estimate (odds ratio 3.12; Section 3.7, Figure 3); it is not a valid analysis. Multiplying its exposure betas and their standard errors by λ = 4.425, or equivalently dividing the MR log-odds estimate and its standard error by λ, gives the arithmetic scale conversion reported in Section 3.7 (odds ratio 1.29, P unchanged).', None, 'Reproduction',
    'IVW with multiplicative random effects: β_IVW = Σ(b_Xj·b_Yj / se_Yj²) ÷ Σ(b_Xj² / se_Yj²); SE(β_IVW) = √[1 / Σ(b_Xj² / se_Yj²)] × √max(1, Q/(K − 1)), with Q = Σ(b_Yj − β_IVW·b_Xj)² / se_Yj², where b_Xj and b_Yj are the exposure and outcome betas of variant j, se_Yj the outcome standard error and K the number of variants. Recomputed from these sheets, the six main estimates equal those reported in the main text and in Online Resource 1 (Table S5).', None, 'Where a release does not report sample size per variant, the N columns give the analysis sample size from Table 1 of the main text.']
INDEX_NOTE = 'Instrument strength (mean and minimum F) and the inverse-variance-weighted estimates for the six main analyses are given in Online Resource 1 (Tables S8 and S5). The estimate from the original mixed-scale analysis (Table S2-8) and its arithmetic scale conversion are reported in Section 3.7 and Figure 3 of the main text.'
INDEX_COLUMNS = ["Table", "Sheet", "Exposure", "Exposure ID", "Outcome", "Outcome ID", "Variants",
                 "Palindromic", "Outcome proxies"]


def fill_n(values, release_id):
    if values.notna().all():
        return values.astype(int).astype(object)
    assert values.isna().all(), "partially missing sample sizes are not handled"
    fill = {"finn-b-J10_ASTHMA": N_FINNGEN_ASTHMA, "ieu-b-102": N_IEU_B_102}[release_id]
    return pd.Series([fill] * len(values), index=values.index, dtype=object)


def sheet(stem):
    d = pd.read_csv(os.path.join(ROOT, "data", "mr", stem + "_harmonised_kept.csv"))
    d = d[d["mr_keep"] == True].copy()  # noqa: E712
    d = d.sort_values(["chr", "pos"], kind="mergesort")
    eid, oid = d["id.exposure"].iloc[0], d["id.outcome"].iloc[0]
    if "proxy.outcome" in d:
        proxy = d["proxy.outcome"].map(lambda v: v is True or str(v).upper() == "TRUE")
        proxy_snp = d["proxy_snp.outcome"].astype(object)
    else:
        proxy = pd.Series(False, index=d.index)
        proxy_snp = pd.Series(None, index=d.index, dtype=object)
    out = pd.DataFrame({
        "SNP": d["SNP"],
        "Chromosome": d["chr"].astype(int),
        "Position": d["pos"].astype(int),
        "Effect allele": d["effect_allele.exposure"],
        "Other allele": d["other_allele.exposure"],
        "EAF (exposure)": d["eaf.exposure"].astype(object).where(d["eaf.exposure"].notna(), "NR"),
        "EAF (outcome)": d["eaf.outcome"].astype(object).where(d["eaf.outcome"].notna(), "NR"),
        "Beta (exposure)": d["beta.exposure"],
        "SE (exposure)": d["se.exposure"],
        "P (exposure)": d["pval.exposure"],
        "N (exposure)": fill_n(d["samplesize.exposure"], eid),
        "Beta (outcome, log odds)": d["beta.outcome"],
        "SE (outcome)": d["se.outcome"],
        "P (outcome)": d["pval.outcome"],
        "N (outcome)": fill_n(d["samplesize.outcome"], oid),
        "F statistic": d["F_stat"],  # as stored in the harmonised file, (beta/se)^2
        "Palindromic": np.where(d["palindromic"].astype(bool), "Yes", "No"),
        "Outcome proxy used": np.where(proxy, "Yes", "No"),
        "Outcome proxy SNP": np.where(proxy, proxy_snp, "\u2014"),
    })[COLUMNS].reset_index(drop=True)
    return out, eid, oid


def ivw_or(x):
    bx, by, se = x["Beta (exposure)"], x["Beta (outcome, log odds)"], x["SE (outcome)"]
    return float(np.exp(np.sum(bx * by / se ** 2) / np.sum(bx ** 2 / se ** 2)))


def main():
    sheets, index = {}, []
    for table, name, exposure, stem, outcome in ANALYSES:
        x, eid, oid = sheet(stem)
        sheets[name] = x
        index.append([table, name, exposure, eid, outcome, oid, len(x),
                      int((x["Palindromic"] == "Yes").sum()), int((x["Outcome proxy used"] == "Yes").sum())])
    ref = pd.read_csv(os.path.join(ROOT, "results", "MR_rerun_five_estimators.csv"))
    ref = ref[ref["method"].str.startswith("Inverse variance weighted")]
    for (table, name, _e, stem, _o), row in zip(ANALYSES, index):
        if "superseded" in stem:
            continue
        r = ref[(ref["id.exposure"] == row[3]) & (ref["id.outcome"] == row[5])]
        assert len(r) == 1 and int(r["nsnp"].iloc[0]) == row[6], (table, len(r))
        assert abs(ivw_or(sheets[name]) - float(r["OR"].iloc[0])) < 1e-9, table
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with pd.ExcelWriter(OUT, engine="openpyxl") as xw:
        pd.DataFrame({"x": README_ROWS}).to_excel(xw, sheet_name="README", index=False, header=False)
        idx = pd.DataFrame(index, columns=INDEX_COLUMNS)
        idx = pd.concat([idx, pd.DataFrame([[None] * 9, [INDEX_NOTE] + [None] * 8], columns=INDEX_COLUMNS)])
        idx.to_excel(xw, sheet_name="Index", index=False)
        for name, x in sheets.items():
            x.to_excel(xw, sheet_name=name, index=False)
    print(f"wrote {OUT}: {len(sheets)} analyses; six IVW odds ratios equal results/MR_rerun_five_estimators.csv")


if __name__ == "__main__":
    main()
