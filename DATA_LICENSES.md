# Licences by file category

This repository mixes code written for this study, outputs and documentation written for this study,
SNP-level subsets derived from third-party GWAS summary statistics, and one third-party source file.
They are not under one licence. Each category is listed below; where a file falls under more than one
description, the more restrictive entry applies.

| Category | Files | Terms |
|---|---|---|
| Code written for this study | `*.R`, `*.py`, `*.sh` and `.github/workflows/*`, except the third-party file below | MIT licence, `LICENSE` |
| Aggregate outputs, tables, figures and documentation written for this study | `results/`, `audit/three_signature/*.csv`, `case_study/pinned_mrpresso/*.csv`, `case_study/leave_one_out/*.csv`, `case_study/direct_diagnostic/direct_diagnostic_summary.csv`, `sim/**/*.csv`, `sim/**/*.md`, `sim/**/*.txt`, `figures/submitted/` image files and plotted values, `scheme_c/results/`, `scheme_c/source_version_check/Scheme_C_source_version_reverse_check.csv`, the Markdown files at the top level | Creative Commons Attribution 4.0 International (CC BY 4.0), <https://creativecommons.org/licenses/by/4.0/> |
| SNP-level subsets derived from third-party GWAS summary statistics | `data/mr/`, `data/calibration/`, `case_study/direct_diagnostic/*_variants.csv`, `scheme_c/source_version_check/retained_instrument_effects.csv` | **Not relicensed.** The terms of each original data provider apply (UK Biobank, PGC, FinnGen, NHGRI-EBI GWAS Catalog, IEU OpenGWAS; see `DATASET_SOURCE_MANIFEST.csv` and `DATA_PROVENANCE_AND_REUSE.md`). These files are supplied so that the reported analyses can be verified; further redistribution is subject to the original providers' terms. |
| Repository metadata derived from the OpenGWAS index | `data/repository/` | **Not relicensed.** OpenGWAS terms apply, <https://opengwas.io/policies/eula>. |
| Record-level simulation outputs (archive assets, not in Git) | listed in `release_assets/RELEASE_ASSETS.csv` | CC BY 4.0. These are simulated and contain no real data. |
| Third-party source file | `sim/3_posthoc_comparator/comparator_source/CheckSumStats_make_plot_predlnor_ba67076.R` | MIT licence, copyright 2021 Philip Haycock (CheckSumStats); notice in the same directory |
| Software environment records | `*sessionInfo*.txt`, `rerun_session_info.txt`, `ENVIRONMENT.txt`, `FREEZE.txt` | Distributed as generated, as provenance, with the local software-environment path replaced by `<LOCAL_CONDA_ENV>` (`SESSION_RECORD_REDACTION.csv`). |

Copyright holder for the MIT and CC BY 4.0 entries: the authors (see `CITATION.cff`).
