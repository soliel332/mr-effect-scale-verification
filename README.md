# Effect-scale verification for summary-data Mendelian randomization

This repository contains the code, derived inputs, machine-readable results and provenance records
supporting the study *Effect-scale metadata discordance in binary-trait GWAS summary statistics evades conventional Mendelian randomization sensitivity diagnostics*. It reproduces the invariance checks, calibration analyses, simulations,
repository audit and worked depression–asthma case study reported in the manuscript and its Online
Resources 1 and 2.

Version 1.0 accompanies the manuscript submitted to the *European Journal of Epidemiology*.

`MANUSCRIPT_TO_REPOSITORY_MAP.csv` links each table, figure and result block of the main text and of
Online Resources 1 and 2 to its evidence files. `results/claim_evidence_manifest.csv` lists the
numerical claims, each with its evidence file and manuscript location; `tests/validate_claims.py`
recomputes every one of them.

## Interpretation corrections

Two historical records are kept unchanged, because one is an output of the frozen core script and the other
carries a recorded hash. Read them with the linked notes:

- `results/repository_scale_summary.csv`: the rows `median_classification_margin` and
  `minimum_classification_margin` are margins from the core script's two-signature comparison. For
  interpretation they are superseded by the three-signature geometry (Online Resource 1, Table S3b); every
  assignment is the same. See [`results/repository_scale_summary_NOTE.md`](results/repository_scale_summary_NOTE.md).
- `sim/2_gate_analysis/VALIDATION_GATE.md`: rows of its gate-verdict table are revised by errata E2–E4 later in
  the same document; where they differ, the errata apply. See
  [`sim/2_gate_analysis/READ_FIRST_ERRATA.md`](sim/2_gate_analysis/READ_FIRST_ERRATA.md).

## Repository structure

| Folder or file | Contents | Manuscript |
|---|---|---|
| `effect_scale_mr_analysis.R` | Core analysis: invariance checks, calibration, repository audit, corrected case-study MR | Methods and Results on invariance, calibration, repository audit and the case study |
| `data/` | Derived inputs: calibration pairs (`calibration/`), harmonised case-study instruments (`mr/`), repository catalogue and audit sample (`repository/`) | Tables 1–4 |
| `results/` | Machine-readable outputs of the core analysis and the claim–evidence manifest | Tables 2–4; Online Resource 1 |
| `sim/` | Pre-specified simulation (frozen protocol and code), its analysis, and the post hoc CheckSumStats comparison | Methods and Results on simulated operating characteristics; Online Resource 1, Tables S16–S20 |
| `audit/` | Three-signature implementation of the audit rule and its sensitivity to reported sample size | Methods on the conditional scale-attribution rule; Results on the repository audit; Figure 2; Online Resource 1, Tables S3b, S14 |
| `case_study/` | Direct scale diagnostic, MR-PRESSO in the pinned environment, leave-one-out, post hoc bridge comparison (`bridge/`) | Results on the effect-scale assessment, replacement-release reanalysis and corrected estimates; Online Resource 1, Tables S9, S21, S22 |
| `prior_studies/` | Post hoc provenance check of the cited depression–asthma MR studies and re-expression of one published estimate | Discussion |
| `scheme_c/` | Exploratory phenotype-architecture analysis (LDSC, Genomic SEM) | Online Resource 1, Section S9 |
| `online_resource_2/` | Export of the instrument-level workbook (Online Resource 2) | Online Resource 2 |
| `figures/` | Plotting code, submitted figure files and plotted values | Figures 1–3 |
| `tests/` | Checks that recompute the reported values; run by `.github/workflows/reproduce-core.yml` | — |
| `release_assets/` | Names and hashes of the three record-level simulation files attached to release `v1.0`, and a download script | — |

Section numbers are not used here because they may change in revision; `MANUSCRIPT_TO_REPOSITORY_MAP.csv` gives the mapping at the level of individual tables and figures.

## Quick start

Requirements: R 4.4.1 with the packages listed under *Software environments*, and Python 3.11 with
`pandas`, `numpy`, `scipy` and `openpyxl`. No OpenGWAS token is needed; all analyses read the derived
files in `data/`.

```sh
Rscript effect_scale_mr_analysis.R         # core analyses; writes results/
python3 tests/validate_claims.py           # recomputes every numerical claim in the manuscript
```

`Running the checks` below lists the remaining checks.

## What the repository reproduces

- **Invariance.** Numerical checks under uniform positive rescaling of the exposure coefficients
  (`results/invariance_*_unit_test.csv`). The domain claimed in the manuscript is exposure F-statistics;
  ratio z-statistics; the IVW estimate, equivariant with invariant z and P; Cochran's Q and I²; the
  MR-Egger intercept and its test; the MR-Egger slope, equivariant with invariant z and P; and the
  MR-PRESSO global statistic, its decision and the flagged outliers. The estimator unit test also
  reports the median and mode estimators; the claim is not extended to them.
- **Calibration.** Three paired-release calibrations of the scale signatures
  (`results/calibration_reanalysis.csv`; Main Table 2, Figure 1).
- **Pre-specified simulation** (`sim/1_frozen_prespecified/`, `sim/2_gate_analysis/`). Protocol,
  simulator, runner and analysis script were frozen, with the hashes in `FREEZE.txt`, before any
  simulated release was classified. In the core grid, 168,000 releases with correct sample size and
  case counts received no wrong scale attribution.
- **Post hoc comparison with CheckSumStats** (`sim/3_posthoc_comparator/`). **Not pre-specified.**
  CheckSumStats `R/make_plot_predlnor.R` at commit `ba670767d42a9f6041b65e4a155b64a3ac096434` was
  reimplemented and applied to releases regenerated from the frozen seeds (86 test cells, 200 replicates
  each; 0 label mismatches against the frozen records).
- **Repository audit.** The family-stratified sample of 98 records, 41 classifiable, 11 declaring
  standardised units that their statistics do not support (`results/`, `data/repository/`; Main Table 3).
  The audit rule as reported compares three signatures (`audit/three_signature/`). The retained audit
  output, `results/repository_scale_reanalysis.csv`, was produced by an earlier implementation that
  compared two signatures; the three-signature implementation reproduces all 41 assignments, and
  perturbing the reported sample size by a factor of 0.8 to 1.2 changes none of them.
- **Case study.** The direct scale diagnostic of the two case-study releases from their variant-level
  inputs (`case_study/direct_diagnostic/`); corrected univariable MR against UK Biobank and FinnGen
  asthma with five estimators, heterogeneity, MR-Egger intercepts and instrument strength (`results/`);
  MR-PRESSO for the six analyses in the pinned environment (`case_study/pinned_mrpresso/`);
  leave-one-out estimates (`case_study/leave_one_out/`); and the three distinct quantities of Figure 3
  — the superseded original mixed-scale estimate (OR 3.12), its arithmetic rescaling (OR 1.29) and the
  replacement-GWAS reanalysis (OR 1.07). These are separate quantities, not one corrected sequence.
- **Exploratory phenotype-architecture analysis** (`scheme_c/`; Online Resource 1, Section S9),
  including the source-version check of the public input files against the retained instrument-level
  effects (`scheme_c/source_version_check/`).
- **Online Resource 2.** `online_resource_2/export_or2.py` writes the instrument-level workbook from
  `data/mr/`; `online_resource_2/check_or2_docx.py` checks the Word rendering cell by cell.
- **Figures 1–3** (`figures/`): plotting code, submitted files, plotted values and export records.

## What is not included

- Raw genome-wide summary statistics and LD-score reference files. They are retrieved from their
  sources by the scripts in `scheme_c/code/`.
- Steiger directionality results. Under the quantitative-trait approximation used by TwoSampleMR they
  are not interpretable for these binary phenotypes, and the case and control counts needed to
  pre-calculate the correlations are not in the retained harmonised files.

## Running the checks

From the repository root:

```sh
Rscript effect_scale_mr_analysis.R              # core MR, calibration and audit outputs (pinned R 4.4.1)
Rscript tests/validate_outputs.R                # headline values of the core outputs
python3 tests/validate_claims.py                # every claim in results/claim_evidence_manifest.csv
python3 tests/validate_sync_layers.py           # simulation, comparator, audit, case-study and figure layers
python3 tests/validate_prior_studies.py         # prior-study provenance and re-expression (v1.1)
python3 tests/validate_bridge.py                # post hoc bridge comparison (v1.1)
python3 tests/validate_scheme_c.py              # phenotype-architecture outputs and source-version check
python3 tests/validate_repository_map.py        # manuscript map, source manifest and FILE_MANIFEST.csv
python3 tests/compare_regenerated_outputs.py    # after the R rerun: regenerated tables against the committed ones
```

The R analysis uses the derived CSV files only and does not need an OpenGWAS token. Stochastic
procedures use the seed `20260909`; the simulation uses `numpy.random.SeedSequence(20260924)`.

### Record-level simulation outputs

Three record-level files (about 50 MB in total) are attached to release `v1.0` of this repository
rather than stored as Git objects. Their names, sizes and SHA-256 hashes are in
`release_assets/RELEASE_ASSETS.csv`.

```sh
python3 release_assets/fetch_release_assets.py --base-url https://github.com/soliel332/mr-effect-scale-verification/releases/download/v1.0    # or --from-dir <folder>
python3 tests/validate_sync_layers.py      # now also regenerates the comparator tables from the records
```

Without the assets, every other check runs, and Online Resource 1 Tables S16–S20 are rebuilt from the
distributed layer tables.

## Software environments

Each analysis records the environment in which it was run. The records are distributed as generated, except that the local path of the software environment is replaced by `<LOCAL_CONDA_ENV>`; `SESSION_RECORD_REDACTION.csv` lists the nine affected files with their SHA-256 before and after.

| Component | Environment | Record |
|---|---|---|
| Core MR, calibration, audit | R 4.4.1, TwoSampleMR 0.6.22, MRPRESSO 1.0, ggplot2 4.0.0 | `results/sessionInfo.txt` |
| MR-PRESSO, six analyses | R 4.4.1, MRPRESSO 1.0 (commit `3e3c92d`) | `case_study/pinned_mrpresso/mrpresso_pinned_sessionInfo.txt` |
| Leave-one-out | R 4.5.3, TwoSampleMR 0.7.9; reproduced independently in Python by `check_leave_one_out.py` | `case_study/leave_one_out/rerun_session_info.txt` |
| Simulation and comparator | Python 3.11.15, NumPy 2.4.6, SciPy 1.17.1, pandas 2.3.3 | `sim/1_frozen_prespecified/FREEZE.txt`, `sim/ENVIRONMENT.txt` |
| Figures 1 and 3, submitted export | R 4.5.3, ggplot2 4.0.3 | `figures/submitted/fig13_sessionInfo.txt` (original record) |
| Figure 2, submitted export | Matplotlib 3.11.0; Python version not recorded | header of `figures/submitted/fig2_margins.eps` |
| Phenotype architecture | R 4.4.1, GenomicSEM 0.0.5, LDSC commit `6c67395` | `scheme_c/provenance/` |

Figure checks made on 17 October 2026 in the current environment (`figures/verification/`): Figures 1
and 3 regenerated from the repository data with plotted values identical to the distributed plotted
values; Figure 2 regenerated under Python 3.11.15 and Matplotlib 3.11.0 with a pixel-identical PNG and
an EPS identical apart from its creation date. The submitted figure files were not replaced, and the
original export records were not rewritten.

### Cross-platform reproduction

The committed result tables were produced on macOS (arm64, R reference BLAS). The workflow in
`.github/workflows/reproduce-core.yml` reruns `effect_scale_mr_analysis.R` on Linux (ubuntu-24.04)
and then requires:

- every tracked file to be byte-identical, except the nine R result tables below and
  `results/sessionInfo.txt`, whose environment-specific difference is printed but not asserted;
- no new untracked file under `results/` or `data/`;
- for the nine tables (`MR_rerun_egger_intercept`, `MR_rerun_five_estimators`,
  `MR_rerun_heterogeneity`, `calibration_reanalysis`, and the `invariance_*_unit_test` and
  `repository_scale_reanalysis` tables), checked by `tests/compare_regenerated_outputs.py`: identical
  columns, rows, row order, text cells, missing-value positions and integer columns; floating-point
  cells within |Δ| ≤ 2×10⁻¹³ + 2×10⁻¹⁴ × max(|old|, |new|) and identical at four decimal places and four
  significant figures; residual columns that are zero in exact arithmetic (`absolute_*_difference`)
  at most 2×10⁻¹³ on both platforms;
- all 35 claims recomputed from the regenerated tables (`tests/validate_claims.py`).

These tables cannot be required to be byte-identical across platforms: in the first Linux run
(17 October 2026) 92 floating-point cells differed in the last one or two of 15 significant digits
(largest absolute difference 1.07×10⁻¹³, in Cochran's Q; largest difference-to-tolerance ratio 0.098),
with no change in structure, text, integers or missing values. The tolerance was set from that run.
A change of runner image (BLAS, compiler) could exceed it; the workflow then fails, and the tolerance
is not widened without examining the differences.

## Interpretation boundaries

- The 1,007-record list contains records with a standard-deviation unit label; it is not a list of
  1,007 confirmed errors.
- The 11 discordant records come from two UK Biobank identifier families and are not 11 independent
  pipelines. Discordance is a statistical inconsistency between declared and indicated scale, not a
  confirmed error.
- `ebi-a-GCST005902` is a targeted example outside the formal 98-record sample.
- The adequacy threshold τ = 0.1633 comes from the simulation; its use on the audit records is
  descriptive and post hoc.
- The MR-PRESSO outlier-corrected estimate is one sensitivity result among several.
- Genome-wide genetic correlation describes shared liability across the genome; it does not exclude
  locus-specific architecture within an MR instrument set. The phenotype-architecture analysis is
  exploratory.
- The paired calibration excludes a "seen a psychiatrist" pair, which regressed coefficients from UK
  Biobank field 2100 on field 2090 and was therefore not a same-field comparison.

## Data provenance and licences

The sampling rule and input checksums are in `data/repository/repository_sampling_provenance.md`.
Dataset identifiers and access pages are in `DATASET_SOURCE_MANIFEST.csv`; reuse conditions are in
`DATA_PROVENANCE_AND_REUSE.md`. Licences differ by file category and are listed in `DATA_LICENSES.md`:
code is under the MIT licence (`LICENSE`); derived SNP-level subsets of third-party GWAS are not
relicensed. No credential is included. Local software-environment paths in the environment records
have been replaced by `<LOCAL_CONDA_ENV>`; the original and redacted hashes are recorded in
`SESSION_RECORD_REDACTION.csv`. `FILE_MANIFEST.csv` gives the size and SHA-256 of every file.

## Citation and contact

Please cite the article and this repository; `CITATION.cff` gives the details. Questions can be raised
as a GitHub issue or sent to the corresponding authors of the article.
