# Phenotype-architecture analysis

This directory accompanies the exploratory analysis in Online Resource 1, Section S9, of the genome-wide genetic architecture of three depression-related exposure definitions. It is a conceptual complement to effect-scale verification, not a newly developed MR estimator and not evidence that genome-wide correlation can exclude instrument-local architecture.

## Contents

- `code/`: source manifest, download and preparation scripts, LDSC/GenomicSEM analysis, figure generation, and the source-version check.
- `results/genomic/`: univariate LDSC QC, the complete seven-trait/21-pair genetic-correlation results, covariance and intercept matrices, and 12 planned tier contrasts.
- `source_version_check/`: the retained instrument-level effects used by the source-version check (SNP, alleles, beta and P value only) and the check output reported in Online Resource 1, Section S9.
- `results/figures/`: publication-source PNG and PDF figures.
- `provenance/`: software versions, input hashes and completed-run logs.

## Reproduction boundary

The full genome-wide workflow is described in `code/README.md` and starts with `code/download_scheme_c_inputs.sh`. Raw or processed genome-wide summary statistics and LD-score reference files are not committed. The scripts retrieve the permitted public source files and fail closed when required resources are absent.

The completed analysis used CBIIT LDSC commit `6c673952cee74bd5c57aef1555a03b1c015399a0`, the MRC IEU European LD-score reference commit `6e03f24973d69ba15ecc814e11e5d6ae896187c9`, R 4.4.1 and GenomicSEM 0.0.5. Exact details and input hashes are in `provenance/SOFTWARE_AND_INPUT_PROVENANCE.md`.

GitHub Actions validates the structure and headline values of the retained results. It does not redownload multi-gigabyte GWAS inputs on every commit.
