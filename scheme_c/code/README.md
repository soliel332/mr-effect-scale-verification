# Scheme C reproduction code

This folder is the public-safe code layer for the seven-trait phenotype-architecture extension. Raw genome-wide summary statistics are deliberately excluded.

## Inputs

`ldsc_input_manifest.tsv` records the seven traits, public source URLs, filenames, reported sample sizes, ancestry, genome build, and redistribution restrictions. `download_scheme_c_inputs.sh` retrieves those sources into a local `private_raw_inputs` folder. The download script is sequential and fails immediately if any source transfer fails.

## Configuration

The scripts use `SCHEME_C_ROOT` as the analysis root. If it is unset, they default to their own folder. External resources can be provided through:

- `HM3_SNP_LIST`: HapMap3 SNP-list path;
- `LDSC_SOFTWARE_DIR`: LDSC software directory;
- `LDSC_PYDEPS_DIR`: compatible Python dependency directory;
- `LD_SCORE_DIR`: European LD-score directory;
- `PYTHON_BIN`: Python executable;

## Run order

1. `download_scheme_c_inputs.sh`
2. `prepare_ldsc_inputs.py`
3. LDSC munging for each prepared trait file
4. `run_scheme_c_full_ldsc.sh`
5. `run_genomicsem_scheme_c.R`
6. `make_scheme_c_figures.py`

`validate_ldsc_sources_against_frozen_instruments.py` is the source-version check reported in Online Resource 1 (Section S9). It compares the prepared public files with the retained instrument-level effects in `../source_version_check/retained_instrument_effects.csv` and writes `../source_version_check/Scheme_C_source_version_reverse_check.csv`. It is not required to reproduce the genetic-correlation results.

## External software snapshots used in the completed run

- CBIIT LDSC commit: `6c673952cee74bd5c57aef1555a03b1c015399a0`
- MRC IEU European LD-score reference commit: `6e03f24973d69ba15ecc814e11e5d6ae896187c9`
- R 4.4.1
- GenomicSEM 0.0.5

Do not commit `private_raw_inputs`, processed genome-wide summary statistics, authentication tokens, or restricted retained objects.
