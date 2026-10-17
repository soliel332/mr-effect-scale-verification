# Software and input provenance

## Software

- CBIIT LDSC commit: `6c673952cee74bd5c57aef1555a03b1c015399a0`
- MRC IEU European LD-score reference commit: `6e03f24973d69ba15ecc814e11e5d6ae896187c9`
- R 4.4.1 on aarch64-apple-darwin20
- GenomicSEM 0.0.5
- GenomicSEM used 200 jackknife blocks.

Complete R package details are preserved in `Scheme_C_GenomicSEM_sessionInfo.txt` and `sessionInfo.txt`.

## Raw-input SHA-256 values (files not redistributed)

The hashes below permit exact local verification without redistributing the files.

- `adhd_PGC2019.gz`: `1d7f3a359f599d320f9d2405f33e633283ab47e4d4fb013e095734c2bbcc632d`
- `asd_PGC2019.gz`: `fa30c2767d6f32a3482aa6bf3f223832f62a28a618a6b40a60475ffdde48086a`
- `asthma_UKB_EBI_h.tsv.gz`: `90e64bec1ec6e7fdca04397269282958655be208360a9253a7878ecc61ae3859`
- `mdd_PGC2019.txt`: `082a1602e405fa5e8b40917521f063444131db68e2f435e7aacfd80084852110`
- `neuroticism_EBI_h.tsv.gz`: `21951930d33a4cc5d3e5780f99a541dded842e3576f1da6fafd37aeae7d20031`
- `t1_mood_swings_EBI_h.tsv.gz`: `2115ded37ff23b38f35510a962bb24384705344c27ccf6bc5e635f4e389312b1`
- `t2_helpseeking_corrected_EBI_h.tsv.gz`: `97f31be7398ad55030052e1c989c9798a65b3526bd766f1cb522c406ec47ae28`

## Source-version boundary

The completed result is valid for the public source files recorded in `ldsc_input_manifest.tsv`. Tier 1, corrected Tier 2, neuroticism, and asthma matched retained instrument effects exactly or to reporting precision. For major depression, all 55 overlapping HapMap3 instrument SNPs had aligned alleles and 53 of 55 effects agreed within 1e-4. The public PGC ADHD file matched SNP and allele identities but was numerically different from the older OpenGWAS `ieu-a-1183` instrument snapshot. The ADHD results therefore describe the recorded public file, not an exact re-download of that historical snapshot (Online Resource 1, Section S9). The comparison is reproduced by `code/validate_ldsc_sources_against_frozen_instruments.py`; its output is `source_version_check/Scheme_C_source_version_reverse_check.csv`.
