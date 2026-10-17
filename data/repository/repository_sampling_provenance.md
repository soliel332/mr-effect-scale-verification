# Provenance of the 98-dataset repository sample

The sampling procedure was recovered from the archived local analysis record that generated the original repository audit. The private archive is not part of this release.

The operative original rule was:

```python
samp = POOL.groupby("pref", group_keys=False).apply(
    lambda d: d.sample(min(len(d), 14), random_state=0)
)
```

`POOL` was the frozen set of 2,658 eligible European-ancestry binary datasets. It contained seven identifier families: `ebi-a`, `finn-b`, `ieu-a`, `ieu-b`, `ukb-a`, `ukb-b`, and `ukb-d`. The rule sampled 14 datasets from each family, giving 98 records. An adjacent `default_rng(0)` declaration in the recovered notebook was not used by the sampling expression; the operative reproducibility parameter was `random_state=0`.

The reconstruction was validated as follows:

- 98 unique sampled IDs were recovered.
- Each of the seven identifier families contributed exactly 14 IDs.
- All 41 IDs in `scale_audit_catalog.csv` are members of the reconstructed 98-ID sample.
- The classifiable counts by family are identical to the original run: ebi-a 8, finn-b 3, ieu-a 8, ieu-b 10, ukb-a 9, ukb-b 2, and ukb-d 1.

The original end-to-end API retrieval and failure log was not retained in a release-ready form. This repository therefore reproduces the sampling operation from the frozen eligible catalogue and verifies the 41 preserved classifiable IDs against that sample. It does not claim to reconstruct dataset-level attrition reasons or to compare native and reference allele-frequency inputs.

Input SHA-256 checksums:

- `eligible_binary_european_catalogue.csv`: `8bd9f4770a5a7fa4b05e33feee253d3baec59db00bbf0b5f800745e34077a922`
- `scale_audit_catalog.csv`: `93776b2dd484ee73e6dc5cebb2958c238d31e8c6f6c92e53ba7fa0252cd3d754`

Reconstructed output SHA-256:

- `repository_sample_98_reconstructed.csv`: `dcee1fe680fcb773ed4a7b65356d54f871bf3dd1e08186493bd32a49fb4ad441`
