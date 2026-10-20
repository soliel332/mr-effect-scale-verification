# Note on `repository_scale_summary.csv`

`repository_scale_summary.csv` is written by the core analysis script, `MR_scale_submission_analysis.R`, whose
audit step compares only the linear-probability and log-odds signatures. Its two margin rows,
`median_classification_margin` (2.08) and `minimum_classification_margin` (1.338), are margins from that
two-signature comparison. They are kept unchanged because the file is a regenerated output of the frozen core
script and is checked by CI against that script.

The manuscript and Online Resource 1 report the three-signature geometry instead (fit, separation and attainable
ceiling; Online Resource 1, Table S3b), produced by `audit/three_signature/signature_geometry.py`. Every assignment
is the same under both implementations. All other rows of the summary are unaffected.
