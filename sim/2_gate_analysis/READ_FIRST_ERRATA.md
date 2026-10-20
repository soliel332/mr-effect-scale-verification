# Read with the errata

`VALIDATION_GATE.md` is kept as written at the time of the gate analysis, and its hash is recorded in
`sim/LAYER_MANIFEST.csv`. Its errata section revises statements made earlier in the same document, including rows of
the gate-verdict table:

- **E2** restricts the statement, in the verdict table and §7, that substituting population prevalence for the
  case fraction (arm A4b) "misattributed 100%" of releases: the standardised pipeline (P2) was assigned correctly in
  every release, because the standardised signature does not depend on the case fraction.
- **E3** withdraws the "eleven coincidences" argument: the 11 discordant records come from two release families that
  may share one convention and are not independent tests.
- **E4** withdraws the verdict-row claim of robustness to "variant-level sample-size heterogeneity including the
  CheckSumStats slope-0.23 configuration": the A7 arm did not reproduce that mechanism, so robustness was shown only
  for heterogeneity as simulated.

Where the verdict table and the errata differ, the errata apply. The manuscript and Online Resource 1 report the
corrected statements.
