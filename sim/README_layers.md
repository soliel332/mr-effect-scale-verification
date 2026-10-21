# Simulation and comparator layers — reproducibility manifest

Three layers, kept physically separate. Nothing in layer 1 has been modified since it was frozen;
the six hashes recorded in `1_frozen_prespecified/FREEZE.txt` were re-verified against the packaged
copies on 25 September 2026 (all match).

| Layer | Status | Contents |
|---|---|---|
| `1_frozen_prespecified/` | **Protocol-frozen** (designed after the repository audit; pre-specified with respect to its own outcomes only; not a preregistration). Frozen 24 Sep 2026 15:54:05 UTC, root seed `SeedSequence(20260924)`, before any simulated release was classified. | protocol, simulator, runner, analysis script (hashed at 15:54:58 UTC, before results existed), pre-execution checks, run log, record-level output (425,600 releases; 1,196 chunks consolidated), τ |
| `2_gate_analysis/` | Output of the frozen analysis script plus integrity checks and the gate report. The gate report carries appended errata E1–E3 and an A7 erratum; the frozen protocol was not rewritten. `gate_audit_coverage_posthoc.csv` is labelled post hoc. | |
| `3_posthoc_comparator/` | **POST HOC — not pre-specified.** CheckSumStats implementation benchmark, commit `ba670767d42a9f6041b65e4a155b64a3ac096434`, reimplemented from `R/make_plot_predlnor.R`. | scripts, record-level outputs for both frames, execution logs, regeneration checks, derived tables, adjudication |

References to "Reviewer A" and "Reviewer B" in the frozen protocol are to two critical reviews of an
earlier draft obtained before submission, not to journal peer review. The protocol is distributed
unchanged because its hash is recorded in `FREEZE.txt`.

## Post hoc comparator — key facts

- **Test cells: 86** (the plan has 110 cells; the other 24 are tuning cells and were not used), first 200 test replicates each.
- Regeneration against the frozen records: declared-log-odds frame 123,200 records, declared-standardised frame 120,400; 0 label mismatches; maximum |Δσ̂| = 0.
- The three-variant minimum for computing a CheckSumStats slope is **our choice**; the package specifies none.
- 1000 Genomes reference-variant exclusion (rsid-based) was not reproduced.

## Regenerating every reported comparator number

Copy the three layer directories' files into one working directory, then:

    python posthoc_comparator_tables.py        # the three paired CSVs, factor quantiles, scalar summaries
    python build_or_sim_tables.py <out_dir>     # Online Resource Tables S16-S20

`3_posthoc_comparator/REGENERATION_TEST_stdout.txt` records a run in a clean temporary directory
on 25 September 2026. All five comparator CSVs and all five Online Resource tables were
byte-identical to the packaged copies.

Regenerating the record-level files themselves (`posthoc_checksumstats.py`, `posthoc_sd.py`) re-simulates
from the frozen seeds and needs the frozen modules `sim_scale.py` and `run_sim.py` on the path. Runtime was not recorded.

## Environment

See `ENVIRONMENT.txt`. `MANIFEST.csv` lists every file with SHA-256 and size. The three chunk
directories are represented by a directory digest rather than copied.
