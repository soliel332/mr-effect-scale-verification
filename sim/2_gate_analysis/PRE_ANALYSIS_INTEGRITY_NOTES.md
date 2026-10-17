# Pre-analysis integrity notes (recorded while the primary run was executing; no simulation output read)

## Dependency versions at relaunch
Runtime after installing pyarrow: numpy 2.4.6 | scipy 1.17.1 | pandas 2.3.3 | pyarrow 25.0.0.
Frozen record: numpy 2.4.6 | scipy 1.17.1 | pandas 2.3.3. **Identical.** pyarrow is used only to
write output files and enters no computation.

## First launch
Failed at the first file write for want of pyarrow. `sim/out/` contained 0 files before the
relaunch (checked). It contributes no observations.

## Code hashes at relaunch
sim_scale.py 2b1bcb04…, run_sim.py b0b7a736…, analyse_sim.py a0e5d614… — identical to FREEZE.txt.

## Deployment rule of the repository audit (validation-gate point 2)
The frozen reproduction package does NOT document the audit's instrument-selection rule; the only
statement in the frozen layer is the manuscript's own sentence in Methods 2.2. The rule was
recovered from an archived candidate script (`MR_scale_full_rerun.R`, untrusted candidate evidence,
outside the frozen package; its folder is marked as containing a credential to be removed, and
only parameter lines were read):

- instruments: `TwoSampleMR::extract_instruments(p1 = 5e-8, clump = TRUE, r2 = 0.001, kb = 10000)`
- missing effect-allele frequencies filled from reference associations (`fill_missing_eaf`)
- `diagnose_scale(..., minimum_snps = 5)`: median implied SD, argmin log distance, three signatures
- **no adequacy or unknown class** — nearest-signature assignment only
- releases with no instruments recorded with `n_usable = 0`, not dropped

Agreement with the simulation: P < 5e-8 on each release's own test — yes; minimum of five — yes;
median and argmin — yes; no-instrument releases retained as abstentions — yes.
Differences: (i) the audit uses ALL clumped significant variants, the simulation the first K
passing; (ii) the audit clumps, the simulation has no LD so needs none; (iii) the audit fills
missing frequencies from a reference, which only arm A6 represents.

## Haycock / CheckSumStats benchmark (validation-gate point 3)
Saved outputs hold per-release summaries only (sigma_hat, distances, separation), not per-variant
beta, SE and frequency, so the benchmark cannot be computed from them. Seeds are deterministic, so
identical data sets can be regenerated and the benchmark applied in a separate post-hoc run that
does not touch the primary outputs. Haycock et al. do not give the binary-trait expected-log-OR
formula in the main text (a simulation method, their reference 40, detailed in supplementary
material), so the benchmark must implement CheckSumStats `predict_lnor_sh` from its source code,
not from the paper's description.

## Adequacy rule status (validation-gate point 1)
Pre-specified in protocol §3.6 as a SECONDARY analysis, with tau set on the tuning set only, and
implemented in analyse_sim.py before any result was inspected (hash recorded in FREEZE.txt). The
PRIMARY rule, which is also the rule the audit deployed, has no unknown class.
