# VALIDATION_GATE — three-signature scale rule, simulation of 24 September 2026

No manuscript file was changed. The primary run was not stopped, modified, restarted or tuned.
Every quantity below comes from the frozen run (`sim/out/`, 1,196 chunks) or from the frozen
analysis script, except where a section is explicitly labelled **post hoc**. Performance is
reported as an unconditional decomposition — correct, wrong by class, abstained by reason — over
all simulated releases, never as a single accuracy among classified releases.

---

## Gate verdict

| Question | Answer |
|---|---|
| Run integrity (point 4) | **Passed.** Details in §4. |
| Can the primary rule abstain when a release fits no signature? | **No.** The primary rule — which is also the rule the audit deployed — is forced nearest-signature classification with abstention only on data sufficiency. An adequacy abstention was pre-specified as a **secondary** analysis (protocol §3.6, τ set on the tuning set before any result was seen). |
| Under correct release-level metadata, does the rule attribute scale correctly? | **Yes, in every simulated condition.** 0 wrong labels in 168,000 core releases, across seven pipelines including Firth and two saddlepoint-based stylisations, K = 5, 10, 30, and case fractions 0.01–0.75. All losses are sufficiency abstentions. |
| Does the rule abstain rather than misattribute when the error is not a scale error? | **Depends on the error.** Robust — no misattribution — to allele-frequency error, allele flips, variant-level sample-size heterogeneity including the CheckSumStats slope-0.23 configuration, and imputation INFO 0.8. **Not robust to release-level misreporting of *N*:** a five-fold overstatement at case fractions 0.35 and 0.50 is misattributed in 100% of P1 and P2 releases, and the secondary adequacy rule does **not** catch it, because the result fits the wrong signature closely. Supplying population prevalence in place of case fraction is misattributed 100% by the primary rule and caught 100% by the adequacy rule. |
| Recommended status | **Conditional attribution rule**, valid given correct release-level *N* and analysis-sample case fraction, and within the simulated pipeline and instrument-count range. Not an unconditional "validated attribution tool". Details in §7. |

---

## 1. Attribution versus forced classification

### 1.1 Core grid, primary rule, unconditional, pooled over case fraction

| K | Pipeline | Correct | Wrong (any class) | Abstain: < 5 GWS variants |
|---|---|---|---|---|
| 5 | P1 / P2 / P7 | 0.987 | 0.000 | 0.013 |
| 5 | P3 / P4 | 0.981 | 0.000 | 0.019 |
| 5 | P5 / P6 | 0.979 | 0.000 | 0.021 |
| 10 | P1 / P2 / P7 | 0.999 | 0.000 | 0.001 |
| 10 | P3 / P4 | 0.998 | 0.000 | 0.002 |
| 10 | P5 / P6 | 0.996 | 0.000 | 0.003–0.004 |
| 30 | P1 / P2 / P7 | 0.937 | 0.000 | 0.063 |
| 30 | P3 | 0.918 | 0.000 | 0.082 |
| 30 | P4 | 0.924 | 0.000 | 0.076 |
| 30 | P5 / P6 | 0.893 | 0.000 | 0.107 |

**Exact count of wrong labels in the core grid: 0 of 168,000 releases.** Full per-cell table:
`gate_unconditional_decomposition.csv`.

The lower "correct" at K = 30 is almost entirely one cell: case fraction 0.01, where 82.8% of P5/P6
releases (and fewer for other pipelines) abstain for want of five genome-wide-significant
variants. This is a property of the simulated architecture, not of the rule: 90 candidate variants
each calibrated to non-centrality 6 at 1% case fraction imply per-variant odds ratios near 2 and a
very large total genetic effect, which dilutes each variant's marginal association. It is reported,
not corrected.

### 1.2 Absolute fit and relative separation

Per-cell median and 95th percentile of the distance to the assigned signature, and median and 5th
percentile of relative separation, are in `gate_fit_selection.csv`. The adequacy threshold set on
the tuning set is **τ = 0.1633**.

### 1.3 Secondary (adequacy) rule — cost under correct metadata

The adequacy rule abstains on correctly scaled releases at these rates (core grid, pooled over *p*):

| K | P1 | P2 | P3 | P4 | P5 | P6 | P7 |
|---|---|---|---|---|---|---|---|
| 5 | 0.000 | 0.000 | 0.032 | 0.035 | 0.013 | 0.000 | 0.000 |
| 10 | 0.000 | 0.000 | 0.031 | 0.034 | 0.007 | 0.000 | 0.000 |
| 30 | 0.000 | 0.000 | 0.035 | 0.042 | 0.005 | 0.000 | 0.000 |

The cost falls on log-odds pipelines because the log-odds signature is an approximation, and τ was
set as a 99th percentile pooled over P1–P3.

### 1.4 Non-scale errors — the decisive quantity

All adversarial arms were run at **K = 10 only**. Nothing below may be extrapolated to K = 30 or to
the 100+ instruments of some audit records.

| Arm | Level | Max wrong, primary | Max wrong, adequacy | Max adequacy abstention | Reading |
|---|---|---|---|---|---|
| A4 | case fraction supplied | 0.000 | 0.000 | 0.000 | correct |
| A4 | **prevalence supplied** | **1.000** | 0.000 | 1.000 | primary misattributes; adequacy abstains |
| A5 | N × 0.5 | 0.740 | 0.000 | 1.000 | primary misattributes at *p* = 0.5; adequacy abstains |
| A5 | N × 2 | 0.240 | 0.000 | 1.000 | primary misattributes at *p* = 0.5; adequacy abstains |
| A5 | **N × 5** | **1.000** | **1.000** | 1.000 | **both rules misattribute at *p* = 0.35 and 0.50** |
| A6 | EAF noise SD 0.2 / 0.5 / 20% flipped | 0.000 | 0.000 | 0.032 / 0.054 / 0.020 | robust |
| A7 | CheckSumStats-matched heterogeneity | 0.000 | 0.000 | 0.008 | robust |
| A7 | severe (half of variants at N/4) | 0.000 | 0.000 | 0.042 | robust |
| A8 | INFO 0.8 | 0.000 | 0.000 | 0.006 | robust |
| A8 | INFO 0.5 | 0.022 | 0.000 | 0.092 | see note |

**Why N × 5 defeats both rules.** σ̂ is proportional to √*N*ᵣ. A five-fold overstatement multiplies
σ̂ by √5 = 2.24. At *p* = 0.5 the linear-probability signature 0.5 becomes 1.12 and the
standardised signature 1 becomes 2.24 — each within τ of the next signature up. Every P1 release was
assigned *standardised* and every P2 release *log odds*, with good fit. This is the "stable wrong
answer" failure mode: a systematic release-level error lands on another signature, and neither more
variants nor an adequacy threshold can detect it. At *p* = 0.10 the same error lands between
signatures and the adequacy rule abstains.

**INFO 0.5, reported both ways.** The frozen analysis script reported a misattribution of 0.314 for
P1 and P2 at *p* = 0.5, computed among classified releases. Unconditionally, 93.0% of those releases
abstained for want of five significant variants, and 11 of the 35 that were classified were wrong:
an unconditional misattribution of 0.022. Both figures are correct; the first alone would mislead.

**Direction relevant to the audit.** The audit's discordance type is "declared standardised,
indicated linear probability". Arm A5 shows that a release-level *N* error can move a release
between those two signatures. For the 11 audit discordances to be truly standardised releases with
misreported *N*, the reported *N* would have to **understate** the true *N* by factors of 4.17,
4.19, 4.27, 4.47, 4.74, 5.07, 5.69, 23.77, 34.56, 49.65 and 65.43, and in each record the factor
would have to equal 1/(*p*(1 − *p*)) to within the observed distances of 0.0015–0.022 log units,
across case fractions from 0.016 to 0.666. **This is an interpretive argument, not a simulation
result**, and it is the correct place to answer the alternative: the simulation shows the confusion
is possible; the audit data show it would require a case-fraction-specific coincidence in eleven
records. Effective sample size, the commonest alternative *N* convention, gives *N*ᵣ/*N* = 4*p*(1 − *p*),
which places a standardised release at twice the linear signature rather than on it, and so does not
produce this pattern.

---

## 2. Instrument-selection conditioning

### 2.1 Releases reaching five genome-wide-significant variants

Minimum over pipelines, core grid:

| K \ case fraction | 0.01 | 0.02 | 0.05 | 0.10 | 0.20 | 0.35 | 0.50 | 0.75 |
|---|---|---|---|---|---|---|---|---|
| 5 | 0.890 | 0.968 | 0.990 | 0.995 | 0.996 | 0.999 | 0.996 | 0.998 |
| 10 | 0.974 | 0.998 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| 30 | 0.172 | 0.974 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |

Every release that failed is retained as a sufficiency abstention (20,820 release records across
all sets); none was dropped. **All operating characteristics in this report are conditional on
instruments selected at *P* < 5 × 10⁻⁸ under each release's own test.**

### 2.2 Alignment with the audit's deployment rule — three levels

| Level | Status |
|---|---|
| **Screening rule** — *P* < 5 × 10⁻⁸ on the release's own test; at least five usable variants | **Aligned** |
| **Computation rule** — median implied SD; minimum log distance over three signatures; no unknown class | **Aligned** |
| **Deployment distribution** — number of instruments per release | **Not aligned.** The audit uses all clumped significant variants; the simulation uses the first K = 5, 10 or 30 passing. |

Audit instrument counts, all 41 classifiable records:
5, 5, 5, 6, 7, 8, 8, 9, 9, 10, 10, 10, 11, 12, 14, 14, 15, 16, 18, 18, 20, 23, 24, 25, 28, 34, 37,
39, 43, 43, 44, 49, 52, 53, 67, 81, 92, 118, 123, 138, 156.

- 5–9: 9 records; 10–29: 16 records; exactly 30: 0; **above 30: 16 records**.
- **25 of 41 fall within the simulated range (5–30); 16 lie above it.**
- Of the 11 discordant records, instrument counts are 5, 7, 8, 9, 10, 18, 18, 23, 24, 43, 43: nine
  within the simulated range, two above.

**No claim is made that performance improves with K beyond 30.** More variants make the median
more stable against random noise, but a systematic error — mixed scales within a release,
structured variant-level *N*, a biased reference frequency — would make the median converge more
stably on a wrong answer, and would raise the apparent separation while doing so. The core grid
shows K = 30 no less accurate than K = 10 under correct metadata; that is all it shows. Any
all-instrument analysis would be **post hoc** and is not part of this report.

---

## 3. Provenance of the audit's parameters

| Parameter | Source | Status |
|---|---|---|
| *P* < 5 × 10⁻⁸, clump r² = 0.001, 10,000 kb | archived candidate script `MR_scale_full_rerun.R` (outside the frozen package) | **found in candidate code; not proven to have generated the frozen outputs** |
| minimum five usable variants | same script (`minimum_snps = 5`) | found in candidate code; **compatible** with frozen outputs: minimum `nSNP` is 5, three records sit at exactly 5 |
| usable-variant count per classified record | frozen `repository_scale_reanalysis.csv` vs frozen 98-record ledger | **compatible**: 41 of 41 equal |
| no-instrument releases recorded as `n_usable = 0` | candidate code | **not matched**: the frozen ledger holds null, not 0, for all 57 unclassified records |
| instrument-selection rule | frozen package | **absent** — documented only in the manuscript's Methods sentence |

**Conclusion: compatible, not provenance-confirmed.** The frozen package contains no file linking
the candidate script to the frozen audit outputs, and re-extraction to test it is blocked by the
metadata service's expired certificate. This is an open provenance gap and is not closed by this
report.

---

## 4. Run integrity

| Check | Result |
|---|---|
| Chunks expected / present / missing / unexpected / leftover temp | 1,196 / 1,196 / 0 / 0 / 0 |
| Chunks with wrong row count | 0 |
| Release records | 425,600; duplicates on (set, cell, replicate, pipeline, supplied-p) = 0 |
| Replicate seeds / unique | 59,800 / 59,800 |
| Data sets: tuning / test | 4,800 / 55,000 |
| Deterministic regeneration | one replicate from each of 110 cells regenerated and compared: 784 records, 0 label / instrument-count / finiteness mismatches, maximum absolute difference in continuous fields 0.0 |
| Code hashes (protocol, simulator, runner, analysis) | all match `FREEZE.txt` |
| Dependency versions at run time | numpy 2.4.6, scipy 1.17.1, pandas 2.3.3 — identical to the frozen record; pyarrow 25.0.0 added, used only for file output |
| First, failed launch | wrote 0 files (checked before relaunch); contributes no observations |
| Logistic non-convergence | 175 data sets, all in three cells: A2 low-MAF *p* = 0.02 (170), A7 matched *p* = 0.02 (4), A7 severe *p* = 0.02 (1) — quasi-separation in rare-variant, rare-case data |
| Firth non-convergence / SPA fallback | 0 / 0 |
| Non-finite σ̂ outside sufficiency abstentions | 0 |

The analysis kernel could not read Parquet after pyarrow was installed into it mid-session (a pandas
extension-type registration error); all analysis was therefore run in fresh processes, whose
scripts are saved (`integrity_check.py`, `gate_tables*.py`).

---

## 5. Haycock / CheckSumStats benchmark

**Not run. It is post hoc and was not requested at this stage.**

- Saved outputs hold per-release summaries only, so the benchmark requires regenerating the data.
  That prerequisite is now verified (§4, deterministic regeneration: 0 mismatches).
- Haycock et al. do not give the binary-trait expected-effect formula in the main text; they cite
  a simulation method and refer to supplementary material. The comparator must therefore be
  described as the **CheckSumStats `predict_lnor_sh` implementation benchmark**, implemented from
  the package source at a recorded commit, not from the paper.
- When run, it should report separately (a) anomaly detection — slope outside 0.8–1.2; (b) correct
  scale attribution — which it does not attempt; and (c) false attribution under non-scale errors,
  on the same regenerated releases as §1.4.

---

## 6. Pre-specified consequences (protocol §9), as evaluated by the frozen script

| Rule | Triggered | Qualification added by this gate |
|---|---|---|
| R1 broad support | yes — core 24/24, A1 8/8, A2 4/4, A3 8/8, A4(a) 2/2 | conditional on correct release-level *N* and case fraction, K ≤ 30, and GWS-selected instruments |
| R2 bounded support | no failing cell | — |
| R3 saddlepoint pipelines fail | no | the P5/P6 stylisations are not REGENIE or SAIGE |
| R4 misattribution > 0.05 | yes — A5 at N × 0.5, 2, 5; A8 INFO 0.5 | INFO 0.5 is 0.022 unconditionally; N × 5 defeats the adequacy rule as well |
| R5 prevalence fails where case fraction passes | yes — 2/2 vs 0/2 | direct evidence that *p* must be the analysis-sample case fraction |

---

## 7. What the method can be called

**Not** a "validated attribution tool" without qualification. The evidence supports:

> Given correct release-level sample size and analysis-sample case fraction, the rule assigned the
> reporting scale correctly in every simulated release that had at least five genome-wide-
> significant instruments — 0 errors in 168,000 — across linear, standardised, logistic, Firth and
> two saddlepoint-based pipelines, with 5 to 30 instruments and case fractions from 0.01 to 0.75.
> It was unaffected by allele-frequency error, variant-level sample-size heterogeneity and moderate
> imputation uncertainty. It is not robust to release-level misreporting of sample size, which can
> move a release onto an adjacent signature; nor to substitution of population prevalence for case
> fraction, which a pre-specified adequacy threshold detects.

That is a **conditional attribution rule** — stronger than a compatibility screen, because under
its stated conditions it attributes rather than merely flags, and weaker than a validated tool,
because one of its conditions, correct release-level *N*, is the very kind of metadata the paper
says cannot be taken on trust.

### Coverage gaps this simulation does not close

1. **K above 30.** 16 of 41 audit records have more instruments than any simulated release.
2. **Meta-analytic log-odds releases.** Post hoc, applying the frozen τ to the 41 audit records:
   5 exceed it — ebi-a-GCST90018915 (0.264), ieu-a-8 (0.246), ieu-a-93 (0.224), ieu-a-91 (0.203),
   ebi-a-GCST90018910 (0.186) — all assigned log odds, none discordant. Real meta-analysed log-odds
   releases therefore sit farther from the log-odds signature than any simulated pipeline, most
   plausibly through heterogeneity and genomic-control inflation of standard errors, which were not
   simulated. The 11 discordant records sit at 0.0015–0.022, far inside τ.
3. **Real REGENIE and SAIGE**, relatedness, structure, leave-one-chromosome-out predictors.
4. **Adversarial arms at K ≠ 10.**

---

## 8. Files

`sim_operating_characteristics.csv`, `sim_claim_rules.csv`, `sim_separation_calibration.csv`,
`sim_tau.txt` — frozen analysis. `gate_unconditional_decomposition.csv`,
`gate_fit_selection.csv`, `gate_adversarial_primary_vs_adequacy.csv`,
`gate_logistic_nonconvergence.csv`, `integrity_check_results.csv` — this gate.
`gate_audit_coverage_posthoc.csv` — **post hoc** application of τ to the audit records.
`PRE_ANALYSIS_INTEGRITY_NOTES.md` — checks made while the run was executing; its line reading
"Agreement with the simulation … yes" is superseded by the three-level distinction in §2.2.

---

## Erratum (appended 25 September 2026, before any comparator result was inspected)

**E1. R4 uses the unconditional proportion.** Protocol §7 defines the misattribution rate as the
"proportion assigned a scale other than the true one", without the "non-abstained" qualifier that
it attaches to accuracy. The frozen analysis script instead computed it as one minus accuracy among
non-abstained releases. For the R4 decision the protocol definition governs: **INFO 0.5 gives
0.022, below 0.05, and does not trigger R4.** R4 remains triggered by arm A5 (N × 0.5, 2 and 5),
whose sufficiency abstention is zero, so the conditional and unconditional values coincide there.
The 0.314 figure is retained only as a secondary, conditional description.

**E2. Prevalence substitution affects only the case-fraction-dependent pipelines.** Supplying
population prevalence 0.01 in place of the case fraction (arm A4b) misassigned every
linear-probability and log-odds release — P1, P3, P4, P5, P6 and P7, all to *standardised* — while
**P2, the standardised pipeline, was assigned correctly in every release**, because the
standardised signature does not depend on *p*. The gate-verdict row and §7 wording "misattributed
100%" are to be read with this restriction.

**E3. The "eleven coincidences" argument in §1.4 is withdrawn.** The 11 discordant records lie in
two release families — nine in `ukb-a`, two in `ukb-b` — and may share a single systematic
convention, so they are not independent tests. The *N* factors that a truly standardised release
would require (4.17 to 65.43, each equal to 1/(*p*(1 − *p*)) within 0.022 log units) are retained
as a descriptive statement only. What they show is that the *N*-misreporting alternative requires
the misreporting factor to track each record's own case fraction; they do not supply eleven
independent pieces of evidence against it, and a family-level convention that produced exactly
this dependence cannot be excluded from the data at hand.

**E4. The A7 "matched" arm did not reproduce the CheckSumStats slope-0.23 mechanism** (found in
the post hoc comparator, `COMPARATOR_ADJUDICATION.md` §4). A7 matched the *study-count
distribution of candidate variants* to the CheckSumStats colorectal example and set variant-level
*N* proportional to study count, i.e. equal-sized studies. After genome-wide-significant selection,
38–41% of selected variants had all 15 studies (median 14, minimum 4), and the CheckSumStats slope
computed on the selected releases had median 0.96–1.03, with at most 3% of releases below 0.8 and
at most 1% at or below 0.30. The 0.23 anomaly was therefore not reproduced, and the gate's
statement that the rule is robust to "variant-level sample-size heterogeneity including the
CheckSumStats slope-0.23 configuration" is withdrawn. What was shown is robustness to the
heterogeneity *as simulated*; robustness to a heterogeneity large enough to produce a 0.23 slope is
untested.
