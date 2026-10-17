# Simulation protocol — operating characteristics of the three-signature scale check

**Status: frozen before any simulated record was classified.** The SHA-256 of this file and of
`sim_scale.py` are recorded in `FREEZE.txt` before execution. The only computations run before
freezing are the code-validation checks in §11, which compare estimator implementations against
reference implementations and measure runtime. None of them computes σ̂, a signature distance or
a classification.

---

## 1. Question

The manuscript applies a rule that assigns a reporting scale — linear probability, standardised,
or log odds — to a binary-trait summary-statistics release from the reported standard errors,
effect-allele frequencies, sample size and case fraction. It has been applied to one direct case
and 41 audit records, and calibrated against three paired releases. It has never been evaluated
on releases whose true scale is known by construction.

Both external reviewers asked for this (Reviewer A, major 2; Reviewer B, ②). The simulation
answers four questions:

1. **Accuracy.** When the true scale is known, how often does the rule assign it, and how often
   does it assign each wrong scale?
2. **Pipeline coverage.** Does the rule hold for estimators outside the three it was derived
   for — Firth-penalised logistic regression and pipelines whose standard error is back-computed
   from a saddlepoint-approximated P value, which are the defaults in current biobank software?
3. **Case fraction versus prevalence.** In case–control samples, is the rule correct when given
   the analysis-sample case fraction, and how wrong is it when given population prevalence
   (Reviewer B, ①)?
4. **Misattribution.** How often does a metadata error that is *not* a scale error — a wrong
   reported *N*, a wrong allele frequency, variant-level sample-size heterogeneity, imputation
   uncertainty — cause the rule to report a scale change? This is the scenario in the
   CheckSumStats documentation in which an expected-versus-reported slope of 0.23 arose from
   variant-level sample-size heterogeneity rather than from scale. A rule that misreads that
   situation as a scale error would be worse than an alarm that does not attribute a cause.

## 2. What the simulation can and cannot establish

It can establish error rates for the rule **under the data-generating and analysis models
specified here**. It cannot establish that real software behaves as the stylised pipelines of §5
do. Relatedness, population structure, leave-one-chromosome-out polygenic predictors, genomic
control and study-specific quality control are not simulated. The results bound the rule's
behaviour under known mechanisms. They are not a validation against every real pipeline, and the
manuscript must not describe them as one.

## 3. The rule under test — fixed exactly as stated in manuscript §2.2

For a release with *K* usable variants, reported sample size *N*ᵣ, reported effect-allele
frequencies *f*ⱼ and supplied case fraction *p*:

1. σ̂ⱼ = SEⱼ × √(2 *N*ᵣ *f*ⱼ(1 − *f*ⱼ)); σ̂ = median over *j*.
2. Expected signatures: σ_lin = √(*p*(1 − *p*)), σ_std = 1, σ_lo = 1/√(*p*(1 − *p*)).
3. Distances d_s = |log(σ̂ / σ_s)|; assign the scale with minimum distance.
4. Separation = second-smallest minus smallest distance; relative separation = separation ÷
   |½ log *p*(1 − *p*)|.
5. Abstention on data sufficiency: fewer than five usable variants → no assignment.
6. **Abstention on adequacy (new; evaluated as a secondary analysis):** abstain if the minimum
   distance exceeds τ. τ is set once, on the tuning set only (§8), as the 99th percentile of the
   distance to the *true* signature among correctly specified P1–P3 records pooled across the
   core grid. It is then frozen and applied unchanged to every test cell.

**The primary analysis applies rules 1–5 only**, because that is the rule the manuscript applied
to the 41 audit records. The adequacy rule is reported alongside, never instead.

## 4. Data-generating model

Individual-level data are generated and analysed by actual regression fits; no summary statistic
is generated from the signature formulas.

- Genotypes: *G*ⱼ ~ Binomial(2, *f*ⱼ), independent across variants (no linkage disequilibrium);
  *f*ⱼ ~ Uniform(0.05, 0.50) in the core grid.
- Covariate: *C* ~ Normal(0, 1), with log-odds effect γ (γ = 0 in the core grid).
- Phenotype: logit P(*Y* = 1) = α + Σⱼ βⱼ*G*ⱼ + γ*C*, with α solved numerically so that the
  expected case fraction equals the target *p*. The realised case fraction ncase/*N* is what the
  analyst reports, as in practice.
- Effect sizes: calibrated to a target non-centrality λ_z. For each variant,
  |βⱼ| = λ_z / √(*N* · 2*f*ⱼ(1 − *f*ⱼ) · *p*(1 − *p*)), with random sign. The core value
  λ_z = 6 gives variants just beyond genome-wide significance. **Rationale:** the rule reads only
  standard errors, but two of the pipelines (P5, P6) back-compute the standard error as |β|/*z*,
  which is unstable when *z* is small. Real instruments are genome-wide significant, so a fixed
  small odds ratio would impose a noise on these pipelines that real releases do not have. The
  consequence is that larger β accompanies rarer case fractions, which is the realistic direction.
- *N* = 20,000 per release in all cells.

## 5. Analysis pipelines

Each simulated data set is analysed by all seven pipelines, yielding seven releases with known
true scale. Each fits one variant at a time, with *C* included whenever γ ≠ 0.

| Code | Pipeline | True scale | Standard error |
|---|---|---|---|
| P1 | Least squares on 0/1 outcome (linear probability; what BOLT-LMM reports in unrelated samples) | linear probability | model-based |
| P2 | Least squares on the standardised outcome (*Y* − *Ȳ*)/sd(*Y*) | standardised | model-based |
| P3 | Logistic regression, Wald | log odds | inverse information |
| P4 | Firth-penalised logistic regression, Wald | log odds | inverse information at the penalised estimate |
| P5 | Logistic β with SE back-computed from a saddlepoint-approximated score-test P value (REGENIE-like stylisation) | log odds | \|β\| / Φ⁻¹(1 − *P*_SPA/2) |
| P6 | Score-based one-step β = *U*/*V* with SE back-computed from the SPA P value (SAIGE-like stylisation) | log odds | \|β\| / Φ⁻¹(1 − *P*_SPA/2) |
| P7 | P1 converted to log odds by dividing β and SE by *p*(1 − *p*) | log odds | converted |

Saddlepoint approximation follows the standard cumulant-generating-function construction for the
covariate-adjusted score, applied when |*z*| > 2 and otherwise replaced by the normal
approximation, as in SAIGE. P5 and P6 are stylisations of the mechanism by which those programs
report standard errors. They are not the programs.

## 6. Design

**Core grid** (24 cells): *p* ∈ {0.01, 0.02, 0.05, 0.10, 0.20, 0.35, 0.50, 0.75} × *K* ∈
{5, 10, 30}; λ_z = 6; *f* ~ U(0.05, 0.50); γ = 0; correct metadata. **1,000 data sets per cell.**

**Arms** (each departs from the core in one factor; *K* = 10; *p* ∈ {0.02, 0.10, 0.35, 0.50};
**500 data sets per arm cell**):

| Arm | Factor | Levels | Question |
|---|---|---|---|
| A1 | Instrument strength λ_z | 10, 20 | Does the log-odds approximation hold for larger effects? |
| A2 | Low allele frequency | *f* ~ U(0.01, 0.05) | Does the rare-variant regime shift the logistic signature? |
| A3 | Covariate adjustment | γ = log 2, log 4 | Does residual-variance reduction (linear) and non-collapsibility (logistic) move the signatures? |
| A4 | Case–control ascertainment | population prevalence 0.01, sampled case fraction 0.5 and 0.2; rule given (a) case fraction, (b) prevalence | Which quantity must the rule be given? |
| A5 | Misreported *N* | reported = true × {0.5, 2, 5} | How large a sample-size error is misread as a scale change? |
| A6 | Allele-frequency error | reference *f* with logit-scale noise SD 0.2, 0.5; and 20% of variants allele-flipped (*f* → 1 − *f*) | Does frequency substitution error shift the assignment? |
| A7 | Variant-level sample-size heterogeneity | (i) matched to the CheckSumStats example: contributing studies per variant with minimum 1, median 14, mean ≈ 11.7, maximum 15, *N*ⱼ = *N* × studiesⱼ/15, reported *N* = *N*; (ii) severe: half of variants at *N*/4 | Is the slope-0.23 mechanism misread as a scale error by a median-based rule? |
| A8 | Imputation uncertainty | INFO = 0.8, 0.5, applied as dosage with squared correlation INFO to the true genotype | How much does SE inflation from imperfect imputation move σ̂? |

For A4 the ascertained case fraction 0.5 and 0.2 replace the four-level *p* set. For A5–A8 the
truth is that the *scale* is correct and only metadata or data quality departs, so any assignment
other than the true scale is a **misattribution**.

## 7. Estimands

For every pipeline × cell:

- **Accuracy** — proportion of non-abstained releases assigned the true scale;
- **Confusion** — proportion assigned each wrong scale;
- **Abstention rate** — overall and by reason;
- **False-discordance rate** — probability that a release whose declared unit equals its true
  scale is flagged as discordant (= 1 − accuracy among non-abstained, for that pipeline);
- **Detection rate for the audit's discordance type** — probability that a linear-probability
  release (P1) declared as standardised is flagged (= proportion of P1 releases not assigned
  standardised);
- **Misattribution rate** (A5–A8) — proportion assigned a scale other than the true one;
- **Separation calibration** — accuracy by relative-separation decile, pooled over the core grid.

Monte Carlo standard error for a proportion near 0.95 is 0.007 at 1,000 replicates and 0.010 at
500; near 0.99, 0.003 and 0.004.

## 8. Tuning–test separation

Seeds for the tuning set (200 data sets per core cell; used only to set τ in §3.6) and for every
test cell are drawn from disjoint streams of a single `numpy.random.SeedSequence(20260924)` via
`spawn`, with the tuning stream spawned first. No quantity is tuned on the 41 audit records, on
the three calibration pairs, or on any test cell.

## 9. Pre-specified consequences for the manuscript

A cell **meets the performance criterion** if, for every pipeline, accuracy among non-abstained
releases is ≥ 0.95. These rules are fixed now, before any result is seen:

- **R1 — broad support.** If every core cell with *K* ≥ 5 meets the criterion, together with A1,
  A2, A3 and A4(a), the manuscript may state that in simulation the rule identified the reporting
  scale correctly across the stated conditions, including Firth and saddlepoint-based pipelines,
  and scale attribution may be listed as a contribution.
- **R2 — bounded support.** If failures are confined to identifiable regions, the manuscript
  states the domain in which the rule is supported, identifies any of the 41 audit records that
  fall outside it, and re-reports the audit result restricted to records within it.
- **R3 — pipeline exclusion.** If P5 or P6 fails in any core cell, the Methods must state that
  the rule is not supported for saddlepoint-based pipelines in that region.
- **R4 — misattribution.** Any arm in A5–A8 with a misattribution rate above 0.05 is reported in
  the Discussion as a failure mode, with the magnitude of the metadata error that produces it.
- **R5 — ascertainment.** If A4(b) fails where A4(a) passes, that is reported as direct evidence
  that the rule requires the analysis-sample case fraction.

No outcome of the simulation alters the invariance result, the calibration pairs or the
case-study estimates, and the manuscript must not imply otherwise.

## 10. Known limitations

No linkage disequilibrium, relatedness, population structure or polygenic background; no
selection of instruments on significance, and so no winner's curse; one sample size; stylised
rather than actual REGENIE and SAIGE; logistic generative model only; a single covariate. The
standardised pipeline is one plausible construction of "standardised effects" for a binary trait,
not the only one.

## 11. Pre-execution checks (the only computations before freezing)

1. Least-squares and logistic estimates and SEs against `statsmodels` on one data set.
2. Firth estimate against direct numerical maximisation of the penalised log-likelihood.
3. Saddlepoint P value against the normal approximation for a common variant at *p* = 0.5
   (should agree) and against a 20,000-permutation P value for a rare-case-fraction variant
   (should be closer than the normal approximation).
4. Runtime per data set in the largest cell.

## 11a. Amendment 1 — made before freezing, before any classification was computed

**Change.** Instruments are now selected on significance, as in practice. Each data set
generates 3*K* candidate variants. For each release, only variants reaching *P* < 5 × 10⁻⁸ under
that pipeline's own test are retained (Wald for P1–P4 and P7; the saddlepoint P value for P5 and
P6), and the first *K* passing variants in random order are used. A release with fewer than five
passing variants abstains on sufficiency, and the number passing is recorded.

**Reason.** Pre-execution check V3 showed that in the rare regime — case fraction 0.01 and minor
allele frequency below 0.02 — the score statistic is strongly skewed. Protective variants cannot
reach large |*z*| because the minor-allele count in cases is bounded below by zero, and they gave
|*z*| ≈ 1–2 while risk variants of the same target non-centrality gave |*z*| up to 14. Without
selection, P5 and P6, whose standard error is back-computed as |β|/*z*, would be scored on
variants that could never serve as instruments, manufacturing a failure that no real release
would show. Selection also removes the limitation in §10 that no winner's curse was represented.

**What this does not change.** The rule under test, the estimands, the cells, the replicate
counts, the seeds and the pre-specified consequences in §9 are unchanged.

**Pre-execution check results** (full output in `precheck_results.csv`): least-squares and
logistic estimates and standard errors matched `statsmodels` to ≤ 1 × 10⁻¹⁴, including with a
covariate and a variant-level mask; Firth estimates matched direct maximisation of the penalised
likelihood to 1.4 × 10⁻⁷; the saddlepoint P value agreed with the normal approximation for common
variants at *p* = 0.5 (relative difference in log *P* ≤ 5 × 10⁻⁴), and for six rare-regime
variants with |*z*| between 2.8 and 3.5 it was within 3 Monte Carlo standard errors of a
100,000-replicate Monte Carlo P value (mean absolute error 1.37 SE), where the normal
approximation erred by a mean of 8.69 SE. A first version of this check, at 20,000 replicates, was
uninformative because two of its three variants fell below the Monte Carlo resolution floor; it
is recorded, not discarded.

## 12. Outputs

`sim_results_records.parquet` (one row per release), `sim_operating_characteristics.csv`
(estimands of §7), `sim_claim_rules.csv` (evaluation of R1–R5), `sim_session_info.txt`,
`FREEZE.txt`.
