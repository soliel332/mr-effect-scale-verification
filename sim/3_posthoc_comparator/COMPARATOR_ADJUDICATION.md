# Comparator adjudication (corrected) — three-signature rule versus CheckSumStats

**Status: POST HOC benchmark**, outside the frozen protocol. The frozen simulation and its outputs were
not modified, and no manuscript file was changed. This version supersedes the first comparator
adjudication. Its corrections are listed in §0.

## 0. Corrections to the first version

1. **Regeneration covered 86 test cells, not 110.** The 110 of the plan include 24 tuning cells, which
   were not part of the benchmark. The log-odds benchmark regenerated 123,200 release records: 86 cells ×
   200 replicates × 7 pipelines, plus the prevalence-supplied duplicates of the two A4 cells. The
   declared-SD benchmark regenerated 120,400 records. Both matched the frozen records exactly: 0 label
   mismatches, maximum absolute difference in σ̂ = 0.
2. **Comparisons are now paired by target.** For *anomaly detection*, the CheckSumStats slope alarm is
   compared with the pre-specified three-signature **adequacy** alarm. For *scale attribution*, the
   slope-implied scale is compared with the three-signature **primary** assignment. The first version
   compared the slope alarm with the primary rule, and on that basis called the A5 result "decisive"
   for CheckSumStats. That description is withdrawn: the adequacy rule detects the same sample-size
   errors (§2.1).
3. **The prevalence result is an input-interface distinction, not a statistical immunity** (§2.3).
4. **The Haycock paper does state the continuous-trait expected-effect formula** in its main text,
   β = z / √(2p(1 − p)(n + z²)). The first version said only that the binary-trait formula was absent.
   That statement remains true, but it omitted this one. The formula is the per-variant term of
   CheckSumStats' `predict_beta_sd()` and is central to §3.

## 1. Sources

- CheckSumStats `github.com/MRCIEU/CheckSumStats`, commit `ba670767d42a9f6041b65e4a155b64a3ac096434`
  (2 May 2023), file `R/make_plot_predlnor.R`: `predict_lnor_sh()`, `format_data_predlnor_sh()`,
  `predict_beta_sd()`, `make_plot_pred_effect()`.
- Haycock et al. 2023, main text (the supplement was not read).
- The CheckSumStats README of 23 September 2026.
- The alarm rule "slope > 1.20 or < 0.8" is from the paper; the code applies no threshold.
- A slope requires at least three variants surviving the package filters. This choice was made here; the
  package sets none.

Denominators are **all releases** in every table. *Not assessable* means that CheckSumStats had fewer
than three usable variants, or that the three-signature rule had fewer than five genome-wide-significant
variants.

---

## 2. Declared-log-odds frame (`predict_lnor_sh`)

### 2.1 Anomaly detection — slope alarm versus adequacy alarm

| Setting | Releases | CheckSumStats alarm (not assessable) | 3-sig adequacy alarm (not assessable) |
|---|---|---|---|
| Core K = 5 | log odds, no error | 0.110 (0.022) | 0.012 (0.019) |
| Core K = 10 | log odds, no error | 0.056 (0.006) | 0.014 (0.002) |
| Core K = 30 | log odds, no error | 0.068 (0.085) | 0.013 (0.092) |
| Core K = 10 | P1 / P2, scale error | 0.997 (0.003) | 0.999 (0.001) |
| A2 low MAF | log odds, no error | 0.067 (0.394) | **0.286** (0.000) |
| A2 low MAF | P1 / P2, scale error | 0.621 (0.379) | 1.000 (0.000) |
| A5 N × 0.5 / × 2 / × 5 | log odds, N error | 0.997 / 0.946 / 0.994 | 1.000 / 1.000 / 1.000 |
| A5 N × 5 | P2, scale + N error | 0.500 | 0.500 |
| A6 frequency error | log odds | 0.034–0.048 | 0.002–0.006 |
| A7 N heterogeneity (as simulated) | log odds | 0.044–0.072 | 0.000–0.004 |
| A8 INFO 0.8 | log odds | 0.050 (0.007) | 0.001 (0.024) |
| A8 INFO 0.5 | log odds | 0.300 (0.662) | 0.066 (0.934) |

When the target is detection of any anomaly, the two alarms have **essentially the same sensitivity to
sample-size errors**. Two things separate them:

- **Baseline rate.** The adequacy alarm has a lower baseline alarm rate with 5–30 instruments.
- **Low-MAF behaviour.** At low MAF, CheckSumStats cannot assess 39% of releases because of its
  minor-allele-count filter, while the adequacy rule assesses them all but raises a false alarm on 29% of
  correctly scaled log-odds releases. Neither behaviour dominates the other.

The combination "standardised effects with a five-fold overstated N at case fraction near 0.5" is missed
by both methods at the same rate.

### 2.2 Scale attribution — slope-implied scale versus primary assignment

The slope-implied scale is the scale among {1/(p(1 − p)), 1/√(p(1 − p)), 1} nearest the CheckSumStats
slope. It is a **post hoc derived quantity, not a package feature**.

| Setting | Slope-implied correct / wrong / none | 3-sig primary correct / wrong / abstain |
|---|---|---|
| Core K = 10, P1 / P2 | 0.996 / 0.000 / 0.004 | 0.999 / 0.000 / 0.001 |
| Core K = 10, log odds | 0.988 / 0.006 / 0.006 | 0.998 / 0.000 / 0.002 |
| Core K = 5, log odds | 0.964 / 0.012 / 0.024 | 0.981 / 0.000 / 0.019 |
| Core K = 30, log odds | 0.895 / 0.018 / 0.086 | 0.908 / 0.000 / 0.092 |
| A2 low MAF, all | 0.594–0.621 / ≤ 0.012 / 0.379–0.394 | 1.000 / 0.000 / 0.000 |
| A5 N × 0.5, P2 | 0.700 / 0.300 / 0 | 0.816 / 0.184 / 0 |
| A5 N × 0.5, log odds | 0.702 / 0.298 / 0 | 0.925 / 0.075 / 0 |
| A5 N × 5, P1 / P2 | 0.250 / 0.750 / 0 | 0.250 / 0.750 / 0 |

Both methods carry nearly the same scale information. The three-signature assignment is slightly more
accurate and never produced a wrong label in the core grid, against 0.6–1.8% for the slope. Under
sample-size errors both methods mislabel, and neither should be used for attribution unless N has been
verified.

### 2.3 Prevalence versus case fraction

`predict_lnor_sh` takes case and control counts; the three-signature rule as simulated takes *p*
directly. **With correct counts, the three-signature implementation derives the analysis-sample case
fraction as ncase / (ncase + ncontrol)** — exactly how the audit catalogue's case fraction was computed,
which was verified to equal ncase / N for every record. Supplied that way, the A4b error cannot arise.
Conversely, CheckSumStats given counts that do not describe the analysed sample would be misled in the
same way. The A4b result is therefore **a property of the input interface, conditional on correct
counts**, not a statistical robustness of CheckSumStats. The practical recommendation is to supply
counts, not *p*.

---

## 3. Declared-SD frame (`predict_beta_sd`) — the audit's own question

### 3.1 The estimators are the same up to a negligible factor

At commit `ba670767`:

- `predict_beta_sd()` computes β_sd = z / √(2f(1 − f)(n + z²)) per variant — the formula printed in the
  Haycock main text.
- It then sets `sd_est = median(β / β_sd)`. Since β / z = SE:

  **sd_est = median[ SE × √(2f(1 − f)(n + z²)) ]**, where *n* is `dat[, sample_size]` and **the default is
  `sample_size = "ncontrol"`**.
- Manuscript equation (1) is **σ̂ = median[ SE × √(2Nf(1 − f)) ]**.

The per-variant ratio is √(1 + z²/N) when *n* = N. In the 120,400 regenerated releases, with *n* set
explicitly to the total analysis sample size:

| Quantity | Minimum | Median | 99th percentile | Maximum |
|---|---|---|---|---|
| Per-release median of √(1 + z²/N) | 1.00016 | 1.00096 | 1.00484 | 1.00548 |
| Largest per-variant factor | 1.00020 | 1.00141 | 1.00612 | 1.00839 |
| Release ratio sd_est / σ̂ | 1.00015 | 1.00099 | 1.00492 | 1.00624 |

Applying the three signatures to CheckSumStats' own sd_est reproduced the three-signature label in
**99.67%** of releases. All 395 disagreements lay at a discriminant boundary: 393 in A5 N × 0.5 and × 2 at
*p* = 0.5, and 2 in A8 INFO 0.5.

**Relying on the default is not harmless.** With `n = ncontrol`, the release ratio ranged from 0.493 to
0.998 (median 0.896), and label accuracy in the core grid fell from 1.000 to 0.801 at *p* = 0.5 and 0.143
at *p* = 0.75. The total analysis sample size must be passed explicitly.

**Conclusion: the manuscript's implied-SD estimator is already implemented in CheckSumStats**, as the
median-ratio `sd_est`, and its per-variant term is published in the Haycock main text. The manuscript
omits the z²/N term. With 5–30 genome-wide-significant variants and N = 20,000, this changes σ̂ by at
most 0.6% and changes no assignment away from a boundary.

### 3.2 Declared-SD comparison

The documented use for already-standardised effects is `make_plot_pred_effect(pred_beta = "beta_sd",
beta = "beta", sd_est = NULL)`. With `sd_est = NULL`, the reported β is compared directly with β_sd.

| Setting | CheckSumStats declared-SD alarm (not assessable) | 3-sig primary alarm | 3-sig adequacy alarm |
|---|---|---|---|
| Core K = 10, P1 linear (the audit's discordance type) | 1.000 (0) | 0.999 | 0.999 |
| Core K = 10, P2 standardised | 0.000 (0) | 0.000 | 0.000 |
| Core K = 10, log odds | 1.000 (0) | 0.998 | 0.998 |
| Core K = 30, P1 | 0.980 (0.020) | 0.929 | 0.929 |
| A5 N × 0.5 / × 2 / × 5, P2 (standardised, N error) | 1.000 / 1.000 / 1.000 | 0.184 / 0.056 / 0.750 | 1.000 / 1.000 / 1.000 |
| A5 N × 5, P1 | 0.500 | 0.250 | 0.500 |
| A6 EAF SD 0.5 / A7 matched / A7 severe, P2 | 0.026 / 0.015 / 0.034 | 0.000 | 0.002 / 0.000 / 0.000 |
| A8 INFO 0.5, P1 | 0.394 (0.606) | 0.069 (0.929) | 0.071 |

- In this setting the CheckSumStats slope is **the reciprocal of σ̂**: the correlation of log slope with
  −log σ̂ is 0.998, and the median of slope × σ̂ is 0.998.
- **CheckSumStats already detects the audit's discordance type** — declared SD, linear-probability
  effects — when used as documented for standardised effects. It does so with slightly better coverage,
  because it needs three variants rather than five.
- **The default, self-normalising call** (sd_est estimated from the same data) raised an alarm on 0.000
  of linear-probability releases. It cannot detect a scale error, because it rescales the reported β by
  the very quantity that carries the error.

### 3.3 What the three-signature framework adds beyond CheckSumStats

1. **Closed-form expectations for binary traits.** CheckSumStats returns sd_est, or a slope, and an
   alarm. It does not state what value sd_est should take under each binary reporting convention. The
   signatures √(p(1 − p)), 1 and 1/√(p(1 − p)) turn an alarm into a named alternative that can be
   compared with a declared unit. This was not found in the Haycock main text, the README or the package
   source; the supplement was not read.
2. **An explicit decision rule and error rates** — a three-class assignment with separation, a
   pre-specified adequacy abstention and discriminant boundaries (p(1 − p))^±¼, validated by simulation
   across seven pipelines. The sources read report no operating characteristics.
3. **Slightly lower wrong-label and false-alarm rates with few instruments**, and assessment at low MAF
   in the log-odds frame. These are modest, and are offset by the adequacy rule's false alarms at low MAF.
4. **Nothing at the estimator level.** Equation (1) is CheckSumStats' sd_est without the z²/N term.

---

## 4. A7

Unchanged from the first version, and reflected in gate erratum E4:

- 38–41% of the selected variants had all 15 studies; the median was 14 studies and the minimum 4.
- The CheckSumStats slope had a median of 0.96–1.03.
- **The slope-0.23 mechanism was not reproduced.**

## 5. Proposed claim boundary

1. **Cite CheckSumStats for the estimator.** Write equation (1) as the release-level median of the
   per-variant ratio used by `predict_beta_sd()` [Haycock et al.], and state that the z²/N term is omitted
   because it changes σ̂ by less than 1% with genome-wide-significant instruments.
2. **The declared-SD check itself is prior capability.** The paper's contribution is the binary-trait
   signatures, the decision rule, its simulated operating characteristics, and the empirical application.
3. **Terminology.**
   - Simulation with verified N and counts: *conditional scale-attribution rule*.
   - Real records: *statistical discordance between the declared unit and the indicated scale*.
4. **Pass counts and total N explicitly.** Never pass *p* or rely on `ncontrol`.
5. **Run both checks, as complements.** An alarm with a primary assignment equal to the declared unit
   points to a non-scale error. A primary assignment to a different, named signature points to a scale
   error — conditional on N being correct.

Files: `posthoc_logodds_anomaly_paired.csv`, `posthoc_logodds_attribution_paired.csv`,
`posthoc_declaredSD_comparison.csv`, `posthoc_sd_regeneration_check.csv`, `posthoc_sd.py`,
`NOVELTY_DELTA.md`.
