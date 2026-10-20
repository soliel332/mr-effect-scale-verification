#!/usr/bin/env Rscript
# Direct verification of the post hoc bridge comparison in the pinned R environment.
#
# usage: Rscript case_study/bridge/verify_bridge_pinned.R      (run from the repository root, after the
#        pinned packages have been installed; run by .github/workflows/reproduce-core.yml)
#
# Recomputes the four rows of case_study/bridge/bridge_estimates.csv with TwoSampleMR::mr_ivw_mre, the
# first-order IVW multiplicative-random-effects estimator used by the core analysis, from the same inputs as
# bridge_analysis.py: the 35 variants of the original Tier 2 analysis, their UK Biobank asthma associations and
# standard errors, and, for the bridge row only, the replacement release's exposure coefficients for the same
# variants. It then checks b, SE, OR, 95% CI, P, Q and I2 against the committed table, and b and SE of the three
# non-bridge rows against the outputs of the core analysis. Nothing is written. Any difference beyond the
# tolerances below stops the script and prints every value compared.
#
# This check supplements, and does not replace, the Python implementation and its cross-validation
# (tests/validate_bridge.py).

args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- normalizePath(file.path(dirname(script), "..", ".."))
f <- function(...) file.path(root, ...)

TOL_ABS <- 1e-8   # b, SE, I2
TOL_REL <- 1e-8   # OR, CI limits, P, Q
Z <- 1.959964     # CI multiplier used in bridge_estimates.csv

orig <- read.csv(f("data", "mr", "Tier_2_superseded_mixed_scale__UK_Biobank_asthma_harmonised_kept.csv"),
                 stringsAsFactors = FALSE)
repl <- read.csv(f("data", "mr", "Tier_2_help_seeking_corrected__UK_Biobank_asthma_harmonised_kept.csv"),
                 stringsAsFactors = FALSE)
cal <- read.csv(f("data", "calibration", "GP_visit_MRC_IEU_SNP_data.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(orig) == 35, !anyDuplicated(orig$SNP), all(cal$id.outcome == "ebi-a-GCST90013959"),
          all(orig$SNP %in% cal$SNP))
cal <- cal[match(orig$SNP, cal$SNP), ]
stopifnot(all(cal$effect_allele.outcome == orig$effect_allele.exposure),
          all(cal$other_allele.outcome == orig$other_allele.exposure))

ivw <- function(bx, sx, by, sy) {
  m <- TwoSampleMR::mr_ivw_mre(bx, by, sx, sy)
  c(b = m$b, se = m$se, P = m$pval, Q = m$Q, I2 = max(0, (m$Q - m$Q_df) / m$Q), nsnp = m$nsnp)
}
lambda <- 1 / (0.3450365567 * (1 - 0.3450365567))   # same expression as effect_scale_mr_analysis.R
res <- rbind(
  original    = ivw(orig$beta.exposure, orig$se.exposure, orig$beta.outcome, orig$se.outcome),
  rescaled    = ivw(orig$beta.exposure * lambda, orig$se.exposure * lambda, orig$beta.outcome, orig$se.outcome),
  bridge      = ivw(cal$beta.outcome, cal$se.outcome, orig$beta.outcome, orig$se.outcome),
  replacement = ivw(repl$beta.exposure, repl$se.exposure, repl$beta.outcome, repl$se.outcome)
)
res <- cbind(res, OR = exp(res[, "b"]), CI_low = exp(res[, "b"] - Z * res[, "se"]),
             CI_high = exp(res[, "b"] + Z * res[, "se"]))

ref <- read.csv(f("case_study", "bridge", "bridge_estimates.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(ref) == 4, startsWith(ref$analysis[3], "Bridge"))
ref_b <- log(ref$OR); ref_se <- (log(ref$CI_high) - log(ref$OR)) / Z

checks <- list()
add <- function(row, quantity, r_value, committed, kind) {
  d <- abs(r_value - committed); if (kind == "rel") d <- d / abs(committed)
  tol <- if (kind == "rel") TOL_REL else TOL_ABS
  checks[[length(checks) + 1]] <<- data.frame(row = row, quantity = quantity, pinned_R = r_value,
                                              committed = committed, difference = d, type = kind,
                                              tolerance = tol, pass = d <= tol)
}
for (i in 1:4) {
  rn <- rownames(res)[i]
  stopifnot(res[i, "nsnp"] == ref$n_variants[i])
  add(rn, "b", res[i, "b"], ref_b[i], "abs");          add(rn, "SE", res[i, "se"], ref_se[i], "abs")
  add(rn, "OR", res[i, "OR"], ref$OR[i], "rel");       add(rn, "CI_low", res[i, "CI_low"], ref$CI_low[i], "rel")
  add(rn, "CI_high", res[i, "CI_high"], ref$CI_high[i], "rel")
  add(rn, "P", res[i, "P"], ref$P[i], "rel");          add(rn, "Q", res[i, "Q"], ref$Q[i], "rel")
  add(rn, "I2", res[i, "I2"], ref$I2[i], "abs")
}
# the three non-bridge rows must also equal the core analysis outputs regenerated in this environment
inv <- read.csv(f("results", "invariance_estimator_unit_test.csv"), stringsAsFactors = FALSE)
inv <- inv[grepl("multiplicative random effects", inv$method), ]
fe <- read.csv(f("results", "MR_rerun_five_estimators.csv"), stringsAsFactors = FALSE)
fe <- fe[grepl("multiplicative random effects", fe$method) &
         fe$model == "Help-seeking (corrected) | UK Biobank asthma", ]
stopifnot(nrow(inv) == 1, nrow(fe) == 1)
add("original", "b vs core output", res["original", "b"], inv$b_original, "abs")
add("original", "SE vs core output", res["original", "se"], inv$se_original, "abs")
add("rescaled", "b vs core output", res["rescaled", "b"], inv$b_rescaled, "abs")
add("rescaled", "SE vs core output", res["rescaled", "se"], inv$se_rescaled, "abs")
add("replacement", "b vs core output", res["replacement", "b"], fe$b, "abs")
add("replacement", "SE vs core output", res["replacement", "se"], fe$se, "abs")

out <- do.call(rbind, checks)
options(width = 160)
print(format(out, digits = 12), row.names = FALSE)
cat("TwoSampleMR", as.character(utils::packageVersion("TwoSampleMR")), "|", R.version.string, "\n")
if (!all(out$pass)) stop(sum(!out$pass), " of ", nrow(out), " comparisons exceed tolerance (see table above)")
cat("bridge comparison verified in the pinned R environment:", nrow(out), "comparisons within tolerance\n")
