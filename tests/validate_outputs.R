#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- dirname(dirname(normalizePath(script_file)))
results <- file.path(root, "results")

assert_equal <- function(actual, expected, label) {
  if (!identical(as.character(actual), as.character(expected))) {
    stop(label, ": expected ", expected, ", got ", actual)
  }
}

assert_close <- function(actual, expected, tolerance, label) {
  if (length(actual) != 1L || !is.finite(actual) || abs(actual - expected) > tolerance) {
    stop(label, ": expected ", expected, " +/- ", tolerance, ", got ", actual)
  }
}

summary <- read.csv(file.path(results, "repository_scale_summary.csv"))
value_for <- function(metric) summary$value[match(metric, summary$metric)]

assert_equal(value_for("reconstructed_stratified_sample_size"), 98, "sample size")
assert_equal(value_for("classifiable_datasets_in_preserved_output"), 41, "classifiable count")
assert_equal(value_for("declared_SD_and_inferred_linear_probability"), 11, "discordant SD-labelled count")
assert_equal(value_for("identifier_family_ukb_a"), 9, "ukb-a count")
assert_equal(value_for("identifier_family_ukb_b"), 2, "ukb-b count")
assert_equal(value_for("identifier_family_ukb_d"), 0, "ukb-d count")
assert_equal(value_for("classification_changes_under_N_perturbation"), 0, "N-perturbation changes")
assert_equal(value_for("SD_label_carriers_in_frozen_catalogue"), 1007, "SD-label carrier count")
assert_equal(value_for("eligible_binary_European_datasets_in_frozen_catalogue"), 2658, "eligible catalogue count")

cal <- read.csv(file.path(results, "calibration_reanalysis.csv"))
assert_equal(nrow(cal), 3, "same-field calibration count")
assert_equal(paste(cal$shared_snps, collapse = ","), "35,15,54", "calibration SNP counts")
for (i in seq_along(c(4.49327980089621, 4.28617287964186, 4.06922121854665))) {
  assert_close(cal$observed_slope[i], c(4.49327980089621, 4.28617287964186, 4.06922121854665)[i], 1e-10,
               paste("calibration slope", i))
}

mr <- read.csv(file.path(results, "MR_rerun_five_estimators.csv"))
ivw <- mr[grepl("Inverse variance weighted", mr$method), ]
tier2_ukb <- ivw[ivw$model == "Help-seeking (corrected) | UK Biobank asthma", ]
assert_equal(tier2_ukb$nsnp, 26, "corrected Tier 2 UKB SNP count")
assert_close(tier2_ukb$OR, 1.06794887832226, 1e-10, "corrected Tier 2 UKB OR")
assert_close(tier2_ukb$pval, 0.507343581367476, 1e-10, "corrected Tier 2 UKB P value")

manifest <- read.csv(file.path(results, "claim_evidence_manifest.csv"))
claim_value <- function(claim) manifest$value[match(claim, manifest$claim)]
assert_close(claim_value("Tier 2 mixed-scale UKB OR"), 3.1217485705749, 1e-10, "mixed-scale OR")
assert_close(claim_value("Tier 2 arithmetic-rescaling OR"), 1.29338307178027, 1e-10, "rescaled OR")
assert_close(claim_value("Tier 2 replacement-release UKB OR"), 1.06794887832226, 1e-10, "replacement OR")

cat("All claim-level validation checks passed.\n")
