#!/usr/bin/env Rscript

# Offline reproduction of the submission-stage results.
# This script uses only the accompanying derived CSV files and requires no API token.

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(MRPRESSO)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)
set.seed(20260909)

script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)])
root <- if (length(script_file)) dirname(normalizePath(script_file)) else getwd()
data_dir <- file.path(root, "data")
result_dir <- file.path(root, "results")
figure_dir <- file.path(result_dir, "figures")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

write_csv <- function(x, filename) {
  write.csv(x, file.path(result_dir, filename), row.names = FALSE, na = "")
}

fmt_p <- function(p) {
  ifelse(p < 0.001, format(p, scientific = TRUE, digits = 2), sprintf("%.3f", p))
}

# -----------------------------------------------------------------------------
# 1. Repository diagnostic re-analysis from the preserved classifiable set
# -----------------------------------------------------------------------------

audit <- read.csv(file.path(data_dir, "repository", "scale_audit_catalog.csv"), check.names = FALSE)
sample98 <- read.csv(
  file.path(data_dir, "repository", "repository_sample_98_reconstructed.csv"),
  check.names = FALSE
)
sample98$classifiable_in_preserved_output <-
  tolower(as.character(sample98$classifiable_in_preserved_output)) == "true"
sample_family_counts <- table(sample98$identifier_family)
stopifnot(
  nrow(sample98) == 98,
  length(sample_family_counts) == 7,
  all(sample_family_counts == 14),
  all(audit$id %in% sample98$id),
  sum(sample98$classifiable_in_preserved_output) == nrow(audit)
)
audit$distance_linear <- abs(log(audit$sigma_obs / audit$exp_linear))
audit$distance_log_odds <- abs(log(audit$sigma_obs / audit$exp_logOR))
audit$inferred_recomputed <- ifelse(
  audit$distance_linear <= audit$distance_log_odds,
  "linear_probability",
  "log_odds"
)
audit$best_distance <- pmin(audit$distance_linear, audit$distance_log_odds)
audit$second_distance <- pmax(audit$distance_linear, audit$distance_log_odds)
audit$classification_margin <- audit$second_distance - audit$best_distance
audit$declared_normalized <- trimws(tolower(audit$unit_declared))
audit$metadata_discordant <- (
  audit$declared_normalized == "sd" & audit$inferred_recomputed == "linear_probability"
) | (
  audit$declared_normalized %in% c("logor", "log odds") &
    audit$inferred_recomputed == "linear_probability"
)

perturb_factors <- c(0.80, 0.90, 0.95, 1.00, 1.05, 1.10, 1.20)
perturb <- do.call(rbind, lapply(perturb_factors, function(factor_n) {
  sigma_p <- audit$sigma_obs * sqrt(factor_n)
  d_linear <- abs(log(sigma_p / audit$exp_linear))
  d_log <- abs(log(sigma_p / audit$exp_logOR))
  inferred <- ifelse(d_linear <= d_log, "linear_probability", "log_odds")
  data.frame(
    id = audit$id,
    sample_size_factor = factor_n,
    inferred_scale = inferred,
    changed_from_primary = inferred != audit$inferred_recomputed
  )
}))

sd_rows <- audit[which(audit$declared_normalized == "sd"), ]
repository_summary <- data.frame(
  metric = c(
    "reconstructed_stratified_sample_size",
    "identifier_families_in_sample",
    "datasets_sampled_per_identifier_family",
    "classifiable_datasets_in_preserved_output",
    "declared_SD_and_inferred_linear_probability",
    "identifier_family_ukb_a",
    "identifier_family_ukb_b",
    "identifier_family_ukb_d",
    "median_classification_margin",
    "minimum_classification_margin",
    "classification_changes_under_N_perturbation",
    "SD_label_carriers_in_frozen_catalogue",
    "eligible_binary_European_datasets_in_frozen_catalogue"
  ),
  value = c(
    nrow(sample98),
    length(sample_family_counts),
    unique(as.numeric(sample_family_counts)),
    nrow(audit),
    nrow(sd_rows[sd_rows$inferred_recomputed == "linear_probability", ]),
    sum(sd_rows$pref == "ukb-a", na.rm = TRUE),
    sum(sd_rows$pref == "ukb-b", na.rm = TRUE),
    sum(sd_rows$pref == "ukb-d", na.rm = TRUE),
    median(audit$classification_margin),
    min(audit$classification_margin),
    sum(perturb$changed_from_primary),
    nrow(read.csv(file.path(data_dir, "repository", "SD_label_carriers.csv"))),
    nrow(read.csv(file.path(data_dir, "repository", "eligible_binary_european_catalogue.csv")))
  )
)

write_csv(audit, "repository_scale_reanalysis.csv")
write_csv(perturb, "N_perturbation_sensitivity.csv")
write_csv(repository_summary, "repository_scale_summary.csv")
write_csv(data.frame(
  metric = c(
    "sample_rows",
    "identifier_families",
    "minimum_rows_per_family",
    "maximum_rows_per_family",
    "classifiable_rows",
    "all_classified_ids_in_sample"
  ),
  value = c(
    nrow(sample98),
    length(sample_family_counts),
    min(sample_family_counts),
    max(sample_family_counts),
    sum(sample98$classifiable_in_preserved_output),
    all(audit$id %in% sample98$id)
  )
), "repository_sample_reconstruction_check.csv")

p_margin <- ggplot(audit, aes(x = reorder(id, classification_margin), y = classification_margin,
                              colour = inferred_recomputed)) +
  geom_point(size = 2.1) +
  coord_flip() +
  scale_colour_manual(values = c(linear_probability = "#2A6F97", log_odds = "#C65D3A")) +
  labs(x = NULL, y = "Classification margin (log-distance units)", colour = "Empirical scale") +
  theme_classic(base_size = 10) +
  theme(legend.position = "top", axis.text.y = element_text(size = 6.5))
ggsave(file.path(figure_dir, "repository_classification_margins.png"), p_margin,
       width = 7.2, height = 8.2, dpi = 320)
ggsave(file.path(figure_dir, "repository_classification_margins.pdf"), p_margin,
       width = 7.2, height = 8.2)

# -----------------------------------------------------------------------------
# 2. Empirical calibration from stored SNP-level associations
# -----------------------------------------------------------------------------

calibration_spec <- data.frame(
  pair = c("Help-seeking: MRC-IEU linear vs logistic",
           "Help-seeking: Neale linear vs logistic",
           "Mood swings: MRC-IEU linear vs logistic"),
  file = c("GP_visit_MRC_IEU_SNP_data.csv",
           "GP_visit_Neale_SNP_data.csv",
           "Mood_swings_SNP_data.csv"),
  prevalence = c(0.3450365567, 0.3441517, 0.4526205),
  placement = c("Main", "Main", "Main")
)

calibration <- do.call(rbind, lapply(seq_len(nrow(calibration_spec)), function(i) {
  x <- calibration_spec[i, ]
  d <- read.csv(file.path(data_dir, "calibration", x$file), check.names = FALSE)
  d <- d[isTRUE(d$mr_keep) | d$mr_keep %in% TRUE, ]
  fit <- lm(beta.outcome ~ 0 + beta.exposure, data = d)
  data.frame(
    pair = x$pair,
    shared_snps = nrow(d),
    prevalence = x$prevalence,
    observed_slope = unname(coef(fit)[1]),
    slope_se = unname(summary(fit)$coefficients[1, 2]),
    predicted_lambda = 1 / (x$prevalence * (1 - x$prevalence)),
    correlation = cor(d$beta.exposure, d$beta.outcome),
    placement = x$placement
  )
}))
write_csv(calibration, "calibration_reanalysis.csv")

cal_points <- do.call(rbind, lapply(seq_len(nrow(calibration_spec)), function(i) {
  x <- calibration_spec[i, ]
  d <- read.csv(file.path(data_dir, "calibration", x$file), check.names = FALSE)
  d$pair <- x$pair
  d$placement <- x$placement
  d
}))
cal_points <- cal_points[cal_points$placement == "Main", ]

p_cal <- ggplot(cal_points, aes(beta.exposure, beta.outcome)) +
  geom_point(size = 1.7, alpha = 0.80, colour = "#2A6F97") +
  geom_smooth(method = "lm", formula = y ~ 0 + x, se = FALSE,
              colour = "#C65D3A", linewidth = 0.8) +
  facet_wrap(~pair, scales = "free", ncol = 3) +
  labs(x = "Linear-probability effect", y = "Log-odds effect") +
  theme_classic(base_size = 9) +
  theme(strip.text = element_text(size = 8))
ggsave(file.path(figure_dir, "calibration_pairs.png"), p_cal,
       width = 8.2, height = 3.5, dpi = 320)
ggsave(file.path(figure_dir, "calibration_pairs.pdf"), p_cal,
       width = 8.2, height = 3.5)

# -----------------------------------------------------------------------------
# 3. Offline UVMR re-analysis from stored harmonised SNP-level inputs
# -----------------------------------------------------------------------------

mr_files <- list.files(file.path(data_dir, "mr"), pattern = "harmonised_kept\\.csv$", full.names = TRUE)
primary_files <- mr_files[!grepl("ADHD|superseded", basename(mr_files))]

model_label <- function(filename) {
  x <- sub("_harmonised_kept\\.csv$", "", basename(filename))
  x <- sub("Tier_1_affective_instability", "Affective instability", x)
  x <- sub("Tier_2_help_seeking_corrected", "Help-seeking (corrected)", x)
  x <- sub("Tier_3_major_depression_meta_analysis", "Major depression meta-analysis", x)
  x <- sub("__UK_Biobank_asthma", " | UK Biobank asthma", x)
  x <- sub("__FinnGen_asthma", " | FinnGen asthma", x)
  x
}

mr_methods <- c("mr_ivw_mre", "mr_weighted_median", "mr_weighted_mode",
                "mr_simple_mode", "mr_egger_regression")

mr_all <- list()
het_all <- list()
egger_all <- list()
strength_all <- list()

for (f in primary_files) {
  d <- read.csv(f, check.names = FALSE)
  d <- d[d$mr_keep %in% TRUE, ]
  label <- model_label(f)
  set.seed(20260909)
  mr_res <- TwoSampleMR::mr(d, method_list = mr_methods)
  mr_res$model <- label
  mr_res$OR <- exp(mr_res$b)
  mr_res$OR_low <- exp(mr_res$b - 1.96 * mr_res$se)
  mr_res$OR_high <- exp(mr_res$b + 1.96 * mr_res$se)
  mr_all[[label]] <- mr_res

  het <- TwoSampleMR::mr_heterogeneity(d, method_list = c("mr_ivw", "mr_egger_regression"))
  het$model <- label
  het$I2 <- pmax(0, (het$Q - het$Q_df) / het$Q)
  het_all[[label]] <- het

  egger <- TwoSampleMR::mr_pleiotropy_test(d)
  egger$model <- label
  egger_all[[label]] <- egger

  strength_all[[label]] <- data.frame(
    model = label,
    harmonised_snps = nrow(d),
    mean_F = mean((d$beta.exposure / d$se.exposure)^2),
    minimum_F = min((d$beta.exposure / d$se.exposure)^2)
  )
}

mr_results <- do.call(rbind, mr_all)
heterogeneity <- do.call(rbind, het_all)
egger_results <- do.call(rbind, egger_all)
instrument_strength <- do.call(rbind, strength_all)
rownames(mr_results) <- rownames(heterogeneity) <- rownames(egger_results) <- NULL

write_csv(mr_results, "MR_rerun_five_estimators.csv")
write_csv(heterogeneity, "MR_rerun_heterogeneity.csv")
write_csv(egger_results, "MR_rerun_egger_intercept.csv")
write_csv(instrument_strength, "MR_rerun_instrument_strength.csv")

ivw <- mr_results[grepl("Inverse variance weighted", mr_results$method), ]
ivw$exposure_group <- sub(" \\|.*", "", ivw$model)
ivw$outcome_group <- ifelse(grepl("UK Biobank", ivw$model), "UK Biobank asthma", "FinnGen asthma")
ivw$exposure_group <- factor(ivw$exposure_group,
  levels = c("Affective instability", "Help-seeking (corrected)", "Major depression meta-analysis"))

p_forest <- ggplot(ivw, aes(x = OR, y = exposure_group, colour = outcome_group)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = OR_low, xmax = OR_high), width = 0.12,
                orientation = "y", position = position_dodge(width = 0.45)) +
  geom_point(size = 2.2, position = position_dodge(width = 0.45)) +
  scale_colour_manual(values = c("UK Biobank asthma" = "#2A6F97", "FinnGen asthma" = "#C65D3A")) +
  labs(x = "Odds ratio for asthma", y = NULL, colour = NULL) +
  theme_classic(base_size = 10) +
  theme(legend.position = "top")
ggsave(file.path(figure_dir, "corrected_MR_forest.png"), p_forest,
       width = 7.2, height = 3.4, dpi = 320)
ggsave(file.path(figure_dir, "corrected_MR_forest.pdf"), p_forest,
       width = 7.2, height = 3.4)

# -----------------------------------------------------------------------------
# 4. Numerical invariance/equivariance test using the superseded original Tier 2 analysis (ukb-b-6991)
# -----------------------------------------------------------------------------

superseded_file <- mr_files[grepl("superseded_mixed_scale.*UK_Biobank", mr_files)]
superseded <- read.csv(superseded_file, check.names = FALSE)
superseded <- superseded[superseded$mr_keep %in% TRUE, ]
lambda <- 1 / (0.3450365567 * (1 - 0.3450365567))
rescaled <- superseded
rescaled$beta.exposure <- rescaled$beta.exposure * lambda
rescaled$se.exposure <- rescaled$se.exposure * lambda
rescaled$exposure <- paste0(rescaled$exposure, " [rescaled]")
rescaled$id.exposure <- paste0(rescaled$id.exposure, "_rescaled")

run_mr <- function(d) {
  set.seed(20260909)
  x <- TwoSampleMR::mr(d, method_list = mr_methods)
  x[, c("method", "nsnp", "b", "se", "pval")]
}

orig_mr <- run_mr(superseded)
scaled_mr <- run_mr(rescaled)
unit_mr <- merge(orig_mr, scaled_mr, by = c("method", "nsnp"), suffixes = c("_original", "_rescaled"))
unit_mr$estimate_ratio_original_to_rescaled <- unit_mr$b_original / unit_mr$b_rescaled
unit_mr$se_ratio_original_to_rescaled <- unit_mr$se_original / unit_mr$se_rescaled
unit_mr$absolute_z_difference <- abs(unit_mr$b_original / unit_mr$se_original -
                                      unit_mr$b_rescaled / unit_mr$se_rescaled)
unit_mr$absolute_p_difference <- abs(unit_mr$pval_original - unit_mr$pval_rescaled)
write_csv(unit_mr, "invariance_estimator_unit_test.csv")

orig_het <- TwoSampleMR::mr_heterogeneity(superseded, method_list = c("mr_ivw", "mr_egger_regression"))
scaled_het <- TwoSampleMR::mr_heterogeneity(rescaled, method_list = c("mr_ivw", "mr_egger_regression"))
unit_het <- merge(orig_het[, c("method", "Q", "Q_df", "Q_pval")],
                  scaled_het[, c("method", "Q", "Q_df", "Q_pval")],
                  by = c("method", "Q_df"), suffixes = c("_original", "_rescaled"))
unit_het$absolute_Q_difference <- abs(unit_het$Q_original - unit_het$Q_rescaled)
unit_het$absolute_p_difference <- abs(unit_het$Q_pval_original - unit_het$Q_pval_rescaled)

orig_egger <- TwoSampleMR::mr_pleiotropy_test(superseded)
scaled_egger <- TwoSampleMR::mr_pleiotropy_test(rescaled)
unit_egger <- data.frame(
  statistic = c("Egger intercept", "Egger intercept SE", "Egger intercept p"),
  original = c(orig_egger$egger_intercept, orig_egger$se, orig_egger$pval),
  rescaled = c(scaled_egger$egger_intercept, scaled_egger$se, scaled_egger$pval)
)
unit_egger$absolute_difference <- abs(unit_egger$original - unit_egger$rescaled)

f_original <- (superseded$beta.exposure / superseded$se.exposure)^2
f_rescaled <- (rescaled$beta.exposure / rescaled$se.exposure)^2
unit_f <- data.frame(
  statistic = c("Mean F", "Maximum SNP-level F difference"),
  original = c(mean(f_original), max(f_original)),
  rescaled = c(mean(f_rescaled), max(f_rescaled)),
  absolute_difference = c(abs(mean(f_original) - mean(f_rescaled)), max(abs(f_original - f_rescaled)))
)

run_presso <- function(d) {
  set.seed(20260909)
  MRPRESSO::mr_presso(
    BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
    SdOutcome = "se.outcome", SdExposure = "se.exposure",
    OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = d,
    NbDistribution = 1000, SignifThreshold = 0.05
  )
}

presso_original <- run_presso(superseded)
presso_rescaled <- run_presso(rescaled)
extract_presso <- function(x) {
  global <- x[["MR-PRESSO results"]][["Global Test"]]
  main <- x[["Main MR results"]]
  data.frame(
    global_RSSobs = as.numeric(global$RSSobs),
    global_p = as.character(global$Pvalue),
    raw_estimate = main$`Causal Estimate`[main$`MR Analysis` == "Raw"],
    raw_se = main$Sd[main$`MR Analysis` == "Raw"],
    raw_p = main$`P-value`[main$`MR Analysis` == "Raw"],
    outliers = paste(x[["MR-PRESSO results"]][["Distortion Test"]][["Outliers Indices"]], collapse = ";")
  )
}
presso_comparison <- rbind(
  cbind(scale = "Original", extract_presso(presso_original)),
  cbind(scale = "Rescaled", extract_presso(presso_rescaled))
)

write_csv(unit_het, "invariance_heterogeneity_unit_test.csv")
write_csv(unit_egger, "invariance_egger_unit_test.csv")
write_csv(unit_f, "invariance_F_unit_test.csv")
write_csv(presso_comparison, "invariance_MRPRESSO_unit_test.csv")

correction <- read.csv(file.path(data_dir, "repository", "Tier2_scale_correction_source_data.csv"))
correction$analysis <- factor(correction$analysis, levels = correction$analysis)
p_correction <- ggplot(correction, aes(x = OR, y = analysis)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = low, xmax = high), width = 0.12,
                orientation = "y", colour = "#2A6F97") +
  geom_point(size = 2.5, colour = "#C65D3A") +
  labs(x = "Odds ratio for asthma", y = NULL) +
  theme_classic(base_size = 10)
ggsave(file.path(figure_dir, "Tier2_scale_correction.png"), p_correction,
       width = 7.2, height = 3.0, dpi = 320)
ggsave(file.path(figure_dir, "Tier2_scale_correction.pdf"), p_correction,
       width = 7.2, height = 3.0)

# -----------------------------------------------------------------------------
# 5. Claim-evidence manifest
# -----------------------------------------------------------------------------
# results/claim_evidence_manifest.csv is maintained by hand (claim, value, evidence file, manuscript
# location) and is not written here; tests/validate_claims.py recomputes every claim from the
# result files this script writes.


session_lines <- capture.output(print(sessionInfo()))
session_lines <- sub("[ \\t]+$", "", session_lines)
writeLines(session_lines, file.path(result_dir, "sessionInfo.txt"))

cat("Offline submission analysis completed.\n")
cat("Outputs:", normalizePath(result_dir), "\n")
