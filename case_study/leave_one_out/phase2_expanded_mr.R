#!/usr/bin/env Rscript
# Phase 2 -- expanded conventional MR outputs for the six corrected tier x outcome analyses.
#
# Inputs : the frozen harmonised-kept files in <repo>/data/mr/ ONLY.
# Outputs: <outdir>/results/*.csv and <outdir>/figures/*.png. Of these, only the per-model
#          leave-one-out tables are distributed; the Steiger table (quantitative-trait
#          approximation) is not interpretable for these binary phenotypes and is not reported.
#
# Definitions, thresholds and estimators follow the frozen repository's own
# effect_scale_mr_analysis.R: TwoSampleMR defaults, MRE-IVW as the primary
# estimator, five-estimator panel, MRPRESSO with NbDistribution = 10000 and
# SignifThreshold = 0.05, and the same seed.
#
# NON-DEFAULT PARAMETERS / DEVIATIONS -- all recorded in phase2_deviations.csv:
#   * "harmonised-all" data cannot be regenerated: the frozen package ships only
#     *_harmonised_kept.csv. Pre-filter rows are not recoverable from it.
#   * The superseded original Tier 2 analysis (ukb-b-6991) is not part of this panel.
#   * Steiger directionality needs r.exposure / r.outcome. Where the inputs to
#     derive them are absent (no outcome sample size on the FinnGen releases, no
#     case/control split in the harmonised files) the test is reported as
#     not applicable rather than run on imputed inputs.

suppressPackageStartupMessages({
  library(TwoSampleMR); library(MRPRESSO); library(ggplot2)
})
options(stringsAsFactors = FALSE)
SEED <- 20260909
set.seed(SEED)

args   <- commandArgs(trailingOnly = TRUE)
repo   <- if (length(args) >= 1) args[1] else "."
outdir <- if (length(args) >= 2) args[2] else "phase2_output"
resdir <- file.path(outdir, "results"); figdir <- file.path(outdir, "figures")
dir.create(resdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

MODELS <- c(
  "Tier_1_affective_instability__UK_Biobank_asthma",
  "Tier_1_affective_instability__FinnGen_asthma",
  "Tier_2_help_seeking_corrected__UK_Biobank_asthma",
  "Tier_2_help_seeking_corrected__FinnGen_asthma",
  "Tier_3_major_depression_meta_analysis__UK_Biobank_asthma",
  "Tier_3_major_depression_meta_analysis__FinnGen_asthma"
)
# Optional third argument restricts the run to one model, so the six can be run in
# parallel; a final invocation with no filter combines whatever per-model files exist.
ONLY <- if (length(args) >= 3 && nzchar(args[3])) args[3] else NA
if (!is.na(ONLY)) MODELS <- MODELS[MODELS == ONLY]

deviations <- list()
note <- function(model, item, detail) {
  deviations[[length(deviations) + 1]] <<-
    data.frame(model = model, item = item, detail = detail)
}

acc <- function(x) if (length(x)) do.call(rbind, x) else NULL
five <- list(); het <- list(); egg <- list(); ssnp <- list(); loo <- list()
stg  <- list(); prs <- list(); fstat <- list(); instr <- list()

# Per-model outputs are written as soon as each model finishes, and a model whose
# marker already exists is skipped, so an interrupted run resumes where it stopped.
perdir <- file.path(resdir, "per_model"); dir.create(perdir, recursive = TRUE, showWarnings = FALSE)

for (m in MODELS) {
  marker <- file.path(perdir, paste0(m, "__mrpresso.csv"))
  if (file.exists(marker)) { cat("skip (already done):", m, "\n"); next }
  f <- file.path(repo, "data", "mr", paste0(m, "_harmonised_kept.csv"))
  d <- read.csv(f, check.names = FALSE)
  d$mr_keep <- TRUE
  set.seed(SEED)

  instr[[m]] <- data.frame(
    model = m, SNP = d$SNP, chr = d$chr, pos = d$pos,
    effect_allele = d$effect_allele.exposure, other_allele = d$other_allele.exposure,
    beta.exposure = d$beta.exposure, se.exposure = d$se.exposure,
    pval.exposure = d$pval.exposure, eaf.exposure = d$eaf.exposure,
    beta.outcome = d$beta.outcome, se.outcome = d$se.outcome,
    pval.outcome = d$pval.outcome, F_stat = d$F_stat)

  # ---- five estimators -------------------------------------------------------
  r5 <- mr(d, method_list = c("mr_ivw_mre", "mr_weighted_median",
                              "mr_weighted_mode", "mr_simple_mode", "mr_egger_regression"))
  r5$model <- m; five[[m]] <- r5

  # ---- heterogeneity, Egger intercept ---------------------------------------
  h <- mr_heterogeneity(d); h$model <- m
  h$I2 <- ifelse(h$Q_df > 0, pmax(0, (h$Q - h$Q_df) / h$Q) * 100, NA); het[[m]] <- h
  e <- mr_pleiotropy_test(d); e$model <- m; egg[[m]] <- e

  # ---- single-SNP and leave-one-out -----------------------------------------
  s <- mr_singlesnp(d); s$model <- m; ssnp[[m]] <- s
  l <- mr_leaveoneout(d); l$model <- m; loo[[m]] <- l

  # ---- instrument strength ---------------------------------------------------
  fstat[[m]] <- data.frame(model = m, harmonised_snps = nrow(d),
                           mean_F = mean(d$F_stat), minimum_F = min(d$F_stat))

  # ---- Steiger ---------------------------------------------------------------
  st <- try(directionality_test(d), silent = TRUE)
  if (inherits(st, "try-error") || is.null(st)) {
    note(m, "Steiger directionality",
         "not applicable: r.exposure/r.outcome could not be derived from the harmonised-kept fields")
  } else { st$model <- m; stg[[m]] <- st }

  # ---- MR-PRESSO -------------------------------------------------------------
  set.seed(SEED)
  pr <- try(mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                      SdOutcome = "se.outcome", SdExposure = "se.exposure",
                      OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = as.data.frame(d),
                      NbDistribution = 10000, SignifThreshold = 0.05, seed = SEED),
            silent = TRUE)
  if (inherits(pr, "try-error")) {
    note(m, "MR-PRESSO", paste("failed:", conditionMessage(attr(pr, "condition"))))
  } else {
    mt <- pr$`Main MR results`
    gt <- pr$`MR-PRESSO results`$`Global Test`
    dt <- pr$`MR-PRESSO results`$`Distortion Test`
    out_idx <- pr$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
    prs[[m]] <- data.frame(
      model = m,
      raw_estimate        = mt$`Causal Estimate`[1],
      raw_se              = mt$Sd[1],
      raw_p               = mt$`P-value`[1],
      corrected_estimate  = if (nrow(mt) > 1) mt$`Causal Estimate`[2] else NA,
      corrected_se        = if (nrow(mt) > 1) mt$Sd[2] else NA,
      corrected_p         = if (nrow(mt) > 1) mt$`P-value`[2] else NA,
      global_RSSobs       = gt$RSSobs,
      global_p            = as.character(gt$Pvalue),
      n_outliers          = if (is.null(out_idx)) 0 else length(out_idx),
      outlier_indices     = if (is.null(out_idx)) "" else paste(out_idx, collapse = ";"),
      outlier_SNPs        = if (is.null(out_idx)) "" else paste(d$SNP[out_idx], collapse = ";"),
      distortion_coef     = if (is.null(dt)) NA else dt$`Distortion Coefficient`,
      distortion_p        = if (is.null(dt)) NA else dt$Pvalue,
      NbDistribution      = 10000, SignifThreshold = 0.05, seed = SEED)
    if (is.null(dt))
      note(m, "MR-PRESSO distortion test",
           "not returned: no significant outliers detected, so no outlier-corrected estimate exists")
  }

  # ---- plots -----------------------------------------------------------------
  p1 <- mr_scatter_plot(r5, d)[[1]]
  ggsave(file.path(figdir, paste0(m, "_scatter.png")), plot = p1, width = 7, height = 6, dpi = 150)
  p2 <- mr_forest_plot(s)[[1]]
  ggsave(file.path(figdir, paste0(m, "_forest.png")), plot = p2, width = 7,
         height = max(5, 0.16 * nrow(d)), dpi = 150, limitsize = FALSE)
  p3 <- mr_funnel_plot(s)[[1]]
  ggsave(file.path(figdir, paste0(m, "_funnel.png")), plot = p3, width = 7, height = 6, dpi = 150)
  p4 <- mr_leaveoneout_plot(l)[[1]]
  ggsave(file.path(figdir, paste0(m, "_leaveoneout.png")), plot = p4, width = 7,
         height = max(5, 0.16 * nrow(d)), dpi = 150, limitsize = FALSE)

  pw <- function(x, tag) if (!is.null(x))
    write.csv(x, file.path(perdir, paste0(m, "__", tag, ".csv")), row.names = FALSE, na = "")
  pw(instr[[m]], "instruments"); pw(r5, "five_estimators"); pw(h, "heterogeneity")
  pw(e, "egger_intercept");      pw(s, "singlesnp");        pw(l, "leaveoneout")
  pw(fstat[[m]], "instrument_strength")
  if (!is.null(stg[[m]])) pw(stg[[m]], "steiger")
  pw(if (is.null(prs[[m]])) data.frame(model = m, status = "failed") else prs[[m]], "mrpresso")
  cat("done:", m, "\n")
}

# Combine whatever per-model files exist (this run plus any earlier resumed run).
combine <- function(tag, outfile) {
  fs <- list.files(perdir, pattern = paste0("__", tag, "\\.csv$"), full.names = TRUE)
  if (!length(fs)) return(invisible(NULL))
  x <- do.call(rbind, lapply(fs, function(p) read.csv(p, check.names = FALSE)))
  write.csv(x, file.path(resdir, outfile), row.names = FALSE, na = "")
}
combine("instruments",         "phase2_instruments.csv")
combine("five_estimators",     "phase2_five_estimators.csv")
combine("heterogeneity",       "phase2_heterogeneity.csv")
combine("egger_intercept",     "phase2_egger_intercept.csv")
combine("singlesnp",           "phase2_singlesnp.csv")
combine("leaveoneout",         "phase2_leaveoneout.csv")
combine("steiger",             "phase2_steiger.csv")
combine("mrpresso",            "phase2_mrpresso.csv")
combine("instrument_strength", "phase2_instrument_strength.csv")
w <- function(x, f) if (!is.null(x)) write.csv(x, file.path(resdir, f), row.names = FALSE, na = "")
note("(all)", "harmonised-all data",
     "not regenerable: the frozen package ships only *_harmonised_kept.csv; pre-filter rows are not recoverable")
w(acc(deviations), "phase2_deviations.csv")

writeLines(capture.output(sessionInfo()), file.path(outdir, "rerun_session_info.txt"))
cat("\nseed:", SEED, "\n")
