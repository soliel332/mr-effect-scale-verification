#!/usr/bin/env Rscript
# MR-PRESSO for the six corrected analyses, under the PINNED environment.
#
# mr_presso() depends on neither TwoSampleMR nor RadialMR -- it takes a plain data
# frame -- so this component runs under the exact pinned stack even where the full
# TwoSampleMR dependency tree is unavailable:
#
#   R 4.4.1
#   MRPRESSO 1.0, source commit 3e3c92d7eda6dce0d1d66077373ec0f7ff4f7e87
#
# Parameters follow the frozen repository: NbDistribution = 10000,
# SignifThreshold = 0.05, seed 20260909, OUTLIERtest and DISTORTIONtest enabled.
#
# PARSER (this is the corrected behaviour)
#   MR-PRESSO returns the CHARACTER sentinel "No significant outliers" in
#   `Outliers Indices` when none are detected. Treating that as an index vector
#   makes length() return 1 and fabricates a phantom outlier. Only a numeric,
#   non-empty vector is a real index set. When there is none, the outlier-corrected
#   estimate and the distortion test are NOT APPLICABLE -- MR-PRESSO does not define
#   them -- which is different from their being missing or having failed.
#
# CONCURRENCY: each model writes only its own prefixed files, each written to a
# .tmp path and then renamed (atomic within a filesystem), finishing with a .done
# marker written last. A model whose marker exists is skipped, so an interrupted
# run resumes without duplicating or half-writing rows.

suppressPackageStartupMessages(library(MRPRESSO))
options(stringsAsFactors = FALSE)
SEED <- 20260909

args   <- commandArgs(trailingOnly = TRUE)
repo   <- if (length(args) >= 1) args[1] else "repo/effect-scale-mr-qc"
outdir <- if (length(args) >= 2) args[2] else "staging/phase2_pinned"
ONLY   <- if (length(args) >= 3 && nzchar(args[3])) args[3] else NA

perdir <- file.path(outdir, "per_model")
dir.create(perdir, recursive = TRUE, showWarnings = FALSE)

ALL_MODELS <- c(
  "Tier_1_affective_instability__UK_Biobank_asthma",
  "Tier_1_affective_instability__FinnGen_asthma",
  "Tier_2_help_seeking_corrected__UK_Biobank_asthma",
  "Tier_2_help_seeking_corrected__FinnGen_asthma",
  "Tier_3_major_depression_meta_analysis__UK_Biobank_asthma",
  "Tier_3_major_depression_meta_analysis__FinnGen_asthma"
)
MODELS <- if (is.na(ONLY)) ALL_MODELS else ALL_MODELS[ALL_MODELS == ONLY]

atomic_csv <- function(x, path) {
  tmp <- paste0(path, ".tmp"); write.csv(x, tmp, row.names = FALSE, na = ""); file.rename(tmp, path)
}
atomic_lines <- function(txt, path) {
  tmp <- paste0(path, ".tmp"); writeLines(txt, tmp); file.rename(tmp, path)
}

for (m in MODELS) {
  marker <- file.path(perdir, paste0(m, ".done"))
  if (file.exists(marker)) { cat("skip (marker present):", m, "\n"); next }

  d <- read.csv(file.path(repo, "data", "mr", paste0(m, "_harmonised_kept.csv")), check.names = FALSE)
  set.seed(SEED)
  pr <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                  SdOutcome = "se.outcome", SdExposure = "se.exposure",
                  OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                  data = as.data.frame(d), NbDistribution = 10000,
                  SignifThreshold = 0.05, seed = SEED)

  mt <- pr$`Main MR results`
  gt <- pr$`MR-PRESSO results`$`Global Test`
  dt <- pr$`MR-PRESSO results`$`Distortion Test`
  oi <- if (is.null(dt)) NULL else dt$`Outliers Indices`

  has_out <- !is.null(oi) && is.numeric(oi) && length(oi) > 0
  idx     <- if (has_out) oi else integer(0)
  ocrow   <- if (nrow(mt) > 1) 2 else NA

  atomic_csv(data.frame(
    model = m,
    outlier_status = if (has_out) "outliers detected" else "no significant outliers",
    n_outliers = length(idx),
    outlier_indices = if (has_out) paste(idx, collapse = ";") else "",
    outlier_SNPs    = if (has_out) paste(d$SNP[idx], collapse = ";") else "",
    raw_estimate = mt$`Causal Estimate`[1], raw_se = mt$Sd[1], raw_p = mt$`P-value`[1],
    corrected_estimate = if (has_out && !is.na(ocrow)) mt$`Causal Estimate`[ocrow] else NA,
    corrected_se       = if (has_out && !is.na(ocrow)) mt$Sd[ocrow] else NA,
    corrected_p        = if (has_out && !is.na(ocrow)) mt$`P-value`[ocrow] else NA,
    corrected_applicable  = has_out,
    global_RSSobs = gt$RSSobs, global_p = as.character(gt$Pvalue),
    distortion_coef = if (has_out) dt$`Distortion Coefficient` else NA,
    distortion_p    = if (has_out) dt$Pvalue else NA,
    distortion_applicable = has_out,
    NbDistribution = 10000, SignifThreshold = 0.05, seed = SEED,
    R_version = as.character(getRversion()),
    MRPRESSO_version = as.character(packageVersion("MRPRESSO"))),
    file.path(perdir, paste0(m, "__mrpresso.csv")))

  atomic_lines(capture.output(sessionInfo()), file.path(perdir, paste0(m, "__sessionInfo.txt")))
  atomic_lines(format(Sys.time(), "%Y-%m-%dT%H:%M:%S"), marker)
  cat("done:", m, "\n")
}

if (is.na(ONLY)) {
  missing <- ALL_MODELS[!file.exists(file.path(perdir, paste0(ALL_MODELS, ".done")))]
  if (length(missing)) {
    cat("NOT COMBINING -- markers absent for:\n"); cat(paste0("  ", missing, collapse = "\n"), "\n")
    quit(status = 1)
  }
  fs <- file.path(perdir, paste0(ALL_MODELS, "__mrpresso.csv"))
  atomic_csv(do.call(rbind, lapply(fs, read.csv, check.names = FALSE)),
             file.path(outdir, "phase2_mrpresso_pinned.csv"))
  atomic_lines(capture.output(sessionInfo()), file.path(outdir, "mrpresso_pinned_sessionInfo.txt"))
  cat("combined.\n")
}
