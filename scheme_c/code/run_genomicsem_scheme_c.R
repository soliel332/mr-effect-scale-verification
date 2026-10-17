suppressPackageStartupMessages(library(GenomicSEM))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
default_root <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
} else {
  getwd()
}
root <- Sys.getenv("SCHEME_C_ROOT", unset = default_root)
trait_names <- c("T1", "T2", "T3", "ADHD", "ASD", "NEURO", "ASTHMA")
traits <- file.path(root, "ldsc_sumstats", paste0(trait_names, ".sumstats.gz"))
stopifnot(all(file.exists(traits)))

ld_path <- Sys.getenv(
  "LD_SCORE_DIR",
  unset = file.path(root, "ldsc_reference_repo", "GRCh37", "results", "ldscores", "EUR", "")
)
result_dir <- file.path(root, "results")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

# All traits are analysed on their observed GWAS scale. Genetic correlations
# are invariant to positive rescaling; liability-scale prevalence assumptions
# are therefore unnecessary for the targeted correlation contrasts.
covstruc <- ldsc(
  traits = traits,
  sample.prev = rep(NA_real_, length(traits)),
  population.prev = rep(NA_real_, length(traits)),
  ld = ld_path,
  wld = ld_path,
  trait.names = trait_names,
  n.blocks = 200,
  ldsc.log = file.path(root, "ldsc_results", "GenomicSEM_scheme_C"),
  stand = TRUE
)

saveRDS(covstruc, file.path(result_dir, "Scheme_C_GenomicSEM_covariance_structure.rds"))
write.csv(covstruc$S, file.path(result_dir, "Scheme_C_genetic_covariance_matrix.csv"), quote = FALSE)
write.csv(covstruc$S_Stand, file.path(result_dir, "Scheme_C_genetic_correlation_matrix.csv"), quote = FALSE)
write.csv(covstruc$I, file.path(result_dir, "Scheme_C_cross_trait_intercept_matrix.csv"), quote = FALSE)

p <- length(trait_names)
se_vec <- sqrt(diag(covstruc$V_Stand))
se_mat <- matrix(NA_real_, p, p, dimnames = list(trait_names, trait_names))
se_mat[lower.tri(se_mat, diag = TRUE)] <- se_vec
se_mat[upper.tri(se_mat)] <- t(se_mat)[upper.tri(se_mat)]
write.csv(se_mat, file.path(result_dir, "Scheme_C_genetic_correlation_SE_matrix.csv"), quote = FALSE)

lower_index <- function(row, col, size = p) {
  if (row < col) {
    tmp <- row
    row <- col
    col <- tmp
  }
  before <- if (col == 1L) 0L else sum(size - seq_len(col - 1L) + 1L)
  before + (row - col + 1L)
}

s_vec <- gdata::lowerTriangle(covstruc$S_Stand, diag = TRUE)
pair_rows <- list()
k <- 1L
for (col in seq_len(p - 1L)) {
  for (row in (col + 1L):p) {
    idx <- lower_index(row, col)
    estimate <- s_vec[idx]
    se <- sqrt(covstruc$V_Stand[idx, idx])
    z <- estimate / se
    pair_rows[[k]] <- data.frame(
      trait_1 = trait_names[col],
      trait_2 = trait_names[row],
      rg = estimate,
      se = se,
      ci_low = estimate - qnorm(0.975) * se,
      ci_high = estimate + qnorm(0.975) * se,
      z = z,
      p_two_sided = 2 * pnorm(-abs(z)),
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
}
pairwise <- do.call(rbind, pair_rows)
write.csv(pairwise, file.path(result_dir, "Scheme_C_GenomicSEM_pairwise_rg.csv"), row.names = FALSE, quote = FALSE)

contrast <- function(tier_a, tier_b, comparator) {
  ia <- lower_index(match(tier_a, trait_names), match(comparator, trait_names))
  ib <- lower_index(match(tier_b, trait_names), match(comparator, trait_names))
  weights <- numeric(length(s_vec))
  weights[ia] <- 1
  weights[ib] <- -1
  estimate <- as.numeric(crossprod(weights, s_vec))
  variance <- as.numeric(t(weights) %*% covstruc$V_Stand %*% weights)
  se <- sqrt(variance)
  z <- estimate / se
  data.frame(
    comparator = comparator,
    contrast = paste0("rg(", tier_a, ",", comparator, ") - rg(", tier_b, ",", comparator, ")"),
    tier_a = tier_a,
    tier_b = tier_b,
    difference = estimate,
    se = se,
    ci_low = estimate - qnorm(0.975) * se,
    ci_high = estimate + qnorm(0.975) * se,
    z = z,
    p_two_sided = 2 * pnorm(-abs(z)),
    stringsAsFactors = FALSE
  )
}

tier_pairs <- list(c("T2", "T1"), c("T3", "T1"), c("T2", "T3"))
comparators <- c("ADHD", "ASD", "NEURO", "ASTHMA")
contrast_rows <- list()
k <- 1L
for (comparator in comparators) {
  for (tiers in tier_pairs) {
    contrast_rows[[k]] <- contrast(tiers[1], tiers[2], comparator)
    k <- k + 1L
  }
}
contrasts <- do.call(rbind, contrast_rows)
contrasts$p_bonferroni_12 <- pmin(1, contrasts$p_two_sided * nrow(contrasts))
contrasts$significant_bonferroni_12 <- contrasts$p_bonferroni_12 < 0.05
write.csv(contrasts, file.path(result_dir, "Scheme_C_formal_tier_rg_contrasts.csv"), row.names = FALSE, quote = TRUE)

capture.output(sessionInfo(), file = file.path(result_dir, "Scheme_C_GenomicSEM_sessionInfo.txt"))
cat("GenomicSEM Scheme C complete:", nrow(pairwise), "pairwise rg estimates and", nrow(contrasts), "formal contrasts\n")
