#!/usr/bin/env Rscript
# Online Resource 1, Fig. S1 (from release v1.2): scatter and leave-one-out plots of the original and the
# arithmetically rescaled Tier 2 analysis. The script also writes the same four-panel plots for the six
# analyses of Table 4 (figS2-figS7), which are not part of the submitted manuscript.
# No new analysis. Estimator lines are drawn from the committed results files
# (results/MR_rerun_five_estimators.csv, results/MR_rerun_egger_intercept.csv,
#  results/invariance_estimator_unit_test.csv, results/invariance_egger_unit_test.csv);
# single-variant Wald ratios and leave-one-out IVW estimates are computed from the harmonised files in
# data/mr/ and the leave-one-out values are checked against case_study/leave_one_out/.
#   Rscript figures/code/make_mr_plots.R . figures/online_resource_1     (run from the repository root)
# Checked by tests/validate_figures.py.
suppressPackageStartupMessages({ library(ggplot2); library(cowplot) })
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[1] else "."
out  <- if (length(args) >= 2) args[2] else "."
dir.create(out, showWarnings = FALSE, recursive = TRUE)
W_IN <- 174 / 25.4
LAMBDA <- 1 / (0.3450365567 * (1 - 0.3450365567))     # as in effect_scale_mr_analysis.R, section 4
COLS <- c("IVW (multiplicative random effects)" = "#C65D3A", "MR-Egger" = "#2A6F97",
          "Weighted median" = "#4D9078", "Weighted mode" = "#8E6C8A", "Simple mode" = "#B8A04A")
MNAME <- c("Inverse variance weighted (multiplicative random effects)" = "IVW (multiplicative random effects)",
           "MR Egger" = "MR-Egger", "Weighted median" = "Weighted median",
           "Weighted mode" = "Weighted mode", "Simple mode" = "Simple mode")
th <- theme_classic(base_size = 7.5) +
  theme(axis.text = element_text(colour = "black", size = 6.5), plot.title = element_text(size = 8, face = "bold"),
        legend.key.height = unit(3, "mm"), legend.text = element_text(size = 6.5), legend.title = element_blank())

read_h <- function(f) {
  d <- read.csv(file.path(root, "data", "mr", f), check.names = FALSE)
  d <- d[d$mr_keep %in% TRUE, ]
  s <- sign(d$beta.exposure); d$bx <- d$beta.exposure * s; d$by <- d$beta.outcome * s   # orient to positive exposure effect
  d$sx <- d$se.exposure; d$sy <- d$se.outcome; d
}
ivw <- function(bx, by, sy) {
  w <- 1 / sy^2; b <- sum(w * bx * by) / sum(w * bx^2); k <- length(bx)
  sig <- sqrt(sum(w * (by - b * bx)^2) / (k - 1)); se <- sqrt(1 / sum(w * bx^2)) * max(1, sig)
  c(b = b, se = se, p = 2 * pnorm(-abs(b / se)))
}
loo <- function(d) do.call(rbind, lapply(seq_len(nrow(d)), function(i) {
  r <- ivw(d$bx[-i], d$by[-i], d$sy[-i]); data.frame(SNP = d$SNP[i], b = r[["b"]], se = r[["se"]], p = r[["p"]]) }))

panels <- function(d, est, egg_int, xlab_unit, tag) {
  d$wr <- d$by / d$bx; d$wr_se <- d$sy / d$bx
  lines <- data.frame(method = MNAME[est$method], slope = est$b,
                      intercept = ifelse(est$method == "MR Egger", egg_int, 0))
  lines$method <- factor(lines$method, levels = names(COLS))
  ivw_b <- est$b[est$method == "Inverse variance weighted (multiplicative random effects)"]
  ivw_se <- est$se[est$method == "Inverse variance weighted (multiplicative random effects)"]
  egg_b <- est$b[est$method == "MR Egger"]
  pa <- ggplot(d, aes(bx, by)) +
    geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
    geom_errorbar(aes(ymin = by - 1.96 * sy, ymax = by + 1.96 * sy), colour = "grey65", linewidth = 0.25, width = 0) +
    geom_errorbarh(aes(xmin = bx - 1.96 * sx, xmax = bx + 1.96 * sx), colour = "grey65", linewidth = 0.25, height = 0) +
    geom_point(size = 0.9) +
    geom_abline(data = lines, aes(slope = slope, intercept = intercept, colour = method), linewidth = 0.5) +
    scale_colour_manual(values = COLS, drop = FALSE) +
    labs(x = paste0("Variant effect on exposure (", xlab_unit, ")"), y = "Variant effect on asthma (log odds)",
         title = "a  Scatter plot") + expand_limits(x = 0) + th + theme(legend.position = "none")
  pb <- ggplot(d, aes(wr, 1 / wr_se)) + geom_point(size = 0.9) +
    geom_vline(xintercept = ivw_b, colour = COLS[1], linewidth = 0.5) +
    geom_vline(xintercept = egg_b, colour = COLS[2], linewidth = 0.5, linetype = 2) +
    labs(x = "Single-variant Wald ratio", y = "1 / SE of Wald ratio", title = "b  Funnel plot") + th
  ord <- order(d$wr); f <- data.frame(SNP = c(d$SNP[ord], "IVW"), b = c(d$wr[ord], ivw_b),
                                      se = c(d$wr_se[ord], ivw_se), pooled = c(rep(FALSE, nrow(d)), TRUE))
  f$SNP <- factor(f$SNP, levels = rev(f$SNP))
  pc <- ggplot(f, aes(b, SNP, colour = pooled)) + geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_errorbarh(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se), height = 0, linewidth = 0.3) +
    geom_point(aes(shape = pooled), size = 1) + scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = COLS[[1]]), guide = "none") +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18), guide = "none") +
    labs(x = "Wald ratio (95% CI)", y = NULL, title = "c  Single-variant estimates") + th +
    theme(axis.text.y = element_text(size = 5.5))
  l <- loo(d); l <- l[order(l$b), ]; l <- rbind(l, data.frame(SNP = "All variants", b = ivw_b, se = ivw_se, p = NA))
  l$all <- l$SNP == "All variants"; l$SNP <- factor(l$SNP, levels = rev(l$SNP))
  pd <- ggplot(l, aes(b, SNP, colour = all)) + geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_errorbarh(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se), height = 0, linewidth = 0.3) +
    geom_point(aes(shape = all), size = 1) + scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = COLS[[1]]), guide = "none") +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18), guide = "none") +
    labs(x = "IVW estimate with variant omitted (95% CI)", y = NULL, title = "d  Leave-one-out") + th +
    theme(axis.text.y = element_text(size = 5.5))
  leg <- get_legend(pa + theme(legend.position = "bottom", legend.direction = "horizontal") +
                    guides(colour = guide_legend(nrow = 1, override.aes = list(linewidth = 0.8))))
  list(pa = pa, pb = pb, pc = pc, pd = pd, loo = l, d = d, leg = leg)
}
save_fig <- function(p, stem, h_mm) {
  ggsave(file.path(out, paste0(stem, ".pdf")), p, width = W_IN, height = h_mm / 25.4, device = cairo_pdf)
  ggsave(file.path(out, paste0(stem, ".png")), p, width = W_IN, height = h_mm / 25.4, dpi = 600, bg = "white")
}

five <- read.csv(file.path(root, "results", "MR_rerun_five_estimators.csv"))
egg <- read.csv(file.path(root, "results", "MR_rerun_egger_intercept.csv"))
an <- data.frame(
  file = c("Tier_2_help_seeking_corrected__UK_Biobank_asthma", "Tier_2_help_seeking_corrected__FinnGen_asthma",
           "Tier_1_affective_instability__UK_Biobank_asthma", "Tier_1_affective_instability__FinnGen_asthma",
           "Tier_3_major_depression_meta_analysis__UK_Biobank_asthma", "Tier_3_major_depression_meta_analysis__FinnGen_asthma"),
  model = c("Help-seeking (corrected) | UK Biobank asthma", "Help-seeking (corrected) | FinnGen asthma",
            "Affective instability | UK Biobank asthma", "Affective instability | FinnGen asthma",
            "Major depression meta-analysis | UK Biobank asthma", "Major depression meta-analysis | FinnGen asthma"),
  fig = paste0("figS", 2:7), stringsAsFactors = FALSE)
chk <- list()
for (i in seq_len(nrow(an))) {
  d <- read_h(paste0(an$file[i], "_harmonised_kept.csv"))
  est <- five[five$model == an$model[i], ]; stopifnot(nrow(est) == 5, all(est$nsnp == nrow(d)))
  P <- panels(d, est, egg$egger_intercept[egg$model == an$model[i]], "log odds", an$fig[i])
  ref <- read.csv(file.path(root, "case_study", "leave_one_out", paste0(an$file[i], "__leaveoneout.csv")))
  ref <- ref[ref$SNP != "All", ]; m <- merge(P$loo[!P$loo$all, ], ref, by = "SNP")
  iv <- ivw(d$bx, d$by, d$sy); ivr <- est[est$method == "Inverse variance weighted (multiplicative random effects)", ]
  chk[[i]] <- data.frame(figure = an$fig[i], analysis = an$model[i], K = nrow(d),
                         ivw_matches_results = isTRUE(all.equal(unname(iv[c("b","se")]), c(ivr$b, ivr$se), tolerance = 1e-8)),
                         loo_rows_matched = nrow(m), loo_max_abs_diff_b = max(abs(m$b.x - m$b.y)), loo_max_abs_diff_se = max(abs(m$se.x - m$se.y)))
  h <- 62 + 7 + 2.9 * (nrow(d) + 1) + 14
  g <- plot_grid(plot_grid(P$pa, P$pb, ncol = 2), P$leg, plot_grid(P$pc, P$pd, ncol = 2), ncol = 1,
                 rel_heights = c(62, 7, h - 62))
  save_fig(g, an$fig[i], h)
}
# ---- Figure S1: invariance, original mixed-scale vs arithmetically rescaled (Tier 2 x UK Biobank)
sup <- read_h("Tier_2_superseded_mixed_scale__UK_Biobank_asthma_harmonised_kept.csv")
resc <- sup; resc$bx <- resc$bx * LAMBDA; resc$sx <- resc$sx * LAMBDA
inv <- read.csv(file.path(root, "results", "invariance_estimator_unit_test.csv"))
ie <- read.csv(file.path(root, "results", "invariance_egger_unit_test.csv"))
e_o <- data.frame(method = inv$method, b = inv$b_original, se = inv$se_original)
e_r <- data.frame(method = inv$method, b = inv$b_rescaled, se = inv$se_rescaled)
int <- ie$original[ie$statistic == "Egger intercept"]
Po <- panels(sup, e_o, int, "linear probability, read as log odds", "o")
Pr <- panels(resc, e_r, int, "log odds, rescaled by \u03bb", "r")
lab <- function(e) { r <- e[e$method == "Inverse variance weighted (multiplicative random effects)", ]
  sprintf("IVW OR %.2f, P = %.2g", exp(r$b), 2 * pnorm(-abs(r$b / r$se))) }
Po$pa <- Po$pa + ggtitle("a  Original mixed-scale analysis") + annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.6, size = 2.3, label = lab(e_o))
Pr$pa <- Pr$pa + ggtitle("b  Arithmetically rescaled") + annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.6, size = 2.3, label = lab(e_r))
Po$pd <- Po$pd + ggtitle("c  Leave-one-out, original"); Pr$pd <- Pr$pd + ggtitle("d  Leave-one-out, rescaled")
lo <- Po$loo[!Po$loo$all, ]; lr <- Pr$loo[!Pr$loo$all, ]; mm <- merge(lo, lr, by = "SNP")
chk[[length(chk) + 1]] <- data.frame(figure = "figS1", analysis = "Tier 2 original vs rescaled, UK Biobank", K = nrow(sup),
  ivw_matches_results = isTRUE(all.equal(unname(ivw(sup$bx, sup$by, sup$sy)[c("b","se")]), c(e_o$b[1], e_o$se[1]), tolerance = 1e-8)) &&
                        isTRUE(all.equal(unname(ivw(resc$bx, resc$by, resc$sy)[c("b","se")]), c(e_r$b[1], e_r$se[1]), tolerance = 1e-8)),
  loo_rows_matched = nrow(mm), loo_max_abs_diff_b = max(abs(mm$b.x / mm$b.y - LAMBDA)), loo_max_abs_diff_se = max(abs(mm$p.x - mm$p.y)))
h1 <- 62 + 7 + 2.9 * (nrow(sup) + 1) + 14
g1 <- plot_grid(plot_grid(Po$pa, Pr$pa, ncol = 2), Po$leg, plot_grid(Po$pd, Pr$pd, ncol = 2), ncol = 1, rel_heights = c(62, 7, h1 - 62))
save_fig(g1, "figS1", h1)
C <- do.call(rbind, chk); write.csv(C, file.path(out, "mr_plots_check.csv"), row.names = FALSE); print(C)
