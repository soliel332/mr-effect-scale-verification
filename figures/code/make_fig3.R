#!/usr/bin/env Rscript
# Main-text Figure 3 (from release v1.2): the four case-study quantities, including the post hoc bridge comparison.
# Replaces the three-row Figure 3 of make_fig1_fig3.R; the v1.1 files are kept in figures/superseded_v1.1/.
# Values are read from case_study/bridge/bridge_estimates.csv (written by bridge_analysis.py and
# verified in the pinned R environment by verify_bridge_pinned.R); rows 1, 2 and 4 are checked
# against data/repository/Tier2_scale_correction_source_data.csv, the source of the earlier figure.
#   Rscript figures/code/make_fig3.R . figures/submitted     (run from the repository root)
suppressPackageStartupMessages({ library(ggplot2) })
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[1] else "."
out  <- if (length(args) >= 2) args[2] else "."
dir.create(out, showWarnings = FALSE, recursive = TRUE)
b <- read.csv(file.path(root, "case_study", "bridge", "bridge_estimates.csv"))
stopifnot(nrow(b) == 4)
old <- read.csv(file.path(root, "data", "repository", "Tier2_scale_correction_source_data.csv"))
stopifnot(isTRUE(all.equal(sort(round(old$OR, 6)), sort(round(b$OR[c(1, 2, 4)], 6)))))
lab <- c("Original mixed-scale analysis\n(35 variants)",
         "Arithmetic rescaling\n(35 variants)",
         "Bridge: replacement-release\nexposure coefficients (35 variants)",
         "Replacement release, own\ninstrument set (26 variants)")
d <- data.frame(analysis = factor(lab, levels = rev(lab)), OR = b$OR, low = b$CI_low, high = b$CI_high,
                post_hoc = c(FALSE, FALSE, TRUE, FALSE))
d$txt <- sprintf("%.2f (%.2f to %.2f)", d$OR, d$low, d$high)
p <- ggplot(d, aes(x = OR, y = analysis)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = low, xmax = high), width = 0.15, orientation = "y", colour = "#2A6F97") +
  geom_point(aes(shape = post_hoc), size = 2.6, colour = "#C65D3A", fill = "white", stroke = 1) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 21), guide = "none") +
  geom_text(aes(x = 6.6, label = txt), hjust = 0, size = 8.5 / .pt) +
  annotate("text", x = 6.6, y = 4.55, label = "OR (95% CI)", hjust = 0, size = 8.5 / .pt, fontface = "bold") +
  coord_cartesian(xlim = c(0.8, 6), clip = "off") +
  scale_x_log10(breaks = c(0.8, 1, 1.5, 2, 3, 5),
                labels = c("0.8", "1", "1.5", "2", "3", "5")) +
  labs(x = "Odds ratio for asthma (log scale)", y = NULL) +
  theme_classic(base_size = 10) +
  theme(axis.text.y = element_text(size = 8.5, lineheight = 0.95, colour = "black"),
        axis.text.x = element_text(size = 8.5, colour = "black"),
        plot.margin = margin(14, 92, 4, 4))
W_IN <- 174 / 25.4; H_IN <- 3.0
ggsave(file.path(out, "fig3_case.pdf"), p, width = W_IN, height = H_IN, device = cairo_pdf)
ggsave(file.path(out, "fig3_case.eps"), p, width = W_IN, height = H_IN, device = cairo_ps, fallback_resolution = 1200)
ggsave(file.path(out, "fig3_case.tiff"), p, width = W_IN, height = H_IN, dpi = 1200, compression = "lzw", bg = "white")
ggsave(file.path(out, "fig3_case.png"), p, width = W_IN, height = H_IN, dpi = 600, bg = "white")
g <- ggplot_build(p); pt <- g$data[[3]]
ok <- isTRUE(all.equal(sort(10^pt$x), sort(d$OR))) && isTRUE(all.equal(sort(10^g$data[[2]]$xmin), sort(d$low))) && isTRUE(all.equal(sort(10^g$data[[2]]$xmax), sort(d$high)))
write.csv(d[, c("analysis", "OR", "low", "high", "post_hoc")], file.path(out, "fig3_plotted_values.csv"), row.names = FALSE)
cat("plotted = source:", ok, "\n"); print(d[, c("txt", "post_hoc")])
