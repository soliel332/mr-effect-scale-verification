#!/usr/bin/env Rscript
# Publication exports of main-text Figures 1 and 3 from their genuine sources.
#
# Source: effect_scale_mr_analysis.R in the frozen repository -- section 2 (p_cal, lines 156-205)
# for Figure 1 and the Tier 2 correction block (p_correction, lines 382-394) for Figure 3. The data
# preparation and every geom, colour and mapping below are copied verbatim from that script. Only the
# export changes:
#   * width 174 mm (journal maximum; the frozen exports were 8.2 in and 7.2 in wide and would be
#     reduced on placement, taking axis text below 7 pt);
#   * Figure 1 strip labels wrapped onto two lines and x axis given fewer breaks, so that
#     nothing is clipped or overlapping at 174 mm (label wording unchanged);
#   * text sizes raised so that no lettering is below 8 pt at that width
#     (theme_classic base_size 10, which gives 8 pt axis text; facet strips 9 pt);
#   * formats: vector PDF (cairo) and EPS, 1,200-dpi LZW TIFF (line art) and a 600-dpi PNG for the DOCX.
# No value is read from any image. The data actually drawn are written out and checked against the
# source files (fig13_plotted_data_check.csv).
#
#   Rscript make_fig1_fig3.R <repository_root> <out_dir>
suppressPackageStartupMessages({ library(ggplot2) })
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[1] else "effect-scale-mr-qc"
out  <- if (length(args) >= 2) args[2] else "."
dir.create(out, showWarnings = FALSE, recursive = TRUE)
data_dir <- file.path(root, "data")
W_IN <- 174 / 25.4

# ---------------------------------------------------------------- Figure 1 (verbatim data prep)
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
  facet_wrap(~pair, scales = "free", ncol = 3,
             labeller = label_wrap_gen(width = 26)) +          # export: strip text wrapped, wording unchanged
  scale_x_continuous(n.breaks = 4) +                          # export: fewer breaks, no tick-label collision
  labs(x = "Linear-probability effect", y = "Log-odds effect") +
  theme_classic(base_size = 10) +                       # was 9
  theme(strip.text = element_text(size = 9))             # was 8

# ---------------------------------------------------------------- Figure 3 (verbatim)
correction <- read.csv(file.path(data_dir, "repository", "Tier2_scale_correction_source_data.csv"))
correction$analysis <- factor(correction$analysis, levels = correction$analysis)
p_correction <- ggplot(correction, aes(x = OR, y = analysis)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = low, xmax = high), width = 0.12,
                orientation = "y", colour = "#2A6F97") +
  geom_point(size = 2.5, colour = "#C65D3A") +
  labs(x = "Odds ratio for asthma", y = NULL) +
  theme_classic(base_size = 10)

# ---------------------------------------------------------------- plotted-data check
b1 <- ggplot_build(p_cal)
pts <- b1$data[[1]]
chk1 <- all.equal(sort(pts$x), sort(cal_points$beta.exposure)) && all.equal(sort(pts$y), sort(cal_points$beta.outcome))
b3 <- ggplot_build(p_correction)
eb <- b3$data[[2]]; pt <- b3$data[[3]]
chk3 <- isTRUE(all.equal(sort(pt$x), sort(correction$OR))) &&
        isTRUE(all.equal(sort(eb$xmin), sort(correction$low))) &&
        isTRUE(all.equal(sort(eb$xmax), sort(correction$high)))
# smallest rendered text size (pt) from the complete themes
min_pt <- function(p) {
  th <- ggplot2:::plot_theme(p)
  el <- c("axis.text.x", "axis.text.y", "axis.title.x", "axis.title.y", "strip.text")
  s <- sapply(el, function(e) { v <- calc_element(e, th); if (inherits(v, "element_text")) v$size else NA })
  min(s, na.rm = TRUE)
}
write.csv(data.frame(figure = c("Figure 1", "Figure 3"),
                     points_drawn = c(nrow(pts), nrow(pt)),
                     source_rows = c(nrow(cal_points), nrow(correction)),
                     plotted_values_equal_source = c(isTRUE(chk1), isTRUE(chk3)),
                     min_text_pt = c(min_pt(p_cal), min_pt(p_correction)),
                     width_mm = 174),
          file.path(out, "fig13_plotted_data_check.csv"), row.names = FALSE)
write.csv(cal_points[, c("pair", "SNP", "beta.exposure", "beta.outcome")],
          file.path(out, "fig1_plotted_values.csv"), row.names = FALSE)
write.csv(correction, file.path(out, "fig3_plotted_values.csv"), row.names = FALSE)

# ---------------------------------------------------------------- exports
save_all <- function(p, stem, h_in) {
  ggsave(file.path(out, paste0(stem, ".pdf")), p, width = W_IN, height = h_in, device = cairo_pdf)
  ggsave(file.path(out, paste0(stem, ".eps")), p, width = W_IN, height = h_in, device = cairo_ps,
         fallback_resolution = 1200)
  ggsave(file.path(out, paste0(stem, ".tiff")), p, width = W_IN, height = h_in, dpi = 1200,
         compression = "lzw", bg = "white")
  ggsave(file.path(out, paste0(stem, ".png")), p, width = W_IN, height = h_in, dpi = 600, bg = "white")
}
save_all(p_cal, "fig1_calibration", 3.3)
save_all(p_correction, "fig3_case", 2.6)
cat("Figure 1 plotted = source:", isTRUE(chk1), "| Figure 3 plotted = source:", isTRUE(chk3),
    "| min text pt:", min_pt(p_cal), min_pt(p_correction), "\n")
writeLines(capture.output(sessionInfo()), file.path(out, "fig13_sessionInfo.txt"))
