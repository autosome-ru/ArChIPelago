# =============================================================================
# Figure_S5_saturation.R  --  Fig. S5: performance vs the
# number of PWMs used, drawn in R in the style of the main and supplementary
# figures of the paper (Arial, theme_minimal(base_size = 21, base_family = "Arial"), Set2 palette,
# dotted zero lines, cairo_pdf output).
#
# Input : Figures/source_data/Figure_S5_source_data.csv   (per-TF curves)
#         Figures/source_data/Figure_S5_median_curves.csv (median curves)
#         (written by scripts/make_figures_S2_S5_S6.py)
# Usage : Rscript Figure_S5_saturation.R [source_data_dir] [output_dir]
#         (defaults: ../Figures/source_data and ../Figures/panels)
# Mouse test set = mouse chr1/8/19; reference = best single monoPWM.
# =============================================================================

suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2); library(patchwork)})

script_dir <- local({
  f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})
args     <- commandArgs(trailingOnly = TRUE)
pkg_dir <- normalizePath(file.path(script_dir, ".."))   # Manuscript_analysis/
src_dir <- if (length(args) >= 1) args[1] else file.path(pkg_dir, "Figures", "source_data")
out_dir <- if (length(args) >= 2) args[2] else file.path(pkg_dir, "Figures", "panels")
per_tf_path <- file.path(src_dir, "Figure_S5_source_data.csv")
median_path <- file.path(src_dir, "Figure_S5_median_curves.csv")
for (p in c(per_tf_path, median_path)) if (!file.exists(p)) stop("Input not found: ", p)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat("Source data :", src_dir, "\nOutput dir  :", out_dir, "\n")
if (!interactive()) pdf(NULL)

# ---- house style ------------------------------------------------------------
theme_set(theme_get() + theme(text = element_text(family = "Arial")))
SET2_GREEN  <- "#66c2a5"   # Set2[1], as in Figure 2 (monoPWM reference line)
SET2_ORANGE <- "#fc8d62"   # Set2[2], as in Figure 2 (diPWM reference line)
GREY        <- "gray80"    # as in Figure 4 (per-TF points)

LAB_RANDOM <- "median over TFs, random subsets of k PWMs"
LAB_TOPK   <- "median over TFs, k best PWMs by training auROC"
LAB_TF     <- "single TF (mean over replicates), ending at its own number of PWMs P"
COLS <- setNames(c(SET2_GREEN, SET2_ORANGE, GREY), c(LAB_RANDOM, LAB_TOPK, LAB_TF))
LTYS <- setNames(c("solid", "dashed", "solid"),    c(LAB_RANDOM, LAB_TOPK, LAB_TF))

K_GRID  <- c(1, 2, 4, 8, 16, 32, 64, 128)
X_ALL   <- 9              # x position (log2 scale) of the "all (P)" tick
METRICS <- list(c("delta_auroc_H", "auROC, human test set (vs best monoPWM)", "ΔauROC"),
                c("delta_auprc_H", "auPRC, human test set (vs best monoPWM)", "ΔauPRC"),
                c("delta_auroc_M", "auROC, mouse test set (vs best monoPWM)", "ΔauROC"),
                c("delta_auprc_M", "auPRC, mouse test set (vs best monoPWM)", "ΔauPRC"))

# ---- data -------------------------------------------------------------------
per_tf <- read.csv(per_tf_path, stringsAsFactors = FALSE) %>%
  filter(design == "random") %>%
  filter(k_requested == "P" | suppressWarnings(as.numeric(k_requested)) < P) %>%
  mutate(x = log2(k))

med <- read.csv(median_path, stringsAsFactors = FALSE, check.names = FALSE) %>%
  rename_with(~ sub("^d_", "delta_", .x)) %>%
  mutate(x = ifelse(k == "P", X_ALL, log2(suppressWarnings(as.numeric(k)))),
         series = ifelse(design == "random", LAB_RANDOM, LAB_TOPK))

panel <- function(metric, subtitle, ylab, letter, show_x) {
  p <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_line(data = per_tf, aes(x = x, y = .data[[metric]], group = TF, color = LAB_TF),
              linewidth = 0.4, alpha = 0.9) +
    geom_line(data = med, aes(x = x, y = .data[[metric]], color = series, linetype = series),
              linewidth = 1.3) +
    geom_point(data = filter(med, k == "P"), aes(x = x, y = .data[[metric]], color = series),
               size = 3.5, show.legend = FALSE) +
    scale_color_manual(values = COLS, breaks = c(LAB_RANDOM, LAB_TOPK, LAB_TF), name = NULL) +
    scale_linetype_manual(values = LTYS, breaks = c(LAB_RANDOM, LAB_TOPK, LAB_TF), name = NULL) +
    scale_x_continuous(breaks = c(log2(K_GRID), X_ALL),
                       labels = c(as.character(K_GRID), "all\n(P)")) +
    labs(title = letter, subtitle = subtitle, y = ylab,
         x = if (show_x) "number of PWMs used (k)" else NULL) +
    theme_minimal(base_size = 21, base_family = "Arial") +
    theme(legend.position = "none",
          axis.text = element_text(color = "black"),
          plot.title = element_text(face = "bold", hjust = 0, size = 25),
          plot.subtitle = element_text(size = 16, hjust = 0),
          plot.margin = margin(6, 10, 6, 10))
  if (!show_x) p <- p + theme(axis.text.x = element_blank())
  p
}

panels <- Map(function(m, letter, show_x) panel(m[1], m[2], m[3], letter, show_x),
              METRICS, c("A", "B", "C", "D"), c(FALSE, FALSE, TRUE, TRUE))

legend_df <- data.frame(x = 0, y = 0,
                        series = factor(c(LAB_RANDOM, LAB_TOPK, LAB_TF),
                                        levels = c(LAB_RANDOM, LAB_TOPK, LAB_TF)))
legend_plot <- ggplot(legend_df, aes(x = x, y = y, color = series, linetype = series)) +
  geom_line(aes(linewidth = series)) +
  scale_color_manual(values = COLS, name = NULL) +
  scale_linetype_manual(values = LTYS, name = NULL) +
  scale_linewidth_manual(values = setNames(c(1.3, 1.3, 0.6), names(COLS)), name = NULL) +
  theme_minimal(base_size = 21, base_family = "Arial") +
  theme(legend.position = "bottom", legend.direction = "vertical",
        legend.text = element_text(size = 16))

combined <- (panels[[1]] | panels[[2]]) / (panels[[3]] | panels[[4]]) /
  cowplot::get_legend(legend_plot) + plot_layout(heights = c(1, 1, 0.22))

pdf_out <- file.path(out_dir, "Figure_S5_saturation.pdf")
ggsave(pdf_out, combined, device = grDevices::cairo_pdf,
       width = 320, height = 260, units = "mm", dpi = 600)

# ---- stdout diagnostics -----------------------------------------------------
cat("\n==== Fig. S5: median delta vs k (random subsets) ====\n")
print(as.data.frame(med %>% filter(design == "random") %>%
        select(k, delta_auroc_H, delta_auprc_H, delta_auroc_M, delta_auprc_M)), digits = 3, row.names = FALSE)
full <- med %>% filter(design == "random", k == "P")
half <- med %>% filter(design == "random", k %in% as.character(K_GRID)) %>%
  mutate(frac_auroc_H = delta_auroc_H / full$delta_auroc_H,
         frac_auprc_H = delta_auprc_H / full$delta_auprc_H)
cat("fraction of the full gain reached (human test set):\n")
print(as.data.frame(half %>% select(k, frac_auroc_H, frac_auprc_H)), digits = 2, row.names = FALSE)
cat("TFs plotted:", length(unique(per_tf$TF)), "\n")
cat("Written:", pdf_out, "\n")
