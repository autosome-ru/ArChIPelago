# Figure_S6_cross_species.R - Fig. S6: cross-species
# transfer and the mouse-trained control, drawn in R in the style of the main
# and supplementary figures of the paper (Arial, theme_minimal(base_size = 21, base_family = "Arial"),
# Set2 palette, dotted zero lines, ggrepel labels, cairo_pdf output).
#
# Input : Figures/source_data/Figure_S6_source_data.csv (extract of Sup. Table 5,
#         written by make_figure_source_data.py)
# Usage : Rscript Figure_S6.R [source_data_dir] [output_dir]
#         (defaults: ../Figures/source_data and ../Figures/panels)
# Mouse test set = mouse chr1/8/19; baseline = best single monoPWM. The TFs
# below the baseline after transfer are drawn last and labelled; jitter in panel D is seeded.

suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2)
                                library(ggrepel); library(patchwork)})

script_dir <- local({
  f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})
args    <- commandArgs(trailingOnly = TRUE)
source(file.path(script_dir, "..", "..", "style", "archi_style.R"))
pkg_dir <- ARCHI_PKG   # released package, read only
src_dir <- if (length(args) >= 1) args[1] else file.path(pkg_dir, "Figures", "source_data")
out_dir <- if (length(args) >= 2) args[2] else script_dir
src_path <- file.path(src_dir, "Figure_S6_source_data.csv")
if (!file.exists(src_path)) stop("Input not found: ", src_path)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat("Source data :", src_path, "\nOutput dir  :", out_dir, "\n")
if (!interactive()) pdf(NULL)

# style
theme_set(theme_get() + theme(text = element_text(family = "Arial")))
# Restyled: palette from ../archi_style.R. Human-trained Random Forest points = RF green, TFs below the
# baseline after transfer = red (the colour of the baseline, best single monoPWM); delta auROC / auPRC = dark / light grey.
# Opaque marks (no alpha). Data handling and layout unchanged.
SET2_GREEN  <- "#5A5A5A"
SET2_ORANGE <- "white"   # below-baseline TFs: black open diamond (red is reserved for the best monoPWM)
SET2_BLUE   <- RF_COL

LAB_OK   <- "all other TFs"
LAB_FAIL <- "below the baseline after transfer"
FILLS <- setNames(c(SET2_BLUE, SET2_ORANGE), c(LAB_OK, LAB_FAIL))
SHAPES <- setNames(c(21, 23), c(LAB_OK, LAB_FAIL))
MET_COLS <- setNames(c("#555555", "#FFFFFF"), c("ΔauROC", "ΔauPRC"))   # Set2[1], Set2[4]; blue is reserved for "all other TFs" in A-C

base_theme <- function() {
  theme_minimal(base_size = 21, base_family = "Arial") +
    theme(panel.grid = element_blank(), text = element_text(colour = "black"), axis.text = element_text(color = "black"),
          axis.line = element_line(colour = "black", linewidth = 0.6), axis.ticks = element_line(colour = "black"),
          plot.background = element_rect(fill = "white", colour = NA),
          plot.title = element_text(face = "bold", hjust = 0, size = 25),
          plot.subtitle = element_text(size = 16, hjust = 0),
          plot.margin = margin(6, 12, 6, 12),
          legend.text = element_text(size = 16))
}

# data
d <- read.csv(src_path, stringsAsFactors = FALSE, check.names = FALSE)
col <- function(nm) { if (!nm %in% names(d)) stop("CHECK FAILED Fig. S6: column missing in the source data: ", nm); d[[nm]] }
# The human-trained models on the mouse test set as the mean of Random Forest seeds 0 and 1
HM <- function(x) paste0("H>M: ", x, " (mean of 2 seeds)")
d$class <- factor(ifelse(as.character(col(HM("below baseline on >=1 metric"))) %in% c("True", "TRUE"),
                         LAB_FAIL, LAB_OK), levels = c(LAB_OK, LAB_FAIL))
fail_tfs <- d$TF[d$class == LAB_FAIL]
SIM_COL <- "Similarity of the baseline human monoPWM to the nearest mouse monoPWM (Pearson r of aligned columns)"

# label placement only (data untouched): nudges (dx, dy) in data units that move labels into empty space
# With the two-seed mean, TF65 (auROC) and PRGR (auPRC) are the TFs below the baseline; labels beside their diamonds
NUDGE <- list(A = list(TF65 = c(-0.012, 0.010), PRGR = c(-0.030, 0.007)),
              B = list(PRGR = c(-0.035, 0)),   # PRGR off the points
              C = list(TF65 = c(-0.012, -0.001), PRGR = c(-0.007, 0.008)))
# Extra room left of x = 0 in panel B, so the PRGR label sits beside its diamond, clear of the points
X_ROOM <- list(A = 0, B = 0.04)
nudge_of <- function(panel, tf, i) vapply(tf, function(t) { v <- NUDGE[[panel]][[t]]; if (is.null(v)) 0 else v[i] }, 0)

# panels A / B: human-trained vs mouse-trained gain on the mouse test set
panel_ab <- function(metric, letter) {
  pn <- if (metric == "ROC") "A" else "B"
  x <- col(HM(paste0("dau", metric))); y <- col(paste0("M>M: dau", metric))
  lim <- c(min(x, y) - 0.015, max(x, y) + 0.015)
  dd <- data.frame(TF = d$TF, x = x, y = y, class = d$class)
  dd <- dd[order(dd$class == LAB_FAIL), ]      # highlighted TFs last in the data = drawn on top
  ggplot(dd, aes(x = x, y = y)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_vline(xintercept = 0, linetype = "dotted") +
    geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
    geom_point(aes(fill = class, shape = class), color = "black", size = 4, stroke = 0.9) +
    geom_text_repel(aes(label = ifelse(class == LAB_FAIL, TF, "")),
                    size = 5, family = "Arial", point.padding = 0.5, box.padding = 0.3, force = 0.5,
                    min.segment.length = 0, max.overlaps = Inf, seed = 1,
                    nudge_x = nudge_of(pn, dd$TF, 1), nudge_y = nudge_of(pn, dd$TF, 2),
                    max.time = 1, max.iter = 1e5, segment.size = 0.4) +
    scale_fill_manual(values = FILLS, name = NULL, drop = FALSE) + scale_shape_manual(values = SHAPES, name = NULL, drop = FALSE) +
    coord_cartesian(xlim = c(lim[1] - X_ROOM[[pn]], lim[2]), ylim = lim) +
    labs(title = letter,
         subtitle = "Mouse test set",
         x = paste0("Δ", ifelse(metric == "ROC", "auROC", "auPRC"), ", human-trained"),
         y = paste0("Δ", ifelse(metric == "ROC", "auROC", "auPRC"), ", mouse-trained")) +
    base_theme()
}

# panel C: motif similarity vs cross-species gain
sim <- col(SIM_COL); gain <- col(HM("dauROC"))
ct <- suppressWarnings(cor.test(sim, gain, method = "spearman"))
dd_c <- data.frame(TF = d$TF, x = sim, y = gain, class = d$class)
xlim_c <- range(sim) + c(-1, 1) * 0.05 * diff(range(sim))
ylim_c <- range(gain) + c(-1, 1) * 0.05 * diff(range(gain))
panel_c <- ggplot(dd_c, aes(x = x, y = y)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_point(aes(fill = class, shape = class), color = "black", size = 4, stroke = 0.9) +
  geom_text_repel(aes(label = ifelse(class == LAB_FAIL, TF, "")),
                  size = 5, family = "Arial", point.padding = 0.5, box.padding = 0.3, force = 0.5,
                  min.segment.length = 0, max.overlaps = Inf, seed = 1,
                  nudge_x = nudge_of("C", dd_c$TF, 1), nudge_y = nudge_of("C", dd_c$TF, 2),
                  max.time = 1, max.iter = 1e5, segment.size = 0.4) +
  scale_fill_manual(values = FILLS, name = NULL, drop = FALSE) + scale_shape_manual(values = SHAPES, name = NULL, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = 0.07)) +
  scale_y_continuous(expand = expansion(mult = 0.07)) +
  labs(title = "C",
       subtitle = sprintf("Mouse test set, Spearman ρ = %.2f, P = %.2f", unname(ct$estimate), ct$p.value),
       x = "Similarity of the baseline human monoPWM\nto the nearest mouse monoPWM (Pearson r)",
       y = "ΔauROC, human-trained") +
  base_theme()

# panel D: the gain in the three settings
settings <- c("human-trained,\nhuman test", "human-trained,\nmouse test", "mouse-trained,\nmouse test")
dd_d <- bind_rows(
  data.frame(setting = settings[1], `ΔauROC` = col("H>H: dauROC (reference)"),
             `ΔauPRC` = col("H>H: dauPRC (reference)"), check.names = FALSE),
  data.frame(setting = settings[2], `ΔauROC` = col(HM("dauROC")),
             `ΔauPRC` = col(HM("dauPRC")), check.names = FALSE),
  data.frame(setting = settings[3], `ΔauROC` = col("M>M: dauROC"),
             `ΔauPRC` = col("M>M: dauPRC"), check.names = FALSE)) %>%
  pivot_longer(-setting, names_to = "metric", values_to = "value") %>%
  mutate(setting = factor(setting, levels = settings),
         metric  = factor(metric, levels = names(MET_COLS)))

panel_d <- ggplot(dd_d, aes(x = setting, y = value, fill = metric)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.8, seed = 1),
              shape = 21, colour = "black", stroke = 0.6, size = 3) +
  geom_boxplot(aes(group = interaction(setting, metric)), fill = NA, color = "black",
               linewidth = 1, outlier.shape = NA, position = position_dodge(0.8)) +
  scale_fill_manual(values = MET_COLS, name = NULL) +
  labs(title = "D", subtitle = "Gain over best monoPWM",
       x = NULL, y = "Δ (model minus best monoPWM)") +
  base_theme() + theme(axis.text.x = element_text(size = 15), axis.title.y = element_text(size = 17))

combined <- (panel_ab("ROC", "A") | panel_ab("PRC", "B")) / (panel_c | panel_d) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.box = "horizontal")

pdf_out <- file.path(out_dir, "Figure_S6.pdf")
save_both(combined, file.path(out_dir, "Figure_S6"), 320 / 25.4, 280 / 25.4)

# stdout diagnostics
cat("\n==== Fig. S6 ====\n")
cat("TFs:", nrow(d), "| below the baseline after transfer:", paste(fail_tfs, collapse = ", "), "\n")
cat(sprintf("Spearman (motif similarity vs H>M dauROC): rho = %.3f, P = %.3f\n",
            unname(ct$estimate), ct$p.value))
print(as.data.frame(dd_d %>% group_by(setting, metric) %>%
        summarise(n = n(), median = median(value), min = min(value), .groups = "drop")), digits = 3)
cat("Written:", pdf_out, "\n")
