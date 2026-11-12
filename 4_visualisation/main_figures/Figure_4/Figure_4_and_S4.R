# Figure_4_and_S4_restyled.R - Fig. 4 (human test set) and Fig. S4 (mouse test set), restyled.
# Based on the released Figure_4_and_S4.R; same data, same model order and names. Presentation only:
#   * x axis zoomed to the 2nd-98th percentile of all deltas of the panel (padded by 8 %, always containing 0); points beyond the
#     range are not drawn but marked at the edge with an arrow and their count per model
#   * boxes with a light fill and black outline, thick median, whiskers 1.5 x IQR (as before): quartiles are visible
#   * solid black zero line; per model, on the right: "TFs above 0 / n" and the exact two-sided Wilcoxon signed-rank P;
#     the stats table keeps both (p one-sided, p_two_sided)
#   * Slim / LSlim / diChIPMunk (single-motif models) = grey triangles, ArChIPelago RF models = green circles (RF colour of Fig. 2)
#   * the zero-line explanation is in the figure itself: delta = model minus the best single monoPWM of the same TF
# Input : <ARCHI_PKG>/Figures/source_data/Figure_4_source_data.csv, Figure_S4_source_data.csv
# Usage : Rscript Figure_4_and_S4_restyled.R [output_dir]
# Writes: Figure_4_restyled.{pdf,png}, Figure_S4_restyled.{pdf,png}, Figure_4_stats.csv, Figure_S4_stats.csv
# n TFs is taken from the data and must equal the 35 reference TFs (check_tfs);
# ARCHI_DRAW=H draws Fig. 4 only; the stats table also gives the two-sided Wilcoxon P (quoted in the Results).
# P53-free input: ARCHI_PKG=3_analysis/integration/inputs (made by 3_analysis/integration/make_inputs.py).
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(ggbeeswarm); library(patchwork)})
script_dir <- local({
  f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})
source(file.path(script_dir, "..", "archi_style.R"))
source(file.path(script_dir, "..", "..", "../3_analysis/integration", "p53drop.R"))   # check_tfs(), ARCHI_DRAW
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else script_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
if (!interactive()) pdf(NULL)
src_dir <- file.path(ARCHI_PKG, "Figures", "source_data")

ORDER <- rev(c("Slim m=0", "Slim m=1", "LSlim m=-5", "diChIPMunk", "ArChIPelago RF2f",
               "ArChIPelago RF2f + diChIPMunk", "ArChIPelago RF2f + Slim m=1", "ArChIPelago RF2f + LSlim m=-5",
               "ArChIPelago RF2f + Slim m=1, LSlim m=-5 + diChIPMunk", "ArChIPelago RF on all PWMs"))
COMPETITORS <- c("Slim m=0", "Slim m=1", "LSlim m=-5", "diChIPMunk")
PCT <- c(0.02, 0.98); PAD <- 0.08

row_label <- function(x) {   # display names (data keys unchanged)
  x <- gsub(" m=", " m = ", x, fixed = TRUE); x <- sub("ArChIPelago RF on all PWMs", "ArChIPelago RF (all PWMs)", x, fixed = TRUE); x
}
make_panel <- function(data, title, x_label, show_y_labels, letter, ntf) {
  if (anyNA(data$value)) stop("missing values in the source data")
  data <- data %>% mutate(cls = ifelse(as.character(variable) %in% COMPETITORS, "single", "archi"))
  stats <- data %>% group_by(variable) %>%
    summarise(n = n(), n_pos = sum(value > 0), med = median(value), p = wilcox.test(value, alternative = "greater")$p.value,
              p_two_sided = wilcox.test(value)$p.value, .groups = "drop")
  if (any(stats$n != ntf)) stop(sprintf("every model must have n = %d TFs", ntf))
  q <- quantile(data$value, PCT); span <- diff(q)
  lo <- min(q[1] - PAD * span, -0.1 * span); hi <- max(q[2] + PAD * span, 0.1 * span)
  stats <- stats %>% left_join(data %>% group_by(variable) %>%
             summarise(n_lo = sum(value < lo), n_hi = sum(value > hi), .groups = "drop"), by = "variable") %>%
    mutate(lab_n = sprintf("%d/%d", n_pos, ntf), lab_p = vapply(p_two_sided, fmt_p, ""), y = as.numeric(variable))   # two-sided P printed
  inside <- data %>% filter(value >= lo, value <= hi)
  bx <- data %>% group_by(variable) %>% summarise(xlower = quantile(value, 0.25), xmiddle = median(value), xupper = quantile(value, 0.75),
                                                  xmin = min(value[value >= xlower - 1.5 * (xupper - xlower)]),
                                                  xmax = max(value[value <= xupper + 1.5 * (xupper - xlower)]), .groups = "drop") %>%
    mutate(across(c(xmin, xlower, xmiddle, xupper, xmax), ~ pmin(pmax(.x, lo), hi)))
  w <- hi - lo
  p <- ggplot(data, aes(x = value, y = variable)) +
    geom_rect(data = stats, aes(xmin = -Inf, xmax = Inf, ymin = y - 0.5, ymax = y + 0.5), inherit.aes = FALSE,
              fill = ifelse(stats$y %% 2 == 0, "#F2F2F2", "white")) +
    geom_vline(xintercept = 0, colour = "black", linewidth = 0.6, linetype = "dotted") +
    # box statistics computed here and whiskers clamped to the drawn range [lo, hi]; with clip = "off" a whisker
    # beyond the range ran into the "TFs above 0" column (Fig. S4B, five-feature model)
    geom_boxplot(data = bx, aes(y = variable, xmin = xmin, xlower = xlower, xmiddle = xmiddle, xupper = xupper, xmax = xmax, group = variable),
                 inherit.aes = FALSE, stat = "identity", orientation = "y", fill = NA, colour = "black", linewidth = 0.8,
                 median.linewidth = 1.6, width = 0.62, staplewidth = 0.4) +
    geom_quasirandom(data = inside, aes(fill = cls, shape = cls), colour = "black", stroke = 0.35, size = 1.9, groupOnX = FALSE,
                     width = 0.28) +
    geom_segment(data = stats %>% filter(n_lo > 0), aes(x = lo + 0.02 * w, xend = lo, y = y, yend = y),
                 arrow = grid::arrow(length = grid::unit(0.16, "cm"), type = "closed"), linewidth = 0.8, colour = "black") +
    geom_text(data = stats %>% filter(n_lo > 0), aes(x = lo + 0.03 * w, y = y, label = n_lo), hjust = 0, vjust = -0.7,
              size = 3.8, colour = "black", family = "Arial") +
    geom_segment(data = stats %>% filter(n_hi > 0), aes(x = hi - 0.02 * w, xend = hi, y = y, yend = y),
                 arrow = grid::arrow(length = grid::unit(0.16, "cm"), type = "closed"), linewidth = 0.8, colour = "black") +
    geom_text(data = stats %>% filter(n_hi > 0), aes(x = hi - 0.03 * w, y = y, label = n_hi), hjust = 1, vjust = -0.7,
              size = 3.8, colour = "black", family = "Arial") +
    geom_text(data = stats, aes(x = hi + 0.20 * w, y = y, label = lab_n), hjust = 1, size = 4.1, colour = "black", family = "Arial Unicode MS") +
    geom_text(data = stats, aes(x = hi + 0.24 * w, y = y, label = lab_p), hjust = 0, size = 4.1, colour = "black", family = "Arial Unicode MS") +
    scale_fill_manual(values = c(single = SINGLE_COL, archi = RF_COL), guide = "none") +
    # boxes use a light tint of the series colour so that the quartiles stay readable under the points
    scale_shape_manual(values = c(single = 24, archi = 21), guide = "none") +
    scale_x_continuous(breaks = pretty(c(lo, hi), n = 4)) +
    scale_y_discrete(labels = row_label) +
    coord_cartesian(xlim = c(lo, hi), clip = "off") +
    labs(title = title, x = x_label, tag = letter) +
    theme_archi(14) + theme(plot.margin = margin(6, 150, 6, 6), plot.title = element_text(face = "plain", size = 14, hjust = 0.5),
                            plot.tag = element_text(face = "bold", size = 20), plot.tag.position = c(0.005, 1))
  if (show_y_labels) p <- p + labs(y = NULL) else
    p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank()) + labs(y = NULL)
  list(plot = p, stats = stats %>% mutate(lo = lo, hi = hi))
}

draw <- function(csv, fig_name, label, test_lab) {
  d <- read.csv(file.path(src_dir, csv), stringsAsFactors = FALSE) %>% mutate(variable = factor(model, levels = ORDER))
  if (anyNA(d$variable)) stop("unknown model label in ", csv)
  ntf <- check_tfs(d$TF, csv)
  if (nrow(d) != ntf * length(ORDER)) stop(csv, ": unexpected number of rows")
  pa <- make_panel(d %>% transmute(TF, variable, value = d_auROC), test_lab, "ΔauROC (model minus best monoPWM)", TRUE, "A", ntf)
  pb <- make_panel(d %>% transmute(TF, variable, value = d_auPRC), test_lab, "ΔauPRC (model minus best monoPWM)", FALSE, "B", ntf)
  both <- bind_rows(pa$stats %>% mutate(metric = "auROC"), pb$stats %>% mutate(metric = "auPRC")) %>%
    mutate(model = as.character(variable)) %>% select(metric, model, n, n_pos, med, p, p_two_sided, n_lo, n_hi, lo, hi)
  write.csv(both, file.path(out_dir, paste0(fig_name, "_stats.csv")), row.names = FALSE)
  cat("\n====", fig_name, label, "====\n"); print(as.data.frame(both), digits = 4)
  combined <- (pa$plot | pb$plot) + plot_layout(widths = c(1.45, 1))
  save_both(combined, file.path(out_dir, paste0(fig_name, "_restyled")), 16, 6.2)
}
if ("H" %in% ARCHI_DRAW) draw("Figure_4_source_data.csv", "Figure_4", "(human test set)", "Human test set")
if ("M" %in% ARCHI_DRAW) draw("Figure_S4_source_data.csv", "Figure_S4", "(mouse test set)", "Mouse test set")
