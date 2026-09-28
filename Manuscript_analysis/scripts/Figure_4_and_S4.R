# =============================================================================
# Figure_4_and_S4.R -- Fig. 4 (human test set) and Fig. S4 (mouse test set, mouse chr1/8/19).
#
# Quasirandom points + boxplots per model, dashed zero line, one-sided Wilcoxon stars,
# theme_classic; the reference is the best single monoPWM (the *_PWM columns of the
# results table); every row holds all 36 TFs (a missing value stops the script).
#
# Input : Figures/source_data/Figure_4_source_data.csv, Figure_S4_source_data.csv
#         (written by analysis/fig4_refit/assemble_fig4.py)
# Usage : Rscript Figure_4_and_S4.R [source_data_dir] [output_dir]
#         (defaults: ../Figures/source_data and ../Figures/panels)
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(ggbeeswarm); library(patchwork)})

script_dir <- local({
  f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})
args    <- commandArgs(trailingOnly = TRUE)
pkg_dir <- normalizePath(file.path(script_dir, ".."))   # Manuscript_analysis/
src_dir <- if (length(args) >= 1) args[1] else file.path(pkg_dir, "Figures", "source_data")
out_dir <- if (length(args) >= 2) args[2] else file.path(pkg_dir, "Figures", "panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!interactive()) pdf(NULL)

ORDER <- rev(c("Slim m=0", "Slim m=1", "LSlim m=-5", "diChIPMunk", "ArChIPelago RF2f",
               "ArChIPelago RF2f + diChIPMunk", "ArChIPelago RF2f + Slim m=1", "ArChIPelago RF2f + LSlim m=-5",
               "ArChIPelago RF2f + Slim m=1, LSlim m=-5 + diChIPMunk", "ArChIPelago RF on all PWMs"))

create_panel <- function(data, title, x_label, show_y_labels = TRUE) {
  if (anyNA(data$value)) stop("missing values in the source data (every model must have all 36 TFs)")
  stats <- data %>% group_by(variable) %>%
    summarise(n = n(), p = wilcox.test(value, alternative = "greater")$p.value, .groups = "drop") %>%
    mutate(lbl = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ ""))
  if (any(stats$n != 36)) stop("every model must have n = 36 TFs; found: ", paste(stats$n, collapse = ", "))
  p <- ggplot(data, aes(x = value, y = variable)) +
    geom_quasirandom(color = "gray80", size = 1, alpha = 0.8, groupOnX = FALSE) +
    geom_boxplot(fill = NA, color = "black", outlier.shape = NA) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "darkred", alpha = 1) +
    geom_text(data = stats, aes(x = -0.1, label = lbl), color = "red", size = 5, hjust = 1, vjust = 0.3,
              family = "Arial") +
    labs(title = title, x = x_label) +
    theme_classic(base_family = "Arial") +
    theme(text = element_text(color = "black"), axis.text = element_text(color = "black"),
          axis.title = element_text(color = "black"), plot.title = element_text(color = "black", face = "bold"),
          axis.line = element_line(color = "black"), axis.ticks = element_line(color = "black"))
  if (show_y_labels) p <- p + labs(y = "Model") else
    p <- p + theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
                   axis.ticks.y = element_blank(), axis.line.y = element_blank()) + labs(y = NULL)
  p
}

draw <- function(csv, out_pdf, label) {
  d <- read.csv(file.path(src_dir, csv), stringsAsFactors = FALSE) %>%
    mutate(variable = factor(model, levels = ORDER))
  if (anyNA(d$variable)) stop("unknown model label in ", csv)
  if (nrow(d) != 36 * length(ORDER)) stop(csv, ": expected ", 36 * length(ORDER), " rows, found ", nrow(d))
  roc <- d %>% transmute(TF, variable, value = d_auROC)
  prc <- d %>% transmute(TF, variable, value = d_auPRC)
  pa <- create_panel(roc, "A", expression(Delta ~ "auROC"))
  pb <- create_panel(prc, "B", expression(Delta ~ "auPRC"), show_y_labels = FALSE)
  combined <- pa + pb + plot_layout(widths = c(1.5, 1))
  ggsave(out_pdf, combined, device = grDevices::cairo_pdf, width = 8, height = 6)
  cat("\n====", label, ": median delta vs the best monoPWM, TFs above 0, one-sided Wilcoxon P ====\n")
  print(as.data.frame(d %>% group_by(model) %>%
          summarise(n = n(), med_dROC = round(median(d_auROC), 4), pos_ROC = sum(d_auROC > 0),
                    p_ROC = signif(wilcox.test(d_auROC, alternative = "greater")$p.value, 2),
                    med_dPRC = round(median(d_auPRC), 4), pos_PRC = sum(d_auPRC > 0),
                    p_PRC = signif(wilcox.test(d_auPRC, alternative = "greater")$p.value, 2),
                    .groups = "drop") %>% arrange(match(model, rev(ORDER)))), row.names = FALSE)
  cat("Written:", out_pdf, "\n")
}

draw("Figure_4_source_data.csv", file.path(out_dir, "Figure_4_human_test.pdf"), "Figure 4 (human test set)")
draw("Figure_S4_source_data.csv", file.path(out_dir, "Figure_S4_mouse_test.pdf"), "Figure S4 (mouse test set, chr1/8/19)")
