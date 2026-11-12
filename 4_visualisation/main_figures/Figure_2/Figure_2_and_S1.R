# Figure_2_and_S1.R - Fig. 2 (human test set) and Fig. S1 (mouse test set), restyled.
# Based on the released Figure_2_human_test.R / Figure_S1_mouse_test.R (results/round1_package/Manuscript_analysis/scripts/); the data handling
# (join with training_positives_per_TF, Seq_count > 100, one row per TF x model x PWM set) is unchanged. Presentation only:
#   * palette and shapes from ../../style/archi_style.R (single PWMs = triangles, ArChIPelago models = circles, opaque points)
#   * one facet per feature set (monoPWMs | diPWMs | monoPWMs + diPWMs) with 6 / 6 / 5 columns, short labels under the columns
#   * number above every box: TFs (of 36) where the model is above the best single monoPWM (the *_PWM column), same input table
#   * dashed lines = median of the best single monoPWM (red) and best single diPWM (orange)
# Input : <ARCHI_PKG>/results_table.csv, <ARCHI_PKG>/analysis/inputs/training_positives_per_TF.csv
# Usage : Rscript Figure_2_and_S1.R [output_dir]     (default: Fig. 2 here, Fig. S1 in ../../supplementary_figures/Figure_S1)
# Writes: Figure_2.{pdf,png}, Figure_S1.{pdf,png}, Figure_2_plotdata.csv, Figure_S1_plotdata.csv
# The TF set must equal the 35 reference TFs (check_tfs); ARCHI_DRAW=H draws Fig. 2 only.
# P53-free input: ARCHI_PKG=3_analysis/integration/inputs (made by 3_analysis/integration/make_inputs.py).
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(patchwork)})
script_dir <- local({
  f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})
source(file.path(script_dir, "..", "..", "style", "archi_style.R"))
source(file.path(script_dir, "..", "..", "../3_analysis/integration", "p53drop.R"))   # check_tfs(), ARCHI_DRAW
# display names (labels only; order, colours and shapes unchanged)
MODEL_LABELS <- c("best monoPWM", "best diPWM", "Logistic Regression (LogReg)", "Random Forest (RF)", "XGBoost",
                  "Bagging (XGBoost)", "Bagging (LogReg)")
MODEL_SHORT  <- c("best\nmonoPWM", "best\ndiPWM", "LogReg", "RF", "XGBoost", "Bag-XGB", "Bag-LogReg")
MODEL_COL <- setNames(unname(MODEL_COL), MODEL_LABELS); MODEL_SHAPE <- setNames(unname(MODEL_SHAPE), MODEL_LABELS)
args <- commandArgs(trailingOnly = TRUE)
# output folder: Fig. 2 next to this script, Fig. S1 in supplementary_figures/Figure_S1; an argument overrides both
out_dir_of <- function(fig_name) {
  d <- if (length(args) >= 1) args[1] else if (fig_name == "Figure_2") script_dir else
    file.path(script_dir, "..", "..", "supplementary_figures", fig_name)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  d
}
if (!interactive()) pdf(NULL)

df <- read.delim(file.path(ARCHI_PKG, "results_table.csv"), sep = "\t")
ddf <- read.delim(file.path(ARCHI_PKG, "analysis", "inputs", "training_positives_per_TF.csv"), sep = "\t")
df <- df %>% left_join(ddf, by = c("TF_name" = "Names")) %>% filter(Seq_count.y > 100) %>%
  distinct(TF_name, Model, PWM, .keep_all = TRUE)
df$Model <- factor(MODEL_LABELS[match(df$Model, MODEL_KEYS)], levels = MODEL_LABELS)
df$PWM <- factor(PWMSET_LABELS[match(df$PWM, PWMSET_KEYS)], levels = PWMSET_LABELS)
if (anyNA(df$Model) || anyNA(df$PWM)) stop("unknown model or PWM-set label in the results table")
SHORT <- setNames(MODEL_SHORT, MODEL_LABELS)

draw <- function(sp, fig_name, test_label) {
  out_dir <- out_dir_of(fig_name)
  set.seed(1)   # fixed jitter, so that two renders of the same data are identical
  roc <- paste0("roc_auc_test_", sp); prc <- paste0("pr_auc_test_", sp)
  roc_b <- paste0(roc, "_PWM"); prc_b <- paste0(prc, "_PWM")        # best single monoPWM of the same TF
  d <- df %>% transmute(TF = TF_name, Model, PWM, auROC = .data[[roc]], auPRC = .data[[prc]],
                        base_auROC = .data[[roc_b]], base_auPRC = .data[[prc_b]])
  ntf <- check_tfs(d$TF, paste(fig_name, "plot data"))
  if (any(table(d$Model, d$PWM)[table(d$Model, d$PWM) > 0] != ntf)) stop("CHECK FAILED ", fig_name, ": a box without all ", ntf, " TFs")
  med <- function(m, col) median(d[[col]][d$Model == m], na.rm = TRUE)
  ref <- list(roc_mono = med("best monoPWM", "auROC"), roc_di = med("best diPWM", "auROC"),
              prc_mono = med("best monoPWM", "auPRC"), prc_di = med("best diPWM", "auPRC"))
  cat(sprintf("\n==== %s (%s test set): medians of the reference series ====\n", fig_name, test_label))
  cat(sprintf("auROC mono %.4f di %.4f | auPRC mono %.4f di %.4f\n", ref$roc_mono, ref$roc_di, ref$prc_mono, ref$prc_di))
  wins <- d %>% group_by(PWM, Model) %>%
    summarise(n = n(), win_auROC = sum(auROC > base_auROC), win_auPRC = sum(auPRC > base_auPRC),
              med_auROC = median(auROC), med_auPRC = median(auPRC), .groups = "drop")
  print(as.data.frame(wins), digits = 4)
  # two-line check on one TF: the count is built from exactly these comparisons
  chk <- d %>% filter(TF == "ANDR", Model == "Random Forest (RF)", PWM == "monoPWMs + diPWMs")
  cat(sprintf("check ANDR, RF on monoPWMs + diPWMs: auROC %.4f vs best monoPWM %.4f -> %s | auPRC %.4f vs %.4f -> %s\n",
              chk$auROC, chk$base_auROC, chk$auROC > chk$base_auROC, chk$auPRC, chk$base_auPRC, chk$auPRC > chk$base_auPRC))
  tot <- d %>% filter(Model == "Random Forest (RF)", PWM == "monoPWMs + diPWMs")
  cat(sprintf("check: RF on monoPWMs + diPWMs, auROC above the best monoPWM in %d of %d TFs (sum of the per-TF comparison)\n",
              sum(tot$auROC > tot$base_auROC), nrow(tot)))
  write.csv(d, file.path(out_dir, paste0(fig_name, "_plotdata.csv")), row.names = FALSE)

  panel <- function(metric, base, hl, ylab, letter, ylim, breaks, ytxt, bottom) {
    wn <- wins %>% mutate(lab = ifelse(Model == "best monoPWM", "ref.",
                                       sprintf("%d/%d", .data[[ifelse(metric == "auROC", "win_auROC", "win_auPRC")]], n)),
                          y = ytxt)
    p <- ggplot(d, aes(x = Model, y = .data[[metric]])) +
      geom_hline(yintercept = hl[1], colour = MODEL_COL[["best monoPWM"]], linetype = "dashed", linewidth = 0.9) +
      geom_hline(yintercept = hl[2], colour = MODEL_COL[["best diPWM"]], linetype = "22", linewidth = 0.9) +
      geom_jitter(aes(fill = Model, shape = Model), colour = "black", stroke = 0.45, size = 2.7, width = 0.2, height = 0,
                  alpha = 1) +
      geom_boxplot(fill = NA, colour = "black", linewidth = 0.8, width = 0.62, outlier.shape = NA) +
      geom_text(data = wn, aes(y = y, label = lab), family = "Arial", size = 4.6, colour = "black") +
      facet_grid(. ~ PWM, scales = "free_x", space = "free_x") +
      scale_fill_manual(values = MODEL_COL, name = NULL, drop = FALSE) +
      scale_shape_manual(values = MODEL_SHAPE, name = NULL, drop = FALSE) +
      scale_y_continuous(breaks = breaks) +
      scale_x_discrete(labels = SHORT) +
      coord_cartesian(ylim = ylim, clip = "off") +
      labs(x = NULL, y = ylab, title = trimws(paste(test_label, "test set", if (sp == "M") "(human-trained)" else "")), tag = letter) +
      theme_archi(16) +
      theme(legend.position = "bottom", panel.spacing.x = unit(0.5, "cm"),
            plot.title = element_text(face = "plain", size = 16, hjust = 0.5), plot.tag = element_text(face = "bold", size = 22),
            plot.tag.position = c(0.005, 1),
            axis.text.x = element_text(size = 9.5, lineheight = 0.85),
            axis.ticks.x = element_blank())
    if (!bottom) p <- p + theme(axis.text.x = element_blank(), axis.line.x = element_blank())
    else p <- p + theme(strip.text = element_blank(), strip.background = element_blank())
    p
  }
  # lower auROC limit from the data (mouse points below a fixed limit would run into the panel B title)
  roc_lo <- min(0.72, floor((min(d$auROC, na.rm = TRUE) - 0.01) * 50) / 50)
  p1 <- panel("auROC", "base_auROC", c(ref$roc_mono, ref$roc_di), "auROC", "A", c(roc_lo, 1.06), seq(0.6, 1, 0.1), 1.035, FALSE)
  p2 <- panel("auPRC", "base_auPRC", c(ref$prc_mono, ref$prc_di), "auPRC", "B", c(0, 1.09), seq(0, 1, 0.2), 1.045, TRUE)
  combined <- (p1 / p2) + plot_layout(guides = "collect", heights = c(1, 1)) +
    plot_annotation(theme = theme_archi(16)) &
    theme(legend.position = "bottom", legend.text = element_text(size = 14))
  combined <- combined & guides(fill = guide_legend(nrow = 2, override.aes = list(size = 3.6)),
                                shape = guide_legend(nrow = 2))
  save_both(combined, file.path(out_dir, fig_name), 14.6, 11.5)
}
if ("H" %in% ARCHI_DRAW) draw("H", "Figure_2", "Human")
if ("M" %in% ARCHI_DRAW) draw("M", "Figure_S1", "Mouse")
