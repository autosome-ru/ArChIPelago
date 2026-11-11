# Figure_3_and_S3_restyled.R - Fig. 3 A-C (human test set) and Fig. S3 A-C (mouse test set), restyled.
# Based on the released Figure_3ABC_human_test.R / Figure_S3ABC_mouse_test.R; data handling unchanged. Presentation only:
#   * A, B (with marginal histograms drawn as outlines, so that overlapping sets stay visible): Random Forest on monoPWMs / diPWMs / monoPWMs + diPWMs against the best single monoPWM; the feature sets are three
#     lightness steps of the Random Forest green (../archi_style.R: PWMSET_COL), opaque circles with black outline, with a legend
#   * C: bubbles opaque (no alpha), viridis option D instead of plasma (plasma uses the red/orange of the PWM reference series)
#   * panels assembled with cowplot so that a complete figure PDF/PNG exists (previously assembled in Illustrator; panel D = logos,
#     not drawn here)
# Input : <ARCHI_PKG>/results_table.csv     Usage: Rscript Figure_3_and_S3_restyled.R [output_dir]
# Writes: Figure_3_ABC_restyled.{pdf,png}, Figure_S3_ABC_restyled.{pdf,png}, *_plotdata.csv
# The TF set must equal the 35 reference TFs (check_tfs); ARCHI_DRAW=H draws Fig. 3 only.
# P53-free input: ARCHI_PKG=3_analysis/integration/inputs (made by 3_analysis/integration/make_inputs.py).
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(cowplot); library(ggrepel); library(ggExtra); library(viridis)})
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
tab <- read.csv(file.path(ARCHI_PKG, "results_table.csv"), sep = "\t")
F3C_ROC_GAP <- 0.075   # = ROC_GAP of 3_analysis/transfer_drivers/B5_01_drivers_and_ceiling.py (checked by a test)
F3C_PRC_GAP <- 0.15    # = PRC_GAP of the same script

draw <- function(sp, fig_name, label) {
  test_lab <- if (sp == "H") "Human test set" else "Mouse test set\n(human-trained)"
  roc <- paste0("roc_auc_test_", sp); prc <- paste0("pr_auc_test_", sp)
  dA <- tab %>% filter(Model == "RandomForestClassifier") %>%
    transmute(TF_name, PWM, base_roc = .data[[paste0(roc, "_PWM")]], roc = .data[[roc]],
              base_prc = .data[[paste0(prc, "_PWM")]], prc = .data[[prc]], Count, Seq_count)
  if (sp == "M") {
    # the human-trained monoPWM + diPWM model on the mouse test set as the mean of Random Forest seeds
    # 0 and 1 (hm_mm copy of the package); the monoPWM-only and diPWM-only models have seed 0 only. Seed 0 and baselines asserted.
    hm <- read.csv(file.path(ARCHI_PKG, "analysis", "mouse_transfer", "mouse_transfer_results.csv"))
    k <- dA$PWM == "mono+di"; i <- match(dA$TF_name[k], hm$TF)
    if (anyNA(i) || max(abs(dA$roc[k] - hm$HM_rf_test_M_auroc_s0[i])) > 1e-6 || max(abs(dA$prc[k] - hm$HM_rf_test_M_auprc_s0[i])) > 1e-6 ||   # results table: ~8 digits
        max(abs(dA$base_roc[k] - hm$HM_broc_test_M_auroc[i])) > 1e-9 || max(abs(dA$base_prc[k] - hm$HM_bprc_test_M_auprc[i])) > 1e-9)
      stop("CHECK FAILED ", fig_name, ": results table mono+di (mouse test) is not hm_mm seed 0 / its baselines")
    dA$roc[k] <- hm$HM_rf_test_M_auroc_mean[i]; dA$prc[k] <- hm$HM_rf_test_M_auprc_mean[i]
    dA$RF_seed <- ifelse(k, "mean of 0 and 1", "0")
  }
  dA$PWMset <- factor(PWMSET_LABELS[match(dA$PWM, PWMSET_KEYS)], levels = PWMSET_LABELS)
  ntf <- check_tfs(dA$TF_name, paste(fig_name, "A-B"))
  if (any(table(dA$PWMset) != ntf)) stop("CHECK FAILED ", fig_name, ": a feature set without all ", ntf, " TFs")
  dC <- dA %>% group_by(TF_name) %>% mutate(Count_PWM = sum(Count[PWM == "mono"]) + sum(Count[PWM == "di"])) %>% ungroup() %>% filter(PWM == "mono+di") %>%
    mutate(dROC = roc - base_roc, dPRC = prc - base_prc) %>% group_by(TF_name) %>%
    summarise(dROC = mean(dROC, na.rm = TRUE), dPRC = mean(dPRC, na.rm = TRUE), Count_PWM = first(Count_PWM),
              Seq_count_min = min(Seq_count, na.rm = TRUE), .groups = "drop")
  check_tfs(dC$TF_name, paste(fig_name, "C"))
  # PWM count = monoPWMs + diPWMs (one RF mono row and one RF di row per TF); sum(unique()) collapsed equal counts
  if (any(dA %>% count(TF_name, PWM) %>% filter(PWM %in% c("mono", "di")) %>% pull(n) != 1)) stop("CHECK FAILED ", fig_name, " C: not one mono and one di row per TF")
  write.csv(dA, file.path(out_dir, paste0(fig_name, "_AB_plotdata.csv")), row.names = FALSE)
  write.csv(dC, file.path(out_dir, paste0(fig_name, "_C_plotdata.csv")), row.names = FALSE)
  cat(sprintf("\n==== %s (%s): Random Forest medians per feature set ====\n", fig_name, label))
  print(as.data.frame(dA %>% group_by(PWMset) %>% summarise(n = n(), med_roc = median(roc), med_base_roc = median(base_roc),
        med_prc = median(prc), med_base_prc = median(base_prc), .groups = "drop")), digits = 4)
  cat(sprintf("panel C medians: dROC %.4f dPRC %.4f (n = %d TFs)\n", median(dC$dROC), median(dC$dPRC), nrow(dC)))

  sc <- function(x, y, xl, yl) {
    p <- ggplot(dA %>% arrange(desc(Seq_count)), aes(x = .data[[x]], y = .data[[y]], fill = PWMset, colour = PWMset)) +
      geom_abline(intercept = 0, linetype = "dotted", linewidth = 0.8) +
      geom_point(shape = 21, size = 3.6, stroke = 0.5, colour = "black") +
      scale_fill_manual(values = PWMSET_COL, name = NULL) + scale_colour_manual(values = PWMSET_COL, name = NULL) +
      annotate("text", x = 1, y = 0.02, hjust = 1, vjust = 0, label = test_lab, size = if (sp == "H") 5.2 else 4.4, lineheight = 0.85, family = "Arial", colour = "black") +
      scale_x_continuous(breaks = c(0, 0.5, 1)) + scale_y_continuous(breaks = c(0, 0.5, 1)) +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
      labs(x = xl, y = yl) +
      theme_archi(18) + theme(legend.position = "none", plot.margin = margin(8, 8, 8, 30))
    hx <- ggplot(dA, aes(x = .data[[x]], colour = PWMset)) + geom_histogram(aes(y = after_stat(count)), bins = 25, fill = NA,
            position = "identity", linewidth = 0.7, boundary = 0) +
      scale_colour_manual(values = c("#7FBFA6", "#009E73", "#00503B")) + coord_cartesian(xlim = c(0, 1)) + theme_void() +
      theme(legend.position = "none")
    hy <- ggplot(dA, aes(x = .data[[y]], colour = PWMset)) + geom_histogram(bins = 25, fill = NA, position = "identity",
            linewidth = 0.7, boundary = 0) +
      scale_colour_manual(values = c("#7FBFA6", "#009E73", "#00503B")) + coord_flip(xlim = c(0, 1)) + theme_void() +
      theme(legend.position = "none")
    ggdraw(insert_yaxis_grob(insert_xaxis_grob(p, hx, grid::unit(1.6, "cm"), position = "top"), hy, grid::unit(1.6, "cm"),
                             position = "right"))
  }
  a <- sc("base_roc", "roc", "auROC best monoPWM", "auROC")
  b <- sc("base_prc", "prc", "auPRC best monoPWM", "auPRC")
  leg <- cowplot::get_legend(ggplot(dA, aes(roc, prc, fill = PWMset)) + geom_point(shape = 21, colour = "black", size = 4) +
    scale_fill_manual(values = PWMSET_COL, name = "Random Forest on") + theme_archi(18) +
    theme(legend.position = "bottom", legend.key.width = unit(1.2, "cm"), legend.key.spacing.x = unit(1.0, "cm"), legend.title = element_text(size = 15), legend.text = element_text(size = 15)) +
    guides(fill = guide_legend(nrow = 1)))
  lab_df <- subset(dC, dROC < 0 | dROC > 0.05 | dPRC < 0 | dPRC > 0.2 | TF_name == "MAFK")
  # Fig. S3C legend says the smallest gains are labelled too: add the three smallest auROC gains (mouse test set)
  if (sp == "M") lab_df <- unique(rbind(lab_df, dC %>% arrange(dROC) %>% head(3)))
  # label placement only (data untouched): push labels from crowded bubbles into empty space
  NUDGE <- if (sp == "H") list(SRF = c(0.03, 0.004), RXRA = c(0.045, 0.005), IRF4 = c(-0.03, 0.004), ANDR = c(-0.05, 0), MAFK = c(0.035, 0.006), E2F4 = c(-0.01, -0.006), STAT3 = c(0.012, 0.045)) else list(ERG = c(0.035, 0.011), MAFK = c(0.03, 0), PRGR = c(-0.05, 0.004), E2F4 = c(-0.012, 0.004), STA5A = c(0.09, -0.0005), USF2 = c(0.05, -0.0045), CTCF = c(0.03, -0.006), P53 = c(0.015, 0.008), IRF4 = c(-0.012, 0.014), SRF = c(0.012, -0.013), TF65 = c(0.05, -0.001))   # TF65 is below the baseline in auROC (two-seed mean), label right of its diamond; STAT3 (human) label off the bubbles
  lab_df$nx <- vapply(lab_df$TF_name, function(t) if (!is.null(NUDGE[[t]])) NUDGE[[t]][1] else 0, 0)
  lab_df$ny <- vapply(lab_df$TF_name, function(t) if (!is.null(NUDGE[[t]])) NUDGE[[t]][2] else 0, 0)
  # breaks that the two-decimal labels show exactly (steps 0.01 / 0.02 / 0.05 / 0.1 / 0.2)
  nice_breaks <- function(lim) {
    step <- c(0.01, 0.02, 0.05, 0.1, 0.2)[which(diff(lim) / c(0.01, 0.02, 0.05, 0.1, 0.2) <= 7)[1]]
    round(seq(ceiling(lim[1] / step - 1e-9) * step, floor(lim[2] / step + 1e-9) * step, by = step), 10) + 0
  }
  dC$below <- dC$dROC < 0 | dC$dPRC < 0   # below the baseline: black open diamond (spec rule 1)
  ct <- ggplot(dC %>% arrange(desc(Seq_count_min)), aes(x = dPRC, y = dROC, size = Seq_count_min, fill = Count_PWM)) +
    geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.8) + geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.8) +
    # Fig. 3C group thresholds of 3_analysis/transfer_drivers/B5_01_drivers_and_ceiling.py (ROC_GAP, PRC_GAP), set post hoc
    { if (sp == "H") list(geom_hline(yintercept = F3C_ROC_GAP, linetype = "dashed", linewidth = 0.6, colour = "grey45"),
                          geom_vline(xintercept = F3C_PRC_GAP, linetype = "dashed", linewidth = 0.6, colour = "grey45"),
                          annotate("text", x = Inf, y = 0.033, hjust = 1.02, vjust = 0.5, size = 4.4, family = "Arial", colour = "black", lineheight = 0.9,
                                   label = "dashed lines:\ngroup thresholds,\nset post hoc")) else list() } +
    geom_point(data = function(d) d[!d$below, ], shape = 21, colour = "black", stroke = 0.6) +
    geom_point(data = function(d) d[d$below, ], shape = 23, fill = "white", colour = "black", stroke = 1.1) +
    scale_size(range = c(1.5, 13), breaks = c(500, 1000, 5000, 10000, 50000), name = "Training positives") +
    scale_fill_viridis_c(option = "D", name = "Number of PWMs", guide = guide_colourbar(raster = FALSE)) +   # vector colour bar
    geom_text_repel(data = lab_df, aes(label = TF_name), nudge_x = lab_df$nx, nudge_y = lab_df$ny, size = 5.5,
                    family = "Arial", colour = "black", point.padding = 0.3, box.padding = 0.6, force = 2, force_pull = 0.5, max.overlaps = Inf, min.segment.length = 0, segment.colour = "black", segment.size = 0.5, segment.alpha = 1, max.time = 2, max.iter = 1e5,
                    seed = 1, show.legend = FALSE) +
    # mouse corner text moved to the empty top right, room below zero for the smallest-gain labels (USF2, CTCF, STA5A)
    { if (sp == "H") annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -2.2, label = test_lab, size = 6, family = "Arial", colour = "black")
      else list(annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.2, label = test_lab, size = 6, family = "Arial", colour = "black"),
                # room left of x = 0 so the PRGR label sits beside its diamond, not over the bubbles
                expand_limits(y = -0.009, x = -0.075)) } +
    scale_x_continuous(breaks = nice_breaks, labels = function(x) sprintf("%.2f", x)) + scale_y_continuous(breaks = nice_breaks, labels = function(x) sprintf("%.2f", x)) +
    guides(size = guide_legend(override.aes = list(shape = 21, fill = "white", colour = "black"))) +
    labs(x = "ΔauPRC", y = "ΔauROC") + theme_archi(18) +
    theme(legend.position = "right", plot.margin = margin(8, 8, 8, 40), panel.grid = element_blank(), axis.title.y = element_text(margin = margin(r = 12)))
  left <- plot_grid(a, b, nrow = 2, labels = c("A", "B"), label_size = 26, label_fontfamily = "Arial", label_colour = "black")
  top <- plot_grid(left, ct, ncol = 2, rel_widths = c(1, 2.2), labels = c("", "C"), label_size = 26,
                   label_fontfamily = "Arial", label_colour = "black")
  full <- plot_grid(top, leg, ncol = 1, rel_heights = c(1, 0.07))
  save_both(function() print(full), file.path(out_dir, paste0(fig_name, "_ABC_restyled")), 310 / 25.4, 200 / 25.4)
}
if ("H" %in% ARCHI_DRAW) draw("H", "Figure_3", "human test set")
if ("M" %in% ARCHI_DRAW) draw("M", "Figure_S3", "mouse test set")
