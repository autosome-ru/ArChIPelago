# archi_style.R - shared colours, shapes and theme of all ArChIPelago figures.
# Sourced by every producer in 4_visualisation/figures/. Python twin: archi_style.py (same hex values, checked by check_palette.py).
#
# Series (Fig. 2, S1, S2 and wherever models are compared):
#   single best monoPWM  red     #D62728  triangle   (reference series, dashed line = median over TFs)
#   single best diPWM    orange  #FF7F0E  triangle
#   ArChIPelago models   circles; Logistic Regression #0072B2, Random Forest #009E73, XGBoost #56B4E9,
#                        Bagging (Gradient Boosting) #CC79A7, Bagging (Logistic Regression) #332288
# Okabe-Ito based; minimum CIEDE2000 between any two of the 7 colours: 22.3 (normal vision), 15.7 (deuteranopia),
# 12.2 (protanopia), 10.2 (tritanopia) - see check_palette.py.
# PWM feature sets of one model (Fig. 3, S2, S3: Random Forest): three lightness steps of the Random Forest green,
# black outline, so the Random Forest keeps its colour and the sets are told apart by lightness.
# Everything opaque (no alpha: Illustrator fails on transparency), all text black, Arial.

ARCHI_PKG <- Sys.getenv("ARCHI_PKG", "../../results/round1_package/Manuscript_analysis")   # released package, read only

MODEL_KEYS <- c("Single best mono PWM", "Single best di PWM", "LogisticRegression", "RandomForestClassifier",
                "XGBClassifier", "BaggingClassifier_XGBClassifier", "BaggingClassifier_LogisticRegression")
MODEL_LABELS <- c("Single best monoPWM", "Single best diPWM", "Logistic Regression", "Random Forest (RF)",
                  "Gradient Boosting (XGBoost)", "Bagging (Gradient Boosting)", "Bagging (Logistic Regression)")
MODEL_SHORT  <- c("best\nmonoPWM", "best\ndiPWM", "LogReg", "RF", "XGBoost", "Bag-\nGB", "Bag-\nLogReg")
MODEL_COL <- setNames(c("#D62728", "#FF7F0E", "#0072B2", "#009E73", "#56B4E9", "#CC79A7", "#332288"), MODEL_LABELS)
MODEL_SHAPE <- setNames(c(24, 24, 21, 21, 21, 21, 21), MODEL_LABELS)   # 24 = triangle (single PWMs), 21 = circle (ArChIPelago)
RF_COL <- unname(MODEL_COL["Random Forest (RF)"])

PWMSET_KEYS   <- c("mono", "di", "mono+di")
PWMSET_LABELS <- c("monoPWMs", "diPWMs", "monoPWMs + diPWMs")
PWMSET_COL <- setNames(c("#B7E3D3", "#009E73", "#00503B"), PWMSET_LABELS)

GREY_PTS  <- "#BDBDBD"   # per-TF background points / lines
GREY_BOX  <- "#EDEDED"
SINGLE_COL <- "#5A5A5A"  # Slim / diChIPMunk models (single-motif competitors), triangles in Fig. 4 / S4

theme_archi <- function(base_size = 16) {
  ggplot2::theme_classic(base_size = base_size, base_family = "Arial") +
    ggplot2::theme(
      text = ggplot2::element_text(colour = "black"),
      axis.text = ggplot2::element_text(colour = "black"),
      axis.title = ggplot2::element_text(colour = "black"),
      axis.line = ggplot2::element_line(colour = "black", linewidth = 0.7),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.7),
      axis.ticks.length = grid::unit(0.18, "cm"),
      strip.text = ggplot2::element_text(colour = "black", size = base_size),
      strip.background = ggplot2::element_rect(fill = "#EDEDED", colour = NA),
      plot.title = ggplot2::element_text(colour = "black", face = "bold", hjust = 0, size = base_size * 1.4),
      plot.subtitle = ggplot2::element_text(colour = "black", hjust = 0),
      plot.caption = ggplot2::element_text(colour = "black", hjust = 0, size = base_size * 0.8),
      legend.text = ggplot2::element_text(colour = "black"),
      legend.title = ggplot2::element_text(colour = "black"),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA))
}

# Save as cairo PDF (text stays text) and as a 300 dpi PNG rendered from that PDF.
# `path_noext` without extension; `plot` is a ggplot/patchwork or a function that draws (for gtables).
save_both <- function(plot, path_noext, width_in, height_in) {
  pdf_path <- paste0(path_noext, ".pdf")
  grDevices::cairo_pdf(pdf_path, width = width_in, height = height_in, onefile = FALSE)
  if (is.function(plot)) plot() else print(plot)
  grDevices::dev.off()
  png_stub <- path_noext
  st <- system2(Sys.getenv("PDFTOPPM", "pdftoppm"), c("-r", "300", "-png", "-singlefile", shQuote(pdf_path), shQuote(png_stub)))
  if (st != 0) stop("pdftoppm failed for ", pdf_path)
  cat("Written:", pdf_path, "and", paste0(path_noext, ".png"), "\n")
}

# P-value as text with Unicode superscripts: 1.2e-11 -> "P = 1.2×10⁻¹¹"
fmt_p <- function(p) {
  if (is.na(p)) return("P = n/a")
  if (p >= 0.01) return(sprintf("P = %.2f", p))
  e <- floor(log10(p)); m <- p / 10^e
  sup <- c("0"="⁰","1"="¹","2"="²","3"="³","4"="⁴","5"="⁵","6"="⁶","7"="⁷","8"="⁸","9"="⁹","-"="⁻")
  ex <- paste(sup[strsplit(as.character(e), "")[[1]]], collapse = "")
  sprintf("P = %.1f×10%s", m, ex)
}
