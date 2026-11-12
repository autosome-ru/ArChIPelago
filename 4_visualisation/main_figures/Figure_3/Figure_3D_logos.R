# =============================================================================
# Figure_3D_logos.R -- Fig. 3D from the final importances. Columns RXRA, SRF, ANDR; rows the two
# monoPWMs and the two diPWMs with the highest Gini importance (Figure_3D_source_data.csv, build_Figure_3D_source_data_h3.py).
# Logos as in Fig. S7 (ggseqlogo, information content in bits; diPWMs as their mononucleotide marginal; rc = reverse
# complement for display so that all logos of a TF share the strand of its top monoPWM).
# Output: Figure_3D_logos.pdf (cairo_pdf, vector, width of Figure_3_ABC_restyled.pdf) and Figure_3D_logos.png (300 dpi).
# Usage : Rscript Figure_3D_logos.R
# =============================================================================
suppressPackageStartupMessages({library(ggplot2); library(ggseqlogo); library(patchwork)})
script_dir <- local({
  f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})
if (!interactive()) pdf(NULL)
TFS <- c("RXRA", "SRF", "ANDR")
W_IN <- 878 / 72; H_IN <- 5.0          # width of the A-C page (878 pt)
d <- read.csv(file.path(script_dir, "Figure_3D_source_data.csv"), stringsAsFactors = FALSE)
if (!setequal(unique(d$TF), TFS) || any(table(d$TF) != 4)) stop("expected 4 PWMs for each of ", paste(TFS, collapse = ", "))
read_matrix <- function(f) {
  m <- as.matrix(read.delim(file.path(script_dir, f), row.names = 1, check.names = FALSE))
  if (!identical(rownames(m), c("A", "C", "G", "T"))) stop("bad matrix rows in ", f)
  colnames(m) <- NULL
  m
}
x_max <- max(d$motif_length)
logo <- function(r) {
  m <- read_matrix(r$matrix_file)
  ttl <- sprintf("%s %d · %s%s · importance %.3f", ifelse(r$pwm_class == "mono", "mono", "di"), r$rank_within_class, r$dataset,
                 ifelse(r$strand_for_display == "-", " (rc)", ""), r$importance)
  ggplot() + geom_logo(m, method = "bits", seq_type = "dna", font = "helvetica_bold") +
    theme_logo(base_family = "Arial") +
    coord_cartesian(xlim = c(0.5, x_max + 0.5), ylim = c(0, 2), expand = FALSE) +
    scale_y_continuous(breaks = c(0, 1, 2)) +
    labs(title = ttl, x = NULL, y = NULL) +
    theme(plot.title = element_text(size = 9, hjust = 0, colour = "black", margin = margin(b = 1)),
          axis.text.x = element_blank(), axis.text.y = element_text(size = 7, colour = "black"),
          panel.grid = element_blank(), legend.position = "none", plot.margin = margin(2, 8, 2, 2))
}
cols <- lapply(TFS, function(tf) {
  x <- d[d$TF == tf, ]; x <- x[order(x$pwm_class != "mono", x$rank_within_class), ]
  head_txt <- sprintf("%s   ΔauROC %+.3f, ΔauPRC %+.3f", tf, x$delta_auROC[1], x$delta_auPRC[1])
  wrap_plots(lapply(seq_len(nrow(x)), function(i) logo(x[i, ])), ncol = 1) +
    plot_annotation(title = head_txt, theme = theme(plot.title = element_text(family = "Arial", size = 12, hjust = 0.5, colour = "black")))
})
p <- (wrap_elements(cols[[1]]) | wrap_elements(cols[[2]]) | wrap_elements(cols[[3]])) +
  # panel header says "motif variants" (RXRA DR1 and half-site are variants, not HOCOMOCO subtypes)
  plot_annotation(title = "D", subtitle = "The two most important monoPWMs and diPWMs of each TF: motif variants (gains: human test set)",
                  theme = theme(plot.title = element_text(family = "Arial", face = "bold", size = 26, hjust = 0, colour = "black"),
                                plot.subtitle = element_text(family = "Arial", size = 13, hjust = 0.5, colour = "black", margin = margin(t = -22, b = 6))))
for (ext in c("pdf", "png")) {
  f <- file.path(script_dir, paste0("Figure_3D_logos.", ext))
  if (ext == "pdf") ggsave(f, p, width = W_IN, height = H_IN, device = cairo_pdf)
  else ggsave(f, p, width = W_IN, height = H_IN, dpi = 300, bg = "white")
}
cat("Written: Figure_3D_logos.pdf/png;", nrow(d), "logos\n")
