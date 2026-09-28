# =============================================================================
# Figure_S7_motif_subtypes.R -- Fig. S7:
# for the TFs with the largest ArChIPelago gains over the best monoPWM (human test
# set, mono+di RF), the PWMs the Random Forest relies on most. One row per TF:
# the top-4 monoPWMs and the top-2 diPWMs by RF feature importance, drawn as
# information-content logos (ggseqlogo, bits), TF name and dauROC / dauPRC at the left.
#
# Input : Figure_S7_source_data.csv (built by build_source_data.py) + the probability
#         matrices it points to (matrices/<TF>_<feature>.txt, 4 x L, rows A C G T).
#         diPWMs are shown as the mononucleotide marginal of the dinucleotide model
#         (see build_source_data.py); logos on the '-' strand relative to the top
#         monoPWM of the TF are reverse-complemented for display ("rc" in the title).
# Output: ../../Figures/panels/Figure_S7_motif_subtypes.pdf (cairo_pdf, Arial), ../../Figures/previews/
#         Figure_S7_motif_subtypes.png, preview/<TF>.png (150 dpi, one row per TF) for inspection.
# Usage : Rscript Figure_S7_motif_subtypes.R [TF ...]   (edit TFS below, or list TFs on the command line)
# Style : as scripts/Figure_S6_cross_species.R (Arial, theme_minimal,
#         cairo_pdf); logo letters use ggseqlogo's built-in glyphs.
# =============================================================================

suppressPackageStartupMessages({library(ggplot2); library(ggseqlogo); library(patchwork)})

# ---- choose the TFs to draw (row order of the figure); NULL = all candidates -------------
TFS <- c("E2F4", "RXRA", "TAL1", "TFE2")   # the TFs whose top-ranked PWMs show distinct subtypes (tf_ranking.csv)
PREVIEWS <- TRUE                  # also write preview/<TF>.png for every TF drawn
cli <- commandArgs(trailingOnly = TRUE)      # Rscript Figure_S7_motif_subtypes.R TAL1 RXRA  overrides TFS
if (length(cli)) TFS <- cli

script_dir <- local({
  f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})
src_path <- file.path(script_dir, "Figure_S7_source_data.csv")
if (!file.exists(src_path)) stop("Input not found: ", src_path)
if (!interactive()) pdf(NULL)

# ---- house style ------------------------------------------------------------
theme_set(theme_get() + theme(text = element_text(family = "Arial")))
N_MONO <- 4; N_DI <- 2; N_COL <- N_MONO + N_DI
LOGO_FONT <- "helvetica_bold"

d <- read.csv(src_path, stringsAsFactors = FALSE, check.names = FALSE)
all_tfs <- unique(d$TF)                       # ranking order from build_source_data.py
if (is.null(TFS)) TFS <- all_tfs
missing <- setdiff(TFS, all_tfs)
if (length(missing)) stop("TF(s) not in the source data: ", paste(missing, collapse = ", "))
d <- d[d$TF %in% TFS, ]
x_max <- max(d$motif_length)                  # common bp scale for all logos

read_matrix <- function(f) {
  m <- as.matrix(read.delim(file.path(script_dir, f), row.names = 1, check.names = FALSE))
  if (!identical(rownames(m), c("A", "C", "G", "T"))) stop("bad matrix rows in ", f)
  colnames(m) <- NULL
  m
}

short_cell <- function(s, n = 16) ifelse(nchar(s) == 0, "", ifelse(nchar(s) > n, paste0(substr(s, 1, n - 1), "…"), s))

logo_panel <- function(r) {
  m <- read_matrix(r$matrix_file)
  ttl <- sprintf("%s %d · %s.%s%s", ifelse(r$pwm_class == "mono", "mono", "di"),
                 r$rank_within_class, r$dataset, r$length_class,
                 ifelse(r$strand_for_display == "-", " (rc)", ""))
  cell <- short_cell(r$cell_line)
  sub <- sprintf("importance %.3f%s", r$importance, ifelse(cell == "", "", paste0(" · ", cell)))
  ggplot() +
    geom_logo(m, method = "bits", seq_type = "dna", font = LOGO_FONT) +
    theme_logo(base_family = "Arial") +
    coord_cartesian(xlim = c(0.5, max(ncol(m), 12) + 0.5), ylim = c(0, 2), expand = FALSE) +   # own length, at least 12 bp wide
    scale_y_continuous(breaks = c(0, 1, 2)) +
    labs(title = ttl, subtitle = sub, x = NULL, y = NULL) +
    theme(plot.title = element_text(size = 9.5, hjust = 0, margin = margin(b = 1)),
          plot.subtitle = element_text(size = 8.5, hjust = 0, colour = "grey25", margin = margin(b = 2)),
          axis.text.x = element_blank(), axis.text.y = element_text(size = 7.5, colour = "grey40"),
          panel.grid = element_blank(),
          plot.margin = margin(4, 10, 2, 2),
          legend.position = "none")
}

label_panel <- function(tf, r1, letter) {
  txt <- sprintf("%s\nΔauROC %+.3f\nΔauPRC %+.3f", tf, r1$delta_auROC, r1$delta_auPRC)
  ggplot() + annotate("text", x = 0, y = 0.38, label = txt, hjust = 0, vjust = 0.5,
                      size = 4.2, family = "Arial", lineheight = 0.95) +
    annotate("text", x = 0, y = 1.02, label = letter, hjust = 0, vjust = 1, size = 6.5, fontface = "bold", family = "Arial") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") + theme_void() +
    theme(plot.margin = margin(4, 2, 2, 6))
}

tf_row <- function(tf, letter = "") {
  dd <- d[d$TF == tf, ]
  dd <- dd[order(dd$panel), ]
  if (nrow(dd) != N_COL) stop(tf, ": expected ", N_COL, " panels, got ", nrow(dd))
  c(list(label_panel(tf, dd[1, ], letter)), lapply(seq_len(nrow(dd)), function(i) logo_panel(dd[i, ])))
}

assemble <- function(tfs) {
  panels <- unlist(lapply(seq_along(tfs), function(i) tf_row(tfs[i], LETTERS[i])), recursive = FALSE)
  wrap_plots(panels, ncol = N_COL + 1, widths = c(0.62, rep(1, N_COL)), byrow = TRUE)
}

# ---- per-TF previews --------------------------------------------------------
if (PREVIEWS) {
  dir.create(file.path(script_dir, "preview"), showWarnings = FALSE)
  for (tf in TFS) {
    ggsave(file.path(script_dir, "preview", paste0(tf, ".png")), assemble(tf),
           width = 400, height = 34, units = "mm", dpi = 150, bg = "white")
  }
}

# ---- the figure -------------------------------------------------------------
fig <- assemble(TFS) +
  plot_annotation(
    title = "PWMs with the highest Random Forest feature importance (mono+di ArChIPelago, human test set)",
    subtitle = paste0("top-4 monoPWMs and top-2 diPWMs per TF; logos in bits, each on its own length; ",
                      "diPWMs shown as the mononucleotide marginal;\nrc = reverse complement of the file, ",
                      "oriented as monoPWM 1; Δ = gain of ArChIPelago over the best single monoPWM"),
    theme = theme(plot.title = element_text(family = "Arial", face = "bold", size = 15, hjust = 0),
                  plot.subtitle = element_text(family = "Arial", size = 11, hjust = 0, colour = "grey25")))

H <- 20 + 32 * length(TFS)
pdf_out <- file.path(script_dir, "..", "..", "Figures", "panels", "Figure_S7_motif_subtypes.pdf")
png_out <- file.path(script_dir, "..", "..", "Figures", "previews", "Figure_S7_motif_subtypes.png")
dir.create(dirname(png_out), showWarnings = FALSE)
ggsave(pdf_out, fig, device = grDevices::cairo_pdf, width = 400, height = H, units = "mm", limitsize = FALSE)
ggsave(png_out, fig, width = 400, height = H, units = "mm", dpi = 150, bg = "white", limitsize = FALSE)

cat("TFs:", paste(TFS, collapse = ", "), "\n")
cat("Panels:", nrow(d), "| bp scale 1-", x_max, "\n")
cat("Written:", pdf_out, "\n         ", png_out, "\n")
