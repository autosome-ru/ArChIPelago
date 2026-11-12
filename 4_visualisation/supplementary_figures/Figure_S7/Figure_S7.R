# =============================================================================
# Figure_S7.R -- supplementary figure:
# for the TFs with the largest ArChIPelago gains over the best monoPWM (human test
# set, mono+di RF), the PWMs the Random Forest relies on most. One row per TF:
# the top-4 monoPWMs and the top diPWMs (three; builder N_DI) by RF feature importance, drawn as
# information-content logos (ggseqlogo, bits), TF name and dauROC / dauPRC at the left.
#
# Input : Figure_S7_source_data.csv (built by build_source_data.py) + the probability
#         matrices it points to (matrices/<TF>_<feature>.txt, 4 x L, rows A C G T).
#         diPWMs are shown as the mononucleotide marginal of the dinucleotide model
#         (see build_source_data.py); logos on the '-' strand relative to the top
#         monoPWM of the TF are reverse-complemented for display ("rc" in the title).
# Output: Figure_S7.pdf (cairo_pdf, Arial), Figure_S7.png,
#         preview/<TF>.png (150 dpi, one row per TF) for inspection.
# Usage : Rscript Figure_S7.R [TF ...]   (edit TFS below, or list TFs on the command line)
# Style : as results/round1_package/Resubmission/scripts/Figure_S6_cross_species.R (Arial, theme_minimal,
#         cairo_pdf); logo letters use ggseqlogo's built-in glyphs.
# =============================================================================

suppressPackageStartupMessages({library(ggplot2); library(ggseqlogo); library(patchwork)})

# ---- choose the TFs to draw (row order of the figure); NULL = all candidates -------------
# HOCOMOCO v11: subtype labels from the PWM collection, not from logo reading. Rows and roles from
# subtype_verdicts_v11_h3.csv (build_Figure_S7_source_data_h3.py): v11 subtype examples STAT1 and TAL1, then own motif plus
# co-factor motif IRF4 and TFE2. Each tile names its nearest own-TF v11 model and r; a co-factor tile also names the nearest
# v11 model of any TF and its r.
verd <- read.csv(file.path(local({f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
                                  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()}), "subtype_verdicts_v11_h3.csv"),
                 stringsAsFactors = FALSE)
verd <- verd[verd$drawn == "True", ]
TFS <- verd$TF[order(verd$figure_row)]
PREVIEWS <- TRUE                  # also write preview/<TF>.png for every TF drawn
cli <- commandArgs(trailingOnly = TRUE)      # Rscript Figure_S7.R TAL1 RXRA  overrides TFS
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
LOGO_FONT <- "helvetica_bold"

d <- read.csv(src_path, stringsAsFactors = FALSE, check.names = FALSE)
all_tfs <- unique(d$TF)                       # ranking order from build_source_data.py
# panels per TF from the source data (builder parameters N_MONO, N_DI: 4 + 3); every TF must have the same numbers
.n <- table(d$TF, d$pwm_class)
if (length(unique(.n[, "mono"])) != 1 || length(unique(.n[, "di"])) != 1) stop("TFs differ in their number of mono / di panels")
N_MONO <- unname(.n[1, "mono"]); N_DI <- unname(.n[1, "di"]); N_COL <- N_MONO + N_DI
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

# cell type on its own subtitle line; long GTRD descriptions get a short display name here (plotting layer only, data untouched)
CELL_SHORT <- c("megakaryocytes cultured from cord blood CD34-positive cells" = "megakaryocytes (cord blood CD34+)",
                "Mesenchymal stem cells-derived adipocytes" = "adipocytes (from MSCs)",
                "CD133+ expanded umbilical cord blood (UCB) cells" = "CD133+ cord blood cells",
                "HeLa S3" = "HeLa-S3")   # one name for the cell line (GTRD writes both)
short_cell <- function(s, n = 34) {
  s <- ifelse(s %in% names(CELL_SHORT), CELL_SHORT[s], s)
  if (any(nchar(s) > n)) stop("CHECK FAILED Fig. S7: cell name longer than ", n, " characters without a display name: ", s[nchar(s) > n][1])
  unname(s)
}

logo_panel <- function(r) {
  m <- read_matrix(r$matrix_file)
  ttl <- sprintf("%s %d · %s.%s%s", ifelse(r$pwm_class == "mono", "mono", "di"),
                 r$rank_within_class, r$dataset, r$length_class,
                 ifelse(r$strand_for_display == "-", " (rc)", ""))
  cell <- short_cell(r$cell_line)
  # four subtitle lines in every tile (empty where there is no cell type or co-factor) so that all logos share one baseline
  sub <- sprintf("importance %.3f\n%s\nv11 %s, r %.3f\n%s", r$importance, ifelse(cell == "", " ", cell),
                 r$best_own_v11, r$r_best_own,
                 ifelse(r$cofactor == "True", sprintf("any TF: %s, r %.3f", r$best_any_v11, r$r_best_any), " "))
  ggplot() +
    geom_logo(m, method = "bits", seq_type = "dna", font = LOGO_FONT) +
    theme_logo(base_family = "Arial") +
    coord_cartesian(xlim = c(0.5, max(ncol(m), 12) + 0.5), ylim = c(0, 2), expand = FALSE) +   # own length, at least 12 bp wide
    scale_y_continuous(breaks = c(0, 1, 2)) +
    labs(title = ttl, subtitle = sub, x = NULL, y = NULL) +
    theme(plot.title = element_text(size = 9.5, hjust = 0, margin = margin(b = 1)),
          plot.subtitle = element_text(size = 7.5, hjust = 0, colour = "black", lineheight = 1.0, margin = margin(b = 7)),
          axis.text.x = element_blank(), axis.text.y = element_text(size = 7.5, colour = "black"),
          panel.grid = element_blank(),
          plot.margin = margin(4, 16, 2, 2),
          legend.position = "none")
}

label_panel <- function(tf, r1, letter) {
  v <- verd[verd$TF == tf, ]
  role <- if (v$role == "v11 subtype example") "v11 subtypes\n" else "own motif +\nco-factor motif\n"
  txt <- sprintf("%s\nΔauROC %+.3f\nΔauPRC %+.3f\nhuman test set\n%s", tf, r1$delta_auROC, r1$delta_auPRC, role)   # test set named
  ggplot() + annotate("text", x = 0, y = 0.76, label = txt, hjust = 0, vjust = 1,
                      size = 4.0, family = "Arial", lineheight = 0.95) +
    annotate("text", x = 0, y = 1.02, label = letter, hjust = 0, vjust = 1, size = 6.5, fontface = "bold", family = "Arial") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") + theme_void() +
    theme(plot.margin = margin(4, 2, 6, 6))
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
           width = round(400 * (N_COL + 0.62) / 6.62), height = 42, units = "mm", dpi = 300, bg = "white")
  }
}

# ---- the figure -------------------------------------------------------------
fig <- assemble(TFS)   # no title inside the figure (spec rule 8)

H <- 4 + 48 * length(TFS) + 6
fig <- fig + plot_annotation(theme = theme(plot.margin = margin(0, 0, 8, 0, "mm")))   # room for the longer label of the last row
H <- H + 8
W <- round(400 * (N_COL + 0.62) / 6.62)   # 400 mm for six logo columns; wider when N_DI grows, so logos keep their size
pdf_out <- file.path(script_dir, "Figure_S7.pdf")
png_out <- file.path(script_dir, "Figure_S7.png")
ggsave(pdf_out, fig, device = grDevices::cairo_pdf, width = W, height = H, units = "mm", limitsize = FALSE)
ggsave(png_out, fig, width = W, height = H, units = "mm", dpi = 300, bg = "white", limitsize = FALSE)
# ggsave without ragg writes no pHYs chunk; re-save with a 300 dpi tag,
# pixel values unchanged (8-bit RGBA read and written back exactly)
set_png_dpi <- function(f, dpi = 300) png::writePNG(png::readPNG(f), f, dpi = dpi)
set_png_dpi(png_out)

cat("TFs:", paste(TFS, collapse = ", "), "\n")
cat("Panels:", nrow(d), "| bp scale 1-", x_max, "\n")
cat("Written:", pdf_out, "\n         ", png_out, "\n")
