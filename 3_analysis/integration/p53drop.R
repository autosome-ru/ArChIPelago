# p53drop.R - R side of 3_analysis/integration/p53drop.py.
# check_tfs(): stop unless the TF set equals the 35 reference TFs (TF column of results/datasets/peak_ids/mouse_sym_peak_ids.csv).
# ARCHI_DRAW: comma-separated species to draw ("H", "M"; default both), so that the human or the mouse figures can be redrawn alone.
.p53_dir <- local({
  f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  d <- if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
  d
})
REF_TF_FILE <- normalizePath(file.path(.p53_dir, "..", "..", "../results/datasets/peak_ids", "mouse_sym_peak_ids.csv"))
REF_TFS <- sort(unique(read.csv(REF_TF_FILE, stringsAsFactors = FALSE)$TF))
ARCHI_DRAW <- strsplit(Sys.getenv("ARCHI_DRAW", "H,M"), ",")[[1]]
check_tfs <- function(tfs, label) {
  got <- sort(unique(as.character(tfs)))
  if (!identical(got, REF_TFS))
    stop(sprintf("CHECK FAILED %s: %d TFs, extra [%s], missing [%s]", label, length(got),
                 paste(setdiff(got, REF_TFS), collapse = ","), paste(setdiff(REF_TFS, got), collapse = ",")))
  cat(sprintf("check %s: %d TFs, P53 absent\n", label, length(got)))
  length(got)
}
