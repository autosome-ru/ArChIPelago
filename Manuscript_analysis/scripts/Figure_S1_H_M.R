# =============================================================================
# Figure_S1_H_M.R -- Fig. S1: auROC (A) and auPRC (B) of the five algorithms on the monoPWM, diPWM
# and monoPWM+diPWM feature sets, mouse test set (chr1, 8, 19), with the best single monoPWM and the best
# single diPWM (selected on the human training set) as reference series (one dot per TF) and as
# horizontal lines (median over TFs).
#
# Input : ../HUMAN_MOUSE_total_100k.csv (make_results_table.py) and ../analysis/inputs/df_Names_seq_list.csv
# Usage : Rscript Figure_S1_H_M.R [input_table] [output_dir]     (default output: ../Figures/panels)
# Writes: Figure_S1_panels_mouse_test.pdf
# =============================================================================

NAMES_INPUT <- "df_Names_seq_list.csv"
OUT_SUFFIX  <- ""

# ---- locate project dir / parse optional CLI overrides ----------------------
script_dir <- local({
  f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})
PROJECT_DIR <- normalizePath(file.path(script_dir, ".."))   # Manuscript_analysis/

args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1) args[1] else file.path(PROJECT_DIR, "HUMAN_MOUSE_total_100k.csv")
names_path <- file.path(PROJECT_DIR, "analysis", "inputs", NAMES_INPUT)
out_dir    <- if (length(args) >= 2) args[2] else file.path(PROJECT_DIR, "Figures", "panels")
if (!file.exists(input_path)) stop("Input table not found: ", input_path)
if (!file.exists(names_path)) stop("Names table not found: ", names_path)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- function(name) file.path(out_dir, sub("(\\.[^.]+)$", paste0(OUT_SUFFIX, "\\1"), name))
cat("Input table :", input_path, "\n")
cat("Names table :", names_path, "\n")
cat("Output dir  :", out_dir, "\n")
if (!interactive()) pdf(NULL)   # Rscript: swallow implicit print()s so no stray Rplots.pdf is written

# ---- packages (guarded so Rscript works without the full tidyverse) ---------
if (requireNamespace("tidyverse", quietly = TRUE)) {
  suppressPackageStartupMessages(library(tidyverse))
} else {
  suppressPackageStartupMessages(library(dplyr))   # %>%, left_join, filter, distinct, recode_factor
}
library(ggplot2)
library(patchwork)

# Read the main dataframe
df_total <- read.delim(input_path, sep="\t")

# Read the second dataframe
ddf <- read.delim(names_path, sep='\t')

# Merge dataframes
df_total <- df_total %>%
  left_join(ddf, by = c("TF_name" = "Names"))

# Filter by Seq_count.y (from the second dataframe)
df_total <- df_total %>%
  filter(Seq_count.y > 100)

# Remove duplicates; one row per TF x model x PWM set (the single-best reference rows are
# repeated once per algorithm in the table and differ only in their std_* columns)
df_total <- df_total %>%
  distinct(TF_name, Model, PWM, .keep_all = TRUE)

# Define hue order for consistent plotting
hue_order <- c("Single best mono PWM", "Single best di PWM", 
               "LogisticRegression", "RandomForestClassifier", 
               "XGBClassifier", "BaggingClassifier_XGBClassifier", 
               "BaggingClassifier_LogisticRegression")

# Convert Model to factor with specified levels
df_total$Model <- factor(df_total$Model, levels = hue_order)
df_total$PWM <- factor(df_total$PWM, levels = c("mono", "di", "mono+di"))

# Rename Model levels for display
df_total$Model <- recode_factor(df_total$Model,
                                "Single best mono PWM" = "Single best monoPWM",
                                "Single best di PWM" = "Single best diPWM",
                                "LogisticRegression" = "Logistic Regression",
                                "RandomForestClassifier" = "Random Forest (RF)",
                                "XGBClassifier" = "Gradient Boosting (XGBoost)",
                                "BaggingClassifier_XGBClassifier" = "Bagging (Gradient Boosting)",
                                "BaggingClassifier_LogisticRegression" = "Bagging (Logistic Regression)"
)

# Rename PWM levels for display
df_total$PWM <- recode_factor(df_total$PWM,
                              "mono" = "monoPWMs",
                              "di" = "diPWMs",
                              "mono+di" = "monoPWMs+diPWMs"
)

# Calculate medians for horizontal lines - ROC
median_mono_roc <- median(df_total$roc_auc_test_M[df_total$Model == "Single best monoPWM"], na.rm = TRUE)
median_di_roc <- median(df_total$roc_auc_test_M[df_total$Model == "Single best diPWM"], na.rm = TRUE)

# Calculate medians for horizontal lines - PRC
median_mono_prc <- median(df_total$pr_auc_test_M[df_total$Model == "Single best monoPWM"], na.rm = TRUE)
median_di_prc <- median(df_total$pr_auc_test_M[df_total$Model == "Single best diPWM"], na.rm = TRUE)

# ---- stdout diagnostics ------------------------------------------------------
cat("\n==== KEY MEDIANS (mouse test set, Seq_count.y > 100 filter) ====\n")
cat(sprintf("median_mono_roc = %.4f   median_di_roc = %.4f\n", median_mono_roc, median_di_roc))
cat(sprintf("median_mono_prc = %.4f   median_di_prc = %.4f\n", median_mono_prc, median_di_prc))
cat("Per Model x PWM medians (n = TFs):\n")
print(as.data.frame(df_total %>%
  group_by(Model, PWM) %>%
  summarise(n = n(),
            median_auROC = median(roc_auc_test_M, na.rm = TRUE),
            median_auPRC = median(pr_auc_test_M, na.rm = TRUE), .groups = "drop")), digits = 4)
cat("=================================================================\n\n")

# Panel A - auROC (keep all three: monoPWMs, diPWMs, monoPWMs+diPWMs)
p_roc <- ggplot(df_total, aes(x = PWM, y = roc_auc_test_M, color = Model)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              size = 3, alpha = 0.6) +
  geom_boxplot(aes(group = interaction(PWM, Model)), 
               fill = NA, color = "black", linewidth = 1, 
               outlier.shape = NA, position = position_dodge(0.8)) +
  geom_hline(yintercept = median_mono_roc, color = "#66c2a5", linetype = "dotdash", linewidth = 1) +
  geom_hline(yintercept = median_di_roc, color = "#fc8d62", linetype = "dotdash", linewidth = 1) +
  scale_color_brewer(palette = "Set2") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_x_discrete(limits = c("monoPWMs", "diPWMs", "monoPWMs+diPWMs")) +
  labs(x = NULL, y = "auROC", title = "A") +
  theme_classic(base_size = 20) +
  theme(
    axis.ticks.length = unit(0.2, "cm"),
    axis.text = element_text(color = "black", face = "plain"),
    axis.text.x = element_blank(),  # Remove x-axis text for top panel
    axis.title = element_text(face = "plain"),
    legend.position = "none",
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    plot.title = element_text(face = "bold", hjust = 0, size = 24)
  )

# Panel B - auPRC
p_prc <- ggplot(df_total, aes(x = PWM, y = pr_auc_test_M, color = Model)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              size = 3, alpha = 0.6) +
  geom_boxplot(aes(group = interaction(PWM, Model)), 
               fill = NA, color = "black", linewidth = 1, 
               outlier.shape = NA, position = position_dodge(0.8)) +
  geom_hline(yintercept = median_mono_prc, color = "#66c2a5", linetype = "dotdash", linewidth = 1) +
  geom_hline(yintercept = median_di_prc, color = "#fc8d62", linetype = "dotdash", linewidth = 1) +
  scale_color_brewer(palette = "Set2") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  labs(x = "PWM type", y = "auPRC", title = "B") +
  theme_classic(base_size = 20) +
  theme(
    axis.ticks.length = unit(0.2, "cm"),
    axis.text = element_text(color = "black", face = "plain"),
    axis.title = element_text(face = "plain"),
    legend.position = "none",
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    plot.title = element_text(face = "bold", hjust = 0, size = 24)
  )

# Combine panels
combined_plot <- p_roc / p_prc

# Save the combined plot with dpi=600
ggsave(out_file("Figure_S1_panels_mouse_test.pdf"), 
       plot = combined_plot, width = 10, height = 10, dpi = 600, device = cairo_pdf)

# Display the plot
print(combined_plot)

# Extract and save the legend separately
library(cowplot)

# Create a temporary plot WITH legend to extract it
p_with_legend <- ggplot(df_total, aes(x = PWM, y = roc_auc_test_M, color = Model)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              size = 3, alpha = 0.6) +
  scale_color_brewer(palette = "Set2") +
  theme_classic(base_size = 20) +
  theme(
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 16)
  )

# Extract the legend
legend <- get_legend(p_with_legend)

# Save the legend as a separate PDF with dpi=600
ggsave(out_file("Models_comparison_legend_21012026.pdf"), 
       plot = legend, width = 5, height = 4, dpi = 600, device = cairo_pdf)

print("Legend and combined plot saved!")
cat("Written:", out_file("Figure_S1_panels_mouse_test.pdf"), "\n")
cat("Written:", out_file("Models_comparison_legend_21012026.pdf"), "\n")
