# =============================================================================
# Figure_S3ABC_mouse_test.R -- Fig. S3, panels A-C: Random Forest on monoPWMs, diPWMs and monoPWMs+diPWMs against
# the best single monoPWM (selected on the human training set), mouse test set (chr1, 8, 19); C = gain of the
# monoPWM+diPWM model per TF (colour: number of PWMs, size: training positives).
#
# Input : ../results_table.csv (make_results_table.py)
# Usage : Rscript Figure_S3ABC_mouse_test.R [input_table] [output_dir]     (default output: ../Figures/panels)
# Writes: Figure_S3ABC_mouse_test.pdf
# Optional packages (hrbrthemes, gapminder, extrafont, viridis, ggpubr) are used when installed;
# otherwise ggplot2 / cowplot fallbacks are used.
# =============================================================================

OUT_SUFFIX  <- ""

# ---- optional CLI overrides -------------------------------------------------
script_dir <- local({
  f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})
args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1) args[1] else file.path(script_dir, "..", "results_table.csv")
out_dir    <- if (length(args) >= 2) args[2] else file.path(script_dir, "..", "Figures", "panels")
if (!file.exists(input_path)) stop("Input table not found: ", input_path)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- function(name) file.path(out_dir, sub("(\\.[^.]+)$", paste0(OUT_SUFFIX, "\\1"), name))
cat("Input table :", input_path, "\n")
cat("Output dir  :", out_dir, "\n")
if (!interactive()) pdf(NULL)   # Rscript: swallow implicit print()s so no stray Rplots.pdf is written

# Libraries
library(ggplot2)
library(dplyr)
if (requireNamespace("hrbrthemes", quietly = TRUE)) library(hrbrthemes)   # not used by the plots
if (requireNamespace("viridis", quietly = TRUE)) {
  library(viridis)
} else {
  message("viridis not installed: using ggplot2::scale_fill_viridis_c/_d fallback")
  scale_fill_viridis <- function(..., discrete = FALSE)
    if (discrete) ggplot2::scale_fill_viridis_d(...) else ggplot2::scale_fill_viridis_c(...)
}
if (requireNamespace("ggpubr", quietly = TRUE)) {
  library(ggpubr)
} else {
  message("ggpubr not installed: using cowplot::plot_grid fallback for ggarrange()")
  library(cowplot)
  ggarrange <- function(..., nrow = NULL, ncol = NULL, labels = NULL, font.label = list(size = 14),
                        vjust = NULL, widths = NULL, heights = NULL)
    cowplot::plot_grid(..., nrow = nrow, ncol = ncol, labels = labels,
                       label_size = font.label$size, rel_widths = if (is.null(widths)) 1 else widths,
                       rel_heights = if (is.null(heights)) 1 else heights)
  text_grob <- function(label, ...) grid::textGrob(label)
  annotate_figure <- function(p, top = NULL, ...) p
}
library(ggrepel)
if (requireNamespace("gapminder", quietly = TRUE)) library(gapminder)     # not used by the plots
if (requireNamespace("ggExtra", quietly = TRUE)) {
  library(ggExtra)
} else {
  message("ggExtra not installed: marginal histograms (ggMarginal) omitted")
  ggMarginal <- function(p, ...) p
}
if (requireNamespace("extrafont", quietly = TRUE)) {
  library(extrafont)
  # use the fonts already imported by extrafont (font_import() is interactive and slow)
  try(suppressMessages(loadfonts(quiet = TRUE)), silent = TRUE)
} else {
  message("extrafont not installed: skipping font import (device default fonts used)")
}


theme_set(theme_get() + theme(text = element_text(family = 'Arial')))


##### CHS ####

Model_key_list=c("RandomForestClassifier", "LogisticRegression", "XGBClassifier", "BaggingClassifier_XGBClassifier", "BaggingClassifier_LogisticRegression")




my_pal <- function(range = c(1, 6)) {
  force(range)
  function(x) scales::rescale(x, to = range, from = c(0, 1))
}

# 
# 
# data1 = data %>% 
#   dplyr::filter(Model=="RandomForestClassifier") %>% 
#   select(c("TF_name", 
#            "roc_auc_test_M_PWM", "roc_auc_test_M",
#            "pr_auc_test_M_PWM", "pr_auc_test_M", "Count", "Seq_count"))
# 
# data = data1 %>% 
#   mutate(ROC_delta_M = roc_auc_test_M-roc_auc_test_M_PWM, PR_delta_M = pr_auc_test_M-pr_auc_test_M_PWM) %>% 
#   group_by(TF_name) %>% 
#   summarise(ROC_delta_M=mean(ROC_delta_M), PR_delta_M=mean(PR_delta_M), Count_PWM=sum(unique(Count)), Seq_count_min=min(Seq_count)) %>% 
#   ungroup()



data <- read.csv(input_path, sep="\t")
data1 <- data %>% 
  dplyr::filter(Model == "RandomForestClassifier") %>% 
  select(
    TF_name, PWM,
    roc_auc_test_M_PWM, roc_auc_test_M,
    pr_auc_test_M_PWM, pr_auc_test_M,
    Count, Seq_count
  )
data2 <- data1 %>% 
  group_by(TF_name) %>% 
  mutate(Count_PWM = sum(unique(Count))) %>% 
  ungroup()
data3 <- data2 %>% 
  filter(PWM=="mono+di")
data <- data3 %>% 
  mutate(
    ROC_delta_M = roc_auc_test_M - roc_auc_test_M_PWM,
    PR_delta_M  = pr_auc_test_M  - pr_auc_test_M_PWM
  ) %>% 
  group_by(TF_name) %>% 
  summarise(
    ROC_delta_M = mean(ROC_delta_M, na.rm = TRUE),
    PR_delta_M  = mean(PR_delta_M,  na.rm = TRUE),
    Count_PWM  = first(Count_PWM),
    Seq_count_min = min(Seq_count, na.rm = TRUE)
  ) %>% 
  ungroup()


c = data %>%
  arrange(desc(Seq_count_min)) %>%
  mutate(TF_name = factor(TF_name)) %>%
  ggplot(aes(x=PR_delta_M, y=ROC_delta_M, size=Seq_count_min, fill=Count_PWM)) + # 
  geom_point(alpha=0.5, shape=21, color="black") +
  scale_size(range = c(1, 15), breaks= c(500, 1000, 5000, 10000, 50000),
             #labels = c(5*10^2, 10^3, 5*10^3, 10^4, 5*10^4),
             name="Size of \nthe positive \nset") +
  scale_fill_viridis(discrete=F, option="C", name="Number of PWMs") +
  theme_minimal(base_size = 21) +
  theme(legend.position="right",
        legend.box="vertical", 
        legend.margin=margin())+
  ylab("\u0394auROC") +
  xlab("\u0394auPRC") + 
  geom_hline(yintercept = 0, linetype="dotted") + 
  geom_vline(xintercept = 0, linetype="dotted") + 
  geom_text_repel(data = subset(data, (ROC_delta_M < 0)|(ROC_delta_M>0.05)|(PR_delta_M<0)|(PR_delta_M>0.2)),
                  aes(label = TF_name), 
                  point.padding = 1,
                  min.segment.length = 3,
                  max.time = 1, max.iter = 1e5, seed = 1,
                  #box.padding = 0.3, 
                  segment.curvature = -0.1,
                  #segment.ncp = 3,
                  segment.angle = 20,
                  size=5)+
  guides(color = guide_legend(order=1),
         size = guide_legend(order=2),
         shape = "none")



data <- read.csv(input_path, sep="\t")

data = data %>% 
  dplyr::filter(Model=="RandomForestClassifier") %>% 
  select(c("TF_name", 
           "roc_auc_test_M_PWM", "roc_auc_test_M",
           "pr_auc_test_M_PWM", "pr_auc_test_M", "Count", "Seq_count", "PWM"))

a = data %>%
  arrange(desc(Seq_count)) %>%
  ggplot(aes(x=roc_auc_test_M_PWM, y=roc_auc_test_M, color=PWM)) + # 
  geom_point(alpha=0.5, shape=20, size=3)+
  theme_minimal(base_size = 21) +
  theme(legend.position="none",
        axis.text.x=element_text(angle=90, vjust = 0.5))+
  ylab("auROC") +
  xlab("auROC best PWM") + 
  expand_limits(x=0, y=0) +
  geom_abline(intercept = 0, linetype="dotted") +
  guides(color=guide_legend(nrow=3, byrow=TRUE)) +
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1)
  )

a = ggMarginal(a, type = "histogram", 
               #margins = "x",
               #color = "gray",
               fill = "white", 
               alpha=0.5,
               size=8, groupColour = TRUE)

b = data %>%
  arrange(desc(Seq_count))  %>%
  ggplot(aes(x=pr_auc_test_M_PWM, y=pr_auc_test_M, color=PWM)) + # 
  geom_point(alpha=0.5, shape=20, size=3)+
  theme_minimal(base_size = 21) +
  theme(legend.position="none",
        axis.text.x=element_text(angle=90, vjust = 0.5))+
  ylab("auPRC") +
  xlab("auPRC best PWM") + 
  expand_limits(x=0, y=0) +
  geom_abline(intercept = 0, linetype="dotted") +
  guides(color=guide_legend(nrow=3, byrow=TRUE)) +
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1)
  )

b = ggMarginal(b, type = "histogram", 
               #margins = "x",
               #color = "gray",
               fill = "white", 
               alpha=0.5,
               size=8, groupColour = TRUE)

plot = ggarrange(
  
  ggarrange(a, b, nrow = 2, ncol = 1, labels = c("A", "B"),
            font.label=list(size=25), vjust = -0.3 
            #heights = c(1, 1.5)
  ),
  
  ggarrange(c, 
            labels = c("C"),
            ncol = 1, nrow = 1, 
            font.label=list(size=25), vjust = -0.3),  
  widths = c(1, 2.5), ncol = 2, nrow = 1
)

image = annotate_figure(plot, top = text_grob("", 
                                      color = "black", face = "bold", size = 25))

pdf_out <- out_file('Figure_S3ABC_mouse_test.pdf')
if (requireNamespace("Cairo", quietly = TRUE)) {
  library(Cairo)
  Cairo(file=pdf_out, type="pdf", width=310, height=185, units="mm")
} else {
  message("Cairo package not installed: using grDevices::cairo_pdf")
  cairo_pdf(file=pdf_out, width=310/25.4, height=185/25.4)
}
print(image)
dev.off()
cat("Written:", pdf_out, "\n")


data4 <- data3 %>% 
  mutate(
    ROC_delta_M = roc_auc_test_M - roc_auc_test_M_PWM,
    PR_delta_M  = pr_auc_test_M  - pr_auc_test_M_PWM
  )
median_comparison <- data4 %>% 
  summarise(
    median_roc_M      = median(roc_auc_test_M, na.rm = TRUE),
    median_roc_M_PWM  = median(roc_auc_test_M_PWM, na.rm = TRUE),
    median_pr_M       = median(pr_auc_test_M, na.rm = TRUE),
    median_pr_M_PWM   = median(pr_auc_test_M_PWM, na.rm = TRUE)
  )
median_gain <- data4 %>% 
  summarise(
    median_ROC_delta_M = median(ROC_delta_M, na.rm = TRUE),
    median_PR_delta_M  = median(PR_delta_M,  na.rm = TRUE)
  )
median_summary <- bind_cols(median_comparison, median_gain)

# ---- stdout diagnostics ------------------------------------------------------
cat("\n==== KEY MEDIANS: RandomForestClassifier, mono+di, mouse test set (n TFs =", nrow(data4), ") ====\n")
print(as.data.frame(median_summary), digits = 4)
cat("Per PWM-type medians (RandomForestClassifier, mouse test set):\n")
print(as.data.frame(data2 %>%
  group_by(PWM) %>%
  summarise(n = n(),
            median_roc_M = median(roc_auc_test_M, na.rm = TRUE),
            median_roc_M_PWM = median(roc_auc_test_M_PWM, na.rm = TRUE),
            median_pr_M = median(pr_auc_test_M, na.rm = TRUE),
            median_pr_M_PWM = median(pr_auc_test_M_PWM, na.rm = TRUE), .groups = "drop")), digits = 4)
cat("=================================================================================\n")
