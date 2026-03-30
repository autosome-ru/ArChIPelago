# for absolute values


library(tidyverse)
library(ggbeeswarm)
library(patchwork)

# 1. Load the baseline data
df_total <- read_tsv("./HUMAN_MOUSE_total_100k.csv")

# --- 1. Updated Data Processing Function for Absolute Values ---
get_absolute_data <- function(file_path, df_total, metric_type) {
  
  # Define target variables (Test models)
  test_vars <- c(
    "Test SLIM m=0 HUMAN", "Test SLIM m=1 HUMAN", "Test SLIM m=-5 HUMAN", 
    "Test diChIPMunk HUMAN", "Test RandomForest HUMAN", "Test RF+diChIPMunk HUMAN", 
    "Test RF+SLIM m=1 HUMAN", "Test RF+SLIM m=-5 HUMAN", 
    "Test RF+SLIM m=1,-5 +diChIPMunk HUMAN", "Test full RandomForest HUMAN"
  )
  
  # Read file with custom mapping
  first_line <- read_lines(file_path, n_max = 1)
  col_count <- length(str_split(first_line, "\t")[[1]])
  custom_names <- paste0("V", 1:col_count)
  
  custom_names[1]  <- "TF_name"
  custom_names[11] <- "Test SLIM m=0 HUMAN"
  custom_names[12] <- "Test SLIM m=1 HUMAN"
  custom_names[13] <- "Test SLIM m=-5 HUMAN"
  custom_names[14] <- "Test diChIPMunk HUMAN"
  custom_names[25] <- "Test full RandomForest HUMAN"
  custom_names[26] <- "Test RandomForest HUMAN"
  custom_names[27] <- "Test RF+diChIPMunk HUMAN"
  custom_names[28] <- "Test RF+SLIM m=1 HUMAN"
  custom_names[29] <- "Test RF+SLIM m=-5 HUMAN"
  custom_names[30] <- "Test RF+SLIM m=1,-5 +diChIPMunk HUMAN"
  
  slim_df <- read_tsv(file_path, col_names = custom_names, show_col_types = FALSE)
  slim_df <- slim_df %>% mutate(TF_name = str_sub(TF_name, 1, -7))
  
  # Get baseline column
  metric_col <- ifelse(metric_type == "ROC", "roc_auc_test_H_PWM_mono", "pr_auc_test_H_PWM_mono")
  
  baseline_df <- df_total %>%
    filter(PWM == "mono", Model == "Single best mono PWM") %>%
    select(TF_name, `Baseline (Mono PWM)` = !!sym(metric_col)) %>%
    distinct()
  
  # Merge. We keep absolute values (no subtraction)
  # We include the baseline as a variable to plot
  all_vars <- c("Baseline (Mono PWM)", test_vars)
  
  final_df <- slim_df %>%
    inner_join(baseline_df, by = "TF_name") %>%
    select(TF_name, all_of(all_vars)) %>%
    pivot_longer(-TF_name, names_to = "variable", values_to = "value") %>%
    mutate(variable = factor(variable, levels = rev(all_vars)))
  
  return(final_df)
}

# --- 2. Updated Plotting Function ---
create_panel_absolute <- function(data, title, x_label, show_y_labels = TRUE) {
  
  # 1. Prepare Baseline reference for the paired test
  baseline_subset <- data %>%
    filter(variable == "Baseline (Mono PWM)") %>%
    arrange(TF_name)
  
  baseline_values <- baseline_subset$value
  
  # Calculate the central tendency of the baseline to draw the vertical line
  # Using mean here, but you can change to median if preferred
  baseline_intercept <- median(baseline_values, na.rm = TRUE)
  
  # 2. Calculate Paired Wilcoxon p-values
  stats <- data %>%
    filter(variable != "Baseline (Mono PWM)") %>%
    group_by(variable) %>%
    summarise(
      p = {
        current_values <- value[order(TF_name)]
        # Ensure lengths match before testing
        if(length(current_values) == length(baseline_values)) {
          wilcox.test(current_values, baseline_values, paired = TRUE, alternative = "greater")$p.value
        } else {
          NA
        }
      },
      .groups = "drop"
    ) %>%
    mutate(lbl = case_when(
      p < 0.001 ~ "***", 
      p < 0.01 ~ "**", 
      p < 0.05 ~ "*", 
      TRUE ~ ""
    ))
  
  # 3. Plotting
  p <- ggplot(data, aes(x = value, y = variable)) +
    # Add the vertical line representing the average baseline performance
    geom_vline(xintercept = baseline_intercept, linetype = "dashed", color = "darkred", alpha = 1) +
    geom_quasirandom(color = "gray80", size = 1, alpha = 0.8, groupOnX = FALSE) +
    geom_boxplot(fill = NA, color = "black", outlier.shape = NA) +
    # Position stars at the far right
    geom_text(data = stats, aes(x = max(data$value + 0.05, na.rm = TRUE), label = lbl), 
              color = "red", size = 5, hjust = 1.1, vjust = 0.3) +
    labs(title = title, x = x_label) + 
    theme_classic() + theme(
      text = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(color = "black", face = "bold"),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black")
    )
  
  if (show_y_labels) {
    p <- p + labs(y = "Model")
  } else {
    p <- p + theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(), 
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank()
    ) + labs(y = NULL)
  }
  
  return(p)
}

# --- Re-run the processing and plotting ---

df_roc_plot <- get_absolute_data("./HUMAN_MOUSE_SLIM_roc_mono_di_RandomForestClassifier.txt", df_total, "ROC")
df_pr_plot  <- get_absolute_data("./HUMAN_MOUSE_SLIM_pr_mono_di_RandomForestClassifier.txt", df_total, "PR")

panel_a <- create_panel_absolute(df_roc_plot, "A", "auROC")
panel_b <- create_panel_absolute(df_pr_plot, "B", "auPRC", show_y_labels = FALSE)

combined <- panel_a + panel_b + plot_layout(widths = c(1.5, 1))
combined
# Save
ggsave("./Figure 4 abs. values.pdf", 
       combined, 
       device = grDevices::cairo_pdf,
       width = 8, 
       height = 6)

