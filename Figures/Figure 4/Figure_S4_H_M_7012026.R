

library(tidyverse)
library(ggbeeswarm)
library(patchwork)

# 1. Load the baseline data
df_total <- read_tsv("./HUMAN_MOUSE_total_100k.csv")





####



# Load required libraries
library(tidyverse)
library(ggbeeswarm)
library(patchwork)

# --- 1. Data Processing Function ---
get_delta_data <- function(file_path, df_total, metric_type) {
  
  # Define target variables as they appear in the final plot
  target_vars <- c(
    "Test SLIM m=0 MOUSE", "Test SLIM m=1 MOUSE", "Test SLIM m=-5 MOUSE", 
    "Test diChIPMunk MOUSE", "Test RandomForest MOUSE", "Test RF+diChIPMunk MOUSE", 
    "Test RF+SLIM m=1 MOUSE", "Test RF+SLIM m=-5 MOUSE", 
    "Test RF+SLIM m=1,-5 +diChIPMunk MOUSE", "Test full RandomForest MOUSE"
  )
  
  # Read the first line to count columns accurately
  first_line <- read_lines(file_path, n_max = 1)
  col_count <- length(str_split(first_line, "\t")[[1]])
  
  # Create unique name placeholders
  custom_names <- paste0("V", 1:col_count)
  
  # Mapping based on your specific file indices
  custom_names[1]  <- "TF_name"
  custom_names[11] <- "Test SLIM m=0 MOUSE"
  custom_names[12] <- "Test SLIM m=1 MOUSE"
  custom_names[13] <- "Test SLIM m=-5 MOUSE"
  custom_names[14] <- "Test diChIPMunk MOUSE"
  custom_names[25] <- "Test full RandomForest MOUSE"
  custom_names[26] <- "Test RandomForest MOUSE"
  custom_names[27] <- "Test RF+diChIPMunk MOUSE"
  custom_names[28] <- "Test RF+SLIM m=1 MOUSE"
  custom_names[29] <- "Test RF+SLIM m=-5 MOUSE"
  custom_names[30] <- "Test RF+SLIM m=1,-5 +diChIPMunk MOUSE"
  
  # Read file - ensuring tab separator
  slim_df <- read_tsv(file_path, col_names = custom_names, show_col_types = FALSE)
  
  # Clean TF_name (remove last 6 characters)
  slim_df <- slim_df %>%
    mutate(TF_name = str_sub(TF_name, 1, -7))
  
  # Identify baseline column from df_total
  metric_col <- ifelse(metric_type == "ROC", "roc_auc_test_M_PWM_mono", "pr_auc_test_M_PWM_mono")
  
  baseline_df <- df_total %>%
    filter(PWM == "mono", Model == "Single best mono PWM") %>%
    select(TF_name, baseline = !!sym(metric_col)) %>%
    distinct()
  
  # Merge and calculate Delta
  final_df <- slim_df %>%
    inner_join(baseline_df, by = "TF_name") %>%
    mutate(across(all_of(target_vars), ~ . - baseline)) %>%
    select(TF_name, all_of(target_vars)) %>%
    pivot_longer(-TF_name, names_to = "variable", values_to = "value") %>%
    mutate(variable = factor(variable, levels = rev(target_vars)))
  
  return(final_df)
}




###


# --- 1. Updated Plotting Function ---
create_panel <- function(data, title, x_label, show_y_labels = TRUE) {
  # Calculate Wilcoxon p-values for stars
  stats <- data %>%
    group_by(variable) %>%
    summarise(
      p = wilcox.test(value, alternative = "greater")$p.value
    ) %>%
    mutate(lbl = case_when(
      p < 0.001 ~ "***", 
      p < 0.01 ~ "**", 
      p < 0.05 ~ "*", 
      TRUE ~ ""
    ))
  
  p <- ggplot(data, aes(x = value, y = variable)) +
    geom_quasirandom(color = "gray80", size = 1, alpha = 0.8, groupOnX = FALSE) +
    geom_boxplot(fill = NA, color = "black", outlier.shape = NA) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "darkred", alpha = 1) +
    # Asterisks to the left of the 0 line
    geom_text(data = stats, aes(x = -0.1, label = lbl), 
              color = "red", size = 5, hjust = 1, vjust = 0.3) +
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
    # Label for Panel A
    p <- p + labs(y = "Model")
  } else {
    # Complete removal of Y axis elements for Panel B
    p <- p + theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(), 
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank()
    ) + labs(y = NULL)
  }
  
  return(p)
}

# --- 2. Execution with new labels ---

# Process data (assuming df_total and get_delta_data are already defined)
df_roc_plot <- get_delta_data("./HUMAN_MOUSE_SLIM_roc_mono_di_RandomForestClassifier.txt", df_total, "ROC")
df_pr_plot  <- get_delta_data("./HUMAN_MOUSE_SLIM_pr_mono_di_RandomForestClassifier.txt", df_total, "PR")

# Apply the requested panel titles and expression labels
panel_a <- create_panel(df_roc_plot, "A", expression(Delta ~ "auROC"))
panel_b <- create_panel(df_pr_plot, "B", expression(Delta ~ "auPRC"), show_y = FALSE)

combined <- panel_a + panel_b + plot_layout(widths = c(1.5, 1))
combined
# --- 3. Save ---
ggsave("./Figure S4.pdf", 
       combined, 
       device = grDevices::cairo_pdf,
       width = 8, 
       height = 6)
