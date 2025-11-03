# Create Comprehensive Visualizations for ARIMA Study
# This script generates all plots for the research report


# Set seed for reproducibility
set.seed(12345)

get_script_dir <- function() {
  frame_files <- vapply(sys.frames(), function(f) {
    if (!is.null(f$ofile)) f$ofile else NA_character_
  }, character(1), USE.NAMES = FALSE)
  frame_files <- frame_files[!is.na(frame_files)]
  if (length(frame_files) > 0) {
    return(dirname(normalizePath(frame_files[1])))
  }

  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- cmd_args[grepl("^--file=", cmd_args)]
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }

  normalizePath(getwd())
}

script_dir <- get_script_dir()
project_root <- normalizePath(file.path(script_dir, ".."))
original_wd <- getwd()
on.exit(setwd(original_wd), add = TRUE)
setwd(project_root)

library(ggplot2)
library(gridExtra)
library(RColorBrewer)

# Load results
output_dir <- "results"
if (!dir.exists(output_dir)) {
  stop("Results directory not found. Run scripts/comprehensive_study.R first.")
}

results <- readRDS(file.path(output_dir, "all_results.rds"))
fitted_models <- readRDS(file.path(output_dir, "fitted_models.rds"))

# Create plots directory
plots_dir <- file.path(output_dir, "plots")
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}

cat("Creating visualizations...\n\n")

# Color palette
colors <- brewer.pal(3, "Set1")
method_colors <- c("CSS-ML" = colors[1], "PMM2" = colors[2])

# ==================== PLOT 1: AIC Comparison ====================
cat("1. AIC Comparison... ")
p1 <- ggplot(results, aes(x = Model, y = AIC, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = method_colors) +
  labs(title = "AIC Comparison: PMM2 vs CSS-ML",
       subtitle = "Lower AIC indicates better model fit",
       x = "Model Specification",
       y = "AIC") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        legend.position = "top")

ggsave(file.path(plots_dir, "01_aic_comparison.png"), p1,
       width = 10, height = 6, dpi = 300)
cat("✓\n")

# ==================== PLOT 2: BIC Comparison ====================
cat("2. BIC Comparison... ")
p2 <- ggplot(results, aes(x = Model, y = BIC, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = method_colors) +
  labs(title = "BIC Comparison: PMM2 vs CSS-ML",
       subtitle = "Lower BIC indicates better model fit with parsimony",
       x = "Model Specification",
       y = "BIC") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        legend.position = "top")

ggsave(file.path(plots_dir, "02_bic_comparison.png"), p2,
       width = 10, height = 6, dpi = 300)
cat("✓\n")

# ==================== PLOT 3: RMSE Comparison ====================
cat("3. RMSE Comparison... ")
p3 <- ggplot(results, aes(x = Model, y = RMSE, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = method_colors) +
  labs(title = "RMSE Comparison: PMM2 vs CSS-ML",
       subtitle = "Root Mean Square Error of residuals",
       x = "Model Specification",
       y = "RMSE") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        legend.position = "top")

ggsave(file.path(plots_dir, "03_rmse_comparison.png"), p3,
       width = 10, height = 6, dpi = 300)
cat("✓\n")

# ==================== PLOT 4: Computation Time ====================
cat("4. Computation Time... ")
p4 <- ggplot(results, aes(x = Model, y = Time_sec, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = method_colors) +
  labs(title = "Computation Time Comparison",
       subtitle = "Time taken to fit each model (seconds)",
       x = "Model Specification",
       y = "Time (seconds)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        legend.position = "top")

ggsave(file.path(plots_dir, "04_computation_time.png"), p4,
       width = 10, height = 6, dpi = 300)
cat("✓\n")

# ==================== PLOT 5: Residual Kurtosis ====================
cat("5. Residual Kurtosis... ")
p5 <- ggplot(results, aes(x = Model, y = Kurtosis, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = method_colors) +
  labs(title = "Excess Kurtosis of Residuals",
       subtitle = "Measures tail heaviness (0 = normal distribution)",
       x = "Model Specification",
       y = "Excess Kurtosis") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        legend.position = "top")

ggsave(file.path(plots_dir, "05_kurtosis_comparison.png"), p5,
       width = 10, height = 6, dpi = 300)
cat("✓\n")

# ==================== PLOT 6: Residual Skewness ====================
cat("6. Residual Skewness... ")
p6 <- ggplot(results, aes(x = Model, y = Skewness, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = method_colors) +
  labs(title = "Skewness of Residuals",
       subtitle = "Measures asymmetry (0 = symmetric distribution)",
       x = "Model Specification",
       y = "Skewness") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        legend.position = "top")

ggsave(file.path(plots_dir, "06_skewness_comparison.png"), p6,
       width = 10, height = 6, dpi = 300)
cat("✓\n")

# ==================== PLOT 7: Performance Heatmap ====================
cat("7. Performance Heatmap... ")

# Normalize metrics to [0, 1] for comparison
results_norm <- results
for (metric in c("AIC", "BIC", "RMSE", "MAE")) {
  results_norm[[paste0(metric, "_norm")]] <-
    (results[[metric]] - min(results[[metric]])) /
    (max(results[[metric]]) - min(results[[metric]]))
}

# Create long format for heatmap without relying on tidyr
metrics <- c("AIC", "BIC", "RMSE", "MAE")
norm_columns <- paste0(metrics, "_norm")
results_subset <- results_norm[, c("Model", "Method", norm_columns)]

results_long <- do.call(rbind, lapply(seq_along(metrics), function(i) {
  data.frame(
    Model = results_subset$Model,
    Method = results_subset$Method,
    Metric = metrics[i],
    Value = results_subset[[norm_columns[i]]],
    stringsAsFactors = FALSE
  )
}))

results_long$Method_Model <- paste(results_long$Method, results_long$Model, sep = "\n")

p7 <- ggplot(results_long, aes(x = Metric, y = Method_Model, fill = Value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "green", mid = "yellow", high = "red",
                       midpoint = 0.5, name = "Normalized\nValue") +
  labs(title = "Performance Heatmap (Normalized Metrics)",
       subtitle = "Green = Better, Red = Worse",
       x = "Metric",
       y = "Method & Model") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 10))

ggsave(file.path(plots_dir, "07_performance_heatmap.png"), p7,
       width = 8, height = 10, dpi = 300)
cat("✓\n")

# ==================== PLOT 8: Method Difference Plot ====================
cat("8. Method Differences... ")

comparison <- read.csv(file.path(output_dir, "method_comparison.csv"))

# Reshape for plotting without tidyr
diff_metrics <- c("AIC", "BIC", "RMSE")
diff_columns <- paste0(diff_metrics, "_Diff")

comparison_long <- do.call(rbind, lapply(seq_along(diff_metrics), function(i) {
  data.frame(
    Model = comparison$Model,
    Metric = diff_metrics[i],
    Difference = comparison[[diff_columns[i]]],
    stringsAsFactors = FALSE
  )
}))

p8 <- ggplot(comparison_long, aes(x = Model, y = Difference, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "PMM2 vs CSS-ML: Performance Differences",
       subtitle = "Negative values indicate PMM2 is better",
       x = "Model Specification",
       y = "Difference (PMM2 - CSS-ML)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        legend.position = "top")

ggsave(file.path(plots_dir, "08_method_differences.png"), p8,
       width = 10, height = 6, dpi = 300)
cat("✓\n")

# ==================== PLOT 9: Best Model Diagnostics ====================
cat("9. Best Model Diagnostics... ")

# Find best model
best_bic <- results[which.min(results$BIC), ]
best_model_key <- paste0(best_bic$Model, "_", gsub("-", "", best_bic$Method))
best_fit <- fitted_models[[best_model_key]]

if (!is.null(best_fit)) {
  # Extract residuals
  if (best_bic$Method == "CSS-ML") {
    res <- residuals(best_fit)
  } else {
    res <- best_fit@residuals
  }
  res_clean <- res[is.finite(res)]

  # Create diagnostic plots
  png(file.path(plots_dir, "09_best_model_diagnostics.png"),
      width = 12, height = 10, units = "in", res = 300)

  par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))

  # 1. Residuals over time
  plot(res_clean, type = "l", col = "steelblue",
       main = paste("Residuals:", best_bic$Model, "-", best_bic$Method),
       xlab = "Time", ylab = "Residuals")
  abline(h = 0, col = "red", lty = 2)

  # 2. Histogram
  hist(res_clean, breaks = 50, col = "lightblue", border = "white",
       main = "Histogram of Residuals",
       xlab = "Residuals", freq = FALSE)
  curve(dnorm(x, mean = mean(res_clean), sd = sd(res_clean)),
        add = TRUE, col = "red", lwd = 2)

  # 3. Q-Q plot
  qqnorm(res_clean, main = "Normal Q-Q Plot", col = "steelblue")
  qqline(res_clean, col = "red", lwd = 2)

  # 4. ACF
  acf(res_clean, main = "ACF of Residuals", lag.max = 40)

  dev.off()
  cat("✓\n")
}

# ==================== PLOT 10: Summary Statistics ====================
cat("10. Summary Statistics... ")

summary_stats <- aggregate(cbind(AIC, BIC, RMSE, Time_sec) ~ Method,
                           data = results, FUN = mean)
names(summary_stats)[names(summary_stats) == "AIC"] <- "Mean_AIC"
names(summary_stats)[names(summary_stats) == "BIC"] <- "Mean_BIC"
names(summary_stats)[names(summary_stats) == "RMSE"] <- "Mean_RMSE"
names(summary_stats)[names(summary_stats) == "Time_sec"] <- "Mean_Time"

# Create multi-panel summary plot
p10_1 <- ggplot(summary_stats, aes(x = Method, y = Mean_AIC, fill = Method)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = method_colors) +
  labs(title = "Average AIC", y = "Mean AIC") +
  theme_minimal() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

p10_2 <- ggplot(summary_stats, aes(x = Method, y = Mean_BIC, fill = Method)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = method_colors) +
  labs(title = "Average BIC", y = "Mean BIC") +
  theme_minimal() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

p10_3 <- ggplot(summary_stats, aes(x = Method, y = Mean_RMSE, fill = Method)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = method_colors) +
  labs(title = "Average RMSE", y = "Mean RMSE") +
  theme_minimal() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

p10_4 <- ggplot(summary_stats, aes(x = Method, y = Mean_Time, fill = Method)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = method_colors) +
  labs(title = "Average Time", y = "Mean Time (sec)") +
  theme_minimal() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

p10 <- grid.arrange(p10_1, p10_2, p10_3, p10_4, ncol = 2,
                    top = "Summary Statistics by Method")

ggsave(file.path(plots_dir, "10_summary_statistics.png"), p10,
       width = 10, height = 8, dpi = 300)
cat("✓\n")

cat("\n================================================================\n")
cat("All visualizations created successfully!\n")
cat(sprintf("Location: %s\n", plots_dir))
cat("================================================================\n\n")
cat("Files created:\n")
created_files <- list.files(plots_dir, pattern = "\\.png$", full.names = FALSE)
if (length(created_files) > 0) {
  cat(paste0("  - ", created_files), sep = "\n")
}
cat("\n")
