# ==============================================================================
# Visualization of Monte Carlo Results for Seasonal Models
# ==============================================================================

library(ggplot2)
library(gridExtra)

# Load results
if (!file.exists("monte_carlo_seasonal_results.rds")) {
  stop("Results file not found. Please run monte_carlo_seasonal_comparison.R first.")
}

results <- readRDS("monte_carlo_seasonal_results.rds")

cat("==============================================================================\n")
cat("Visualizing Monte Carlo Results\n")
cat("==============================================================================\n")
cat("Data loaded from:", Sys.time(), "\n\n")

# ==============================================================================
# Helper Functions
# ==============================================================================

#' Create comparison data frame
create_comparison_df <- function(model_results, param_names) {
  sample_sizes <- names(model_results)

  df_list <- list()

  for (n_label in sample_sizes) {
    n <- as.numeric(gsub("n", "", n_label))
    res <- model_results[[n_label]]

    if (length(param_names) == 1) {
      # Single parameter (SMA)
      df_list[[n_label]] <- data.frame(
        n = n,
        parameter = param_names,
        pmm2_bias = res$pmm2[[param_names]]$bias,
        pmm2_rmse = res$pmm2[[param_names]]$rmse,
        pmm2_var = res$pmm2[[param_names]]$variance,
        css_bias = res$css[[param_names]]$bias,
        css_rmse = res$css[[param_names]]$rmse,
        css_var = res$css[[param_names]]$variance,
        var_reduction = 1 - res$pmm2[[param_names]]$variance / res$css[[param_names]]$variance,
        mean_g = res$pmm2$mean_g
      )
    } else {
      # Multiple parameters
      for (j in seq_along(param_names)) {
        df_list[[paste0(n_label, "_", j)]] <- data.frame(
          n = n,
          parameter = param_names[j],
          pmm2_bias = res$pmm2$coefficients[[j]]$bias,
          pmm2_rmse = res$pmm2$coefficients[[j]]$rmse,
          pmm2_var = res$pmm2$coefficients[[j]]$variance,
          css_bias = res$css$coefficients[[j]]$bias,
          css_rmse = res$css$coefficients[[j]]$rmse,
          css_var = res$css$coefficients[[j]]$variance,
          var_reduction = 1 - res$pmm2$coefficients[[j]]$variance / res$css$coefficients[[j]]$variance,
          mean_g = res$pmm2$mean_g
        )
      }
    }
  }

  do.call(rbind, df_list)
}

# ==============================================================================
# Prepare Data
# ==============================================================================

# SAR
df_sar <- create_comparison_df(results$sar, c("ar", "sar"))

# SMA
df_sma <- create_comparison_df(results$sma, c("sma"))

# SARMA
df_sarma <- rbind(
  create_comparison_df(results$sarma, c("AR", "SAR", "SMA"))
)

# SARIMA
df_sarima <- rbind(
  create_comparison_df(results$sarima, c("AR", "SAR", "SMA"))
)

# ==============================================================================
# Plot 1: Variance Comparison
# ==============================================================================

cat("Creating variance comparison plots...\n")

# SAR
p1_sar <- ggplot(df_sar, aes(x = factor(n), fill = parameter)) +
  geom_bar(aes(y = pmm2_var), stat = "identity", position = "dodge", alpha = 0.7) +
  geom_bar(aes(y = css_var), stat = "identity", position = "dodge",
           alpha = 0.3, color = "black", linewidth = 0.3) +
  labs(title = "SAR(1,1)_12: Variance Comparison",
       subtitle = "Solid = PMM2, Transparent = CSS",
       x = "Sample Size", y = "Variance", fill = "Parameter") +
  theme_minimal() +
  theme(legend.position = "bottom")

# SMA
p1_sma <- ggplot(df_sma, aes(x = factor(n))) +
  geom_bar(aes(y = pmm2_var, fill = "PMM2"), stat = "identity",
           position = "dodge", alpha = 0.7) +
  geom_bar(aes(y = css_var, fill = "CSS"), stat = "identity",
           position = "dodge", alpha = 0.7) +
  labs(title = "SMA(1)_12: Variance Comparison",
       x = "Sample Size", y = "Variance", fill = "Method") +
  scale_fill_manual(values = c("PMM2" = "#00BFC4", "CSS" = "#F8766D")) +
  theme_minimal() +
  theme(legend.position = "bottom")

# ==============================================================================
# Plot 2: Variance Reduction
# ==============================================================================

cat("Creating variance reduction plots...\n")

p2_sar <- ggplot(df_sar, aes(x = factor(n), y = var_reduction * 100,
                              fill = parameter)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "SAR(1,1)_12: Variance Reduction (PMM2 vs CSS)",
       x = "Sample Size", y = "Variance Reduction (%)", fill = "Parameter") +
  theme_minimal() +
  theme(legend.position = "bottom")

p2_sma <- ggplot(df_sma, aes(x = factor(n), y = var_reduction * 100)) +
  geom_bar(stat = "identity", fill = "#00BFC4") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "SMA(1)_12: Variance Reduction (PMM2 vs CSS)",
       x = "Sample Size", y = "Variance Reduction (%)") +
  theme_minimal()

p2_sarma <- ggplot(df_sarma, aes(x = factor(n), y = var_reduction * 100,
                                  fill = parameter)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "SARMA(1,0)×(1,1)_12: Variance Reduction",
       x = "Sample Size", y = "Variance Reduction (%)", fill = "Parameter") +
  theme_minimal() +
  theme(legend.position = "bottom")

p2_sarima <- ggplot(df_sarima, aes(x = factor(n), y = var_reduction * 100,
                                    fill = parameter)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "SARIMA(1,1,0)×(1,1,1)_12: Variance Reduction",
       x = "Sample Size", y = "Variance Reduction (%)", fill = "Parameter") +
  theme_minimal() +
  theme(legend.position = "bottom")

# ==============================================================================
# Plot 3: RMSE Comparison
# ==============================================================================

cat("Creating RMSE comparison plots...\n")

# Prepare data for RMSE comparison
df_rmse_sar <- df_sar
df_rmse_sar$method_pmm2 <- "PMM2"
df_rmse_sar$method_css <- "CSS"

df_rmse_long_sar <- rbind(
  data.frame(n = df_rmse_sar$n, parameter = df_rmse_sar$parameter,
             method = "PMM2", rmse = df_rmse_sar$pmm2_rmse),
  data.frame(n = df_rmse_sar$n, parameter = df_rmse_sar$parameter,
             method = "CSS", rmse = df_rmse_sar$css_rmse)
)

p3_sar <- ggplot(df_rmse_long_sar, aes(x = factor(n), y = rmse,
                                        fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ parameter, scales = "free_y") +
  labs(title = "SAR(1,1)_12: RMSE Comparison",
       x = "Sample Size", y = "RMSE", fill = "Method") +
  theme_minimal() +
  theme(legend.position = "bottom")

df_rmse_long_sma <- rbind(
  data.frame(n = df_sma$n, method = "PMM2", rmse = df_sma$pmm2_rmse),
  data.frame(n = df_sma$n, method = "CSS", rmse = df_sma$css_rmse)
)

p3_sma <- ggplot(df_rmse_long_sma, aes(x = factor(n), y = rmse, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "SMA(1)_12: RMSE Comparison",
       x = "Sample Size", y = "RMSE", fill = "Method") +
  theme_minimal() +
  theme(legend.position = "bottom")

# ==============================================================================
# Plot 4: Efficiency Factor g
# ==============================================================================

cat("Creating efficiency factor plots...\n")

# Combine all g values
df_g <- rbind(
  data.frame(model = "SAR", n = df_sar$n[!duplicated(df_sar$n)],
             g = sapply(split(df_sar$mean_g, df_sar$n), mean)),
  data.frame(model = "SMA", n = df_sma$n, g = df_sma$mean_g),
  data.frame(model = "SARMA", n = df_sarma$n[!duplicated(df_sarma$n)],
             g = sapply(split(df_sarma$mean_g, df_sarma$n), mean)),
  data.frame(model = "SARIMA", n = df_sarima$n[!duplicated(df_sarima$n)],
             g = sapply(split(df_sarima$mean_g, df_sarima$n), mean))
)

p4 <- ggplot(df_g, aes(x = factor(n), y = g, color = model, group = model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  labs(title = "Efficiency Factor g by Model and Sample Size",
       subtitle = "g = Var[PMM2] / Var[OLS] (lower is better)",
       x = "Sample Size", y = "g", color = "Model") +
  theme_minimal() +
  theme(legend.position = "bottom")

# ==============================================================================
# Save Plots
# ==============================================================================

cat("\nSaving plots...\n")

pdf("monte_carlo_seasonal_plots.pdf", width = 12, height = 8)

# Page 1: Variance comparisons
grid.arrange(p1_sar, p1_sma, ncol = 2,
             top = "Variance Comparison: PMM2 vs CSS")

# Page 2: Variance reduction
grid.arrange(p2_sar, p2_sma, ncol = 2,
             top = "Variance Reduction")

# Page 3: SARMA and SARIMA variance reduction
grid.arrange(p2_sarma, p2_sarima, ncol = 2,
             top = "Variance Reduction: Complex Models")

# Page 4: RMSE comparison
grid.arrange(p3_sar, p3_sma, ncol = 2,
             top = "RMSE Comparison")

# Page 5: Efficiency factor
print(p4)

dev.off()

cat("Plots saved to: monte_carlo_seasonal_plots.pdf\n")

# ==============================================================================
# Create Summary Table
# ==============================================================================

cat("\nCreating summary table...\n")

create_summary_table <- function(df, model_name) {
  df_summary <- aggregate(cbind(pmm2_var, css_var, var_reduction) ~ n + parameter,
                         data = df, FUN = mean)

  df_summary$model <- model_name
  df_summary$var_reduction_pct <- round(df_summary$var_reduction * 100, 2)

  df_summary[, c("model", "n", "parameter", "pmm2_var", "css_var", "var_reduction_pct")]
}

table_sar <- create_summary_table(df_sar, "SAR(1,1)_12")
table_sma <- create_summary_table(df_sma, "SMA(1)_12")
table_sarma <- create_summary_table(df_sarma, "SARMA(1,0)×(1,1)_12")
table_sarima <- create_summary_table(df_sarima, "SARIMA(1,1,0)×(1,1,1)_12")

summary_table <- rbind(table_sar, table_sma, table_sarma, table_sarima)

cat("\n==============================================================================\n")
cat("SUMMARY TABLE\n")
cat("==============================================================================\n\n")
print(summary_table, row.names = FALSE)

# Save summary table
write.csv(summary_table, "monte_carlo_summary_table.csv", row.names = FALSE)
cat("\nSummary table saved to: monte_carlo_summary_table.csv\n")

# ==============================================================================
# Final Statistics
# ==============================================================================

cat("\n==============================================================================\n")
cat("OVERALL STATISTICS\n")
cat("==============================================================================\n\n")

overall_stats <- aggregate(var_reduction_pct ~ model, data = summary_table,
                          FUN = function(x) c(mean = mean(x), min = min(x), max = max(x)))

cat("Average Variance Reduction by Model:\n")
for (i in 1:nrow(overall_stats)) {
  cat(sprintf("%-25s: Mean = %5.1f%%, Range = [%5.1f%%, %5.1f%%]\n",
              overall_stats$model[i],
              overall_stats$var_reduction_pct[i, "mean"],
              overall_stats$var_reduction_pct[i, "min"],
              overall_stats$var_reduction_pct[i, "max"]))
}

cat("\n==============================================================================\n")
cat("Visualization completed!\n")
cat("==============================================================================\n")
