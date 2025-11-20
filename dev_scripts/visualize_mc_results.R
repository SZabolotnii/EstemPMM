#!/usr/bin/env Rscript
# ==============================================================================
# Visualization of Monte Carlo Results
# Generate plots and tables for CRAN submission documentation
# ==============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# Load results
results <- readRDS("mc_results_final_comprehensive/comprehensive_mc_results.rds")
summary_df <- read.csv("mc_results_final_comprehensive/summary_report.csv")

cat("==============================================================================\n")
cat("Monte Carlo Results Visualization\n")
cat("==============================================================================\n\n")

# ==============================================================================
# 1. Variance Reduction by Model and Innovation Type
# ==============================================================================

cat("Creating variance reduction plots...\n")

# Filter out NA values and focus on main parameter
variance_plot_data <- summary_df %>%
  filter(!is.na(Var_Reduction), Var_Reduction > -100) %>%
  group_by(Scenario, Innovation) %>%
  summarise(
    Mean_Var_Reduction = mean(Var_Reduction, na.rm = TRUE),
    SD_Var_Reduction = sd(Var_Reduction, na.rm = TRUE),
    .groups = "drop"
  )

p1 <- ggplot(variance_plot_data, 
       aes(x = Scenario, y = Mean_Var_Reduction, fill = Innovation)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = Mean_Var_Reduction - SD_Var_Reduction,
                    ymax = Mean_Var_Reduction + SD_Var_Reduction),
                position = position_dodge(width = 0.9),
                width = 0.25) +
  labs(title = "PMM2 Variance Reduction by Model and Innovation Type",
       subtitle = "Average across all sample sizes and replication counts",
       x = "Model Type",
       y = "Variance Reduction (%)",
       fill = "Innovation Distribution") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 14),
        legend.position = "bottom") +
  scale_fill_brewer(palette = "Set2")

ggsave("mc_results_final_comprehensive/variance_reduction_plot.png",
       p1, width = 10, height = 6, dpi = 300)

# ==============================================================================
# 2. RMSE Improvement by Sample Size
# ==============================================================================

cat("Creating RMSE improvement plots...\n")

rmse_plot_data <- summary_df %>%
  filter(!is.na(RMSE_Improvement), RMSE_Improvement > -50) %>%
  group_by(Scenario, n) %>%
  summarise(
    Mean_RMSE_Imp = mean(RMSE_Improvement, na.rm = TRUE),
    SD_RMSE_Imp = sd(RMSE_Improvement, na.rm = TRUE),
    .groups = "drop"
  )

p2 <- ggplot(rmse_plot_data, 
       aes(x = factor(n), y = Mean_RMSE_Imp, color = Scenario, group = Scenario)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Mean_RMSE_Imp - SD_RMSE_Imp,
                    ymax = Mean_RMSE_Imp + SD_RMSE_Imp),
                width = 0.2) +
  labs(title = "RMSE Improvement vs Sample Size",
       subtitle = "PMM2 compared to CSS method",
       x = "Sample Size (n)",
       y = "RMSE Improvement (%)",
       color = "Model") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        legend.position = "right") +
  scale_color_brewer(palette = "Set1")

ggsave("mc_results_final_comprehensive/rmse_improvement_plot.png",
       p2, width = 10, height = 6, dpi = 300)

# ==============================================================================
# 3. Convergence Rates Comparison
# ==============================================================================

cat("Creating convergence rate plots...\n")

conv_data <- summary_df %>%
  filter(!is.na(PMM2_Conv)) %>%
  group_by(Scenario, Innovation) %>%
  summarise(
    PMM2_Conv_Rate = mean(PMM2_Conv, na.rm = TRUE) * 100,
    CSS_Conv_Rate = mean(CSS_Conv, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(PMM2_Conv_Rate, CSS_Conv_Rate),
               names_to = "Method",
               values_to = "Convergence_Rate") %>%
  mutate(Method = gsub("_Conv_Rate", "", Method))

p3 <- ggplot(conv_data, 
       aes(x = Scenario, y = Convergence_Rate, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  facet_wrap(~ Innovation, ncol = 2) +
  labs(title = "Convergence Rates: PMM2 vs CSS",
       x = "Model Type",
       y = "Convergence Rate (%)",
       fill = "Method") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 14),
        legend.position = "bottom",
        strip.text = element_text(face = "bold")) +
  scale_fill_manual(values = c("PMM2" = "#2E86AB", "CSS" = "#A23B72"))

ggsave("mc_results_final_comprehensive/convergence_rates_plot.png",
       p3, width = 12, height = 8, dpi = 300)

# ==============================================================================
# 4. g-Factor Distribution
# ==============================================================================

cat("Creating g-factor distribution plots...\n")

g_data <- summary_df %>%
  filter(!is.na(Mean_g), Mean_g > 0)

p4 <- ggplot(g_data, aes(x = Mean_g, fill = Scenario)) +
  geom_density(alpha = 0.6) +
  facet_wrap(~ Innovation, ncol = 2) +
  labs(title = "Distribution of PMM2 g-Factor",
       subtitle = "Lower g-factor indicates higher efficiency gain",
       x = "g-Factor",
       y = "Density",
       fill = "Model") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        legend.position = "bottom",
        strip.text = element_text(face = "bold")) +
  scale_fill_brewer(palette = "Dark2") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8)

ggsave("mc_results_final_comprehensive/g_factor_distribution_plot.png",
       p4, width = 12, height = 8, dpi = 300)

# ==============================================================================
# 5. Sample Size Effect
# ==============================================================================

cat("Creating sample size effect plot...\n")

sample_size_data <- summary_df %>%
  filter(!is.na(Var_Reduction), Scenario %in% c("SAR_1_1_12", "SMA_1_12", "SARMA_1_0_1_1_12")) %>%
  group_by(Scenario, n, Innovation) %>%
  summarise(Mean_Var_Red = mean(Var_Reduction, na.rm = TRUE),
            .groups = "drop")

p5 <- ggplot(sample_size_data,
       aes(x = n, y = Mean_Var_Red, color = Innovation, shape = Scenario)) +
  geom_point(size = 4, alpha = 0.7) +
  geom_line(aes(group = interaction(Scenario, Innovation)), linewidth = 1) +
  labs(title = "Variance Reduction vs Sample Size",
       subtitle = "Effect of sample size on PMM2 efficiency",
       x = "Sample Size (n)",
       y = "Variance Reduction (%)",
       color = "Innovation",
       shape = "Model") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        legend.position = "right") +
  scale_color_brewer(palette = "Set2")

ggsave("mc_results_final_comprehensive/sample_size_effect_plot.png",
       p5, width = 10, height = 6, dpi = 300)

# ==============================================================================
# 6. Summary Table
# ==============================================================================

cat("\nCreating summary tables...\n")

summary_table <- summary_df %>%
  filter(!is.na(Var_Reduction), Scenario != "SARIMA_1_1_0_1_1_1_12") %>%
  group_by(Scenario, Innovation) %>%
  summarise(
    Mean_Var_Reduction = round(mean(Var_Reduction, na.rm = TRUE), 2),
    Mean_RMSE_Imp = round(mean(RMSE_Improvement, na.rm = TRUE), 2),
    Mean_PMM2_Conv = round(mean(PMM2_Conv, na.rm = TRUE) * 100, 1),
    Mean_CSS_Conv = round(mean(CSS_Conv, na.rm = TRUE) * 100, 1),
    Mean_g = round(mean(Mean_g, na.rm = TRUE), 4),
    .groups = "drop"
  ) %>%
  arrange(Scenario, Innovation)

write.csv(summary_table, 
          "mc_results_final_comprehensive/detailed_summary_table.csv",
          row.names = FALSE)

print(summary_table)

# ==============================================================================
# 7. Best Performance Summary
# ==============================================================================

cat("\n==============================================================================\n")
cat("BEST PERFORMANCE SUMMARY\n")
cat("==============================================================================\n\n")

best_by_model <- summary_df %>%
  filter(!is.na(Var_Reduction), Scenario != "SARIMA_1_1_0_1_1_1_12") %>%
  group_by(Scenario) %>%
  slice_max(Var_Reduction, n = 1) %>%
  select(Scenario, Innovation, n, Reps, Var_Reduction, RMSE_Improvement, Mean_g)

cat("Top configurations by variance reduction:\n\n")
print(best_by_model, n = Inf)

# ==============================================================================
# 8. Innovation Type Comparison
# ==============================================================================

cat("\nCreating innovation type comparison...\n")

innov_comparison <- summary_df %>%
  filter(!is.na(Var_Reduction), Scenario != "SARIMA_1_1_0_1_1_1_12") %>%
  group_by(Innovation) %>%
  summarise(
    Mean_Var_Red = round(mean(Var_Reduction, na.rm = TRUE), 2),
    SD_Var_Red = round(sd(Var_Reduction, na.rm = TRUE), 2),
    Mean_RMSE_Imp = round(mean(RMSE_Improvement, na.rm = TRUE), 2),
    Mean_g = round(mean(Mean_g, na.rm = TRUE), 4),
    .groups = "drop"
  ) %>%
  arrange(desc(Mean_Var_Red))

cat("\n")
cat("==============================================================================\n")
cat("INNOVATION TYPE RANKING\n")
cat("==============================================================================\n\n")
print(innov_comparison)

write.csv(innov_comparison,
          "mc_results_final_comprehensive/innovation_comparison.csv",
          row.names = FALSE)

# ==============================================================================
# Create Combined Summary Plot
# ==============================================================================

cat("\nCreating combined summary plot...\n")

# Combine multiple plots
combined_plot <- grid.arrange(
  p1 + theme(legend.position = "none"),
  p2 + theme(legend.position = "none"),
  ncol = 2,
  top = "Monte Carlo Results Summary: PMM2 vs CSS"
)

ggsave("mc_results_final_comprehensive/combined_summary_plot.png",
       combined_plot, width = 16, height = 8, dpi = 300)

# ==============================================================================
# Final Summary
# ==============================================================================

cat("\n")
cat("==============================================================================\n")
cat("VISUALIZATION COMPLETE\n")
cat("==============================================================================\n\n")

cat("Generated files:\n")
cat("  1. variance_reduction_plot.png\n")
cat("  2. rmse_improvement_plot.png\n")
cat("  3. convergence_rates_plot.png\n")
cat("  4. g_factor_distribution_plot.png\n")
cat("  5. sample_size_effect_plot.png\n")
cat("  6. combined_summary_plot.png\n")
cat("  7. detailed_summary_table.csv\n")
cat("  8. innovation_comparison.csv\n\n")

cat("All visualizations saved to: mc_results_final_comprehensive/\n")

cat("\n==============================================================================\n")
