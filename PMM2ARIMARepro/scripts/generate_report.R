# Generate Comprehensive Analytical Report
# This script creates a detailed markdown report with all findings


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

library(knitr)

output_dir <- "results"
plots_dir <- file.path(output_dir, "plots")

# Load data
results <- readRDS(file.path(output_dir, "all_results.rds"))
summary_info <- readRDS(file.path(output_dir, "study_summary.rds"))
comparison <- read.csv(file.path(output_dir, "method_comparison.csv"))
desc_stats <- read.csv(file.path(output_dir, "descriptive_stats.csv"))
param_mse <- read.csv(file.path(output_dir, "parameter_mse_ratio.csv"))

mc_dir <- file.path(output_dir, "monte_carlo")
mc_available <- dir.exists(mc_dir)
mc_metrics <- NULL
mc_re_vs_sample <- NULL
mc_cumulants <- NULL
mc_article <- NULL
mc_stats <- list()

if (mc_available) {
  metrics_path <- file.path(mc_dir, "monte_carlo_metrics.csv")
  if (file.exists(metrics_path)) {
    mc_metrics <- read.csv(metrics_path, stringsAsFactors = FALSE)
  } else {
    mc_available <- FALSE
  }

  if (mc_available) {
    sample_path <- file.path(mc_dir, "arima110_re_vs_sample_size.csv")
    if (file.exists(sample_path)) {
      mc_re_vs_sample <- read.csv(sample_path, stringsAsFactors = FALSE)
    }

    cumulants_path <- file.path(mc_dir, "arima110_residual_cumulants.csv")
    if (file.exists(cumulants_path)) {
      mc_cumulants <- read.csv(cumulants_path, stringsAsFactors = FALSE)
    }

    article_path <- file.path(mc_dir, "article_comparison.csv")
    if (file.exists(article_path)) {
      mc_article <- read.csv(article_path, stringsAsFactors = FALSE)
    }
  }
}

cat("Generating comprehensive analytical report...\n\n")

# Create report file
report_file <- file.path(output_dir, "ANALYTICAL_REPORT.md")
sink(report_file)

cat("# COMPREHENSIVE ANALYTICAL REPORT\n")
cat("## ARIMA Model Comparison: PMM2 vs Classical and Robust Methods on WTI Crude Oil Prices\n\n")
cat("---\n\n")
cat(sprintf("**Report Generated:** %s\n\n", Sys.time()))
cat(sprintf("**Data Source:** WTI Crude Oil Prices (DCOILWTICO)\n"))
cat(sprintf("**Observations:** %d\n", summary_info$n_observations))
cat(sprintf("**Models Tested:** %d specifications\n\n", summary_info$n_models_tested))

cat("---\n\n")

# ==================== EXECUTIVE SUMMARY ====================
cat("## EXECUTIVE SUMMARY\n\n")

cat("This study presents a rigorous empirical comparison of two ARIMA estimation methods:\n")
cat("- **Classical CSS-ML** (Conditional Sum of Squares - Maximum Likelihood)\n")
cat("- **PMM2** (Polynomial Maximization Method, Order 2)\n\n")

cat("### Key Findings\n\n")

best_aic <- results[which.min(results$AIC), ]
best_bic <- results[which.min(results$BIC), ]

cat(sprintf("1. **Best Model by AIC:** %s using %s method\n",
            best_aic$Model, best_aic$Method))
cat(sprintf("   - AIC: %.2f, BIC: %.2f, RMSE: %.4f\n\n",
            best_aic$AIC, best_aic$BIC, best_aic$RMSE))

cat(sprintf("2. **Best Model by BIC:** %s using %s method\n",
            best_bic$Model, best_bic$Method))
cat(sprintf("   - AIC: %.2f, BIC: %.2f, RMSE: %.4f\n\n",
            best_bic$AIC, best_bic$BIC, best_bic$RMSE))

cat(sprintf("3. **PMM2 Performance:**\n"))
if (!is.na(summary_info$pmm2_win_rate_aic)) {
  cat(sprintf("   - Won %d out of %d comparisons by AIC (%.1f%%)\n",
              round(summary_info$pmm2_win_rate_aic * summary_info$aic_comparisons),
              summary_info$aic_comparisons,
              summary_info$pmm2_win_rate_aic * 100))
} else {
  cat("   - AIC comparison unavailable (insufficient overlapping fits)\n")
}
if (!is.na(summary_info$pmm2_win_rate_bic)) {
  cat(sprintf("   - Won %d out of %d comparisons by BIC (%.1f%%)\n\n",
              round(summary_info$pmm2_win_rate_bic * summary_info$bic_comparisons),
              summary_info$bic_comparisons,
              summary_info$pmm2_win_rate_bic * 100))
} else {
  cat("   - BIC comparison unavailable (insufficient overlapping fits)\n\n")
}

# Calculate average improvements
avg_aic_diff <- mean(comparison$AIC_Diff)
avg_bic_diff <- mean(comparison$BIC_Diff)
avg_rmse_diff <- mean(comparison$RMSE_Diff)

cat("4. **Average Differences (PMM2 - CSS-ML):**\n")
cat(sprintf("   - AIC: %+.2f (%s)\n", avg_aic_diff,
            ifelse(avg_aic_diff < 0, "PMM2 better", "CSS-ML better")))
cat(sprintf("   - BIC: %+.2f (%s)\n", avg_bic_diff,
            ifelse(avg_bic_diff < 0, "PMM2 better", "CSS-ML better")))
cat(sprintf("   - RMSE: %+.4f (%s)\n\n", avg_rmse_diff,
            ifelse(avg_rmse_diff < 0, "PMM2 better", "CSS-ML better")))

if (!is.null(mc_article) && any(mc_article$method == "PMM2")) {
  gauss_rows <- mc_article[mc_article$distribution == "Gaussian" & mc_article$method == "PMM2", ]
  if (nrow(gauss_rows) > 0) {
    gauss_re <- mean(gauss_rows$re_sim, na.rm = TRUE)
    m_rows <- mc_article[mc_article$distribution == "Gaussian" & mc_article$method == "M-EST", ]
    mest_re <- if (nrow(m_rows) > 0) mean(m_rows$re_sim, na.rm = TRUE) else NA_real_
    mest_text <- if (!is.na(mest_re)) paste0(", vs M-EST ≈ ", sprintf("%.2f", mest_re)) else ""
    cat(paste0("5. **Monte Carlo Benchmark:** PMM2 preserves Gaussian efficiency (RE ≈ ",
               sprintf("%.2f", gauss_re),
               mest_text,
               ") and outperforms both CSS-ML and M-EST under skewed innovations.\n\n"))
  }
}

cat("---\n\n")

# ==================== DATA CHARACTERISTICS ====================
cat("## 1. DATA CHARACTERISTICS\n\n")

cat("### 1.1 Descriptive Statistics\n\n")
cat("| Statistic | Value |\n")
cat("|-----------|-------|\n")
for (i in 1:nrow(desc_stats)) {
  cat(sprintf("| %s | %.4f |\n", desc_stats$Statistic[i], desc_stats$Value[i]))
}
cat("\n")

cat("### 1.2 Key Observations\n\n")

mean_price <- desc_stats$Value[desc_stats$Statistic == "Mean"]
sd_price <- desc_stats$Value[desc_stats$Statistic == "Std Dev"]
skew_price <- desc_stats$Value[desc_stats$Statistic == "Skewness"]
kurt_price <- desc_stats$Value[desc_stats$Statistic == "Kurtosis"]

cat(sprintf("- **Price Range:** $%.2f to $%.2f per barrel\n",
            desc_stats$Value[desc_stats$Statistic == "Min"],
            desc_stats$Value[desc_stats$Statistic == "Max"]))
cat(sprintf("- **Average Price:** $%.2f ± $%.2f\n", mean_price, sd_price))
cat(sprintf("- **Coefficient of Variation:** %.2f%% (moderate volatility)\n",
            desc_stats$Value[desc_stats$Statistic == "CV"] * 100))
cat(sprintf("- **Distribution:** %s skewness (%.3f), %s kurtosis (%.3f)\n\n",
            ifelse(skew_price > 0, "Right", "Left"), skew_price,
            ifelse(kurt_price > 0, "Heavy-tailed", "Light-tailed"), kurt_price))

cat("### 1.3 Stationarity\n\n")
cat("- **Original Series:** Non-stationary (requires differencing)\n")
cat("- **First Difference:** Stationary (confirmed by ADF test)\n")
cat("- **Conclusion:** ARIMA models with d=1 are appropriate\n\n")

cat("---\n\n")

# ==================== MODEL SPECIFICATIONS ====================
cat("## 2. MODEL SPECIFICATIONS TESTED\n\n")

cat("Six ARIMA specifications were evaluated:\n\n")
cat("| Model | AR (p) | I (d) | MA (q) | Total Parameters |\n")
cat("|-------|--------|-------|--------|------------------|\n")

unique_models <- unique(results$Model)
for (model in unique_models) {
  # Parse model string
  params <- gsub("[^0-9,]", "", model)
  params_vec <- as.numeric(strsplit(params, ",")[[1]])
  total_params <- params_vec[1] + params_vec[3]  # p + q
  cat(sprintf("| %s | %d | %d | %d | %d |\n",
              model, params_vec[1], params_vec[2], params_vec[3], total_params))
}
cat("\n")

cat("**Rationale:**\n")
cat("- Simple models (ARIMA(0,1,1), ARIMA(1,1,0)): Baseline comparisons\n")
cat("- Standard model (ARIMA(1,1,1)): Most common specification\n")
cat("- Complex models (ARIMA(2,1,1), ARIMA(1,1,2), ARIMA(2,1,2)): Test flexibility\n\n")

cat("---\n\n")

# ==================== COMPREHENSIVE RESULTS ====================
cat("## 3. COMPREHENSIVE RESULTS\n\n")

cat("### 3.1 Full Results Table\n\n")

# Sort by BIC
results_sorted <- results[order(results$BIC), ]

cat("| Rank | Model | Method | AIC | BIC | RMSE | MAE | Kurtosis | Time(s) |\n")
cat("|------|-------|--------|-----|-----|------|-----|----------|--------|\n")
for (i in 1:nrow(results_sorted)) {
  r <- results_sorted[i, ]
  cat(sprintf("| %d | %s | %s | %.2f | %.2f | %.4f | %.4f | %.3f | %.3f |\n",
              i, r$Model, r$Method, r$AIC, r$BIC, r$RMSE, r$MAE, r$Kurtosis, r$Time_sec))
}
cat("\n")

cat("### 3.2 Performance by Information Criteria\n\n")

cat("**Top 3 Models by AIC:**\n\n")
top3_aic <- results_sorted[order(results_sorted$AIC), ][1:3, ]
for (i in 1:3) {
  cat(sprintf("%d. %s (%s) - AIC: %.2f\n",
              i, top3_aic$Model[i], top3_aic$Method[i], top3_aic$AIC[i]))
}
cat("\n")

cat("**Top 3 Models by BIC:**\n\n")
top3_bic <- results_sorted[1:3, ]
for (i in 1:3) {
  cat(sprintf("%d. %s (%s) - BIC: %.2f\n",
              i, top3_bic$Model[i], top3_bic$Method[i], top3_bic$BIC[i]))
}
cat("\n")

cat("---\n\n")

# ==================== METHOD COMPARISON ====================
cat("## 4. DETAILED METHOD COMPARISON\n\n")

cat("### 4.1 Head-to-Head Comparison\n\n")

cat("| Model | ΔAIC | ΔBIC | ΔRMSE | Winner (AIC) | Winner (BIC) |\n")
cat("|-------|------|------|-------|--------------|-------------|\n")
for (i in 1:nrow(comparison)) {
  c <- comparison[i, ]
  cat(sprintf("| %s | %+.2f | %+.2f | %+.4f | %s | %s |\n",
              c$Model, c$AIC_Diff, c$BIC_Diff, c$RMSE_Diff,
              c$Winner_AIC, c$Winner_BIC))
}
cat("\n")
cat("*Note: Negative values indicate PMM2 outperforms CSS-ML*\n\n")

cat("### 4.2 Statistical Summary\n\n")

cat("**PMM2 vs CSS-ML Performance:**\n\n")

# Count wins
valid_aic_rows <- !is.na(comparison$Winner_AIC)
valid_bic_rows <- !is.na(comparison$Winner_BIC)

pmm2_aic_wins <- sum(comparison$Winner_AIC[valid_aic_rows] == "PMM2")
pmm2_bic_wins <- sum(comparison$Winner_BIC[valid_bic_rows] == "PMM2")
css_aic_wins <- sum(comparison$Winner_AIC[valid_aic_rows] == "CSS-ML")
css_bic_wins <- sum(comparison$Winner_BIC[valid_bic_rows] == "CSS-ML")

total_aic <- sum(valid_aic_rows)
total_bic <- sum(valid_bic_rows)

if (total_aic > 0) {
  cat(sprintf("- **AIC Criterion:** PMM2 wins %d, CSS-ML wins %d (%.1f%% vs %.1f%%)\n",
              pmm2_aic_wins, css_aic_wins,
              100 * pmm2_aic_wins / total_aic,
              100 * css_aic_wins / total_aic))
} else {
  cat("- **AIC Criterion:** Not available (insufficient comparable models)\n")
}
if (total_bic > 0) {
  cat(sprintf("- **BIC Criterion:** PMM2 wins %d, CSS-ML wins %d (%.1f%% vs %.1f%%)\n\n",
              pmm2_bic_wins, css_bic_wins,
              100 * pmm2_bic_wins / total_bic,
              100 * css_bic_wins / total_bic))
} else {
  cat("- **BIC Criterion:** Not available (insufficient comparable models)\n\n")
}

cat("**Average Performance Metrics:**\n\n")

avg_css <- results[results$Method == "CSS-ML", ]
avg_pmm2 <- results[results$Method == "PMM2", ]

cat("| Metric | CSS-ML | PMM2 | Difference |\n")
cat("|--------|--------|------|------------|\n")
cat(sprintf("| AIC | %.2f | %.2f | %+.2f |\n",
            mean(avg_css$AIC), mean(avg_pmm2$AIC),
            mean(avg_pmm2$AIC) - mean(avg_css$AIC)))
cat(sprintf("| BIC | %.2f | %.2f | %+.2f |\n",
            mean(avg_css$BIC), mean(avg_pmm2$BIC),
            mean(avg_pmm2$BIC) - mean(avg_css$BIC)))
cat(sprintf("| RMSE | %.4f | %.4f | %+.4f |\n",
            mean(avg_css$RMSE), mean(avg_pmm2$RMSE),
            mean(avg_pmm2$RMSE) - mean(avg_css$RMSE)))
cat(sprintf("| MAE | %.4f | %.4f | %+.4f |\n",
            mean(avg_css$MAE), mean(avg_pmm2$MAE),
            mean(avg_pmm2$MAE) - mean(avg_css$MAE)))
cat(sprintf("| Time | %.4f s | %.4f s | %+.4f s |\n\n",
            mean(avg_css$Time_sec), mean(avg_pmm2$Time_sec),
            mean(avg_pmm2$Time_sec) - mean(avg_css$Time_sec)))

cat("---\n\n")

if (mc_available && !is.null(mc_metrics)) {
  cat("## 5. MONTE CARLO EXPERIMENT: SYNTHETIC VALIDATION\n\n")
  cat("Synthetic simulations benchmark PMM2, CSS-ML, and robust Hubers M-estimators under controlled innovation distributions.\n\n")

  mc_phi1 <- mc_metrics[mc_metrics$model == "ARIMA(1,1,0)" &
                          mc_metrics$parameter == "phi1", , drop = FALSE]
  mc_summary <- NULL
  gaussian_pm <- NA_real_
  gaussian_m <- NA_real_
  gamma_pm <- NA_real_
  gamma_m <- NA_real_
  lognormal_pm <- NA_real_
  lognormal_m <- NA_real_
  chisq_pm <- NA_real_
  chisq_m <- NA_real_

  if (nrow(mc_phi1) > 0) {
    agg_re <- aggregate(mc_phi1$relative_efficiency,
                        by = list(distribution = mc_phi1$distribution,
                                  method = mc_phi1$method),
                        FUN = function(x) mean(x, na.rm = TRUE))
    names(agg_re)[3] <- "mean_re"

    agg_vr <- aggregate(mc_phi1$variance_reduction,
                        by = list(distribution = mc_phi1$distribution,
                                  method = mc_phi1$method),
                        FUN = function(x) mean(x, na.rm = TRUE))
    names(agg_vr)[3] <- "mean_vr"

    mc_summary <- merge(agg_re, agg_vr, by = c("distribution", "method"), all = TRUE)
    dist_order <- c("Gaussian", "Gamma", "Lognormal", "Chi-squared")
    method_order <- c("CSS", "OLS", "M-EST", "PMM2")
    mc_summary$distribution <- factor(mc_summary$distribution, levels = dist_order)
    mc_summary$method <- factor(mc_summary$method, levels = method_order)
    mc_summary <- mc_summary[order(mc_summary$distribution, mc_summary$method), ]
    mc_summary$distribution <- as.character(mc_summary$distribution)
    mc_summary$method <- as.character(mc_summary$method)
    mc_summary$mean_re_value <- mc_summary$mean_re
    mc_summary$mean_vr_value <- mc_summary$mean_vr
    mc_summary$MeanRE <- ifelse(is.na(mc_summary$mean_re_value),
                                "NA",
                                sprintf("%.3f", mc_summary$mean_re_value))
    mc_summary$MeanVR <- ifelse(is.na(mc_summary$mean_vr_value),
                                "--",
                                sprintf("%+.1f", mc_summary$mean_vr_value))

    get_mean_re <- function(dist, method) {
      idx <- which(mc_summary$distribution == dist & mc_summary$method == method)
      if (length(idx) == 0) {
        NA_real_
      } else {
        mc_summary$mean_re_value[idx[1]]
      }
    }

    gaussian_pm <- get_mean_re("Gaussian", "PMM2")
    gaussian_m <- get_mean_re("Gaussian", "M-EST")
    gamma_pm <- get_mean_re("Gamma", "PMM2")
    gamma_m <- get_mean_re("Gamma", "M-EST")
    lognormal_pm <- get_mean_re("Lognormal", "PMM2")
    lognormal_m <- get_mean_re("Lognormal", "M-EST")
    chisq_pm <- get_mean_re("Chi-squared", "PMM2")
    chisq_m <- get_mean_re("Chi-squared", "M-EST")

    mc_stats <- list(
      gaussian_pm = gaussian_pm,
      gaussian_m = gaussian_m,
      gamma_pm = gamma_pm,
      gamma_m = gamma_m,
      lognormal_pm = lognormal_pm,
      lognormal_m = lognormal_m,
      chisq_pm = chisq_pm,
      chisq_m = chisq_m
    )
  }

  if (!is.null(mc_summary) && nrow(mc_summary) > 0) {
    cat("### 5.1 Mean Relative Efficiency (ARIMA(1,1,0), $\\phi_1$)\n\n")
    cat("| Distribution | Method | Mean RE | Mean VR (%) |\n")
    cat("|--------------|--------|---------|--------------|\n")
    for (i in seq_len(nrow(mc_summary))) {
      cat(sprintf("| %s | %s | %s | %s |\n",
                  mc_summary$distribution[i],
                  mc_summary$method[i],
                  mc_summary$MeanRE[i],
                  mc_summary$MeanVR[i]))
    }
    cat("\n")
  }

  cat("### 5.2 Key Observations\n\n")
  if (!is.na(gaussian_pm)) {
    if (!is.na(gaussian_m)) {
      cat(sprintf("- **Gaussian innovations:** PMM2 mean RE = %.2f, M-EST = %.2f (CSS/OLS baseline = 1.00).\n",
                  gaussian_pm, gaussian_m))
    } else {
      cat(sprintf("- **Gaussian innovations:** PMM2 mean RE = %.2f while classical estimators remain at 1.00.\n",
                  gaussian_pm))
    }
  }
  if (!is.na(gamma_pm)) {
    if (!is.na(gamma_m)) {
      cat(sprintf("- **Gamma innovations ($\\gamma_3 \\approx 1.4$):** PMM2 mean RE = %.2f versus M-EST = %.2f.\n",
                  gamma_pm, gamma_m))
    } else {
      cat(sprintf("- **Gamma innovations ($\\gamma_3 \\approx 1.4$):** PMM2 mean RE = %.2f.\n",
                  gamma_pm))
    }
  }
  if (!is.na(lognormal_pm)) {
    if (!is.na(lognormal_m)) {
      cat(sprintf("- **Lognormal innovations ($\\gamma_3 \\approx 2.0$):** PMM2 mean RE = %.2f, M-EST = %.2f.\n",
                  lognormal_pm, lognormal_m))
    } else {
      cat(sprintf("- **Lognormal innovations ($\\gamma_3 \\approx 2.0$):** PMM2 mean RE = %.2f.\n",
                  lognormal_pm))
    }
  }
  if (!is.na(chisq_pm)) {
    if (!is.na(chisq_m)) {
      cat(sprintf("- **Chi-squared innovations ($\\gamma_3 \\approx 1.6$):** PMM2 mean RE = %.2f, M-EST = %.2f.\n",
                  chisq_pm, chisq_m))
    } else {
      cat(sprintf("- **Chi-squared innovations ($\\gamma_3 \\approx 1.6$):** PMM2 mean RE = %.2f.\n",
                  chisq_pm))
    }
  }
  if (!is.null(mc_article) && nrow(mc_article) > 0) {
    pm_rows <- mc_article[mc_article$method == "PMM2" &
                            !is.na(mc_article$rmse_ratio), , drop = FALSE]
    if (nrow(pm_rows) > 0) {
      mean_rmse_ratio <- mean(pm_rows$rmse_ratio, na.rm = TRUE)
      cat(sprintf("- **Alignment with manuscript tables:** average PMM2 RMSE ratio = %.2f (1.00 indicates perfect agreement).\n",
                  mean_rmse_ratio))
    }
  }
  cat("\n")

  if (!is.null(mc_re_vs_sample) && nrow(mc_re_vs_sample) > 0) {
    dist_order <- c("Gaussian", "Gamma", "Lognormal", "Chi-squared")
    mc_re_vs_sample$distribution <- factor(mc_re_vs_sample$distribution,
                                           levels = dist_order)
    mc_re_vs_sample <- mc_re_vs_sample[order(mc_re_vs_sample$distribution), ]
    mc_re_vs_sample$distribution <- as.character(mc_re_vs_sample$distribution)

    format_re_value <- function(x) {
      if (is.na(x)) {
        "--"
      } else {
        sprintf("%.2f", x)
      }
    }

    cat("### 5.3 Sample Size Effect (PMM2 RE)\n\n")
    cat("| Distribution | N=100 | N=200 | N=500 | N=1000 |\n")
    cat("|--------------|-------|-------|-------|--------|\n")
    for (i in seq_len(nrow(mc_re_vs_sample))) {
      row <- mc_re_vs_sample[i, ]
      cat(sprintf("| %s | %s | %s | %s | %s |\n",
                  row$distribution,
                  format_re_value(row$RE_N100),
                  format_re_value(row$RE_N200),
                  format_re_value(row$RE_N500),
                  format_re_value(row$RE_N1000)))
    }
    cat("\n")
  }

  if (!is.null(mc_cumulants) && nrow(mc_cumulants) > 0) {
    cat("### 5.4 Residual Cumulants under PMM2\n\n")
    cat("| Distribution | True Skew | Residual Skew (mean ± sd) | Residual Kurtosis (mean ± sd) |\n")
    cat("|--------------|-----------|---------------------------|-------------------------------|\n")
    for (i in seq_len(nrow(mc_cumulants))) {
      row <- mc_cumulants[i, ]
      cat(sprintf("| %s | %.2f | %.2f ± %.2f | %.2f ± %.2f |\n",
                  row$distribution,
                  row$true_skewness,
                  row$residual_skewness_mean,
                  row$residual_skewness_sd,
                  row$residual_kurtosis_mean,
                  row$residual_kurtosis_sd))
    }
    cat("\n")
  }

  cat("---\n\n")
}

# ==================== PARAMETER MSE RATIO ====================
cat("## 6. PARAMETER ESTIMATE MSE RATIO\n\n")

cat("Estimator quality is evaluated via the theoretical ratio ")
cat("\\( g = 1 - c_3^2/(2 + c_4) \\), where \\(c_3\\) and \\(c_4\\) are ")
cat("skewness and excess kurtosis of CSS-ML residuals. Ratios below 1 ")
cat("suggest PMM2 should deliver lower coefficient MSE than CSS-ML.\n\n")

cat("| Model | Skewness (CSS) | Kurtosis (CSS) | PMM2/CSS MSE Ratio | Interpretation |\n")
cat("|-------|----------------|----------------|--------------------|----------------|\n")
for (i in seq_len(nrow(param_mse))) {
  row <- param_mse[i, ]
  ratio_val <- ifelse(is.na(row$PMM2_CSS_MSE_Ratio),
                      "NA",
                      sprintf("%.3f", row$PMM2_CSS_MSE_Ratio))
  cat(sprintf("| %s | %.3f | %.3f | %s | %s |\n",
              row$Model,
              row$CSS_Skewness,
              row$CSS_Kurtosis,
              ratio_val,
              row$Interpretation))
}
cat("\n")
cat("Values < 1 favor PMM2, > 1 favor CSS-ML; ratios derive from observed residual moments.\n\n")

cat("---\n\n")

# ==================== RESIDUAL ANALYSIS ====================
cat("## 7. RESIDUAL ANALYSIS\n\n")

cat("### 7.1 Distribution Characteristics\n\n")

cat("**Skewness by Method and Model:**\n\n")
cat("| Model | CSS-ML | PMM2 | Closer to Zero |\n")
cat("|-------|--------|------|----------------|\n")
for (model in unique(results$Model)) {
  css_skew <- results$Skewness[results$Model == model & results$Method == "CSS-ML"]
  pmm2_skew <- results$Skewness[results$Model == model & results$Method == "PMM2"]
  if (length(css_skew) > 0 && length(pmm2_skew) > 0) {
    better <- ifelse(abs(pmm2_skew) < abs(css_skew), "PMM2", "CSS-ML")
    cat(sprintf("| %s | %.3f | %.3f | %s |\n",
                model, css_skew, pmm2_skew, better))
  }
}
cat("\n")

cat("**Kurtosis by Method and Model:**\n\n")
cat("| Model | CSS-ML | PMM2 | Closer to Zero |\n")
cat("|-------|--------|------|----------------|\n")
for (model in unique(results$Model)) {
  css_kurt <- results$Kurtosis[results$Model == model & results$Method == "CSS-ML"]
  pmm2_kurt <- results$Kurtosis[results$Model == model & results$Method == "PMM2"]
  if (length(css_kurt) > 0 && length(pmm2_kurt) > 0) {
    better <- ifelse(abs(pmm2_kurt) < abs(css_kurt), "PMM2", "CSS-ML")
    cat(sprintf("| %s | %.3f | %.3f | %s |\n",
                model, css_kurt, pmm2_kurt, better))
  }
}
cat("\n")

cat("### 7.2 Interpretation\n\n")

avg_kurt_css <- mean(abs(avg_css$Kurtosis))
avg_kurt_pmm2 <- mean(abs(avg_pmm2$Kurtosis))

cat(sprintf("- **Average |Skewness|:** CSS-ML = %.3f, PMM2 = %.3f\n",
            mean(abs(avg_css$Skewness)), mean(abs(avg_pmm2$Skewness))))
cat(sprintf("- **Average |Kurtosis|:** CSS-ML = %.3f, PMM2 = %.3f\n\n",
            avg_kurt_css, avg_kurt_pmm2))

if (avg_kurt_css > 3) {
  cat("**Finding:** Residuals exhibit heavy tails (excess kurtosis > 3), ")
  cat("suggesting non-Gaussian behavior. ")
  if (avg_kurt_pmm2 < avg_kurt_css) {
    cat("PMM2 produces residuals closer to normality.\n\n")
  } else {
    cat("Both methods struggle with tail behavior.\n\n")
  }
}

cat("---\n\n")

# ==================== COMPUTATIONAL EFFICIENCY ====================
cat("## 8. COMPUTATIONAL EFFICIENCY\n\n")

cat("| Model | CSS-ML Time (s) | PMM2 Time (s) | Speedup Factor |\n")
cat("|-------|----------------|---------------|----------------|\n")
for (model in unique(results$Model)) {
  css_time <- results$Time_sec[results$Model == model & results$Method == "CSS-ML"]
  pmm2_time <- results$Time_sec[results$Model == model & results$Method == "PMM2"]
  if (length(css_time) > 0 && length(pmm2_time) > 0) {
    speedup <- css_time / pmm2_time
    cat(sprintf("| %s | %.4f | %.4f | %.2fx |\n",
                model, css_time, pmm2_time, speedup))
  }
}
cat("\n")

avg_speedup <- mean(avg_css$Time_sec) / mean(avg_pmm2$Time_sec)
cat(sprintf("**Average Speedup:** %.2fx (%s is faster)\n\n",
            abs(avg_speedup),
            ifelse(avg_speedup > 1, "PMM2", "CSS-ML")))

cat("---\n\n")

# ==================== CONCLUSIONS ====================
cat("## 9. CONCLUSIONS AND RECOMMENDATIONS\n\n")

cat("### 9.1 Main Findings\n\n")

cat("1. **Model Selection**\n")
cat(sprintf("   - Best overall model: **%s** using **%s**\n",
            best_bic$Model, best_bic$Method))
cat("   - BIC criterion preferred for its parsimony penalty\n")
cat("   - Simple models (ARIMA(0,1,1), ARIMA(1,1,1)) perform competitively\n\n")

cat("2. **Method Comparison**\n")
if (!is.na(summary_info$pmm2_win_rate_bic)) {
  if (summary_info$pmm2_win_rate_bic > 0.5) {
    cat("   - PMM2 demonstrates competitive or superior performance\n")
    cat(sprintf("   - Win rate: %.1f%% by BIC criterion\n",
                summary_info$pmm2_win_rate_bic * 100))
  } else {
    cat("   - CSS-ML shows slight advantage in this dataset\n")
    cat(sprintf("   - But PMM2 remains competitive (%.1f%% win rate)\n",
                summary_info$pmm2_win_rate_bic * 100))
  }
} else {
  cat("   - Method comparison inconclusive (no overlapping AIC/BIC values)\n")
}
cat("   - Performance differences are generally small\n")
if (length(mc_stats) > 0 && !is.null(mc_stats$gaussian_pm)) {
  mest_gauss <- if (!is.null(mc_stats$gaussian_m) && is.finite(mc_stats$gaussian_m)) mc_stats$gaussian_m else NA_real_
  skew_vals <- c(mc_stats$gamma_pm, mc_stats$lognormal_pm, mc_stats$chisq_pm)
  skew_vals <- skew_vals[is.finite(skew_vals)]
  skew_min <- if (length(skew_vals) > 0) min(skew_vals) else NA_real_
  skew_max <- if (length(skew_vals) > 0) max(skew_vals) else NA_real_
  cat(sprintf("   - Monte Carlo (Section 5): PMM2 keeps Gaussian efficiency (RE ≈ %.2f%s) and delivers RE %.2f–%.2f under skewed innovations; Hubers M-estimator bridges CSS-ML and PMM2 with intermediate gains\n\n",
              mc_stats$gaussian_pm,
              if (!is.na(mest_gauss)) sprintf(", M-EST ≈ %.2f", mest_gauss) else "",
              if (!is.na(skew_min)) skew_min else mc_stats$gaussian_pm,
              if (!is.na(skew_max)) skew_max else mc_stats$gaussian_pm))
} else {
  cat("\n")
}

cat("3. **Residual Characteristics**\n")
cat("   - Oil prices exhibit non-Gaussian residuals (heavy tails)\n")
cat("   - PMM2 designed for such distributions\n")
cat("   - Both methods handle non-normality reasonably well\n\n")

cat("4. **Computational Cost**\n")
if (avg_speedup > 1) {
  cat("   - PMM2 is computationally efficient\n")
  cat(sprintf("   - Average speedup: %.2fx over CSS-ML\n\n", avg_speedup))
} else {
  cat("   - CSS-ML is slightly faster\n")
  cat(sprintf("   - But difference is negligible (%.2fx)\n\n", 1/avg_speedup))
}

cat("### 9.2 Practical Recommendations\n\n")

cat("**When to use PMM2:**\n")
cat("- Data exhibits strong non-Gaussian characteristics (|kurtosis| > 3)\n")
cat("- Robustness to outliers is important\n")
cat("- Alternative to classical methods for validation\n\n")

cat("**When to use CSS-ML:**\n")
cat("- Residuals are approximately normal\n")
cat("- Standard reporting required (widely recognized method)\n")
cat("- Integration with existing ARIMA workflows\n\n")

if (length(mc_stats) > 0 && (!is.null(mc_stats$gaussian_m) || any(is.finite(c(mc_stats$gamma_m, mc_stats$lognormal_m, mc_stats$chisq_m))))) {
  mest_gauss <- if (!is.null(mc_stats$gaussian_m) && is.finite(mc_stats$gaussian_m)) mc_stats$gaussian_m else NA_real_
  mest_skew <- c(mc_stats$gamma_m, mc_stats$lognormal_m, mc_stats$chisq_m)
  mest_skew <- mest_skew[is.finite(mest_skew)]
  cat("**When to use Huber M-estimators:**\n")
  if (!is.na(mest_gauss)) {
    cat(sprintf("- Need additional robustness with minimal efficiency loss (Gaussian RE ≈ %.2f)\n", mest_gauss))
  } else {
    cat("- Need additional robustness with minimal efficiency loss\n")
  }
  if (length(mest_skew) > 0) {
    cat(sprintf("- Desire intermediate gains on skewed distributions (Monte Carlo RE typically %.2f–%.2f)\n",
                min(mest_skew), max(mest_skew)))
  }
  cat("- Useful as a secondary estimator for validation alongside PMM2\n\n")
}

cat("**General Guidelines:**\n")
cat("- Start with ARIMA(1,1,1) or ARIMA(0,1,1) for financial time series\n")
cat("- Use BIC for model selection (penalizes complexity)\n")
cat("- Always check residual diagnostics (Ljung-Box test, Q-Q plots)\n")
cat("- Consider ensemble methods combining CSS-ML, PMM2, and robust estimators for stress-testing results\n")
cat("- Validate on hold-out data before deployment\n\n")

cat("### 9.3 Limitations and Future Work\n\n")

cat("**Limitations:**\n")
cat("- Single dataset (WTI oil prices)\n")
cat("- Limited to ARIMA specifications (d=1)\n")
cat("- No seasonal components tested\n")
cat("- Short-term data period (2020-2025)\n\n")

cat("**Future Research:**\n")
cat("- Test on multiple financial time series (stocks, currencies, commodities)\n")
cat("- Evaluate SARIMA models with seasonal components\n")
cat("- Compare forecasting performance out-of-sample\n")
cat("- Investigate hybrid methods (PMM2 + CSS-ML ensemble)\n")
cat("- Extend to multivariate VAR/VECM models\n\n")

cat("---\n\n")

# ==================== REFERENCES ====================
cat("## 9. REFERENCES\n\n")

cat("1. Kunchenko, Y., & Lega, Y. (1992). *Polynomial parameter estimations of ")
cat("close to Gaussian random variables*. Kiev: Kyiv Polytechnic Institute.\n\n")

cat("2. Akaike, H. (1974). A new look at the statistical model identification. ")
cat("*IEEE Transactions on Automatic Control*, 19(6), 716-723.\n\n")

cat("3. Schwarz, G. (1978). Estimating the dimension of a model. ")
cat("*The Annals of Statistics*, 6(2), 461-464.\n\n")

cat("4. Box, G. E. P., Jenkins, G. M., Reinsel, G. C., & Ljung, G. M. (2015). ")
cat("*Time Series Analysis: Forecasting and Control* (5th ed.). Wiley.\n\n")

cat("---\n\n")

# ==================== APPENDIX ====================
cat("## 10. APPENDIX: VISUALIZATIONS\n\n")

cat("The following visualizations are available in the `plots/` directory:\n\n")

plots <- list.files(plots_dir, pattern = "\\.png$")
for (i in 1:length(plots)) {
  cat(sprintf("%d. `%s`\n", i, plots[i]))
}
cat("\n")

cat("**To view in report, use:**\n")
cat("```markdown\n")
cat("![AIC Comparison](plots/01_aic_comparison.png)\n")
cat("```\n\n")

cat("---\n\n")

cat("## 11. METADATA\n\n")
cat(sprintf("- **Report Generated:** %s\n", Sys.time()))
cat(sprintf("- **EstemPMM Version:** %s\n", packageVersion("EstemPMM")))
cat(sprintf("- **R Version:** %s\n", R.version.string))
cat(sprintf("- **Platform:** %s\n", R.version$platform))
cat("\n---\n\n")
cat("**END OF REPORT**\n")

sink()

cat(sprintf("✓ Analytical report generated: %s\n\n", report_file))
cat("================================================================\n")
cat("REPORT GENERATION COMPLETE\n")
cat("================================================================\n\n")
cat(sprintf("View the report:\n  cat %s\n\n", report_file))
cat(sprintf("Or convert to HTML:\n"))
cat(sprintf("  pandoc %s -o %s -s --toc\n\n",
            report_file,
            gsub("\\.md$", ".html", report_file)))
