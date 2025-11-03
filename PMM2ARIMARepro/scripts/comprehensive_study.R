# Comprehensive ARIMA Study: PMM2 vs Classical Methods
# This script conducts a rigorous comparison on WTI Oil Price Data

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

library(EstemPMM)

compute_pmm2_logLik <- function(model) {
  # Manual log-likelihood under Gaussian assumption for PMM2 fits
  res <- model@residuals
  res_clean <- res[is.finite(res)]
  n <- length(res_clean)
  if (n == 0) {
    return(NA_real_)
  }
  sigma2 <- mean(res_clean^2)
  if (!is.finite(sigma2) || sigma2 <= 0) {
    return(NA_real_)
  }
  -0.5 * n * (log(2 * pi) + 1 + log(sigma2))
}

safe_logLik <- function(model) {
  if (methods::is(model, "TS2fit")) {
    return(compute_pmm2_logLik(model))
  }
  value <- tryCatch({
    logLik(model)
  }, error = function(e) NA_real_)
  as.numeric(value)
}

safe_AIC <- function(model) {
  if (methods::is(model, "TS2fit")) {
    ll <- safe_logLik(model)
    if (!is.finite(ll)) {
      return(NA_real_)
    }
    p <- length(model@coefficients)
    return(-2 * ll + 2 * p)
  }
  tryCatch(as.numeric(AIC(model)), error = function(e) NA_real_)
}

safe_BIC <- function(model) {
  if (methods::is(model, "TS2fit")) {
    ll <- safe_logLik(model)
    if (!is.finite(ll)) {
      return(NA_real_)
    }
    n <- sum(is.finite(model@residuals))
    if (n <= 0) {
      return(NA_real_)
    }
    p <- length(model@coefficients)
    return(-2 * ll + log(n) * p)
  }
  tryCatch(as.numeric(BIC(model)), error = function(e) NA_real_)
}

# Output directory for results
output_dir <- "results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cat("================================================================\n")
cat("COMPREHENSIVE ARIMA COMPARISON STUDY\n")
cat("PMM2 vs Classical Methods on WTI Crude Oil Prices\n")
cat("================================================================\n\n")

# ==================== DATA PREPARATION ====================
cat("1. DATA PREPARATION\n")
cat("----------------------------------------------------------\n")

data_path <- system.file("extdata", "DCOILWTICO.csv", package = "EstemPMM")
if (!file.exists(data_path) || data_path == "") {
  data_path <- file.path("data", "DCOILWTICO.csv")
}

if (!file.exists(data_path)) {
  stop("Data file DCOILWTICO.csv not found")
}

oil_data <- read.csv(data_path, stringsAsFactors = FALSE)
oil_data$observation_date <- as.Date(oil_data$observation_date)
oil_data$DCOILWTICO <- as.numeric(oil_data$DCOILWTICO)

# Clean data
oil_clean <- oil_data[!is.na(oil_data$DCOILWTICO), ]
prices <- oil_clean$DCOILWTICO
dates <- oil_clean$observation_date

cat(sprintf("Date range: %s to %s\n", min(dates), max(dates)))
cat(sprintf("Total observations: %d\n", length(prices)))
cat(sprintf("Missing values: %d (%.2f%%)\n\n",
            nrow(oil_data) - nrow(oil_clean),
            100 * (nrow(oil_data) - nrow(oil_clean)) / nrow(oil_data)))

# ==================== DESCRIPTIVE STATISTICS ====================
cat("2. DESCRIPTIVE STATISTICS\n")
cat("----------------------------------------------------------\n")

desc_stats <- data.frame(
  Statistic = c("Mean", "Median", "Std Dev", "Min", "Max",
                "Q1", "Q3", "Skewness", "Kurtosis", "CV"),
  Value = c(
    mean(prices),
    median(prices),
    sd(prices),
    min(prices),
    max(prices),
    quantile(prices, 0.25),
    quantile(prices, 0.75),
    pmm_skewness(prices),
    pmm_kurtosis(prices),
    sd(prices) / mean(prices)
  )
)

print(desc_stats, digits = 4)

# Save to CSV
write.csv(desc_stats, file.path(output_dir, "descriptive_stats.csv"), row.names = FALSE)

# ==================== STATIONARITY TESTS ====================
cat("\n3. STATIONARITY ANALYSIS\n")
cat("----------------------------------------------------------\n")

# ADF test on original series
adf_original <- tseries::adf.test(prices)
cat(sprintf("ADF Test (Original Series):\n"))
cat(sprintf("  Statistic: %.4f\n", adf_original$statistic))
cat(sprintf("  p-value: %.4f\n", adf_original$p.value))
cat(sprintf("  Conclusion: %s\n",
            ifelse(adf_original$p.value < 0.05, "Stationary", "Non-stationary")))

# ADF test on differenced series
diff_prices <- diff(prices)
adf_diff <- tseries::adf.test(diff_prices)
cat(sprintf("\nADF Test (First Difference):\n"))
cat(sprintf("  Statistic: %.4f\n", adf_diff$statistic))
cat(sprintf("  p-value: %.4f\n", adf_diff$p.value))
cat(sprintf("  Conclusion: %s\n\n",
            ifelse(adf_diff$p.value < 0.05, "Stationary", "Non-stationary")))

# ==================== MODEL SPECIFICATIONS ====================
cat("4. MODEL SPECIFICATIONS TO TEST\n")
cat("----------------------------------------------------------\n")

model_specs <- list(
  "ARIMA(0,1,1)" = c(0, 1, 1),
  "ARIMA(1,1,0)" = c(1, 1, 0),
  "ARIMA(1,1,1)" = c(1, 1, 1),
  "ARIMA(2,1,1)" = c(2, 1, 1),
  "ARIMA(1,1,2)" = c(1, 1, 2),
  "ARIMA(2,1,2)" = c(2, 1, 2)
)

for (name in names(model_specs)) {
  order <- model_specs[[name]]
  cat(sprintf("  %s: p=%d, d=%d, q=%d\n", name, order[1], order[2], order[3]))
}

# ==================== MODEL FITTING ====================
cat("\n5. FITTING MODELS\n")
cat("----------------------------------------------------------\n")

results <- data.frame(
  Model = character(),
  Method = character(),
  AIC = numeric(),
  BIC = numeric(),
  LogLik = numeric(),
  RSS = numeric(),
  MAE = numeric(),
  RMSE = numeric(),
  MAPE = numeric(),
  Skewness = numeric(),
  Kurtosis = numeric(),
  Convergence = logical(),
  Time_sec = numeric(),
  stringsAsFactors = FALSE
)

# Store fitted models for later analysis
fitted_models <- list()

for (model_name in names(model_specs)) {
  order <- model_specs[[model_name]]
  cat(sprintf("\nFitting %s...\n", model_name))

  # CSS-ML Method
  cat("  Method: CSS-ML... ")
  time_css <- system.time({
    fit_css <- tryCatch({
      arima(prices, order = order, method = "CSS-ML", include.mean = FALSE)
    }, error = function(e) {
      cat(sprintf("FAILED (%s)\n", e$message))
      NULL
    })
  })

  if (!is.null(fit_css)) {
    cat("SUCCESS\n")
    res_css <- residuals(fit_css)
    res_css_clean <- res_css[is.finite(res_css)]

    # Calculate metrics
    fitted_vals <- prices - res_css
    fitted_vals_clean <- fitted_vals[is.finite(fitted_vals)]
    actual_vals <- prices[is.finite(fitted_vals)]

    mape_css <- mean(abs((actual_vals - fitted_vals_clean) / actual_vals)) * 100

    results <- rbind(results, data.frame(
      Model = model_name,
      Method = "CSS-ML",
      AIC = safe_AIC(fit_css),
      BIC = safe_BIC(fit_css),
      LogLik = safe_logLik(fit_css),
      RSS = sum(res_css_clean^2),
      MAE = mean(abs(res_css_clean)),
      RMSE = sqrt(mean(res_css_clean^2)),
      MAPE = mape_css,
      Skewness = pmm_skewness(res_css_clean),
      Kurtosis = pmm_kurtosis(res_css_clean),
      Convergence = TRUE,
      Time_sec = time_css[3]
    ))

    fitted_models[[paste0(model_name, "_CSS")]] <- fit_css
  }

  # PMM2 Method
  cat("  Method: PMM2... ")
  time_pmm2 <- system.time({
    fit_pmm2 <- tryCatch({
      arima_pmm2(prices, order = order, include.mean = FALSE, verbose = FALSE)
    }, error = function(e) {
      cat(sprintf("FAILED (%s)\n", e$message))
      NULL
    })
  })

  if (!is.null(fit_pmm2)) {
    cat("SUCCESS\n")
    res_pmm2 <- fit_pmm2@residuals
    res_pmm2_clean <- res_pmm2[is.finite(res_pmm2)]

    # Calculate metrics
    fitted_vals <- prices[1:length(res_pmm2)] - res_pmm2
    fitted_vals_clean <- fitted_vals[is.finite(fitted_vals)]
    actual_vals <- prices[1:length(fitted_vals)][is.finite(fitted_vals)]

    mape_pmm2 <- mean(abs((actual_vals - fitted_vals_clean) / actual_vals)) * 100

    results <- rbind(results, data.frame(
      Model = model_name,
      Method = "PMM2",
      AIC = safe_AIC(fit_pmm2),
      BIC = safe_BIC(fit_pmm2),
      LogLik = safe_logLik(fit_pmm2),
      RSS = sum(res_pmm2_clean^2),
      MAE = mean(abs(res_pmm2_clean)),
      RMSE = sqrt(mean(res_pmm2_clean^2)),
      MAPE = mape_pmm2,
      Skewness = pmm_skewness(res_pmm2_clean),
      Kurtosis = pmm_kurtosis(res_pmm2_clean),
      Convergence = fit_pmm2@convergence,
      Time_sec = time_pmm2[3]
    ))

    fitted_models[[paste0(model_name, "_PMM2")]] <- fit_pmm2
  }
}

# ==================== RESULTS SUMMARY ====================
cat("\n6. COMPREHENSIVE RESULTS\n")
cat("----------------------------------------------------------\n\n")

# Sort by BIC
results_sorted <- results[order(results$BIC), ]
print(results_sorted, digits = 4, row.names = FALSE)

# Save full results
write.csv(results_sorted, file.path(output_dir, "full_results.csv"), row.names = FALSE)

# ==================== BEST MODEL IDENTIFICATION ====================
cat("\n7. BEST MODEL SELECTION\n")
cat("----------------------------------------------------------\n")

best_aic <- results_sorted[which.min(results_sorted$AIC), ]
best_bic <- results_sorted[which.min(results_sorted$BIC), ]

cat(sprintf("\nBest by AIC: %s (%s)\n", best_aic$Model, best_aic$Method))
cat(sprintf("  AIC: %.2f, BIC: %.2f, RMSE: %.4f\n",
            best_aic$AIC, best_aic$BIC, best_aic$RMSE))

cat(sprintf("\nBest by BIC: %s (%s)\n", best_bic$Model, best_bic$Method))
cat(sprintf("  AIC: %.2f, BIC: %.2f, RMSE: %.4f\n\n",
            best_bic$BIC, best_bic$BIC, best_bic$RMSE))

# ==================== METHOD COMPARISON ====================
cat("8. METHOD COMPARISON (PMM2 vs CSS-ML)\n")
cat("----------------------------------------------------------\n\n")

comparison_table <- data.frame(
  Model = character(),
  AIC_Diff = numeric(),
  BIC_Diff = numeric(),
  RMSE_Diff = numeric(),
  Winner_AIC = character(),
  Winner_BIC = character(),
  stringsAsFactors = FALSE
)

for (model_name in names(model_specs)) {
  css_row <- results[results$Model == model_name & results$Method == "CSS-ML", ]
  pmm2_row <- results[results$Model == model_name & results$Method == "PMM2", ]

  if (nrow(css_row) > 0 && nrow(pmm2_row) > 0) {
    aic_diff <- pmm2_row$AIC - css_row$AIC
    bic_diff <- pmm2_row$BIC - css_row$BIC
    rmse_diff <- pmm2_row$RMSE - css_row$RMSE

    comparison_table <- rbind(comparison_table, data.frame(
      Model = model_name,
      AIC_Diff = aic_diff,
      BIC_Diff = bic_diff,
      RMSE_Diff = rmse_diff,
      Winner_AIC = ifelse(aic_diff < 0, "PMM2", "CSS-ML"),
      Winner_BIC = ifelse(bic_diff < 0, "PMM2", "CSS-ML")
    ))
  }
}

print(comparison_table, digits = 4, row.names = FALSE)
write.csv(comparison_table, file.path(output_dir, "method_comparison.csv"), row.names = FALSE)

# Calculate win rates (exclude models without valid AIC/BIC differences)
valid_aic <- !is.na(comparison_table$AIC_Diff)
valid_bic <- !is.na(comparison_table$BIC_Diff)

pmm2_aic_wins <- sum(comparison_table$Winner_AIC[valid_aic] == "PMM2", na.rm = TRUE)
pmm2_bic_wins <- sum(comparison_table$Winner_BIC[valid_bic] == "PMM2", na.rm = TRUE)

total_aic_comparisons <- sum(valid_aic)
total_bic_comparisons <- sum(valid_bic)

pmm2_win_rate_aic <- if (total_aic_comparisons > 0) pmm2_aic_wins / total_aic_comparisons else NA_real_
pmm2_win_rate_bic <- if (total_bic_comparisons > 0) pmm2_bic_wins / total_bic_comparisons else NA_real_

cat(sprintf("\nPMM2 Win Rate:\n"))
if (!is.na(pmm2_win_rate_aic)) {
  cat(sprintf("  AIC: %d/%d (%.1f%%)\n",
              pmm2_aic_wins, total_aic_comparisons, 100 * pmm2_win_rate_aic))
} else {
  cat("  AIC: Not available (insufficient comparable models)\n")
}
if (!is.na(pmm2_win_rate_bic)) {
  cat(sprintf("  BIC: %d/%d (%.1f%%)\n\n",
              pmm2_bic_wins, total_bic_comparisons, 100 * pmm2_win_rate_bic))
} else {
  cat("  BIC: Not available (insufficient comparable models)\n\n")
}

# ==================== PARAMETER MSE RATIO ====================
cat("9. PARAMETER ESTIMATE MSE RATIO (PMM2 vs CSS-ML)\n")
cat("----------------------------------------------------------\n")

param_mse_results <- data.frame(
  Model = character(),
  CSS_Skewness = numeric(),
  CSS_Kurtosis = numeric(),
  PMM2_CSS_MSE_Ratio = numeric(),
  Interpretation = character(),
  stringsAsFactors = FALSE
)

for (model_name in names(model_specs)) {
  css_row <- results[results$Model == model_name & results$Method == "CSS-ML", ]

  if (nrow(css_row) > 0) {
    skew <- css_row$Skewness
    kurt <- css_row$Kurtosis
    denom <- 2 + kurt

    ratio <- if (!is.na(skew) && !is.na(kurt) && abs(denom) > 1e-8) {
      1 - (skew^2) / denom
    } else {
      NA_real_
    }

    interpretation <- if (is.na(ratio)) {
      "Undefined"
    } else if (ratio < 1) {
      "PMM2 lower parameter MSE"
    } else if (ratio > 1) {
      "CSS-ML lower parameter MSE"
    } else {
      "Equal"
    }

    param_mse_results <- rbind(param_mse_results, data.frame(
      Model = model_name,
      CSS_Skewness = skew,
      CSS_Kurtosis = kurt,
      PMM2_CSS_MSE_Ratio = ratio,
      Interpretation = interpretation,
      stringsAsFactors = FALSE
    ))

    if (!is.na(ratio)) {
      cat(sprintf("%s: ratio = %.3f (%s)\n", model_name, ratio, interpretation))
    } else {
      cat(sprintf("%s: ratio not available (insufficient moment information)\n", model_name))
    }
  }
}

write.csv(param_mse_results,
          file.path(output_dir, "parameter_mse_ratio.csv"),
          row.names = FALSE)

cat("\nNote: Ratio < 1 implies PMM2 is expected to deliver lower coefficient MSE under observed skewness/kurtosis.\n\n")

# ==================== RESIDUAL ANALYSIS ====================
cat("10. RESIDUAL ANALYSIS OF BEST MODEL\n")
cat("----------------------------------------------------------\n")

best_model_key <- paste0(best_bic$Model, "_", gsub("-", "", best_bic$Method))
best_fit <- fitted_models[[best_model_key]]

if (!is.null(best_fit)) {
  if (best_bic$Method == "CSS-ML") {
    res <- residuals(best_fit)
  } else {
    res <- best_fit@residuals
  }
  res_clean <- res[is.finite(res)]

  # Ljung-Box test
  lb_test <- Box.test(res_clean, lag = 20, type = "Ljung-Box")
  cat(sprintf("Ljung-Box Test (lag=20):\n"))
  cat(sprintf("  Statistic: %.4f\n", lb_test$statistic))
  cat(sprintf("  p-value: %.4f\n", lb_test$p.value))
  cat(sprintf("  Conclusion: %s\n\n",
              ifelse(lb_test$p.value > 0.05, "No autocorrelation (good)", "Autocorrelation present")))

  # Shapiro-Wilk normality test (on sample if too large)
  res_sample <- if (length(res_clean) > 5000) sample(res_clean, 5000) else res_clean
  sw_test <- shapiro.test(res_sample)
  cat(sprintf("Shapiro-Wilk Normality Test:\n"))
  cat(sprintf("  Statistic: %.4f\n", sw_test$statistic))
  cat(sprintf("  p-value: %.4f\n", sw_test$p.value))
  cat(sprintf("  Conclusion: %s\n\n",
              ifelse(sw_test$p.value > 0.05, "Normal distribution", "Non-normal distribution")))
}

# ==================== SUMMARY STATISTICS BY METHOD ====================
cat("11. SUMMARY STATISTICS BY METHOD\n")
cat("----------------------------------------------------------\n\n")

method_summary <- aggregate(cbind(AIC, BIC, RMSE, MAE, Time_sec) ~ Method,
                             data = results,
                             FUN = function(x) c(mean = mean(x), sd = sd(x)))

cat("Average Performance (Mean ± SD):\n\n")
for (method in unique(results$Method)) {
  cat(sprintf("%s:\n", method))
  method_data <- results[results$Method == method, ]
  cat(sprintf("  AIC:   %.2f ± %.2f\n", mean(method_data$AIC), sd(method_data$AIC)))
  cat(sprintf("  BIC:   %.2f ± %.2f\n", mean(method_data$BIC), sd(method_data$BIC)))
  cat(sprintf("  RMSE:  %.4f ± %.4f\n", mean(method_data$RMSE), sd(method_data$RMSE)))
  cat(sprintf("  Time:  %.4f ± %.4f sec\n\n", mean(method_data$Time_sec), sd(method_data$Time_sec)))
}

# ==================== SAVE SUMMARY ====================
cat("12. SAVING RESULTS\n")
cat("----------------------------------------------------------\n")

summary_info <- list(
  date = Sys.time(),
  data_file = data_path,
  n_observations = length(prices),
  n_models_tested = length(model_specs),
  best_model_aic = paste0(best_aic$Model, " (", best_aic$Method, ")"),
  best_model_bic = paste0(best_bic$Model, " (", best_bic$Method, ")"),
  pmm2_win_rate_aic = pmm2_win_rate_aic,
  pmm2_win_rate_bic = pmm2_win_rate_bic,
  aic_comparisons = total_aic_comparisons,
  bic_comparisons = total_bic_comparisons,
  parameter_mse_ratio = param_mse_results
)

saveRDS(summary_info, file.path(output_dir, "study_summary.rds"))
saveRDS(results, file.path(output_dir, "all_results.rds"))
saveRDS(fitted_models, file.path(output_dir, "fitted_models.rds"))

cat(sprintf("  ✓ Full results saved to: %s\n", file.path(output_dir, "full_results.csv")))
cat(sprintf("  ✓ Comparison table saved to: %s\n", file.path(output_dir, "method_comparison.csv")))
cat(sprintf("  ✓ R objects saved to: %s\n", output_dir))

cat("\n================================================================\n")
cat("STUDY COMPLETED SUCCESSFULLY\n")
cat("================================================================\n")
cat(sprintf("\nGenerate visualizations by running:\n"))
cat(sprintf("  source('create_visualizations.R')\n\n"))
cat(sprintf("Generate final report by running:\n"))
cat(sprintf("  source('generate_report.R')\n\n"))
