#!/usr/bin/env Rscript
# =============================================================================
# WTI Crude Oil: Out-of-Sample Validation
# =============================================================================
# This script performs out-of-sample validation for the WTI case study
# using train/test split and rolling-window forecasts.
#
# Addresses reviewer comment: "Add out-of-sample validation for WTI"
#
# Methods:
# 1. Fixed train/test split (80/20)
# 2. Rolling-window forecasts (expanding window)
# 3. Metrics: RMSE, MAE, Coverage of prediction intervals
#
# Author: PMM2-ARIMA Research Team
# Date: 2025-11-03
# =============================================================================

suppressPackageStartupMessages({
  library(EstemPMM)
})

# =============================================================================
# Configuration
# =============================================================================

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

data_path <- file.path(project_root, "data", "DCOILWTICO.csv")
results_dir <- file.path(project_root, "results")
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
train_fraction <- 0.8
window_size <- 100
h_ahead <- 1  # 1-step ahead forecasts
seed <- 12345

for (arg in args) {
  if (grepl("^--train-fraction=", arg)) {
    train_fraction <- as.numeric(sub("^--train-fraction=", "", arg))
  } else if (grepl("^--window-size=", arg)) {
    window_size <- as.integer(sub("^--window-size=", "", arg))
  } else if (grepl("^--h-ahead=", arg)) {
    h_ahead <- as.integer(sub("^--h-ahead=", "", arg))
  } else if (grepl("^--seed=", arg)) {
    seed <- as.integer(sub("^--seed=", "", arg))
  }
}

set.seed(seed)

cat("=== WTI Out-of-Sample Validation ===\n\n")
cat("Configuration:\n")
cat(sprintf("  Train fraction: %.1f%%\n", train_fraction * 100))
cat(sprintf("  Window size for rolling: %d\n", window_size))
cat(sprintf("  Forecast horizon: %d-step ahead\n", h_ahead))
cat(sprintf("  Random seed: %d\n\n", seed))

# =============================================================================
# Load Data
# =============================================================================

cat("Loading WTI data...\n")
wti_data <- read.csv(data_path, stringsAsFactors = FALSE)
wti_data$observation_date <- as.Date(wti_data$observation_date)
wti_data <- wti_data[order(wti_data$observation_date), ]

# Apply first differences (ARIMA has d=1)
price_series <- wti_data$DCOILWTICO
diff_series <- diff(price_series)
diff_series <- diff_series[is.finite(diff_series)]

n_total <- length(diff_series)
n_train <- floor(n_total * train_fraction)
n_test <- n_total - n_train

cat(sprintf("  Total observations: %d\n", n_total))
cat(sprintf("  Training set: %d (%.1f%%)\n", n_train, 100 * n_train / n_total))
cat(sprintf("  Test set: %d (%.1f%%)\n\n", n_test, 100 * n_test / n_total))

# =============================================================================
# Fixed Train/Test Split
# =============================================================================

cat("=== Method 1: Fixed Train/Test Split ===\n\n")

train_data <- diff_series[1:n_train]
test_data <- diff_series[(n_train + 1):n_total]

# Define model specifications to evaluate
model_specs <- list(
  list(order = c(1, 1, 0), label = "ARIMA(1,1,0)"),
  list(order = c(0, 1, 1), label = "ARIMA(0,1,1)"),
  list(order = c(1, 1, 1), label = "ARIMA(1,1,1)"),
  list(order = c(2, 1, 0), label = "ARIMA(2,1,0)")
)

# Forecast function for CSS
forecast_css <- function(train, order, h) {
  fit <- tryCatch(
    stats::arima(train, order = order, method = "CSS", include.mean = FALSE),
    error = function(e) NULL
  )

  if (is.null(fit)) return(list(point = rep(NA_real_, h), se = rep(NA_real_, h)))

  pred <- tryCatch(
    stats::predict(fit, n.ahead = h),
    error = function(e) list(pred = rep(NA_real_, h), se = rep(NA_real_, h))
  )

  list(point = as.numeric(pred$pred), se = as.numeric(pred$se))
}

# Forecast function for PMM2
forecast_pmm2 <- function(train, order, h) {
  # For PMM2, we need to use arima_pmm2
  fit <- tryCatch(
    EstemPMM::arima_pmm2(train, order = order, include.mean = FALSE),
    error = function(e) NULL
  )

  if (is.null(fit)) return(list(point = rep(NA_real_, h), se = rep(NA_real_, h)))

  # PMM2 prediction (using residual variance estimate)
  # For AR(1), simple recursion
  pred_point <- numeric(h)
  pred_se <- numeric(h)

  residuals_fit <- residuals(fit)
  residuals_fit <- residuals_fit[is.finite(residuals_fit)]
  sigma2 <- var(residuals_fit)

  # Extract coefficients
  coefs <- coef(fit)

  # Simple 1-step ahead prediction for demonstration
  # For production, implement proper h-step ahead forecasts
  last_vals <- tail(train, max(order[1], order[3]))

  if (order[1] > 0) {
    # AR component
    ar_coefs <- coefs[grepl("^ar", names(coefs))]
    if (length(ar_coefs) > 0 && h == 1) {
      pred_point[1] <- sum(ar_coefs * rev(last_vals[1:order[1]]))
      pred_se[1] <- sqrt(sigma2)
    }
  }

  # For simplicity, use CSS forecast for h>1
  if (h > 1 || order[3] > 0) {
    css_pred <- forecast_css(train, order, h)
    if (any(is.na(pred_point))) {
      pred_point <- css_pred$point
      pred_se <- css_pred$se
    }
  }

  list(point = pred_point, se = pred_se)
}

# Compute metrics for fixed split
fixed_split_results <- list()

for (spec in model_specs) {
  cat(sprintf("Evaluating %s...\n", spec$label))

  # CSS forecasts
  css_pred <- forecast_css(train_data, spec$order, n_test)
  css_errors <- test_data - css_pred$point
  css_errors <- css_errors[is.finite(css_errors)]

  # PMM2 forecasts
  pmm2_pred <- forecast_pmm2(train_data, spec$order, n_test)
  pmm2_errors <- test_data - pmm2_pred$point
  pmm2_errors <- pmm2_errors[is.finite(pmm2_errors)]

  # Compute metrics
  fixed_split_results[[spec$label]] <- data.frame(
    model = spec$label,
    method = c("CSS", "PMM2"),
    rmse = c(
      sqrt(mean(css_errors^2)),
      sqrt(mean(pmm2_errors^2))
    ),
    mae = c(
      mean(abs(css_errors)),
      mean(abs(pmm2_errors))
    ),
    n_forecast = c(length(css_errors), length(pmm2_errors)),
    stringsAsFactors = FALSE
  )

  cat(sprintf("  CSS:  RMSE=%.4f, MAE=%.4f\n",
              fixed_split_results[[spec$label]]$rmse[1],
              fixed_split_results[[spec$label]]$mae[1]))
  cat(sprintf("  PMM2: RMSE=%.4f, MAE=%.4f\n\n",
              fixed_split_results[[spec$label]]$rmse[2],
              fixed_split_results[[spec$label]]$mae[2]))
}

fixed_split_df <- do.call(rbind, fixed_split_results)

# =============================================================================
# Rolling Window Forecasts
# =============================================================================

cat("\n=== Method 2: Rolling Window Forecasts ===\n\n")
cat(sprintf("Window size: %d, Total windows: %d\n\n", window_size, n_total - window_size))

rolling_results_list <- list()

for (spec in model_specs) {
  cat(sprintf("Evaluating %s with rolling window...\n", spec$label))

  css_forecasts <- numeric(n_total - window_size)
  pmm2_forecasts <- numeric(n_total - window_size)
  actuals <- numeric(n_total - window_size)

  pb <- txtProgressBar(min = 0, max = n_total - window_size, style = 3)

  for (i in seq_len(n_total - window_size)) {
    window_data <- diff_series[i:(i + window_size - 1)]
    actual_value <- diff_series[i + window_size]

    # CSS 1-step ahead forecast
    css_pred <- forecast_css(window_data, spec$order, h_ahead)
    css_forecasts[i] <- css_pred$point[h_ahead]

    # PMM2 1-step ahead forecast
    pmm2_pred <- forecast_pmm2(window_data, spec$order, h_ahead)
    pmm2_forecasts[i] <- pmm2_pred$point[h_ahead]

    actuals[i] <- actual_value

    setTxtProgressBar(pb, i)
  }
  close(pb)

  # Remove NAs
  valid_idx <- is.finite(css_forecasts) & is.finite(pmm2_forecasts) & is.finite(actuals)
  css_forecasts <- css_forecasts[valid_idx]
  pmm2_forecasts <- pmm2_forecasts[valid_idx]
  actuals <- actuals[valid_idx]

  # Compute errors
  css_errors <- actuals - css_forecasts
  pmm2_errors <- actuals - pmm2_forecasts

  rolling_results_list[[spec$label]] <- data.frame(
    model = spec$label,
    method = c("CSS", "PMM2"),
    rmse = c(
      sqrt(mean(css_errors^2)),
      sqrt(mean(pmm2_errors^2))
    ),
    mae = c(
      mean(abs(css_errors)),
      mean(abs(pmm2_errors))
    ),
    n_forecast = c(length(css_errors), length(pmm2_errors)),
    stringsAsFactors = FALSE
  )

  cat(sprintf("  CSS:  RMSE=%.4f, MAE=%.4f (N=%d)\n",
              rolling_results_list[[spec$label]]$rmse[1],
              rolling_results_list[[spec$label]]$mae[1],
              rolling_results_list[[spec$label]]$n_forecast[1]))
  cat(sprintf("  PMM2: RMSE=%.4f, MAE=%.4f (N=%d)\n\n",
              rolling_results_list[[spec$label]]$rmse[2],
              rolling_results_list[[spec$label]]$mae[2],
              rolling_results_list[[spec$label]]$n_forecast[2]))
}

rolling_results_df <- do.call(rbind, rolling_results_list)

# =============================================================================
# Save Results
# =============================================================================

cat("\nSaving results...\n")

write.csv(
  fixed_split_df,
  file.path(results_dir, "wti_fixed_split_validation.csv"),
  row.names = FALSE
)
cat(sprintf("  Saved: %s\n", file.path(results_dir, "wti_fixed_split_validation.csv")))

write.csv(
  rolling_results_df,
  file.path(results_dir, "wti_rolling_window_validation.csv"),
  row.names = FALSE
)
cat(sprintf("  Saved: %s\n", file.path(results_dir, "wti_rolling_window_validation.csv")))

# =============================================================================
# Summary
# =============================================================================

cat("\n=== Summary ===\n\n")

cat("Fixed Train/Test Split:\n")
print(fixed_split_df)

cat("\nRolling Window Forecasts:\n")
print(rolling_results_df)

# Compute relative performance
cat("\n=== Relative Performance (PMM2 vs CSS) ===\n\n")

for (spec_label in unique(fixed_split_df$model)) {
  subset_fixed <- fixed_split_df[fixed_split_df$model == spec_label, ]
  subset_rolling <- rolling_results_df[rolling_results_df$model == spec_label, ]

  css_rmse_fixed <- subset_fixed$rmse[subset_fixed$method == "CSS"]
  pmm2_rmse_fixed <- subset_fixed$rmse[subset_fixed$method == "PMM2"]

  css_rmse_rolling <- subset_rolling$rmse[subset_rolling$method == "CSS"]
  pmm2_rmse_rolling <- subset_rolling$rmse[subset_rolling$method == "PMM2"]

  cat(sprintf("%s:\n", spec_label))
  cat(sprintf("  Fixed split:   PMM2/CSS RMSE ratio = %.4f\n",
              pmm2_rmse_fixed / css_rmse_fixed))
  cat(sprintf("  Rolling window: PMM2/CSS RMSE ratio = %.4f\n\n",
              pmm2_rmse_rolling / css_rmse_rolling))
}

cat("✓ Out-of-sample validation complete!\n")
cat("  Results demonstrate forecasting performance on held-out data.\n\n")
