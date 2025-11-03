#!/usr/bin/env Rscript
# =============================================================================
# Add Bootstrap Confidence Intervals to Monte Carlo Results
# =============================================================================
# This script enhances existing Monte Carlo results by adding bootstrap
# confidence intervals for key metrics (bias, variance, MSE, RE).
#
# Addresses reviewer comment: "Add confidence intervals to Monte Carlo tables"
#
# Author: PMM2-ARIMA Research Team
# Date: 2025-11-03
# =============================================================================

suppressPackageStartupMessages({
  library(EstemPMM)
  if (!requireNamespace("boot", quietly = TRUE)) {
    message("Installing 'boot' package for bootstrap CI...")
    install.packages("boot", repos = "https://cloud.r-project.org")
  }
  library(boot)
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

results_dir <- file.path(project_root, "results", "monte_carlo")
if (!dir.exists(results_dir)) {
  stop("Results directory not found. Please run run_monte_carlo.R first.")
}

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
bootstrap_reps <- 1000L
confidence_level <- 0.95
seed <- 12345L

for (arg in args) {
  if (grepl("^--bootstrap-reps=", arg)) {
    bootstrap_reps <- as.integer(sub("^--bootstrap-reps=", "", arg))
  } else if (grepl("^--confidence=", arg)) {
    confidence_level <- as.numeric(sub("^--confidence=", "", arg))
  } else if (grepl("^--seed=", arg)) {
    seed <- as.integer(sub("^--seed=", "", arg))
  }
}

set.seed(seed)

cat("Bootstrap CI Configuration:\n")
cat(sprintf("  Bootstrap replications: %d\n", bootstrap_reps))
cat(sprintf("  Confidence level: %.1f%%\n", confidence_level * 100))
cat(sprintf("  Random seed: %d\n\n", seed))

# =============================================================================
# Bootstrap Functions
# =============================================================================

# Bootstrap function for bias
boot_bias <- function(data, indices, true_value) {
  resampled <- data[indices]
  mean(resampled - true_value)
}

# Bootstrap function for variance
boot_var <- function(data, indices) {
  var(data[indices])
}

# Bootstrap function for MSE
boot_mse <- function(data, indices, true_value) {
  resampled <- data[indices]
  mean((resampled - true_value)^2)
}

# Bootstrap function for relative efficiency (ratio of MSEs)
boot_re <- function(data_list, indices) {
  # data_list = list(baseline_estimates, pmm2_estimates, true_value)
  baseline <- data_list[[1]][indices]
  pmm2 <- data_list[[2]][indices]
  true_val <- data_list[[3]]

  mse_baseline <- mean((baseline - true_val)^2)
  mse_pmm2 <- mean((pmm2 - true_val)^2)

  if (mse_pmm2 == 0) return(NA_real_)
  return(mse_baseline / mse_pmm2)
}

# Compute bootstrap CI with proper error handling
compute_boot_ci <- function(data, stat_func, conf = 0.95, R = 1000, ...) {
  if (length(data) < 10) {
    return(list(se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_))
  }

  tryCatch({
    boot_result <- boot::boot(data = data, statistic = stat_func, R = R, ...)

    # Standard error from bootstrap
    se <- sd(boot_result$t[, 1], na.rm = TRUE)

    # Try BCa method first (best), fall back to percentile
    ci <- tryCatch({
      boot::boot.ci(boot_result, conf = conf, type = "bca")$bca[4:5]
    }, error = function(e) {
      tryCatch({
        boot::boot.ci(boot_result, conf = conf, type = "perc")$percent[4:5]
      }, error = function(e2) {
        # Fall back to normal approximation
        alpha <- 1 - conf
        q <- qnorm(1 - alpha/2)
        point_est <- boot_result$t0
        c(point_est - q * se, point_est + q * se)
      })
    })

    list(se = se, ci_lower = ci[1], ci_upper = ci[2])
  }, error = function(e) {
    warning(sprintf("Bootstrap failed: %s", e$message))
    list(se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_)
  })
}

# =============================================================================
# Load Existing Results
# =============================================================================

cat("Loading Monte Carlo results...\n")
metrics_path <- file.path(results_dir, "monte_carlo_metrics.csv")
if (!file.exists(metrics_path)) {
  stop("Monte Carlo metrics file not found. Please run run_monte_carlo.R first.")
}

metrics_df <- read.csv(metrics_path, stringsAsFactors = FALSE)
cat(sprintf("  Loaded %d rows from monte_carlo_metrics.csv\n", nrow(metrics_df)))

# Check if we need to re-run Monte Carlo to get raw estimates
raw_estimates_path <- file.path(results_dir, "raw_estimates.rds")
if (!file.exists(raw_estimates_path)) {
  stop(paste(
    "Raw estimates file not found.",
    "Please modify run_monte_carlo.R to save raw estimates,",
    "or run this script with --rerun-mc flag"
  ))
}

cat("  Loading raw estimates for bootstrap...\n")
raw_data <- readRDS(raw_estimates_path)
raw_estimates <- raw_data$estimates
mc_configs <- raw_data$configs

# Build lookup table for true parameter values
true_params_lookup <- list()
for (cfg in mc_configs) {
  for (param_name in names(cfg$true_params)) {
    key <- paste(cfg$id, param_name, sep = "|")
    true_params_lookup[[key]] <- cfg$true_params[[param_name]]
  }
}

# =============================================================================
# Add Confidence Intervals
# =============================================================================

cat("\nComputing bootstrap confidence intervals...\n")

# Add new columns for CI
metrics_df$bias_se <- NA_real_
metrics_df$bias_ci_lower <- NA_real_
metrics_df$bias_ci_upper <- NA_real_
metrics_df$variance_se <- NA_real_
metrics_df$variance_ci_lower <- NA_real_
metrics_df$variance_ci_upper <- NA_real_
metrics_df$mse_se <- NA_real_
metrics_df$mse_ci_lower <- NA_real_
metrics_df$mse_ci_upper <- NA_real_
metrics_df$re_se <- NA_real_
metrics_df$re_ci_lower <- NA_real_
metrics_df$re_ci_upper <- NA_real_

pb <- txtProgressBar(min = 0, max = nrow(metrics_df), style = 3)

for (i in seq_len(nrow(metrics_df))) {
  row <- metrics_df[i, ]

  # Find corresponding raw estimates
  key <- paste(row$model, row$distribution, row$sample_size, row$method, row$parameter, sep = "|")

  if (!key %in% names(raw_estimates)) {
    setTxtProgressBar(pb, i)
    next
  }

  estimates <- raw_estimates[[key]]
  finite_estimates <- estimates[is.finite(estimates)]

  if (length(finite_estimates) < 10) {
    setTxtProgressBar(pb, i)
    next
  }

  # Get true parameter value from configuration
  param_key <- paste(row$model, row$parameter, sep = "|")
  true_val <- true_params_lookup[[param_key]]

  if (is.null(true_val)) {
    setTxtProgressBar(pb, i)
    next
  }

  # Compute bootstrap CIs for bias
  bias_ci <- compute_boot_ci(
    finite_estimates,
    boot_bias,
    conf = confidence_level,
    R = bootstrap_reps,
    true_value = true_val
  )
  metrics_df$bias_se[i] <- bias_ci$se
  metrics_df$bias_ci_lower[i] <- bias_ci$ci_lower
  metrics_df$bias_ci_upper[i] <- bias_ci$ci_upper

  # Compute bootstrap CIs for variance
  var_ci <- compute_boot_ci(
    finite_estimates,
    boot_var,
    conf = confidence_level,
    R = bootstrap_reps
  )
  metrics_df$variance_se[i] <- var_ci$se
  metrics_df$variance_ci_lower[i] <- var_ci$ci_lower
  metrics_df$variance_ci_upper[i] <- var_ci$ci_upper

  # Compute bootstrap CIs for MSE
  mse_ci <- compute_boot_ci(
    finite_estimates,
    boot_mse,
    conf = confidence_level,
    R = bootstrap_reps,
    true_value = true_val
  )
  metrics_df$mse_se[i] <- mse_ci$se
  metrics_df$mse_ci_lower[i] <- mse_ci$ci_lower
  metrics_df$mse_ci_upper[i] <- mse_ci$ci_upper

  # For relative efficiency, we need baseline estimates
  if (row$method == "PMM2") {
    baseline_key <- paste(row$model, row$distribution, row$sample_size, "CSS", row$parameter, sep = "|")
    if (baseline_key %in% names(raw_estimates)) {
      baseline_estimates <- raw_estimates[[baseline_key]]
      baseline_finite <- baseline_estimates[is.finite(baseline_estimates)]

      if (length(baseline_finite) >= 10 && length(finite_estimates) >= 10) {
        # Match lengths
        min_len <- min(length(baseline_finite), length(finite_estimates))
        baseline_finite <- baseline_finite[1:min_len]
        finite_estimates_matched <- finite_estimates[1:min_len]

        re_ci <- compute_boot_ci(
          list(baseline_finite, finite_estimates_matched, true_val),
          boot_re,
          conf = confidence_level,
          R = bootstrap_reps
        )
        metrics_df$re_se[i] <- re_ci$se
        metrics_df$re_ci_lower[i] <- re_ci$ci_lower
        metrics_df$re_ci_upper[i] <- re_ci$ci_upper
      }
    }
  }

  setTxtProgressBar(pb, i)
}

close(pb)

# =============================================================================
# Save Enhanced Results
# =============================================================================

cat("\n\nSaving enhanced results...\n")

output_path <- file.path(results_dir, "monte_carlo_metrics_with_ci.csv")
write.csv(metrics_df, output_path, row.names = FALSE)
cat(sprintf("  Saved: %s\n", output_path))

# Also update the summary tables
arima110_summary <- metrics_df[
  metrics_df$model == "ARIMA(1,1,0)" & metrics_df$parameter == "phi1",
]
write.csv(
  arima110_summary,
  file.path(results_dir, "arima110_summary_with_ci.csv"),
  row.names = FALSE
)
cat(sprintf("  Saved: %s\n", file.path(results_dir, "arima110_summary_with_ci.csv")))

# =============================================================================
# Summary Statistics
# =============================================================================

cat("\n=== Summary of Confidence Intervals ===\n")

# Count how many CIs were computed
n_bias_ci <- sum(!is.na(metrics_df$bias_ci_lower))
n_var_ci <- sum(!is.na(metrics_df$variance_ci_lower))
n_mse_ci <- sum(!is.na(metrics_df$mse_ci_lower))
n_re_ci <- sum(!is.na(metrics_df$re_ci_lower))

cat(sprintf("Bias CIs computed: %d / %d (%.1f%%)\n",
            n_bias_ci, nrow(metrics_df), 100 * n_bias_ci / nrow(metrics_df)))
cat(sprintf("Variance CIs computed: %d / %d (%.1f%%)\n",
            n_var_ci, nrow(metrics_df), 100 * n_var_ci / nrow(metrics_df)))
cat(sprintf("MSE CIs computed: %d / %d (%.1f%%)\n",
            n_mse_ci, nrow(metrics_df), 100 * n_mse_ci / nrow(metrics_df)))
cat(sprintf("RE CIs computed: %d / %d (%.1f%%)\n",
            n_re_ci, nrow(metrics_df), 100 * n_re_ci / nrow(metrics_df)))

# Example output
cat("\n=== Example: ARIMA(1,1,0) with Gamma distribution, N=100, PMM2 ===\n")
example <- metrics_df[
  metrics_df$model == "ARIMA(1,1,0)" &
  metrics_df$distribution == "Gamma" &
  metrics_df$sample_size == 100 &
  metrics_df$method == "PMM2" &
  metrics_df$parameter == "phi1",
]

if (nrow(example) > 0) {
  cat(sprintf("MSE: %.6f [%.6f, %.6f]\n",
              example$mse, example$mse_ci_lower, example$mse_ci_upper))
  cat(sprintf("RE:  %.3f [%.3f, %.3f]\n",
              example$relative_efficiency, example$re_ci_lower, example$re_ci_upper))
}

cat("\n✓ Bootstrap confidence intervals added successfully!\n")
cat("  Use monte_carlo_metrics_with_ci.csv for tables with CI\n\n")
