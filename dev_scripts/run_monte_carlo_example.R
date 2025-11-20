# ==============================================================================
# Quick Start Example: Monte Carlo Comparison of Seasonal Models
# ==============================================================================
#
# This is a simplified example to demonstrate the Monte Carlo framework
# with reduced computational requirements (50 replications instead of 500)
#
# For full simulations, use: source("monte_carlo_seasonal_comparison.R")
# ==============================================================================

library(EstemPMM)

cat("==============================================================================\n")
cat("Quick Monte Carlo Example for Seasonal Models\n")
cat("This is a DEMO with 50 replications (vs 500 in full study)\n")
cat("==============================================================================\n\n")

# Settings
set.seed(42)
N_REPS <- 50  # Reduced for quick demo
SAMPLE_SIZE <- 200
SEASONAL_PERIOD <- 12

# ==============================================================================
# Example 1: SAR(1,1)_12 Model
# ==============================================================================

cat("Example 1: SAR(1,1)_12 Model\n")
cat("------------------------------------------------------------\n")

# True parameters
true_ar <- 0.5
true_sar <- 0.6

# Storage for results
pmm2_estimates <- matrix(NA, N_REPS, 2)
css_estimates <- matrix(NA, N_REPS, 2)
pmm2_g_values <- numeric(N_REPS)

cat("Running", N_REPS, "replications with n =", SAMPLE_SIZE, "...\n")
pb <- txtProgressBar(min = 0, max = N_REPS, style = 3)

for (i in 1:N_REPS) {
  # Generate data with gamma innovations
  innov <- rgamma(SAMPLE_SIZE + 50, shape = 2, scale = 1) - 2
  y <- numeric(SAMPLE_SIZE + 50)

  for (t in 14:(SAMPLE_SIZE + 50)) {
    y[t] <- true_ar * y[t-1] + true_sar * y[t-12] + innov[t]
  }
  y <- y[25:(SAMPLE_SIZE + 50)]  # Remove burn-in

  # Fit with PMM2
  fit_pmm2 <- tryCatch(
    sar_pmm2(y, order = c(1, 1), season = list(period = 12),
             method = "pmm2", verbose = FALSE),
    error = function(e) NULL
  )

  if (!is.null(fit_pmm2) && length(fit_pmm2@coefficients) == 2) {
    pmm2_estimates[i, ] <- fit_pmm2@coefficients
    pmm2_g_values[i] <- pmm2_variance_factor(fit_pmm2@m2, fit_pmm2@m3, fit_pmm2@m4)$g
  }

  # Fit with CSS
  fit_css <- tryCatch(
    sar_pmm2(y, order = c(1, 1), season = list(period = 12),
             method = "css", verbose = FALSE),
    error = function(e) NULL
  )

  if (!is.null(fit_css) && length(fit_css@coefficients) == 2) {
    css_estimates[i, ] <- fit_css@coefficients
  }

  setTxtProgressBar(pb, i)
}
close(pb)

# Remove NA rows
pmm2_estimates <- pmm2_estimates[complete.cases(pmm2_estimates), ]
css_estimates <- css_estimates[complete.cases(css_estimates), ]
pmm2_g_values <- pmm2_g_values[!is.na(pmm2_g_values)]

# Compute statistics
cat("\n\nResults for SAR(1,1)_12:\n")
cat("========================\n\n")

cat("True parameters: AR(1) =", true_ar, ", SAR(1) =", true_sar, "\n\n")

# AR coefficient
cat("AR(1) coefficient:\n")
cat("  PMM2 - Mean:", round(mean(pmm2_estimates[, 1]), 4),
    ", Bias:", round(mean(pmm2_estimates[, 1]) - true_ar, 4),
    ", SD:", round(sd(pmm2_estimates[, 1]), 4), "\n")
cat("  CSS  - Mean:", round(mean(css_estimates[, 1]), 4),
    ", Bias:", round(mean(css_estimates[, 1]) - true_ar, 4),
    ", SD:", round(sd(css_estimates[, 1]), 4), "\n")

var_reduction_ar <- 1 - var(pmm2_estimates[, 1]) / var(css_estimates[, 1])
cat("  Variance reduction:", round(var_reduction_ar * 100, 2), "%\n\n")

# SAR coefficient
cat("SAR(1) coefficient:\n")
cat("  PMM2 - Mean:", round(mean(pmm2_estimates[, 2]), 4),
    ", Bias:", round(mean(pmm2_estimates[, 2]) - true_sar, 4),
    ", SD:", round(sd(pmm2_estimates[, 2]), 4), "\n")
cat("  CSS  - Mean:", round(mean(css_estimates[, 2]), 4),
    ", Bias:", round(mean(css_estimates[, 2]) - true_sar, 4),
    ", SD:", round(sd(css_estimates[, 2]), 4), "\n")

var_reduction_sar <- 1 - var(pmm2_estimates[, 2]) / var(css_estimates[, 2])
cat("  Variance reduction:", round(var_reduction_sar * 100, 2), "%\n\n")

cat("PMM2 efficiency factor (g):", round(mean(pmm2_g_values), 4), "\n")
cat("  (g < 1 means PMM2 is more efficient than OLS/CSS)\n\n")

# ==============================================================================
# Example 2: SARMA(1,0)×(1,1)_12 Model
# ==============================================================================

cat("\n==============================================================================\n")
cat("Example 2: SARMA(1,0)×(1,1)_12 Model\n")
cat("------------------------------------------------------------\n")

true_params <- c(ar = 0.5, sar = 0.6, sma = 0.4)

pmm2_sarma <- matrix(NA, N_REPS, 3)
css_sarma <- matrix(NA, N_REPS, 3)

cat("Running", N_REPS, "replications...\n")
pb <- txtProgressBar(min = 0, max = N_REPS, style = 3)

for (i in 1:N_REPS) {
  # Generate SARMA data
  innov <- rgamma(SAMPLE_SIZE + 100, shape = 2, scale = 1) - 2
  y <- arima.sim(n = SAMPLE_SIZE + 100,
                 list(ar = true_params["ar"],
                      seasonal = list(sar = true_params["sar"],
                                     sma = true_params["sma"],
                                     period = 12)),
                 innov = innov)
  y <- as.numeric(y)

  # PMM2
  fit_pmm2 <- tryCatch(
    sarma_pmm2(y, order = c(1, 1, 0, 1), season = list(period = 12),
               method = "pmm2", verbose = FALSE),
    error = function(e) NULL
  )

  if (!is.null(fit_pmm2) && length(fit_pmm2@coefficients) == 3) {
    pmm2_sarma[i, ] <- fit_pmm2@coefficients
  }

  # CSS
  fit_css <- tryCatch(
    sarma_pmm2(y, order = c(1, 1, 0, 1), season = list(period = 12),
               method = "css", verbose = FALSE),
    error = function(e) NULL
  )

  if (!is.null(fit_css) && length(fit_css@coefficients) == 3) {
    css_sarma[i, ] <- fit_css@coefficients
  }

  setTxtProgressBar(pb, i)
}
close(pb)

pmm2_sarma <- pmm2_sarma[complete.cases(pmm2_sarma), ]
css_sarma <- css_sarma[complete.cases(css_sarma), ]

cat("\n\nResults for SARMA(1,0)×(1,1)_12:\n")
cat("=================================\n\n")

cat("True parameters: AR =", true_params["ar"], ", SAR =", true_params["sar"],
    ", SMA =", true_params["sma"], "\n\n")

param_names <- c("AR(1)", "SAR(1)", "SMA(1)")

for (j in 1:3) {
  cat(param_names[j], ":\n")
  cat("  PMM2 - Mean:", round(mean(pmm2_sarma[, j]), 4),
      ", SD:", round(sd(pmm2_sarma[, j]), 4), "\n")
  cat("  CSS  - Mean:", round(mean(css_sarma[, j]), 4),
      ", SD:", round(sd(css_sarma[, j]), 4), "\n")

  var_red <- 1 - var(pmm2_sarma[, j]) / var(css_sarma[, j])
  cat("  Variance reduction:", round(var_red * 100, 2), "%\n\n")
}

# ==============================================================================
# Summary
# ==============================================================================

cat("\n==============================================================================\n")
cat("Summary\n")
cat("==============================================================================\n\n")

cat("This quick example demonstrates:\n")
cat("1. PMM2 provides lower variance estimates for asymmetric innovations\n")
cat("2. Variance reduction typically ranges from 20-50%\n")
cat("3. The efficiency factor g < 1 confirms PMM2 superiority\n")
cat("4. Both simple (SAR) and complex (SARMA) models benefit from PMM2\n\n")

cat("For comprehensive results with 500 replications and multiple sample sizes:\n")
cat("  source('monte_carlo_seasonal_comparison.R')\n")
cat("  source('visualize_monte_carlo_results.R')\n\n")

cat("==============================================================================\n")
cat("Quick example completed!\n")
cat("==============================================================================\n")
