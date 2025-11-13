#!/usr/bin/env Rscript

# Monte Carlo simulation to compare SMA-PMM2 vs CSS estimation

cat("\n")
cat("╔═══════════════════════════════════════════════════════════════════════╗\n")
cat("║     Monte Carlo Simulation: SMA-PMM2 vs CSS Performance Analysis     ║\n")
cat("╚═══════════════════════════════════════════════════════════════════════╝\n\n")

# Load package
cat("→ Loading EstemPMM package...\n")
suppressMessages({
  devtools::load_all(export_all = FALSE, helpers = FALSE, attach = TRUE)
})
cat("✓ Package loaded\n\n")

# Simulation parameters
cat("→ Simulation Configuration:\n")
cat("──────────────────────────────────────────────────────────────────────\n")

n_sim <- 500           # Number of Monte Carlo replications
n_obs <- 120           # Time series length
s <- 12                # Seasonal period
Q <- 1                 # SMA order
theta_true <- 0.6      # True coefficient
innovation_type <- "gamma"  # "gamma", "lognormal", or "normal"

cat(sprintf("  Monte Carlo replications:  %d\n", n_sim))
cat(sprintf("  Time series length:        %d\n", n_obs))
cat(sprintf("  Seasonal period:           %d\n", s))
cat(sprintf("  SMA order:                 %d\n", Q))
cat(sprintf("  True coefficient:          θ = %.2f\n", theta_true))
cat(sprintf("  Innovation distribution:   %s\n", innovation_type))
cat("──────────────────────────────────────────────────────────────────────\n\n")

# Storage for results
results_css <- numeric(n_sim)
results_pmm2 <- numeric(n_sim)
converged_pmm2 <- logical(n_sim)
iterations_pmm2 <- numeric(n_sim)

# Progress tracking
cat("→ Running Monte Carlo simulation...\n\n")
start_time <- Sys.time()

# Set up progress bar
pb_width <- 50
print_progress <- function(i, total) {
  pct <- i / total
  filled <- round(pb_width * pct)
  bar <- paste0(rep("█", filled), rep("░", pb_width - filled), collapse = "")
  cat(sprintf("\r  Progress: [%s] %3.0f%% (%d/%d)", bar, pct * 100, i, total))
  flush.console()
}

# Run simulation
for (i in 1:n_sim) {
  # Generate innovations based on distribution type
  if (innovation_type == "gamma") {
    # Gamma distribution (right-skewed)
    innov <- rgamma(n_obs, shape = 2, scale = 1) - 2
  } else if (innovation_type == "lognormal") {
    # Log-normal distribution (heavy right tail)
    innov <- rlnorm(n_obs, meanlog = 0, sdlog = 0.5) - exp(0.5^2 / 2)
  } else {
    # Normal distribution (for comparison)
    innov <- rnorm(n_obs, mean = 0, sd = 1)
  }

  # Generate SMA(1)_s process
  y <- numeric(n_obs)
  for (t in 1:n_obs) {
    ma_term <- if (t > s) theta_true * innov[t - s] else 0
    y[t] <- innov[t] + ma_term
  }

  # Fit with CSS
  tryCatch({
    fit_css <- sma_pmm2(y, order = Q, season = list(period = s),
                        method = "css", verbose = FALSE)
    results_css[i] <- coef(fit_css)
  }, error = function(e) {
    results_css[i] <- NA
  })

  # Fit with PMM2
  tryCatch({
    fit_pmm2 <- sma_pmm2(y, order = Q, season = list(period = s),
                         method = "pmm2", verbose = FALSE)
    results_pmm2[i] <- coef(fit_pmm2)
    converged_pmm2[i] <- fit_pmm2@convergence
    iterations_pmm2[i] <- fit_pmm2@iterations
  }, error = function(e) {
    results_pmm2[i] <- NA
    converged_pmm2[i] <- FALSE
    iterations_pmm2[i] <- NA
  })

  # Update progress
  if (i %% 10 == 0 || i == n_sim) {
    print_progress(i, n_sim)
  }
}

cat("\n\n")
elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
cat(sprintf("✓ Simulation completed in %.1f seconds\n\n", elapsed))

# Remove NA values
valid_css <- !is.na(results_css)
valid_pmm2 <- !is.na(results_pmm2)

results_css_clean <- results_css[valid_css]
results_pmm2_clean <- results_pmm2[valid_pmm2]

cat("→ Data Quality Check:\n")
cat(sprintf("  CSS:  %d/%d valid estimates (%.1f%%)\n",
            sum(valid_css), n_sim, 100 * mean(valid_css)))
cat(sprintf("  PMM2: %d/%d valid estimates (%.1f%%)\n",
            sum(valid_pmm2), n_sim, 100 * mean(valid_pmm2)))
cat(sprintf("  PMM2 convergence rate: %.1f%%\n\n",
            100 * mean(converged_pmm2[valid_pmm2])))

# Calculate statistics
cat("\n")
cat("╔═══════════════════════════════════════════════════════════════════════╗\n")
cat("║                         SIMULATION RESULTS                            ║\n")
cat("╚═══════════════════════════════════════════════════════════════════════╝\n\n")

cat("True Parameter:\n")
cat("───────────────────────────────────────────────────────────────────────\n")
cat(sprintf("  θ = %.4f\n\n", theta_true))

cat("CSS Estimation:\n")
cat("───────────────────────────────────────────────────────────────────────\n")
css_mean <- mean(results_css_clean)
css_sd <- sd(results_css_clean)
css_bias <- css_mean - theta_true
css_mse <- mean((results_css_clean - theta_true)^2)
css_rmse <- sqrt(css_mse)

cat(sprintf("  Mean:             %.6f\n", css_mean))
cat(sprintf("  Std. Deviation:   %.6f\n", css_sd))
cat(sprintf("  Bias:             %.6f\n", css_bias))
cat(sprintf("  MSE:              %.6f\n", css_mse))
cat(sprintf("  RMSE:             %.6f\n\n", css_rmse))

cat("PMM2 Estimation:\n")
cat("───────────────────────────────────────────────────────────────────────\n")
pmm2_mean <- mean(results_pmm2_clean)
pmm2_sd <- sd(results_pmm2_clean)
pmm2_bias <- pmm2_mean - theta_true
pmm2_mse <- mean((results_pmm2_clean - theta_true)^2)
pmm2_rmse <- sqrt(pmm2_mse)

cat(sprintf("  Mean:             %.6f\n", pmm2_mean))
cat(sprintf("  Std. Deviation:   %.6f\n", pmm2_sd))
cat(sprintf("  Bias:             %.6f\n", pmm2_bias))
cat(sprintf("  MSE:              %.6f\n", pmm2_mse))
cat(sprintf("  RMSE:             %.6f\n\n", pmm2_rmse))

cat("Performance Comparison:\n")
cat("───────────────────────────────────────────────────────────────────────\n")

# Variance ratio (key metric for PMM2)
var_ratio <- pmm2_sd^2 / css_sd^2
var_reduction <- (1 - var_ratio) * 100

cat(sprintf("  Variance Ratio (PMM2/CSS):      %.4f\n", var_ratio))
if (var_ratio < 1) {
  cat(sprintf("  Variance Reduction:             %.2f%% ✓\n", var_reduction))
} else {
  cat(sprintf("  Variance Increase:              %.2f%%\n", -var_reduction))
}

# MSE ratio
mse_ratio <- pmm2_mse / css_mse
cat(sprintf("  MSE Ratio (PMM2/CSS):           %.4f\n", mse_ratio))
if (mse_ratio < 1) {
  mse_improvement <- (1 - mse_ratio) * 100
  cat(sprintf("  MSE Improvement:                %.2f%% ✓\n\n", mse_improvement))
} else {
  mse_degradation <- (mse_ratio - 1) * 100
  cat(sprintf("  MSE Degradation:                %.2f%%\n\n", mse_degradation))
}

# Algorithm statistics
cat("PMM2 Algorithm Statistics:\n")
cat("───────────────────────────────────────────────────────────────────────\n")
cat(sprintf("  Mean iterations:              %.2f\n", mean(iterations_pmm2, na.rm = TRUE)))
cat(sprintf("  Median iterations:            %.0f\n", median(iterations_pmm2, na.rm = TRUE)))
cat(sprintf("  Max iterations:               %.0f\n\n", max(iterations_pmm2, na.rm = TRUE)))

# Theoretical prediction
cat("Theoretical Prediction:\n")
cat("───────────────────────────────────────────────────────────────────────\n")

# Calculate average moments from one fit
set.seed(123)
innov_test <- rgamma(n_obs, shape = 2, scale = 1) - 2
y_test <- numeric(n_obs)
for (t in 1:n_obs) {
  ma_term <- if (t > s) theta_true * innov_test[t - s] else 0
  y_test[t] <- innov_test[t] + ma_term
}
fit_test <- sma_pmm2(y_test, order = Q, season = list(period = s),
                     method = "pmm2", verbose = FALSE)

c3 <- fit_test@m3 / (fit_test@m2^(3/2))
c4 <- fit_test@m4 / (fit_test@m2^2) - 3
g_theory <- 1 - c3^2 / (2 + c4)

cat(sprintf("  Skewness coefficient (c₃):    %.4f\n", c3))
cat(sprintf("  Excess kurtosis (c₄):         %.4f\n", c4))
cat(sprintf("  Theoretical g:                %.4f\n", g_theory))
cat(sprintf("  Expected variance ratio:      %.4f\n", g_theory))
cat(sprintf("  Observed variance ratio:      %.4f\n", var_ratio))
cat(sprintf("  Match:                        %.1f%%\n\n",
            100 * (1 - abs(var_ratio - g_theory) / g_theory)))

# Interpretation
cat("\n")
cat("╔═══════════════════════════════════════════════════════════════════════╗\n")
cat("║                           INTERPRETATION                              ║\n")
cat("╚═══════════════════════════════════════════════════════════════════════╝\n\n")

if (var_ratio < 1 && var_ratio < 0.95) {
  cat("✓ SUCCESS: PMM2 achieves significant variance reduction!\n\n")
  cat("  The PMM2 estimator demonstrates lower variance than CSS, which is the\n")
  cat("  expected behavior for asymmetric innovation distributions. This confirms\n")
  cat("  that PMM2 provides more stable parameter estimates.\n\n")
} else if (var_ratio < 1) {
  cat("✓ MODERATE SUCCESS: PMM2 achieves modest variance reduction.\n\n")
  cat("  The variance reduction is present but modest. This may be due to the\n")
  cat("  specific data characteristics or sample size.\n\n")
} else {
  cat("⚠ UNEXPECTED: PMM2 does not reduce variance in this simulation.\n\n")
  cat("  This could indicate:\n")
  cat("  - Innovations are not sufficiently asymmetric\n")
  cat("  - Sample size is too small\n")
  cat("  - Implementation issues need investigation\n\n")
}

# Save results
cat("→ Saving results...\n")
results_df <- data.frame(
  replication = 1:n_sim,
  css = results_css,
  pmm2 = results_pmm2,
  converged = converged_pmm2,
  iterations = iterations_pmm2
)

if (!dir.exists("test_results")) {
  dir.create("test_results")
}

output_file <- sprintf("test_results/sma_monte_carlo_%s_%dreps.csv",
                       format(Sys.Date(), "%Y%m%d"), n_sim)
write.csv(results_df, output_file, row.names = FALSE)
cat(sprintf("✓ Results saved to: %s\n\n", output_file))

cat("╔═══════════════════════════════════════════════════════════════════════╗\n")
cat("║                      SIMULATION COMPLETED                             ║\n")
cat("╚═══════════════════════════════════════════════════════════════════════╝\n\n")
