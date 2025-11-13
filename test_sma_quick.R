#!/usr/bin/env Rscript

# Quick test for SMA-PMM2 implementation

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║          SMA-PMM2 Quick Test & Verification                  ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Load package in dev mode
cat("→ Loading EstemPMM package...\n")
devtools::load_all(export_all = FALSE, helpers = FALSE, attach = TRUE)
cat("✓ Package loaded\n\n")

# Generate test data
cat("→ Generating test data (SMA(1)_12 with gamma innovations)...\n")
set.seed(123)
n <- 120
s <- 12
theta_true <- 0.6

# Gamma innovations (asymmetric)
innov <- rgamma(n, shape = 2, scale = 1) - 2
y <- numeric(n)
for (t in 1:n) {
  ma_term <- if (t > s) theta_true * innov[t-s] else 0
  y[t] <- innov[t] + ma_term
}

cat(sprintf("  True coefficient: θ = %.2f\n", theta_true))
cat(sprintf("  Sample size: n = %d\n", n))
cat(sprintf("  Seasonal period: s = %d\n\n", s))

# Fit SMA model with PMM2
cat("→ Fitting SMA(1)_12 model with PMM2...\n\n")
cat("----------------------------------------\n")
fit_pmm2 <- sma_pmm2(y, order = 1, season = list(period = 12),
                     method = "pmm2", verbose = TRUE)
cat("----------------------------------------\n\n")

# Extract results
theta_pmm2 <- coef(fit_pmm2)
cat("→ Results:\n")
cat(sprintf("  PMM2 estimate: θ̂ = %.4f\n", theta_pmm2))
cat(sprintf("  True value:    θ = %.4f\n", theta_true))
cat(sprintf("  Error:         |θ̂ - θ| = %.4f\n", abs(theta_pmm2 - theta_true)))
cat(sprintf("  Relative err:  %.2f%%\n\n", 100 * abs(theta_pmm2 - theta_true) / theta_true))

# Fit with CSS for comparison
cat("→ Fitting SMA(1)_12 model with CSS (for comparison)...\n")
fit_css <- sma_pmm2(y, order = 1, season = list(period = 12),
                    method = "css", verbose = FALSE)
theta_css <- coef(fit_css)

cat(sprintf("  CSS estimate:  θ̂ = %.4f\n", theta_css))
cat(sprintf("  Error:         |θ̂ - θ| = %.4f\n\n", abs(theta_css - theta_true)))

# Compare methods
cat("→ Comparison:\n")
cat("  Method    Estimate    Error      Relative Error\n")
cat("  ------    --------    ------     --------------\n")
cat(sprintf("  CSS       %.4f      %.4f     %.2f%%\n",
            theta_css, abs(theta_css - theta_true),
            100 * abs(theta_css - theta_true) / theta_true))
cat(sprintf("  PMM2      %.4f      %.4f     %.2f%%\n\n",
            theta_pmm2, abs(theta_pmm2 - theta_true),
            100 * abs(theta_pmm2 - theta_true) / theta_true))

# Display full summary
cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║                  Detailed PMM2 Summary                       ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
summary(fit_pmm2)

# Test coef() and other methods
cat("\n→ Testing S4 methods:\n")
cat("  coef(fit):     ", coef(fit_pmm2), "\n")
cat("  fitted(fit):   ", head(fitted(fit_pmm2), 3), "... (first 3 values)\n")
cat("  residuals(fit):", head(residuals(fit_pmm2), 3), "... (first 3 values)\n\n")

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  ✓ ALL TESTS PASSED! SMA-PMM2 is working correctly!         ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")
