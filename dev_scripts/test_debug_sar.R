#!/usr/bin/env Rscript
# Debugging script for SAR-PMM2 Monte Carlo issue

library(EstemPMM)

cat("=== Test 1: Generate SAR process ===\n")
set.seed(123)

# Test parameters (same as in run_sar_tests.R scenario1)
n <- 120
ar_coef <- numeric(0)  # This is the issue - should be NULL
sar_coef <- 0.6
period <- 12

cat("Parameters:\n")
cat("  ar_coef:", ar_coef, "(length:", length(ar_coef), ")\n")
cat("  sar_coef:", sar_coef, "(length:", length(sar_coef), ")\n")
cat("  period:", period, "\n\n")

# Generate simple SAR process manually
y <- numeric(n)
y[1:period] <- rnorm(period, sd = 1)
for (t in (period + 1):n) {
  y[t] <- sar_coef * y[t - period] + rnorm(1, sd = 1)
}

cat("Generated series length:", length(y), "\n")
cat("Series range:", range(y), "\n\n")

cat("=== Test 2: Fit with OLS ===\n")
p <- length(ar_coef)
P <- length(sar_coef)
cat("p (AR order):", p, "\n")
cat("P (SAR order):", P, "\n\n")

# Try OLS manually
tryCatch({
  y_centered <- y - mean(y)
  X <- create_sar_matrix(y_centered, p = p, P = P, s = period)
  cat("Design matrix X dimensions:", dim(X), "\n")
  cat("Number of coefficients expected:", ncol(X), "\n")

  max_lag <- max(p, P * period)
  y_dep <- y_centered[(max_lag + 1):length(y_centered)]
  cat("y_dep length:", length(y_dep), "\n")
  cat("X rows:", nrow(X), "\n\n")

  beta_ols <- solve(t(X) %*% X, t(X) %*% y_dep)
  cat("OLS estimate:", as.numeric(beta_ols), "\n")
  cat("True value:", sar_coef, "\n")
  cat("Error:", abs(as.numeric(beta_ols) - sar_coef), "\n\n")

}, error = function(e) {
  cat("ERROR in OLS:", conditionMessage(e), "\n\n")
})

cat("=== Test 3: Fit with sar_pmm2() ===\n")
tryCatch({
  fit <- sar_pmm2(
    y,
    order = c(p, P),
    season = list(period = period),
    method = "pmm2",
    include.mean = TRUE,
    verbose = TRUE
  )

  cat("\nFit object class:", class(fit), "\n")
  cat("Converged:", fit@convergence, "\n")
  cat("Coefficients:", coef(fit), "\n")
  cat("Coefficient names:", names(coef(fit)), "\n\n")

  summary(fit)

}, error = function(e) {
  cat("ERROR in sar_pmm2():", conditionMessage(e), "\n")
  cat("Full error:\n")
  print(e)
})

cat("\n=== Test 4: Check true_coef construction ===\n")
true_coef <- c(ar_coef, sar_coef)
cat("true_coef:", true_coef, "\n")
cat("length(true_coef):", length(true_coef), "\n")

# Test the problematic indexing
cat("\nTesting parameter indexing:\n")
cat("ar_coef[1]:", ar_coef[1], "(is.na:", is.na(ar_coef[1]), ")\n")
cat("numeric(0)[1]:", numeric(0)[1], "\n")

# Proper way to get true value
for (i in 1:length(true_coef)) {
  if (i <= p) {
    true_val <- ar_coef[i]
  } else {
    true_val <- sar_coef[i - p]
  }
  cat(sprintf("Param %d: true_val = %.3f\n", i, true_val))
}

cat("\nDone!\n")
