#!/usr/bin/env Rscript

# Regenerate documentation and reinstall package
cat("Step 1: Regenerating documentation with roxygen2...\n")
roxygen2::roxygenize()

cat("\nStep 2: Installing package...\n")
devtools::install(upgrade = "never")

cat("\nStep 3: Testing sma_pmm2 function...\n")
library(EstemPMM)

# Check if function exists
if (exists("sma_pmm2")) {
  cat("SUCCESS: sma_pmm2 function is available!\n")

  # Test it
  set.seed(123)
  n <- 120
  innov <- rgamma(n, 2, 1) - 2
  y <- numeric(n)
  for (t in 1:n) {
    y[t] <- innov[t] + if(t > 12) 0.6 * innov[t-12] else 0
  }

  cat("\nTesting sma_pmm2()...\n")
  fit <- sma_pmm2(y, order = 1, season = list(period = 12), verbose = TRUE)
  cat("\nFit summary:\n")
  print(summary(fit))

} else {
  cat("ERROR: sma_pmm2 function still not found!\n")
}
