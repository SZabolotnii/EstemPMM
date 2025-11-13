#!/usr/bin/env Rscript

# Complete clean reinstall of EstemPMM package

cat("=" , rep("=", 70), "\n", sep = "")
cat("CLEAN REINSTALL OF EstemPMM PACKAGE\n")
cat("=" , rep("=", 70), "\n\n", sep = "")

# Step 1: Remove installed package
cat("Step 1: Removing installed package...\n")
tryCatch({
  remove.packages("EstemPMM")
  cat("  ✓ Package removed successfully\n\n")
}, error = function(e) {
  cat("  ℹ Package not installed or already removed\n\n")
})

# Step 2: Clean R session
cat("Step 2: Cleaning R session...\n")
rm(list = ls(all.names = TRUE))
gc()
cat("  ✓ Session cleaned\n\n")

# Step 3: Regenerate documentation
cat("Step 3: Regenerating roxygen2 documentation...\n")
tryCatch({
  roxygen2::roxygenize(clean = TRUE)
  cat("  ✓ Documentation regenerated\n\n")
}, error = function(e) {
  cat("  ✗ Error:", conditionMessage(e), "\n\n")
  stop("Failed to regenerate documentation")
})

# Step 4: Install package
cat("Step 4: Installing package from source...\n")
tryCatch({
  devtools::install(upgrade = "never", build_vignettes = FALSE)
  cat("  ✓ Package installed successfully\n\n")
}, error = function(e) {
  cat("  ✗ Error:", conditionMessage(e), "\n\n")
  stop("Failed to install package")
})

# Step 5: Load and test
cat("Step 5: Loading and testing package...\n")
library(EstemPMM)

# Check exported functions
sma_functions <- ls("package:EstemPMM", pattern = "sma")
cat("  SMA-related functions found:", length(sma_functions), "\n")
cat("  Functions:", paste(sma_functions, collapse = ", "), "\n\n")

if ("sma_pmm2" %in% sma_functions) {
  cat("  ✓ sma_pmm2 function is available!\n\n")

  # Run basic test
  cat("Step 6: Running basic test...\n")
  set.seed(123)
  n <- 120
  innov <- rgamma(n, 2, 1) - 2
  y <- numeric(n)
  for (t in 1:n) {
    y[t] <- innov[t] + if(t > 12) 0.6 * innov[t-12] else 0
  }

  fit <- sma_pmm2(y, order = 1, season = list(period = 12), verbose = FALSE)
  cat("  ✓ Test completed successfully\n")
  cat("  Estimated coefficient:", coef(fit), "\n")
  cat("  True coefficient: 0.6\n\n")

  cat("=" , rep("=", 70), "\n", sep = "")
  cat("SUCCESS! Package is ready to use.\n")
  cat("=" , rep("=", 70), "\n", sep = "")

} else {
  cat("  ✗ ERROR: sma_pmm2 function not found!\n")
  cat("  Available functions:\n")
  print(ls("package:EstemPMM"))
}
