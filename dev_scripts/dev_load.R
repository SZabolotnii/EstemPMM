#!/usr/bin/env Rscript

# DEVELOPMENT MODE - Load package without installing
# This bypasses lazy-load database issues

cat("\n")
cat("═══════════════════════════════════════════════════════════\n")
cat("  DEVELOPMENT MODE - Loading EstemPMM without installation\n")
cat("═══════════════════════════════════════════════════════════\n\n")

# Unload if already loaded
if ("package:EstemPMM" %in% search()) {
  detach("package:EstemPMM", unload = TRUE, force = TRUE)
  cat("✓ Unloaded existing EstemPMM\n")
}

# Load using devtools::load_all() - this skips installation
cat("→ Loading package with devtools::load_all()...\n\n")

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Package 'devtools' required. Install: install.packages('devtools')")
}

devtools::load_all(export_all = FALSE, helpers = FALSE, attach = TRUE)

cat("\n✓ Package loaded in development mode\n\n")

# Check functions
cat("→ Checking SMA functions...\n")
if (exists("sma_pmm2")) {
  cat("✓ sma_pmm2 is available!\n\n")

  # Run test
  cat("→ Running test...\n")
  set.seed(123)
  n <- 120
  innov <- rgamma(n, 2, 1) - 2
  y <- numeric(n)
  for (t in 1:n) {
    y[t] <- innov[t] + if(t > 12) 0.6 * innov[t-12] else 0
  }

  fit <- sma_pmm2(y, order = 1, season = list(period = 12), verbose = TRUE)

  cat("\n→ Results:\n")
  cat("  Coefficient:", coef(fit), "\n")
  cat("  True value: 0.6\n")
  cat("  Converged:", fit@convergence, "\n")
  cat("  Iterations:", fit@iterations, "\n\n")

  cat("═══════════════════════════════════════════════════════════\n")
  cat("  ✓ SUCCESS! You can now use all EstemPMM functions\n")
  cat("═══════════════════════════════════════════════════════════\n\n")

  cat("NOTE: This loads the package in development mode.\n")
  cat("      Changes to R files will require re-running this script.\n")
  cat("      For production use, run force_rebuild.R instead.\n\n")

} else {
  cat("✗ ERROR: sma_pmm2 function not found!\n\n")
  cat("Available objects:\n")
  print(ls(envir = asNamespace("EstemPMM"), pattern = "sma"))
}
