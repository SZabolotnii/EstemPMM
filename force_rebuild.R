#!/usr/bin/env Rscript

# AGGRESSIVE PACKAGE REBUILD - Forces complete clean reinstall

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║  FORCE REBUILD EstemPMM - Aggressive Clean Reinstall          ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# Force quit any loaded EstemPMM namespace
cat("→ Step 1: Forcing unload of EstemPMM namespace...\n")
if ("package:EstemPMM" %in% search()) {
  detach("package:EstemPMM", unload = TRUE, force = TRUE)
  cat("  ✓ Package detached\n")
}
try(unloadNamespace("EstemPMM"), silent = TRUE)
cat("  ✓ Namespace unloaded\n\n")

# Remove installed package
cat("→ Step 2: Removing installed package...\n")
lib_path <- .libPaths()[1]
pkg_path <- file.path(lib_path, "EstemPMM")

if (dir.exists(pkg_path)) {
  cat("  Found package at:", pkg_path, "\n")

  # Try normal removal
  try(remove.packages("EstemPMM"), silent = TRUE)

  # If still exists, force delete
  if (dir.exists(pkg_path)) {
    cat("  Normal removal failed, forcing deletion...\n")
    unlink(pkg_path, recursive = TRUE, force = TRUE)
  }

  if (!dir.exists(pkg_path)) {
    cat("  ✓ Package removed successfully\n\n")
  } else {
    cat("  ✗ Could not remove package directory!\n")
    cat("  Please manually delete:", pkg_path, "\n")
    cat("  Then restart R and run this script again.\n\n")
    stop("Manual intervention required")
  }
} else {
  cat("  ℹ Package not installed\n\n")
}

# Clean all temporary build artifacts
cat("→ Step 3: Cleaning build artifacts...\n")
build_files <- c(
  "src/*.o", "src/*.so", "src/*.dll",
  "src-i386", "src-x64",
  ".Rcheck"
)
for (pattern in build_files) {
  files <- Sys.glob(pattern)
  if (length(files) > 0) {
    unlink(files, recursive = TRUE, force = TRUE)
    cat("  Removed:", pattern, "\n")
  }
}
cat("  ✓ Build artifacts cleaned\n\n")

# Force garbage collection
cat("→ Step 4: Clearing R session memory...\n")
rm(list = ls(all.names = TRUE, envir = .GlobalEnv), envir = .GlobalEnv)
gc(verbose = FALSE, full = TRUE)
cat("  ✓ Memory cleared\n\n")

# Regenerate documentation
cat("→ Step 5: Regenerating documentation (roxygen2)...\n")
if (!requireNamespace("roxygen2", quietly = TRUE)) {
  stop("Package 'roxygen2' is required. Install with: install.packages('roxygen2')")
}
roxygen2::roxygenize(clean = TRUE)
cat("  ✓ Documentation regenerated\n\n")

# Build package tarball
cat("→ Step 6: Building package tarball...\n")
if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Package 'devtools' is required. Install with: install.packages('devtools')")
}

tarball <- devtools::build(binary = FALSE, vignettes = FALSE, manual = FALSE)
cat("  ✓ Package built:", basename(tarball), "\n\n")

# Install from tarball
cat("→ Step 7: Installing from tarball...\n")
install.packages(tarball, repos = NULL, type = "source", INSTALL_opts = "--no-multiarch")
cat("  ✓ Package installed\n\n")

# Verify installation
cat("→ Step 8: Verifying installation...\n")
if (!requireNamespace("EstemPMM", quietly = TRUE)) {
  stop("Installation verification failed!")
}

library(EstemPMM)
cat("  ✓ Package loaded successfully\n\n")

# Check for sma_pmm2
sma_funcs <- ls("package:EstemPMM", pattern = "sma")
cat("→ Step 9: Checking SMA functions...\n")
cat("  Found", length(sma_funcs), "SMA-related functions:\n")
cat("  -", paste(sma_funcs, collapse = "\n  - "), "\n\n")

if ("sma_pmm2" %in% sma_funcs) {
  cat("  ✓ sma_pmm2 is available!\n\n")

  # Run quick test
  cat("→ Step 10: Running quick test...\n")
  set.seed(123)
  n <- 120
  innov <- rgamma(n, 2, 1) - 2
  y <- numeric(n)
  for (t in 1:n) {
    y[t] <- innov[t] + if(t > 12) 0.6 * innov[t-12] else 0
  }

  fit <- sma_pmm2(y, order = 1, season = list(period = 12), verbose = FALSE)
  theta_hat <- coef(fit)

  cat("  ✓ Test completed\n")
  cat("  Estimated θ:", round(theta_hat, 4), "\n")
  cat("  True θ:     0.6000\n")
  cat("  Error:      ", round(abs(theta_hat - 0.6), 4), "\n\n")

  cat("╔════════════════════════════════════════════════════════════════╗\n")
  cat("║  ✓ SUCCESS! Package is ready to use.                          ║\n")
  cat("╚════════════════════════════════════════════════════════════════╝\n\n")

} else {
  cat("  ✗ ERROR: sma_pmm2 not found!\n")
  cat("\n  Available functions in EstemPMM:\n")
  all_funcs <- ls("package:EstemPMM")
  cat("  ", paste(head(all_funcs, 20), collapse = "\n   "), "\n")
  if (length(all_funcs) > 20) {
    cat("   ... and", length(all_funcs) - 20, "more\n")
  }
  stop("\nPackage loaded but sma_pmm2 function not exported!")
}
