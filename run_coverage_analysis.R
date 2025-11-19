#!/usr/bin/env Rscript
# Coverage Analysis Script for EstemPMM Package
# This script analyzes test coverage for the EstemPMM package using the covr package

# Install covr if not already installed
if (!requireNamespace("covr", quietly = TRUE)) {
  message("Installing covr package...")
  install.packages("covr")
}

# Load the covr package
library(covr)

# Set working directory to package root
# setwd(dirname(sys.frame(1)$ofile))

cat("========================================\n")
cat("EstemPMM Package Coverage Analysis\n")
cat("========================================\n\n")

# Run package coverage analysis
cat("Running coverage analysis...\n")
coverage <- package_coverage()

# Print coverage summary
cat("\n--- Coverage Summary ---\n")
print(coverage)

# Get detailed coverage by file
cat("\n--- Coverage by File ---\n")
file_coverage <- coverage_to_cobertura(coverage)

# Show percentage covered
cat("\n--- Percentage Covered ---\n")
percent <- percent_coverage(coverage)
cat(sprintf("Overall Coverage: %.2f%%\n", percent))

# Generate HTML report
cat("\n--- Generating HTML Report ---\n")
report_file <- "coverage_report.html"
report(coverage, file = report_file)
cat(sprintf("HTML report saved to: %s\n", report_file))

# Identify areas with low coverage
cat("\n--- Functions with Low Coverage (<80%) ---\n")
zero_cov <- zero_coverage(coverage)
if (length(zero_cov) > 0) {
  print(zero_cov)
} else {
  cat("None! All functions have >80% coverage.\n")
}

# Coverage by file type
cat("\n--- Coverage by Source File ---\n")
file_report <- coverage_to_cobertura(coverage)

# Specific focus on seasonal models
cat("\n--- Seasonal Models Coverage ---\n")
seasonal_files <- c(
  "R/pmm2_ts_main.R",      # Contains sar_pmm2, sma_pmm2, sarma_pmm2, sarima_pmm2
  "R/pmm2_classes.R",       # Contains S4 class definitions
  "R/pmm2_ts_methods.R"     # Contains methods for seasonal classes
)

for (file in seasonal_files) {
  if (file.exists(file)) {
    file_cov <- file_coverage(coverage, file)
    cat(sprintf("\n%s: %.2f%% covered\n", file, attr(file_cov, "coverage")))
  }
}

# Recommendations
cat("\n========================================\n")
cat("Recommendations:\n")
cat("========================================\n")

if (percent < 80) {
  cat("⚠️  Coverage is below 80%. Consider adding more tests.\n")
  cat("   Priority areas:\n")
  cat("   - Edge cases and error handling\n")
  cat("   - S4 method dispatching\n")
  cat("   - Seasonal model convergence scenarios\n")
} else if (percent < 90) {
  cat("✓ Good coverage (>80%), but room for improvement.\n")
  cat("  Consider testing:\n")
  cat("  - Complex seasonal patterns\n")
  cat("  - Boundary conditions\n")
  cat("  - Integration tests\n")
} else {
  cat("✓✓ Excellent coverage (>90%)!\n")
  cat("   Maintain this level with ongoing test additions.\n")
}

cat("\n========================================\n")
cat("Analysis Complete!\n")
cat("========================================\n")

# Save coverage data for CI/CD
saveRDS(coverage, file = "coverage_data.rds")
cat("\nCoverage data saved to: coverage_data.rds\n")

# Return coverage percentage for scripting
invisible(percent)
