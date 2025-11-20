#!/usr/bin/env Rscript
# Regenerate NAMESPACE and documentation using roxygen2

cat("Regenerating package documentation...\n\n")

# Check if roxygen2 is available
if (!requireNamespace("roxygen2", quietly = TRUE)) {
  cat("Installing roxygen2...\n")
  install.packages("roxygen2", repos = "https://cloud.r-project.org")
}

# Load roxygen2
library(roxygen2)

# Regenerate documentation
cat("Running roxygen2::roxygenize()...\n")
roxygen2::roxygenize(".", roclets = c("rd", "collate", "namespace"))

cat("\nDone! NAMESPACE and .Rd files have been regenerated.\n\n")

# Show SAR-related exports
cat("SAR-related exports in NAMESPACE:\n")
namespace_lines <- readLines("NAMESPACE")
sar_lines <- grep("sar|SAR", namespace_lines, value = TRUE, ignore.case = FALSE)
if (length(sar_lines) > 0) {
  cat(paste(sar_lines, collapse = "\n"), "\n")
} else {
  cat("(none found - check for errors)\n")
}
