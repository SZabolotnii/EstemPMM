# Master Script: Run Complete ARIMA Comparison Study
# This script orchestrates the entire research pipeline


# Set seed for reproducibility
set.seed(12345)

get_script_dir <- function() {
  # Try to resolve the path when sourced or executed via Rscript
  frame_files <- vapply(sys.frames(), function(f) {
    if (!is.null(f$ofile)) f$ofile else NA_character_
  }, character(1), USE.NAMES = FALSE)
  frame_files <- frame_files[!is.na(frame_files)]
  if (length(frame_files) > 0) {
    return(dirname(normalizePath(frame_files[1])))
  }

  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- cmd_args[grepl("^--file=", cmd_args)]
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }

  normalizePath(getwd())
}

scripts_dir <- get_script_dir()
project_root <- normalizePath(file.path(scripts_dir, ".."))
original_wd <- getwd()
on.exit(setwd(original_wd), add = TRUE)
setwd(project_root)

cat("================================================================\n")
cat("COMPLETE ARIMA COMPARISON STUDY\n")
cat("PMM2 vs Classical Methods on WTI Oil Prices\n")
cat("================================================================\n\n")

# Check required packages
required_packages <- c("EstemPMM", "ggplot2", "tseries",
                       "gridExtra", "RColorBrewer")

missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
  cat("Missing required packages:\n")
  cat(paste0("  - ", missing_packages, collapse = "\n"), "\n\n")
  cat("Install them with:\n")
  cat(sprintf("  install.packages(c('%s'))\n\n", paste(missing_packages, collapse = "', '")))
  stop("Please install missing packages and try again.")
}

# Load libraries
library(EstemPMM)

cat("All required packages are available.\n\n")

# ==================== STEP 1: DATA ANALYSIS AND MODEL FITTING ====================
cat("STEP 1: Running comprehensive model comparison...\n")
cat("----------------------------------------------------------\n")
cat("This will fit 12 models (6 specifications × 2 methods)\n")
cat("Estimated time: 2-5 minutes\n\n")

start_time <- Sys.time()

tryCatch({
  source(file.path("scripts", "comprehensive_study.R"))
  cat("\n✓ Model fitting completed successfully\n\n")
}, error = function(e) {
  cat(sprintf("\n✗ Error in model fitting: %s\n\n", e$message))
  stop("Study aborted")
})

step1_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

# ==================== STEP 2: CREATE VISUALIZATIONS ====================
cat("STEP 2: Creating visualizations...\n")
cat("----------------------------------------------------------\n")
cat("Generating 10 publication-quality plots\n\n")

step2_start <- Sys.time()

tryCatch({
  source(file.path("scripts", "create_visualizations.R"))
  cat("\n✓ Visualization creation completed successfully\n\n")
}, error = function(e) {
  cat(sprintf("\n✗ Error in visualization: %s\n\n", e$message))
  cat("Continuing to report generation...\n\n")
})

step2_time <- as.numeric(difftime(Sys.time(), step2_start, units = "secs"))

# ==================== STEP 3: GENERATE ANALYTICAL REPORT ====================
cat("STEP 3: Generating analytical report...\n")
cat("----------------------------------------------------------\n")
cat("Creating comprehensive markdown report\n\n")

step3_start <- Sys.time()

tryCatch({
  source(file.path("scripts", "generate_report.R"))
  cat("\n✓ Report generation completed successfully\n\n")
}, error = function(e) {
  cat(sprintf("\n✗ Error in report generation: %s\n\n", e$message))
  stop("Report generation failed")
})

step3_time <- as.numeric(difftime(Sys.time(), step3_start, units = "secs"))

# ==================== COMPLETION SUMMARY ====================
total_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

cat("================================================================\n")
cat("STUDY COMPLETED SUCCESSFULLY\n")
cat("================================================================\n\n")

cat("Timing Summary:\n")
cat(sprintf("  Step 1 (Model Fitting):     %.1f seconds\n", step1_time))
cat(sprintf("  Step 2 (Visualizations):    %.1f seconds\n", step2_time))
cat(sprintf("  Step 3 (Report Generation): %.1f seconds\n", step3_time))
cat(sprintf("  Total Time:                 %.1f seconds (%.1f minutes)\n\n",
            total_time, total_time / 60))

cat("Output Files:\n")
cat("  ├── results/\n")
cat("  │   ├── full_results.csv           - Complete results table\n")
cat("  │   ├── method_comparison.csv      - Head-to-head comparison\n")
cat("  │   ├── descriptive_stats.csv      - Data statistics\n")
cat("  │   ├── parameter_mse_ratio.csv    - Parameter MSE ratios (PMM2 vs CSS)\n")
cat("  │   ├── all_results.rds            - R data objects\n")
cat("  │   ├── fitted_models.rds          - Fitted model objects\n")
cat("  │   ├── study_summary.rds          - Summary information\n")
cat("  │   ├── ANALYTICAL_REPORT.md       - Comprehensive report\n")
cat("  │   └── plots/\n")
cat("  │       ├── 01_aic_comparison.png\n")
cat("  │       ├── 02_bic_comparison.png\n")
cat("  │       ├── 03_rmse_comparison.png\n")
cat("  │       ├── 04_computation_time.png\n")
cat("  │       ├── 05_kurtosis_comparison.png\n")
cat("  │       ├── 06_skewness_comparison.png\n")
cat("  │       ├── 07_performance_heatmap.png\n")
cat("  │       ├── 08_method_differences.png\n")
cat("  │       ├── 09_best_model_diagnostics.png\n")
cat("  │       └── 10_summary_statistics.png\n\n")

cat("Quick Actions:\n")
cat("  # View the analytical report\n")
cat("  cat results/ANALYTICAL_REPORT.md\n\n")

cat("  # View results in R\n")
cat("  results <- readRDS('results/all_results.rds')\n")
cat("  View(results)\n\n")

cat("  # Convert report to HTML\n")
cat("  pandoc results/ANALYTICAL_REPORT.md \\\n")
cat("    -o results/ANALYTICAL_REPORT.html \\\n")
cat("    -s --toc --toc-depth=3 \\\n")
cat("    --metadata title='ARIMA Comparison Study'\n\n")

cat("  # Convert report to PDF (requires LaTeX)\n")
cat("  pandoc results/ANALYTICAL_REPORT.md \\\n")
cat("    -o results/ANALYTICAL_REPORT.pdf \\\n")
cat("    -s --toc --toc-depth=3\n\n")

# Load and display key findings
summary_info <- readRDS("results/study_summary.rds")

cat("================================================================\n")
cat("KEY FINDINGS\n")
cat("================================================================\n\n")
cat(sprintf("Best Model (AIC): %s\n", summary_info$best_model_aic))
cat(sprintf("Best Model (BIC): %s\n\n", summary_info$best_model_bic))
cat(sprintf("PMM2 Win Rate (AIC): %.1f%%\n", summary_info$pmm2_win_rate_aic * 100))
cat(sprintf("PMM2 Win Rate (BIC): %.1f%%\n\n", summary_info$pmm2_win_rate_bic * 100))

cat("For detailed analysis, see: results/ANALYTICAL_REPORT.md\n\n")
