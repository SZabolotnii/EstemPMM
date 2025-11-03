# Quick demo: ARIMA comparison on WTI Oil Price Data
# This is a simplified version that runs quickly without generating full report

library(EstemPMM)

# Set seed for reproducibility
set.seed(12345)

get_script_dir <- function() {
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

script_dir <- get_script_dir()
project_root <- normalizePath(file.path(script_dir, ".."))

cat("=======================================================\n")
cat("ARIMA Model Comparison: PMM2 vs Classical Methods\n")
cat("Data: WTI Crude Oil Prices (DCOILWTICO)\n")
cat("=======================================================\n\n")

# Load data
candidate_paths <- c(
  system.file("extdata", "DCOILWTICO.csv", package = "EstemPMM"),
  file.path(project_root, "data", "DCOILWTICO.csv"),
  file.path(script_dir, "DCOILWTICO.csv")
)

data_path <- ""
for (candidate in candidate_paths) {
  if (nzchar(candidate) && file.exists(candidate)) {
    data_path <- candidate
    break
  }
}

if (!nzchar(data_path) || !file.exists(data_path)) {
  stop("Data file DCOILWTICO.csv not found. Please ensure it exists in the repository data directory.")
}

# Read and clean data
cat("Loading data...\n")
oil_data <- read.csv(data_path, stringsAsFactors = FALSE)
prices <- as.numeric(oil_data$DCOILWTICO)
prices <- prices[!is.na(prices)]

cat(sprintf("Total observations: %d\n\n", length(prices)))

# Test different ARIMA specifications
models_to_test <- list(
  "ARIMA(1,1,1)" = c(1, 1, 1),
  "ARIMA(0,1,1)" = c(0, 1, 1),
  "ARIMA(1,1,0)" = c(1, 1, 0),
  "ARIMA(2,1,1)" = c(2, 1, 1)
)

results <- data.frame(
  Model = character(),
  Method = character(),
  AIC = numeric(),
  BIC = numeric(),
  LogLik = numeric(),
  RSS = numeric(),
  stringsAsFactors = FALSE
)

for (model_name in names(models_to_test)) {
  order <- models_to_test[[model_name]]

  cat(sprintf("Fitting %s...\n", model_name))

  # CSS-ML method
  tryCatch({
    fit_css <- arima(prices, order = order, method = "CSS-ML", include.mean = FALSE)
    res_css <- residuals(fit_css)
    res_css <- res_css[is.finite(res_css)]

    results <- rbind(results, data.frame(
      Model = model_name,
      Method = "CSS-ML",
      AIC = AIC(fit_css),
      BIC = BIC(fit_css),
      LogLik = as.numeric(logLik(fit_css)),
      RSS = sum(res_css^2)
    ))
  }, error = function(e) {
    cat(sprintf("  CSS-ML failed: %s\n", e$message))
  })

  # PMM2 method
  tryCatch({
    fit_pmm2 <- arima_pmm2(prices, order = order, include.mean = FALSE, verbose = FALSE)
    res_pmm2 <- fit_pmm2@residuals
    res_pmm2 <- res_pmm2[is.finite(res_pmm2)]

    results <- rbind(results, data.frame(
      Model = model_name,
      Method = "PMM2",
      AIC = AIC(fit_pmm2),
      BIC = BIC(fit_pmm2),
      LogLik = as.numeric(logLik(fit_pmm2)),
      RSS = sum(res_pmm2^2)
    ))
  }, error = function(e) {
    cat(sprintf("  PMM2 failed: %s\n", e$message))
  })
}

cat("\n=======================================================\n")
cat("RESULTS SUMMARY\n")
cat("=======================================================\n\n")

# Print results table
print(results, digits = 3, row.names = FALSE)

# Find best models
cat("\n-------------------------------------------------------\n")
cat("BEST MODELS BY INFORMATION CRITERIA\n")
cat("-------------------------------------------------------\n")

best_aic <- results[which.min(results$AIC), ]
best_bic <- results[which.min(results$BIC), ]

cat(sprintf("\nBest by AIC: %s (%s)\n", best_aic$Model, best_aic$Method))
cat(sprintf("  AIC = %.2f, BIC = %.2f, LogLik = %.2f\n",
            best_aic$AIC, best_aic$BIC, best_aic$LogLik))

cat(sprintf("\nBest by BIC: %s (%s)\n", best_bic$Model, best_bic$Method))
cat(sprintf("  AIC = %.2f, BIC = %.2f, LogLik = %.2f\n",
            best_bic$AIC, best_bic$BIC, best_bic$LogLik))

# Calculate AIC and BIC improvements
css_ml_rows <- results[results$Method == "CSS-ML", ]
pmm2_rows <- results[results$Method == "PMM2", ]

if (nrow(css_ml_rows) > 0 && nrow(pmm2_rows) > 0) {
  cat("\n-------------------------------------------------------\n")
  cat("PMM2 vs CSS-ML COMPARISON\n")
  cat("-------------------------------------------------------\n\n")

  for (i in 1:nrow(css_ml_rows)) {
    model <- css_ml_rows$Model[i]
    pmm2_row <- pmm2_rows[pmm2_rows$Model == model, ]

    if (nrow(pmm2_row) > 0) {
      aic_diff <- pmm2_row$AIC - css_ml_rows$AIC[i]
      bic_diff <- pmm2_row$BIC - css_ml_rows$BIC[i]

      cat(sprintf("%s:\n", model))
      cat(sprintf("  AIC difference (PMM2 - CSS-ML): %+.2f %s\n",
                  aic_diff,
                  ifelse(aic_diff < 0, "(PMM2 better)", "(CSS-ML better)")))
      cat(sprintf("  BIC difference (PMM2 - CSS-ML): %+.2f %s\n\n",
                  bic_diff,
                  ifelse(bic_diff < 0, "(PMM2 better)", "(CSS-ML better)")))
    }
  }
}

cat("\n=======================================================\n")
cat("For detailed analysis with plots, run:\n")
cat("  source('PMM2-ARIMA/run_full_study.R')\n")
cat("=======================================================\n")
