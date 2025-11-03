#!/usr/bin/env Rscript
# =============================================================================
# WTI Crude Oil: Enhanced Diagnostics
# =============================================================================
# This script generates enhanced diagnostic plots and statistics for the
# WTI case study, including p-values, Q-Q plots, and ACF/PACF.
#
# Addresses reviewer comment: "Add diagnostic plots with p-values"
#
# Outputs:
# - Q-Q plots for residuals (CSS vs PMM2)
# - ACF/PACF plots
# - Diagnostic statistics table with p-values
# - Residual plots
#
# Author: PMM2-ARIMA Research Team
# Date: 2025-11-03
# =============================================================================

suppressPackageStartupMessages({
  library(EstemPMM)
  if (!requireNamespace("tseries", quietly = TRUE)) {
    install.packages("tseries", repos = "https://cloud.r-project.org")
  }
  library(tseries)
})

# =============================================================================
# Configuration
# =============================================================================

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
original_wd <- getwd()
on.exit(setwd(original_wd), add = TRUE)
setwd(project_root)

data_path <- file.path(project_root, "data", "DCOILWTICO.csv")
plots_dir <- file.path(project_root, "results", "plots")
results_dir <- file.path(project_root, "results")

if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
model_order <- c(1, 1, 0)  # Default: ARIMA(1,1,0)

for (arg in args) {
  if (grepl("^--order=", arg)) {
    order_str <- sub("^--order=", "", arg)
    model_order <- as.integer(strsplit(order_str, ",")[[1]])
  }
}

cat("=== WTI Enhanced Diagnostics ===\n\n")
cat(sprintf("Model: ARIMA(%d,%d,%d)\n\n", model_order[1], model_order[2], model_order[3]))

# =============================================================================
# Load Data and Fit Models
# =============================================================================

cat("Loading WTI data...\n")
wti_data <- read.csv(data_path, stringsAsFactors = FALSE)
wti_data$observation_date <- as.Date(wti_data$observation_date)
wti_data <- wti_data[order(wti_data$observation_date), ]

price_series <- wti_data$DCOILWTICO
price_series <- price_series[is.finite(price_series)]

cat(sprintf("  Total observations: %d\n\n", length(price_series)))

cat("Fitting models...\n")

# CSS fit
css_fit <- tryCatch(
  stats::arima(price_series, order = model_order, method = "CSS", include.mean = FALSE),
  error = function(e) {
    cat("Warning: CSS fit failed\n")
    NULL
  }
)

# PMM2 fit
pmm2_fit <- tryCatch(
  EstemPMM::arima_pmm2(price_series, order = model_order, include.mean = FALSE),
  error = function(e) {
    cat("Warning: PMM2 fit failed\n")
    NULL
  }
)

if (is.null(css_fit) || is.null(pmm2_fit)) {
  stop("Model fitting failed. Cannot proceed with diagnostics.")
}

cat("  Models fitted successfully\n\n")

# Extract residuals
css_residuals <- residuals(css_fit)
pmm2_residuals <- residuals(pmm2_fit)

css_residuals <- css_residuals[is.finite(css_residuals)]
pmm2_residuals <- pmm2_residuals[is.finite(pmm2_residuals)]

# =============================================================================
# Diagnostic Statistics
# =============================================================================

cat("Computing diagnostic statistics...\n\n")

# Function to compute all diagnostics
compute_diagnostics <- function(residuals, label) {
  n <- length(residuals)

  # Ljung-Box test (portmanteau test for autocorrelation)
  # Test at lag = sqrt(n) as commonly recommended
  lb_lag <- max(10, ceiling(sqrt(n)))
  lb_test <- tryCatch(
    Box.test(residuals, lag = lb_lag, type = "Ljung-Box"),
    error = function(e) list(statistic = NA, p.value = NA)
  )

  # Jarque-Bera test for normality
  jb_test <- tryCatch(
    tseries::jarque.bera.test(residuals),
    error = function(e) list(statistic = NA, p.value = NA)
  )

  # Shapiro-Wilk test for normality (if n <= 5000)
  sw_test <- if (n <= 5000) {
    tryCatch(
      shapiro.test(residuals),
      error = function(e) list(statistic = NA, p.value = NA)
    )
  } else {
    list(statistic = NA, p.value = NA)
  }

  # ARCH test for heteroskedasticity
  arch_test <- tryCatch({
    # Simplified ARCH test: regress squared residuals on lagged squared residuals
    resid_sq <- residuals^2
    resid_sq_lag <- c(NA, resid_sq[-length(resid_sq)])
    valid_idx <- is.finite(resid_sq) & is.finite(resid_sq_lag)
    if (sum(valid_idx) > 10) {
      lm_fit <- lm(resid_sq[valid_idx] ~ resid_sq_lag[valid_idx])
      lm_summary <- summary(lm_fit)
      f_stat <- lm_summary$fstatistic
      if (!is.null(f_stat)) {
        p_val <- pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE)
        list(statistic = f_stat[1], p.value = p_val)
      } else {
        list(statistic = NA, p.value = NA)
      }
    } else {
      list(statistic = NA, p.value = NA)
    }
  }, error = function(e) list(statistic = NA, p.value = NA))

  # Basic statistics
  list(
    method = label,
    n = n,
    mean = mean(residuals),
    sd = sd(residuals),
    skewness = moments_skewness(residuals),
    kurtosis = moments_kurtosis(residuals),
    ljung_box_stat = lb_test$statistic,
    ljung_box_pval = lb_test$p.value,
    ljung_box_lag = lb_lag,
    jarque_bera_stat = jb_test$statistic,
    jarque_bera_pval = jb_test$p.value,
    shapiro_wilk_stat = sw_test$statistic,
    shapiro_wilk_pval = sw_test$p.value,
    arch_stat = arch_test$statistic,
    arch_pval = arch_test$p.value
  )
}

# Helper functions for moments
moments_skewness <- function(x) {
  n <- length(x)
  if (n < 3) return(NA)
  m <- mean(x)
  m2 <- mean((x - m)^2)
  m3 <- mean((x - m)^3)
  if (m2 <= 0) return(NA)
  m3 / (m2^1.5)
}

moments_kurtosis <- function(x) {
  n <- length(x)
  if (n < 4) return(NA)
  m <- mean(x)
  m2 <- mean((x - m)^2)
  m4 <- mean((x - m)^4)
  if (m2 <= 0) return(NA)
  (m4 / (m2^2)) - 3  # Excess kurtosis
}

# Compute diagnostics
css_diagnostics <- compute_diagnostics(css_residuals, "CSS")
pmm2_diagnostics <- compute_diagnostics(pmm2_residuals, "PMM2")

# Create data frame
diagnostics_df <- data.frame(
  method = c(css_diagnostics$method, pmm2_diagnostics$method),
  n_residuals = c(css_diagnostics$n, pmm2_diagnostics$n),
  mean = c(css_diagnostics$mean, pmm2_diagnostics$mean),
  sd = c(css_diagnostics$sd, pmm2_diagnostics$sd),
  skewness = c(css_diagnostics$skewness, pmm2_diagnostics$skewness),
  kurtosis = c(css_diagnostics$kurtosis, pmm2_diagnostics$kurtosis),
  ljung_box_stat = c(css_diagnostics$ljung_box_stat, pmm2_diagnostics$ljung_box_stat),
  ljung_box_pval = c(css_diagnostics$ljung_box_pval, pmm2_diagnostics$ljung_box_pval),
  ljung_box_lag = c(css_diagnostics$ljung_box_lag, pmm2_diagnostics$ljung_box_lag),
  jarque_bera_stat = c(css_diagnostics$jarque_bera_stat, pmm2_diagnostics$jarque_bera_stat),
  jarque_bera_pval = c(css_diagnostics$jarque_bera_pval, pmm2_diagnostics$jarque_bera_pval),
  shapiro_wilk_stat = c(css_diagnostics$shapiro_wilk_stat, pmm2_diagnostics$shapiro_wilk_stat),
  shapiro_wilk_pval = c(css_diagnostics$shapiro_wilk_pval, pmm2_diagnostics$shapiro_wilk_pval),
  arch_stat = c(css_diagnostics$arch_stat, pmm2_diagnostics$arch_stat),
  arch_pval = c(css_diagnostics$arch_pval, pmm2_diagnostics$arch_pval),
  stringsAsFactors = FALSE
)

# Save diagnostics table
write.csv(
  diagnostics_df,
  file.path(results_dir, "wti_diagnostics_statistics.csv"),
  row.names = FALSE
)

cat("Diagnostic Statistics:\n")
print(diagnostics_df)
cat("\n")

# =============================================================================
# Generate Diagnostic Plots
# =============================================================================

cat("Generating diagnostic plots...\n")

# 1. Q-Q Plots
png(file.path(plots_dir, "wti_qq_plots.png"), width = 1200, height = 600, res = 120)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

qqnorm(css_residuals, main = "Q-Q Plot: CSS Residuals", pch = 20, col = rgb(0, 0, 0, 0.3))
qqline(css_residuals, col = "red", lwd = 2)
grid()

qqnorm(pmm2_residuals, main = "Q-Q Plot: PMM2 Residuals", pch = 20, col = rgb(0, 0, 1, 0.3))
qqline(pmm2_residuals, col = "red", lwd = 2)
grid()

dev.off()
cat("  Saved: wti_qq_plots.png\n")

# 2. ACF/PACF Plots
png(file.path(plots_dir, "wti_acf_pacf.png"), width = 1200, height = 800, res = 120)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

acf(css_residuals, main = "ACF: CSS Residuals", lag.max = 40)
pacf(css_residuals, main = "PACF: CSS Residuals", lag.max = 40)

acf(pmm2_residuals, main = "ACF: PMM2 Residuals", lag.max = 40)
pacf(pmm2_residuals, main = "PACF: PMM2 Residuals", lag.max = 40)

dev.off()
cat("  Saved: wti_acf_pacf.png\n")

# 3. Residual Time Series Plots
png(file.path(plots_dir, "wti_residual_time_series.png"), width = 1200, height = 800, res = 120)
par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))

plot(css_residuals, type = "l", col = "blue",
     main = "CSS Residuals Over Time",
     xlab = "Index", ylab = "Residual")
abline(h = 0, col = "red", lty = 2)
grid()

plot(pmm2_residuals, type = "l", col = "darkgreen",
     main = "PMM2 Residuals Over Time",
     xlab = "Index", ylab = "Residual")
abline(h = 0, col = "red", lty = 2)
grid()

dev.off()
cat("  Saved: wti_residual_time_series.png\n")

# 4. Residual Histograms with Normal Overlay
png(file.path(plots_dir, "wti_residual_histograms.png"), width = 1200, height = 600, res = 120)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

hist(css_residuals, breaks = 30, freq = FALSE, col = "lightblue",
     main = "CSS Residuals Distribution",
     xlab = "Residual", border = "white")
curve(dnorm(x, mean = mean(css_residuals), sd = sd(css_residuals)),
      add = TRUE, col = "red", lwd = 2)
grid()

hist(pmm2_residuals, breaks = 30, freq = FALSE, col = "lightgreen",
     main = "PMM2 Residuals Distribution",
     xlab = "Residual", border = "white")
curve(dnorm(x, mean = mean(pmm2_residuals), sd = sd(pmm2_residuals)),
      add = TRUE, col = "red", lwd = 2)
grid()

dev.off()
cat("  Saved: wti_residual_histograms.png\n")

# =============================================================================
# Summary
# =============================================================================

cat("\n=== Summary of Diagnostic Tests ===\n\n")

cat("Ljung-Box Test (H0: No autocorrelation):\n")
cat(sprintf("  CSS:  Q(%.0f) = %.4f, p-value = %.4f %s\n",
            css_diagnostics$ljung_box_lag,
            css_diagnostics$ljung_box_stat,
            css_diagnostics$ljung_box_pval,
            ifelse(css_diagnostics$ljung_box_pval > 0.05, "[OK]", "[REJECT]")))
cat(sprintf("  PMM2: Q(%.0f) = %.4f, p-value = %.4f %s\n\n",
            pmm2_diagnostics$ljung_box_lag,
            pmm2_diagnostics$ljung_box_stat,
            pmm2_diagnostics$ljung_box_pval,
            ifelse(pmm2_diagnostics$ljung_box_pval > 0.05, "[OK]", "[REJECT]")))

cat("Jarque-Bera Test (H0: Normality):\n")
cat(sprintf("  CSS:  JB = %.4f, p-value = %.4f %s\n",
            css_diagnostics$jarque_bera_stat,
            css_diagnostics$jarque_bera_pval,
            ifelse(css_diagnostics$jarque_bera_pval > 0.05, "[Normal]", "[Non-normal]")))
cat(sprintf("  PMM2: JB = %.4f, p-value = %.4f %s\n\n",
            pmm2_diagnostics$jarque_bera_stat,
            pmm2_diagnostics$jarque_bera_pval,
            ifelse(pmm2_diagnostics$jarque_bera_pval > 0.05, "[Normal]", "[Non-normal]")))

cat("✓ Enhanced diagnostics complete!\n")
cat("  Generated 4 diagnostic plot files in results/plots/\n")
cat("  Diagnostic statistics saved to wti_diagnostics_statistics.csv\n\n")
