# ==============================================================================
# Monte Carlo Comparison of Seasonal Time Series Models
# Comparing PMM2 vs Classical Methods (CSS, ML) for SAR, SMA, SARMA, SARIMA
# Sample sizes: 100, 200, 500
# ==============================================================================

library(EstemPMM)

# Set global parameters
set.seed(12345)
N_REPLICATIONS <- 500  # Number of Monte Carlo replications
SAMPLE_SIZES <- c(100, 200, 500)
SEASONAL_PERIOD <- 12

cat("==============================================================================\n")
cat("Monte Carlo Comparison of Seasonal Models with PMM2\n")
cat("==============================================================================\n")
cat("Replications:", N_REPLICATIONS, "\n")
cat("Sample sizes:", paste(SAMPLE_SIZES, collapse = ", "), "\n")
cat("Seasonal period:", SEASONAL_PERIOD, "\n\n")

# ==============================================================================
# Helper Functions
# ==============================================================================

#' Generate SAR(p,P)_s time series with asymmetric innovations
#'
#' @param n Sample size
#' @param ar Non-seasonal AR coefficients
#' @param sar Seasonal AR coefficients
#' @param s Seasonal period
#' @param innovation_type Distribution type: "gamma", "lognormal", "exponential"
#' @return Time series vector
generate_sar_data <- function(n, ar = NULL, sar = NULL, s = 12,
                              innovation_type = "gamma") {
  p <- length(ar)
  P <- length(sar)

  # Generate innovations
  innov <- switch(innovation_type,
                  "gamma" = rgamma(n, shape = 2, scale = 1) - 2,
                  "lognormal" = rlnorm(n, meanlog = 0, sdlog = 0.5) - exp(0.125),
                  "exponential" = rexp(n, rate = 1) - 1,
                  rnorm(n))

  y <- numeric(n)
  max_lag <- max(p, P * s)

  for (t in (max_lag + 1):n) {
    ar_term <- if (p > 0) sum(ar * y[(t-1):(t-p)]) else 0
    sar_term <- if (P > 0) sum(sar * y[t - s * (1:P)]) else 0
    y[t] <- ar_term + sar_term + innov[t]
  }

  return(y[(max_lag + 10):n])  # Remove burn-in
}

#' Generate SMA(Q)_s time series with asymmetric innovations
generate_sma_data <- function(n, sma = NULL, s = 12, innovation_type = "gamma") {
  Q <- length(sma)

  innov <- switch(innovation_type,
                  "gamma" = rgamma(n, shape = 2, scale = 1) - 2,
                  "lognormal" = rlnorm(n, meanlog = 0, sdlog = 0.5) - exp(0.125),
                  "exponential" = rexp(n, rate = 1) - 1,
                  rnorm(n))

  y <- numeric(n)
  eps <- numeric(n)

  for (t in 1:n) {
    sma_term <- if (t > s && Q > 0) sum(sma * eps[t - s * (1:Q)]) else 0
    eps[t] <- innov[t]
    y[t] <- sma_term + eps[t]
  }

  return(y[(s * Q + 10):n])  # Remove burn-in
}

#' Generate SARMA time series
generate_sarma_data <- function(n, ar = NULL, sar = NULL, ma = NULL, sma = NULL,
                                s = 12, innovation_type = "gamma") {
  # Use arima.sim for SARMA with custom innovations
  innov <- switch(innovation_type,
                  "gamma" = rgamma(n, shape = 2, scale = 1) - 2,
                  "lognormal" = rlnorm(n, meanlog = 0, sdlog = 0.5) - exp(0.125),
                  "exponential" = rexp(n, rate = 1) - 1,
                  rnorm(n))

  model_list <- list()
  if (!is.null(ar)) model_list$ar <- ar
  if (!is.null(ma)) model_list$ma <- ma

  seasonal_list <- list(period = s)
  if (!is.null(sar)) seasonal_list$sar <- sar
  if (!is.null(sma)) seasonal_list$sma <- sma

  if (length(seasonal_list) > 1) {
    model_list$seasonal <- seasonal_list
  }

  y <- arima.sim(n = n, model = model_list, innov = innov)
  return(as.numeric(y))
}

#' Compute estimation metrics
compute_metrics <- function(estimates, true_params) {
  bias <- mean(estimates - true_params)
  rmse <- sqrt(mean((estimates - true_params)^2))
  variance <- var(estimates)
  mae <- mean(abs(estimates - true_params))

  list(bias = bias, rmse = rmse, variance = variance, mae = mae)
}

#' Safe model fitting wrapper
safe_fit <- function(fit_func, ...) {
  tryCatch({
    fit <- fit_func(...)
    # Convert S4 object to list for easier access
    list(
      coefficients = fit@coefficients,
      residuals = fit@residuals,
      convergence = fit@convergence,
      m2 = fit@m2,
      m3 = fit@m3,
      m4 = fit@m4
    )
  },
  error = function(e) {
    list(coefficients = NA, convergence = FALSE,
         residuals = NA, m2 = NA, m3 = NA, m4 = NA)
  })
}

# ==============================================================================
# Scenario 1: SAR(1,1)_12 Model
# ==============================================================================

cat("\n==============================================================================\n")
cat("Scenario 1: SAR(1,1)_12 Model\n")
cat("==============================================================================\n")

true_sar <- list(ar = 0.5, sar = 0.6)

results_sar <- list()

for (n in SAMPLE_SIZES) {
  cat("\nSample size:", n, "\n")

  pmm2_coef1 <- numeric(N_REPLICATIONS)
  pmm2_coef2 <- numeric(N_REPLICATIONS)
  css_coef1 <- numeric(N_REPLICATIONS)
  css_coef2 <- numeric(N_REPLICATIONS)

  pmm2_g <- numeric(N_REPLICATIONS)
  pmm2_converged <- 0
  css_converged <- 0

  pb <- txtProgressBar(min = 0, max = N_REPLICATIONS, style = 3)

  for (i in 1:N_REPLICATIONS) {
    # Generate data
    y <- generate_sar_data(n + 50, ar = true_sar$ar, sar = true_sar$sar,
                          innovation_type = "gamma")

    # PMM2 estimation
    fit_pmm2 <- safe_fit(sar_pmm2, y, order = c(1, 1),
                        season = list(period = SEASONAL_PERIOD),
                        method = "pmm2", verbose = FALSE)

    if (!is.na(fit_pmm2$coefficients[1])) {
      pmm2_coef1[i] <- fit_pmm2$coefficients[1]
      pmm2_coef2[i] <- fit_pmm2$coefficients[2]
      pmm2_g[i] <- pmm2_variance_factor(fit_pmm2$m2, fit_pmm2$m3, fit_pmm2$m4)$g
      if (fit_pmm2$convergence) pmm2_converged <- pmm2_converged + 1
    } else {
      pmm2_coef1[i] <- NA
      pmm2_coef2[i] <- NA
      pmm2_g[i] <- NA
    }

    # CSS estimation
    fit_css <- safe_fit(sar_pmm2, y, order = c(1, 1),
                       season = list(period = SEASONAL_PERIOD),
                       method = "css", verbose = FALSE)

    if (!is.na(fit_css$coefficients[1])) {
      css_coef1[i] <- fit_css$coefficients[1]
      css_coef2[i] <- fit_css$coefficients[2]
      if (fit_css$convergence) css_converged <- css_converged + 1
    } else {
      css_coef1[i] <- NA
      css_coef2[i] <- NA
    }

    setTxtProgressBar(pb, i)
  }
  close(pb)

  # Remove NAs
  valid_pmm2 <- !is.na(pmm2_coef1) & !is.na(pmm2_coef2)
  valid_css <- !is.na(css_coef1) & !is.na(css_coef2)

  # Compute metrics
  metrics_pmm2_ar <- compute_metrics(pmm2_coef1[valid_pmm2], true_sar$ar)
  metrics_pmm2_sar <- compute_metrics(pmm2_coef2[valid_pmm2], true_sar$sar)
  metrics_css_ar <- compute_metrics(css_coef1[valid_css], true_sar$ar)
  metrics_css_sar <- compute_metrics(css_coef2[valid_css], true_sar$sar)

  # Store results
  results_sar[[paste0("n", n)]] <- list(
    pmm2 = list(
      ar = metrics_pmm2_ar,
      sar = metrics_pmm2_sar,
      mean_g = mean(pmm2_g[valid_pmm2], na.rm = TRUE),
      convergence_rate = pmm2_converged / N_REPLICATIONS
    ),
    css = list(
      ar = metrics_css_ar,
      sar = metrics_css_sar,
      convergence_rate = css_converged / N_REPLICATIONS
    )
  )

  # Print summary
  cat("\n\nResults for n =", n, "\n")
  cat("PMM2 - AR(1) coefficient:\n")
  cat("  Bias:", round(metrics_pmm2_ar$bias, 4), "\n")
  cat("  RMSE:", round(metrics_pmm2_ar$rmse, 4), "\n")
  cat("  Variance:", round(metrics_pmm2_ar$variance, 4), "\n")

  cat("PMM2 - SAR(1) coefficient:\n")
  cat("  Bias:", round(metrics_pmm2_sar$bias, 4), "\n")
  cat("  RMSE:", round(metrics_pmm2_sar$rmse, 4), "\n")
  cat("  Variance:", round(metrics_pmm2_sar$variance, 4), "\n")

  cat("PMM2 - Mean g:", round(mean(pmm2_g[valid_pmm2], na.rm = TRUE), 4), "\n")
  cat("PMM2 - Convergence rate:", round(pmm2_converged / N_REPLICATIONS, 3), "\n")

  cat("\nCSS - AR(1) coefficient:\n")
  cat("  Bias:", round(metrics_css_ar$bias, 4), "\n")
  cat("  RMSE:", round(metrics_css_ar$rmse, 4), "\n")
  cat("  Variance:", round(metrics_css_ar$variance, 4), "\n")

  cat("CSS - SAR(1) coefficient:\n")
  cat("  Bias:", round(metrics_css_sar$bias, 4), "\n")
  cat("  RMSE:", round(metrics_css_sar$rmse, 4), "\n")
  cat("  Variance:", round(metrics_css_sar$variance, 4), "\n")

  cat("CSS - Convergence rate:", round(css_converged / N_REPLICATIONS, 3), "\n")

  # Variance reduction
  var_reduction_ar <- 1 - metrics_pmm2_ar$variance / metrics_css_ar$variance
  var_reduction_sar <- 1 - metrics_pmm2_sar$variance / metrics_css_sar$variance

  cat("\nVariance reduction (PMM2 vs CSS):\n")
  cat("  AR coefficient:", round(var_reduction_ar * 100, 2), "%\n")
  cat("  SAR coefficient:", round(var_reduction_sar * 100, 2), "%\n")
}

# ==============================================================================
# Scenario 2: SMA(1)_12 Model
# ==============================================================================

cat("\n\n==============================================================================\n")
cat("Scenario 2: SMA(1)_12 Model\n")
cat("==============================================================================\n")

true_sma <- 0.6

results_sma <- list()

for (n in SAMPLE_SIZES) {
  cat("\nSample size:", n, "\n")

  pmm2_coef <- numeric(N_REPLICATIONS)
  css_coef <- numeric(N_REPLICATIONS)
  pmm2_g <- numeric(N_REPLICATIONS)
  pmm2_converged <- 0
  css_converged <- 0

  pb <- txtProgressBar(min = 0, max = N_REPLICATIONS, style = 3)

  for (i in 1:N_REPLICATIONS) {
    y <- generate_sma_data(n + 50, sma = true_sma, innovation_type = "gamma")

    # PMM2
    fit_pmm2 <- safe_fit(sma_pmm2, y, order = 1,
                        season = list(period = SEASONAL_PERIOD),
                        method = "pmm2", verbose = FALSE)

    if (!is.na(fit_pmm2$coefficients[1])) {
      pmm2_coef[i] <- fit_pmm2$coefficients[1]
      pmm2_g[i] <- pmm2_variance_factor(fit_pmm2$m2, fit_pmm2$m3, fit_pmm2$m4)$g
      if (fit_pmm2$convergence) pmm2_converged <- pmm2_converged + 1
    } else {
      pmm2_coef[i] <- NA
      pmm2_g[i] <- NA
    }

    # CSS
    fit_css <- safe_fit(sma_pmm2, y, order = 1,
                       season = list(period = SEASONAL_PERIOD),
                       method = "css", verbose = FALSE)

    if (!is.na(fit_css$coefficients[1])) {
      css_coef[i] <- fit_css$coefficients[1]
      if (fit_css$convergence) css_converged <- css_converged + 1
    } else {
      css_coef[i] <- NA
    }

    setTxtProgressBar(pb, i)
  }
  close(pb)

  valid_pmm2 <- !is.na(pmm2_coef)
  valid_css <- !is.na(css_coef)

  metrics_pmm2 <- compute_metrics(pmm2_coef[valid_pmm2], true_sma)
  metrics_css <- compute_metrics(css_coef[valid_css], true_sma)

  results_sma[[paste0("n", n)]] <- list(
    pmm2 = list(
      sma = metrics_pmm2,
      mean_g = mean(pmm2_g[valid_pmm2], na.rm = TRUE),
      convergence_rate = pmm2_converged / N_REPLICATIONS
    ),
    css = list(
      sma = metrics_css,
      convergence_rate = css_converged / N_REPLICATIONS
    )
  )

  cat("\n\nResults for n =", n, "\n")
  cat("PMM2 - SMA(1) coefficient:\n")
  cat("  Bias:", round(metrics_pmm2$bias, 4), "\n")
  cat("  RMSE:", round(metrics_pmm2$rmse, 4), "\n")
  cat("  Variance:", round(metrics_pmm2$variance, 4), "\n")
  cat("  Mean g:", round(mean(pmm2_g[valid_pmm2], na.rm = TRUE), 4), "\n")

  cat("\nCSS - SMA(1) coefficient:\n")
  cat("  Bias:", round(metrics_css$bias, 4), "\n")
  cat("  RMSE:", round(metrics_css$rmse, 4), "\n")
  cat("  Variance:", round(metrics_css$variance, 4), "\n")

  var_reduction <- 1 - metrics_pmm2$variance / metrics_css$variance
  cat("\nVariance reduction:", round(var_reduction * 100, 2), "%\n")
}

# ==============================================================================
# Scenario 3: SARMA(1,0)×(1,1)_12 Model
# ==============================================================================

cat("\n\n==============================================================================\n")
cat("Scenario 3: SARMA(1,0)×(1,1)_12 Model\n")
cat("==============================================================================\n")

true_sarma <- list(ar = 0.5, sar = 0.6, sma = 0.4)

results_sarma <- list()

for (n in SAMPLE_SIZES) {
  cat("\nSample size:", n, "\n")

  pmm2_coef <- matrix(NA, N_REPLICATIONS, 3)
  css_coef <- matrix(NA, N_REPLICATIONS, 3)
  pmm2_g <- numeric(N_REPLICATIONS)
  pmm2_converged <- 0
  css_converged <- 0

  pb <- txtProgressBar(min = 0, max = N_REPLICATIONS, style = 3)

  for (i in 1:N_REPLICATIONS) {
    y <- generate_sarma_data(n + 50, ar = true_sarma$ar, sar = true_sarma$sar,
                            sma = true_sarma$sma, innovation_type = "gamma")

    # PMM2
    fit_pmm2 <- safe_fit(sarma_pmm2, y, order = c(1, 1, 0, 1),
                        season = list(period = SEASONAL_PERIOD),
                        method = "pmm2", verbose = FALSE)

    if (length(fit_pmm2$coefficients) == 3 && !any(is.na(fit_pmm2$coefficients))) {
      pmm2_coef[i, ] <- fit_pmm2$coefficients
      pmm2_g[i] <- pmm2_variance_factor(fit_pmm2$m2, fit_pmm2$m3, fit_pmm2$m4)$g
      if (fit_pmm2$convergence) pmm2_converged <- pmm2_converged + 1
    }

    # CSS
    fit_css <- safe_fit(sarma_pmm2, y, order = c(1, 1, 0, 1),
                       season = list(period = SEASONAL_PERIOD),
                       method = "css", verbose = FALSE)

    if (length(fit_css$coefficients) == 3 && !any(is.na(fit_css$coefficients))) {
      css_coef[i, ] <- fit_css$coefficients
      if (fit_css$convergence) css_converged <- css_converged + 1
    }

    setTxtProgressBar(pb, i)
  }
  close(pb)

  valid_pmm2 <- complete.cases(pmm2_coef)
  valid_css <- complete.cases(css_coef)

  true_params <- c(true_sarma$ar, true_sarma$sar, true_sarma$sma)

  metrics_pmm2 <- lapply(1:3, function(j) {
    compute_metrics(pmm2_coef[valid_pmm2, j], true_params[j])
  })

  metrics_css <- lapply(1:3, function(j) {
    compute_metrics(css_coef[valid_css, j], true_params[j])
  })

  results_sarma[[paste0("n", n)]] <- list(
    pmm2 = list(
      coefficients = metrics_pmm2,
      mean_g = mean(pmm2_g[valid_pmm2], na.rm = TRUE),
      convergence_rate = pmm2_converged / N_REPLICATIONS
    ),
    css = list(
      coefficients = metrics_css,
      convergence_rate = css_converged / N_REPLICATIONS
    )
  )

  cat("\n\nResults for n =", n, "\n")
  param_names <- c("AR(1)", "SAR(1)", "SMA(1)")

  for (j in 1:3) {
    cat("\nPMM2 -", param_names[j], ":\n")
    cat("  Bias:", round(metrics_pmm2[[j]]$bias, 4), "\n")
    cat("  RMSE:", round(metrics_pmm2[[j]]$rmse, 4), "\n")
    cat("  Variance:", round(metrics_pmm2[[j]]$variance, 4), "\n")
  }

  cat("\nPMM2 - Mean g:", round(mean(pmm2_g[valid_pmm2], na.rm = TRUE), 4), "\n")

  for (j in 1:3) {
    cat("\nCSS -", param_names[j], ":\n")
    cat("  Bias:", round(metrics_css[[j]]$bias, 4), "\n")
    cat("  RMSE:", round(metrics_css[[j]]$rmse, 4), "\n")
    cat("  Variance:", round(metrics_css[[j]]$variance, 4), "\n")
  }

  cat("\nVariance reduction (PMM2 vs CSS):\n")
  for (j in 1:3) {
    var_red <- 1 - metrics_pmm2[[j]]$variance / metrics_css[[j]]$variance
    cat(" ", param_names[j], ":", round(var_red * 100, 2), "%\n")
  }
}

# ==============================================================================
# Scenario 4: SARIMA(1,1,0)×(1,1,1)_12 Model
# ==============================================================================

cat("\n\n==============================================================================\n")
cat("Scenario 4: SARIMA(1,1,0)×(1,1,1)_12 Model with differencing\n")
cat("==============================================================================\n")

true_sarima <- list(ar = 0.4, sar = 0.5, sma = 0.6)

results_sarima <- list()

for (n in SAMPLE_SIZES) {
  cat("\nSample size:", n, "\n")

  pmm2_coef <- matrix(NA, N_REPLICATIONS, 3)
  css_coef <- matrix(NA, N_REPLICATIONS, 3)
  pmm2_g <- numeric(N_REPLICATIONS)
  pmm2_converged <- 0
  css_converged <- 0

  pb <- txtProgressBar(min = 0, max = N_REPLICATIONS, style = 3)

  for (i in 1:N_REPLICATIONS) {
    # Generate SARIMA data
    innov <- rgamma(n + 100, shape = 2, scale = 1) - 2
    y <- arima.sim(n = n + 100,
                   list(order = c(1, 1, 0),
                        ar = true_sarima$ar,
                        seasonal = list(order = c(1, 1, 1), period = SEASONAL_PERIOD,
                                       sar = true_sarima$sar, sma = true_sarima$sma)),
                   innov = innov)
    y <- as.numeric(y)

    # PMM2
    fit_pmm2 <- safe_fit(sarima_pmm2, y, order = c(1, 1, 0, 1),
                        seasonal = list(order = c(1, 1), period = SEASONAL_PERIOD),
                        method = "pmm2", verbose = FALSE)

    if (length(fit_pmm2$coefficients) == 3 && !any(is.na(fit_pmm2$coefficients))) {
      pmm2_coef[i, ] <- fit_pmm2$coefficients
      pmm2_g[i] <- pmm2_variance_factor(fit_pmm2$m2, fit_pmm2$m3, fit_pmm2$m4)$g
      if (fit_pmm2$convergence) pmm2_converged <- pmm2_converged + 1
    }

    # CSS
    fit_css <- safe_fit(sarima_pmm2, y, order = c(1, 1, 0, 1),
                       seasonal = list(order = c(1, 1), period = SEASONAL_PERIOD),
                       method = "css", verbose = FALSE)

    if (length(fit_css$coefficients) == 3 && !any(is.na(fit_css$coefficients))) {
      css_coef[i, ] <- fit_css$coefficients
      if (fit_css$convergence) css_converged <- css_converged + 1
    }

    setTxtProgressBar(pb, i)
  }
  close(pb)

  valid_pmm2 <- complete.cases(pmm2_coef)
  valid_css <- complete.cases(css_coef)

  true_params <- c(true_sarima$ar, true_sarima$sar, true_sarima$sma)

  metrics_pmm2 <- lapply(1:3, function(j) {
    compute_metrics(pmm2_coef[valid_pmm2, j], true_params[j])
  })

  metrics_css <- lapply(1:3, function(j) {
    compute_metrics(css_coef[valid_css, j], true_params[j])
  })

  results_sarima[[paste0("n", n)]] <- list(
    pmm2 = list(
      coefficients = metrics_pmm2,
      mean_g = mean(pmm2_g[valid_pmm2], na.rm = TRUE),
      convergence_rate = pmm2_converged / N_REPLICATIONS
    ),
    css = list(
      coefficients = metrics_css,
      convergence_rate = css_converged / N_REPLICATIONS
    )
  )

  cat("\n\nResults for n =", n, "\n")
  param_names <- c("AR(1)", "SAR(1)", "SMA(1)")

  for (j in 1:3) {
    cat("\nPMM2 -", param_names[j], ":\n")
    cat("  Bias:", round(metrics_pmm2[[j]]$bias, 4), "\n")
    cat("  RMSE:", round(metrics_pmm2[[j]]$rmse, 4), "\n")
    cat("  Variance:", round(metrics_pmm2[[j]]$variance, 4), "\n")
  }

  cat("\nPMM2 - Mean g:", round(mean(pmm2_g[valid_pmm2], na.rm = TRUE), 4), "\n")

  for (j in 1:3) {
    cat("\nCSS -", param_names[j], ":\n")
    cat("  Bias:", round(metrics_css[[j]]$bias, 4), "\n")
    cat("  RMSE:", round(metrics_css[[j]]$rmse, 4), "\n")
    cat("  Variance:", round(metrics_css[[j]]$variance, 4), "\n")
  }

  cat("\nVariance reduction (PMM2 vs CSS):\n")
  for (j in 1:3) {
    var_red <- 1 - metrics_pmm2[[j]]$variance / metrics_css[[j]]$variance
    cat(" ", param_names[j], ":", round(var_red * 100, 2), "%\n")
  }
}

# ==============================================================================
# Save Results
# ==============================================================================

cat("\n\n==============================================================================\n")
cat("Saving results...\n")
cat("==============================================================================\n")

all_results <- list(
  sar = results_sar,
  sma = results_sma,
  sarma = results_sarma,
  sarima = results_sarima,
  metadata = list(
    n_replications = N_REPLICATIONS,
    sample_sizes = SAMPLE_SIZES,
    seasonal_period = SEASONAL_PERIOD,
    date = Sys.time()
  )
)

saveRDS(all_results, file = "monte_carlo_seasonal_results.rds")
cat("Results saved to: monte_carlo_seasonal_results.rds\n")

# ==============================================================================
# Summary Report
# ==============================================================================

cat("\n\n==============================================================================\n")
cat("SUMMARY REPORT\n")
cat("==============================================================================\n\n")

# SAR Summary
cat("SAR(1,1)_12 Model - Average Variance Reduction:\n")
for (n in SAMPLE_SIZES) {
  res <- results_sar[[paste0("n", n)]]
  var_red_ar <- 1 - res$pmm2$ar$variance / res$css$ar$variance
  var_red_sar <- 1 - res$pmm2$sar$variance / res$css$sar$variance
  cat("  n =", n, ": AR =", round(var_red_ar * 100, 1), "%, SAR =",
      round(var_red_sar * 100, 1), "%\n")
}

# SMA Summary
cat("\nSMA(1)_12 Model - Average Variance Reduction:\n")
for (n in SAMPLE_SIZES) {
  res <- results_sma[[paste0("n", n)]]
  var_red <- 1 - res$pmm2$sma$variance / res$css$sma$variance
  cat("  n =", n, ":", round(var_red * 100, 1), "%\n")
}

# SARMA Summary
cat("\nSARMA(1,0)×(1,1)_12 Model - Average Variance Reduction:\n")
for (n in SAMPLE_SIZES) {
  res <- results_sarma[[paste0("n", n)]]
  avg_var_red <- mean(sapply(1:3, function(j) {
    1 - res$pmm2$coefficients[[j]]$variance / res$css$coefficients[[j]]$variance
  }))
  cat("  n =", n, ":", round(avg_var_red * 100, 1), "%\n")
}

# SARIMA Summary
cat("\nSARIMA(1,1,0)×(1,1,1)_12 Model - Average Variance Reduction:\n")
for (n in SAMPLE_SIZES) {
  res <- results_sarima[[paste0("n", n)]]
  avg_var_red <- mean(sapply(1:3, function(j) {
    1 - res$pmm2$coefficients[[j]]$variance / res$css$coefficients[[j]]$variance
  }))
  cat("  n =", n, ":", round(avg_var_red * 100, 1), "%\n")
}

cat("\n==============================================================================\n")
cat("Monte Carlo comparison completed!\n")
cat("==============================================================================\n")
