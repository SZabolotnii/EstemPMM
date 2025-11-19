# ==============================================================================
# Adaptive MA Estimator: Auto-select PMM vs MLE
# ==============================================================================
#
# Strategy:
#   1. If innovations ~ Gaussian (|γ₃| < threshold) → Use MLE
#   2. If innovations ~ Asymmetric (|γ₃| ≥ threshold) → Use PMM2
#   3. Fallback to MLE if PMM fails
#
# ==============================================================================

source("R/experimental/01_ma_support.R")

#' Adaptive MA estimator with automatic method selection
#'
#' @param x Numeric vector (time series)
#' @param q MA order
#' @param include.mean Include intercept
#' @param asymmetry_threshold Threshold for |γ₃| to use PMM (default 0.3)
#' @param verbose Print diagnostic information
#'
#' @return List with estimation results and selected method
#'
#' @details
#' Decision rule:
#'   - Compute preliminary residuals via MLE
#'   - Estimate skewness γ̂₃
#'   - If |γ̂₃| < threshold: use MLE (optimal for symmetric)
#'   - If |γ̂₃| ≥ threshold: use PMM2 (better for asymmetric)
#'
#' @export
adaptive_ma_estimator <- function(x, q,
                                   include.mean = TRUE,
                                   asymmetry_threshold = 0.3,
                                   verbose = FALSE) {

  # Step 1: Initial MLE fit for diagnostic
  if (verbose) cat("Step 1: Initial MLE fit for diagnostics...\n")

  mle_fit <- tryCatch({
    stats::arima(x, order = c(0, 0, q), method = "ML", include.mean = include.mean)
  }, error = function(e) {
    warning("MLE initialization failed: ", e$message)
    return(NULL)
  })

  if (is.null(mle_fit)) {
    stop("Cannot initialize estimator (MLE failed)")
  }

  # Step 2: Analyze residuals
  residuals <- as.numeric(mle_fit$residuals)
  residuals <- residuals[!is.na(residuals)]

  moments <- compute_moments(residuals)
  gamma3 <- moments$m3 / (moments$m2^(3/2))  # Skewness coefficient
  gamma4 <- moments$m4 / (moments$m2^2) - 3  # Excess kurtosis

  # Theoretical RE
  if (abs(gamma3) < 1e-6) {
    re_theoretical <- 1.0
  } else {
    denominator <- 2 + gamma4 - gamma3^2
    if (denominator <= 0) {
      re_theoretical <- 1.0
    } else {
      re_theoretical <- (2 + gamma4) / denominator
    }
  }

  if (verbose) {
    cat(sprintf("  Residual moments: μ₂=%.4f, γ₃=%.4f, γ₄=%.4f\n",
                moments$m2, gamma3, gamma4))
    cat(sprintf("  Theoretical RE: %.3f\n", re_theoretical))
  }

  # Step 3: Decision rule
  use_pmm <- abs(gamma3) >= asymmetry_threshold

  if (verbose) {
    cat(sprintf("\nStep 2: Method selection (|γ₃| = %.3f, threshold = %.2f)\n",
                abs(gamma3), asymmetry_threshold))
    if (use_pmm) {
      cat("  → Using PMM2 (asymmetry detected)\n")
    } else {
      cat("  → Using MLE (near-symmetric, PMM offers no advantage)\n")
    }
  }

  # Step 4: Estimate with selected method
  if (!use_pmm) {
    # Use MLE (already computed)
    theta_est <- numeric(q)
    for (j in 1:q) {
      coef_name <- paste0("ma", j)
      theta_est[j] <- mle_fit$coef[coef_name]
    }

    mu_est <- if (include.mean && "intercept" %in% names(mle_fit$coef)) {
      mle_fit$coef["intercept"]
    } else {
      0
    }

    result <- list(
      theta_est = theta_est,
      mu_est = mu_est,
      method = "MLE",
      gamma3 = gamma3,
      gamma4 = gamma4,
      re_theoretical = re_theoretical,
      converged = TRUE,
      iterations = NA,
      mle_fit = mle_fit
    )

  } else {
    # Use PMM2 iterative
    if (verbose) cat("\nStep 3: Running PMM2 iterative...\n")

    pmm_result <- tryCatch({
      create_ma_design_matrix(x, q = q,
                              pmm_method = "pmm1",  # Currently only PMM1
                              include.mean = include.mean,
                              verbose = verbose,
                              max_iter = 50,
                              tol = 1e-6)
    }, error = function(e) {
      warning("PMM2 failed: ", e$message, "\nFalling back to MLE.")
      return(NULL)
    })

    if (is.null(pmm_result) || !pmm_result$converged) {
      # Fallback to MLE
      if (verbose) cat("⚠ PMM2 failed or didn't converge. Using MLE instead.\n")

      theta_est <- numeric(q)
      for (j in 1:q) {
        theta_est[j] <- mle_fit$coef[paste0("ma", j)]
      }

      result <- list(
        theta_est = theta_est,
        mu_est = if (include.mean) mle_fit$coef["intercept"] else 0,
        method = "MLE (PMM2 fallback)",
        gamma3 = gamma3,
        gamma4 = gamma4,
        re_theoretical = re_theoretical,
        converged = TRUE,
        iterations = NA,
        mle_fit = mle_fit
      )

    } else {
      # PMM2 success
      result <- list(
        theta_est = pmm_result$theta_est,
        mu_est = pmm_result$mu_est,
        method = "PMM2",
        gamma3 = gamma3,
        gamma4 = gamma4,
        re_theoretical = re_theoretical,
        converged = pmm_result$converged,
        iterations = pmm_result$iterations,
        pmm_result = pmm_result
      )
    }
  }

  if (verbose) {
    cat("\n=== Final Estimates ===\n")
    cat("Method:", result$method, "\n")
    cat("θ̂ = [", paste(round(result$theta_est, 4), collapse = ", "), "]\n")
    if (include.mean) {
      cat("μ̂ =", round(result$mu_est, 4), "\n")
    }
    cat("\n")
  }

  return(result)
}


#' Adaptive seasonal MA estimator
#'
#' @param x Time series
#' @param Q Seasonal MA order
#' @param s Seasonal period
#' @param include.mean Include mean
#' @param asymmetry_threshold Threshold for method selection
#' @param verbose Verbose output
#'
#' @export
adaptive_sma_estimator <- function(x, Q, s,
                                    include.mean = TRUE,
                                    asymmetry_threshold = 0.3,
                                    verbose = FALSE) {

  # Step 1: MLE initialization
  if (verbose) cat("Adaptive SMA(", Q, ")_", s, " estimator\n", sep = "")

  mle_fit <- tryCatch({
    stats::arima(x, order = c(0, 0, 0),
                 seasonal = list(order = c(0, 0, Q), period = s),
                 method = "ML",
                 include.mean = include.mean)
  }, error = function(e) {
    warning("MLE failed for SMA: ", e$message)
    return(NULL)
  })

  if (is.null(mle_fit)) {
    stop("SMA estimation failed")
  }

  # Step 2: Analyze residuals
  residuals <- as.numeric(mle_fit$residuals)
  residuals <- residuals[!is.na(residuals)]

  moments <- compute_moments(residuals)
  gamma3 <- moments$m3 / (moments$m2^(3/2))

  if (verbose) {
    cat(sprintf("Residual skewness: γ₃ = %.3f\n", gamma3))
  }

  # Step 3: Decision
  use_pmm <- abs(gamma3) >= asymmetry_threshold

  if (use_pmm) {
    if (verbose) cat("Using PMM2 approach (not yet fully implemented for SMA)\n")
    # TODO: Implement SMA iterative approach
    # For now, fallback to MLE
    warning("PMM2 for pure SMA not yet implemented. Using MLE.")
  }

  # Extract coefficients
  Theta_est <- numeric(Q)
  for (j in 1:Q) {
    coef_name <- paste0("sma", j)
    if (coef_name %in% names(mle_fit$coef)) {
      Theta_est[j] <- mle_fit$coef[coef_name]
    }
  }

  mu_est <- if (include.mean && "intercept" %in% names(mle_fit$coef)) {
    mle_fit$coef["intercept"]
  } else {
    0
  }

  list(
    Theta_est = Theta_est,
    mu_est = mu_est,
    method = ifelse(use_pmm, "MLE (PMM2 TODO)", "MLE"),
    gamma3 = gamma3,
    mle_fit = mle_fit
  )
}
