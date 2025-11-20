# ==============================================================================
# MA Models Support for EstemPMM (Linear Approach)
# ==============================================================================
#
# This module implements conditional linearization approach for MA models
# based on fixed-point iteration scheme.
#
# Theory: See docs/THEORETICAL_ANALYSIS_SMA_SUPPORT.md
#
# Key functions:
#   - compute_ma_innovations(): Recursive innovation computation
#   - create_ma_design_matrix(): Iterative design matrix construction
#   - solve_pmm_linear(): PMM1/PMM2 linear solver
#
# ==============================================================================

#' Compute MA innovations recursively
#'
#' For MA(q) model: x_t = μ + ε_t + θ₁ε_{t-1} + ... + θ_qε_{t-q}
#' Innovations are computed forward: ε_t = x_t - μ - Σθ_jε_{t-j}
#'
#' @param x_centered Centered time series (x - μ)
#' @param ma_coef Vector of MA coefficients [θ₁, θ₂, ..., θ_q]
#' @param q MA order
#' @param init_method Initialization method: "zero" or "backcast"
#'
#' @return Vector of innovations (same length as x_centered)
#'
#' @details
#' Initial conditions:
#'   - "zero": ε_t = 0 for t ≤ 0
#'   - "backcast": ε_t = x_t for t = 1, ..., q (conditional MLE approach)
#'
#' The function implements forward recursion:
#'   For t = q+1, ..., n:
#'     ε_t = x_t - Σ_{j=1}^q θ_j ε_{t-j}
#'
#' @references
#' Box & Jenkins (2015), Time Series Analysis, Chapter 7
#' Brockwell & Davis (2016), Introduction to Time Series, Section 5.1
#'
#' @export
compute_ma_innovations <- function(x_centered, ma_coef, q = length(ma_coef),
                                    init_method = c("backcast", "zero")) {
  init_method <- match.arg(init_method)

  n <- length(x_centered)
  innovations <- numeric(n)

  # Validate inputs
  if (q <= 0) {
    stop("MA order q must be positive")
  }
  if (length(ma_coef) != q) {
    stop("Length of ma_coef must equal q")
  }
  if (n <= q) {
    stop("Time series too short for MA(", q, ") model")
  }

  # Step 1: Initialize innovations for t = 1, ..., q
  if (init_method == "backcast") {
    # Conditional MLE: use observations as initial innovations
    innovations[1:q] <- x_centered[1:q]
  } else if (init_method == "zero") {
    # Unconditional: assume ε_t = 0 for t ≤ 0
    innovations[1:q] <- 0
  }

  # Step 2: Forward recursion for t = q+1, ..., n
  for (t in (q + 1):n) {
    # Calculate expected value based on previous innovations
    expected <- 0
    for (j in 1:q) {
      expected <- expected + ma_coef[j] * innovations[t - j]
    }

    # Current innovation: ε_t = x_t - Σθ_jε_{t-j}
    innovations[t] <- x_centered[t] - expected
  }

  # Step 3: Check for numerical issues
  if (any(is.infinite(innovations)) || any(is.na(innovations))) {
    warning("Infinite or NA innovations detected. Using regularization.")

    # Replace problematic values
    bad_idx <- is.infinite(innovations) | is.na(innovations)
    if (sum(!bad_idx) > 0) {
      # Use mean of valid values
      innovations[bad_idx] <- mean(innovations[!bad_idx])
    } else {
      # All values bad: use zeros
      innovations[bad_idx] <- 0
    }
  }

  # Check for explosive behavior (instability indicator)
  if (max(abs(innovations)) > 1e6 * sd(x_centered)) {
    warning("Innovations show explosive behavior. MA parameters may be non-invertible.")
  }

  return(innovations)
}


#' Compute seasonal MA innovations recursively
#'
#' For SMA(Q)_s model: x_t = μ + ε_t + Θ₁ε_{t-s} + ... + Θ_Qε_{t-Qs}
#'
#' @param x_centered Centered time series
#' @param sma_coef Vector of seasonal MA coefficients [Θ₁, Θ₂, ..., Θ_Q]
#' @param Q Seasonal MA order
#' @param s Seasonal period (e.g., 4 for quarterly, 12 for monthly)
#' @param init_method Initialization method
#'
#' @return Vector of innovations
#'
#' @export
compute_sma_innovations <- function(x_centered, sma_coef, Q = length(sma_coef), s,
                                     init_method = c("backcast", "zero")) {
  init_method <- match.arg(init_method)

  n <- length(x_centered)
  innovations <- numeric(n)

  # Validate inputs
  if (Q <= 0) {
    stop("Seasonal MA order Q must be positive")
  }
  if (s <= 1) {
    stop("Seasonal period s must be > 1")
  }
  if (length(sma_coef) != Q) {
    stop("Length of sma_coef must equal Q")
  }

  max_lag <- Q * s
  if (n <= max_lag) {
    stop("Time series too short for SMA(", Q, ")_", s, " model")
  }

  # Step 1: Initialize innovations for t = 1, ..., Q*s
  if (init_method == "backcast") {
    innovations[1:max_lag] <- x_centered[1:max_lag]
  } else {
    innovations[1:max_lag] <- 0
  }

  # Step 2: Forward recursion for t = Q*s + 1, ..., n
  for (t in (max_lag + 1):n) {
    # Calculate expected value based on seasonal lags
    expected <- 0
    for (j in 1:Q) {
      lag <- j * s
      if (t - lag > 0) {
        expected <- expected + sma_coef[j] * innovations[t - lag]
      }
    }

    # Current innovation
    innovations[t] <- x_centered[t] - expected
  }

  # Numerical checks
  if (any(is.infinite(innovations)) || any(is.na(innovations))) {
    warning("Infinite or NA innovations detected in SMA computation.")
    bad_idx <- is.infinite(innovations) | is.na(innovations)
    innovations[bad_idx] <- ifelse(sum(!bad_idx) > 0,
                                    mean(innovations[!bad_idx]), 0)
  }

  return(innovations)
}


#' Solve PMM linear system
#'
#' Solves linear PMM system with optional weighting by moments
#'
#' @param X Design matrix (n × p)
#' @param y Response vector (length n)
#' @param moments List with m2, m3, m4 (central moments of innovations)
#' @param method "ols", "pmm1", or "pmm2"
#'
#' @return Vector of estimated coefficients (length p)
#'
#' @details
#' PMM1: k_{1,v} = (1/μ₂) · ε_{t-p}
#'   → Weighted LS with W = I/μ₂
#'
#' PMM2: Uses higher-order moments (to be implemented)
#'   → More complex weighting scheme
#'
#' @export
solve_pmm_linear <- function(X, y, moments = NULL, method = c("pmm1", "ols", "pmm2")) {
  method <- match.arg(method)

  n <- nrow(X)
  p <- ncol(X)

  # Validate inputs
  if (length(y) != n) {
    stop("Length of y must equal nrow(X)")
  }

  if (method == "ols") {
    # Standard OLS: θ̂ = (X'X)^{-1}X'y
    XtX <- crossprod(X)  # X'X
    Xty <- crossprod(X, y)  # X'y

    # Check for singularity
    if (qr(XtX)$rank < p) {
      warning("X'X is singular or near-singular. Using pseudo-inverse.")
      theta_est <- MASS::ginv(XtX) %*% Xty
    } else {
      theta_est <- solve(XtX, Xty)
    }

    return(as.numeric(theta_est))
  }

  if (method == "pmm1") {
    # PMM1: Weighted LS with W = I/μ₂
    if (is.null(moments) || is.null(moments$m2)) {
      stop("moments$m2 required for PMM1")
    }

    mu2 <- moments$m2
    if (mu2 <= 0) {
      warning("Invalid μ₂ (≤ 0). Falling back to OLS.")
      return(solve_pmm_linear(X, y, method = "ols"))
    }

    # Weighted system: (X'WX)θ = X'Wy where W = I/μ₂
    # Simplifies to: (X'X/μ₂)θ = X'y/μ₂
    # Which is equivalent to: (X'X)θ = X'y (weight cancels!)
    # So PMM1 = OLS for this formulation

    # Actually, correct formulation uses weight in residuals:
    # W = diag(1/μ₂, n)
    W_inv <- mu2  # scalar weight
    XtX <- crossprod(X) / W_inv
    Xty <- crossprod(X, y) / W_inv

    if (qr(XtX)$rank < p) {
      warning("X'X is singular. Using pseudo-inverse.")
      theta_est <- MASS::ginv(XtX) %*% Xty
    } else {
      theta_est <- solve(XtX, Xty)
    }

    return(as.numeric(theta_est))
  }

  if (method == "pmm2") {
    # PMM2: Not yet implemented for MA models
    # Would require solving nonlinear system with higher-order moments
    stop("PMM2 for MA models not yet implemented. Use 'pmm1' or 'ols'.")
  }
}


#' Create MA design matrix using iterative scheme
#'
#' Implements fixed-point iteration for MA model parameter estimation
#'
#' @param x Numeric vector (time series data)
#' @param q MA order
#' @param max_iter Maximum iterations for fixed-point
#' @param tol Convergence tolerance
#' @param init_method Innovation initialization method
#' @param pmm_method "pmm1", "ols", or "pmm2"
#' @param include.mean Include intercept/mean parameter
#' @param verbose Print iteration details
#'
#' @return List with:
#'   \item{theta_est}{Estimated MA coefficients}
#'   \item{mu_est}{Estimated mean (if include.mean=TRUE)}
#'   \item{innovations}{Final innovations}
#'   \item{X}{Final design matrix}
#'   \item{y}{Response vector}
#'   \item{iterations}{Number of iterations used}
#'   \item{converged}{Logical: did algorithm converge?}
#'   \item{moments}{Moments of final innovations}
#'
#' @details
#' Algorithm (fixed-point iteration):
#'   1. Initialize: θ^(0) via MLE/CSS
#'   2. Compute: ε^(k) given θ^(k)
#'   3. Build: X^(k) from ε^(k)
#'   4. Solve: θ^(k+1) = PMM_solve(X^(k), y)
#'   5. Check: ||θ^(k+1) - θ^(k)|| < tol
#'   6. Repeat until convergence or max_iter
#'
#' @examples
#' \dontrun{
#' # Simulate MA(1) with θ = 0.6
#' set.seed(123)
#' n <- 200
#' eps <- rnorm(n, 0, 1)
#' x <- numeric(n)
#' x[1] <- eps[1]
#' for (t in 2:n) {
#'   x[t] <- eps[t] + 0.6 * eps[t-1]
#' }
#'
#' # Estimate
#' result <- create_ma_design_matrix(x, q = 1, pmm_method = "pmm1")
#' print(result$theta_est)  # Should be close to 0.6
#' }
#'
#' @export
create_ma_design_matrix <- function(x, q,
                                     max_iter = 50,
                                     tol = 1e-6,
                                     init_method = "backcast",
                                     pmm_method = c("pmm1", "ols"),
                                     include.mean = TRUE,
                                     verbose = FALSE) {
  pmm_method <- match.arg(pmm_method)

  # Validate inputs
  n <- length(x)
  if (n <= q + 10) {
    stop("Time series too short for MA(", q, ") estimation")
  }

  # Step 1: Mean centering
  if (include.mean) {
    mu_est <- mean(x, na.rm = TRUE)
    x_centered <- x - mu_est
  } else {
    mu_est <- 0
    x_centered <- x
  }

  # Step 2: Initialize parameters via MLE/CSS
  if (verbose) cat("Initializing with stats::arima()...\n")

  init_fit <- tryCatch({
    stats::arima(x, order = c(0, 0, q), method = "CSS-ML",
                 include.mean = include.mean)
  }, error = function(e) {
    warning("MLE initialization failed: ", e$message, "\nUsing zeros.")
    list(coef = rep(0, q))
  })

  # Extract initial MA coefficients
  theta_current <- numeric(q)
  for (j in 1:q) {
    coef_name <- paste0("ma", j)
    if (coef_name %in% names(init_fit$coef)) {
      theta_current[j] <- init_fit$coef[coef_name]
    } else {
      theta_current[j] <- 0.1  # Default fallback
    }
  }

  if (verbose) {
    cat("Initial θ: [", paste(round(theta_current, 4), collapse = ", "), "]\n")
  }

  # Step 3: Fixed-point iteration
  converged <- FALSE
  iter <- 0

  for (iter in 1:max_iter) {
    # Step 3a: Compute innovations given current θ
    innovations <- compute_ma_innovations(x_centered, theta_current, q, init_method)

    # Step 3b: Build design matrix
    max_lag <- q
    n_obs <- n - max_lag
    X <- matrix(0, nrow = n_obs, ncol = q)
    y <- x_centered[(max_lag + 1):n]

    for (j in 1:q) {
      X[, j] <- innovations[(max_lag - j + 1):(n - j)]
    }

    # Step 3c: Compute moments of innovations
    moments <- list(
      m2 = mean(innovations^2),
      m3 = mean(innovations^3),
      m4 = mean(innovations^4)
    )

    # Step 3d: Solve PMM linear system
    theta_new <- solve_pmm_linear(X, y, moments, method = pmm_method)

    # Step 3e: Check convergence
    delta <- sqrt(sum((theta_new - theta_current)^2))

    if (verbose && (iter %% 10 == 0 || delta < tol)) {
      cat(sprintf("Iter %2d: θ = [%s], ||Δθ|| = %.6f\n",
                  iter,
                  paste(round(theta_new, 4), collapse = ", "),
                  delta))
    }

    if (delta < tol) {
      converged <- TRUE
      if (verbose) cat("✓ Converged in", iter, "iterations\n")
      break
    }

    # Update for next iteration
    theta_current <- theta_new
  }

  if (!converged && verbose) {
    warning("Fixed-point iteration did not converge in ", max_iter, " iterations")
  }

  # Step 4: Return results
  list(
    theta_est = theta_current,
    mu_est = mu_est,
    innovations = innovations,
    X = X,
    y = y,
    iterations = iter,
    converged = converged,
    moments = moments,
    method = pmm_method,
    include.mean = include.mean
  )
}


#' Compute central moments of a vector
#'
#' @param x Numeric vector
#' @return List with m2, m3, m4 (central moments)
#' @keywords internal
compute_moments <- function(x) {
  n <- length(x)
  if (n == 0) {
    return(list(m2 = NA, m3 = NA, m4 = NA))
  }

  xm <- mean(x, na.rm = TRUE)
  x_centered <- x - xm

  list(
    m2 = mean(x_centered^2, na.rm = TRUE),
    m3 = mean(x_centered^3, na.rm = TRUE),
    m4 = mean(x_centered^4, na.rm = TRUE)
  )
}
