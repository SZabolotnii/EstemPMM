# SAR-PMM2 Концептуальний прототип
# Автор: Теоретичний аналіз 2025-11-13
# Мета: Демонстрація можливої реалізації SAR моделей з PMM2

# =============================================================================
# 1. ФУНКЦІЯ ДЛЯ СТВОРЕННЯ DESIGN MATRIX ДЛЯ SAR
# =============================================================================

#' Create design matrix for seasonal AR model
#'
#' @param x Numeric vector (centered time series)
#' @param p Non-seasonal AR order
#' @param P Seasonal AR order
#' @param s Seasonal period (e.g., 12 for monthly data)
#' @param multiplicative Logical, include cross-terms (default FALSE)
#'
#' @return Design matrix with lagged values
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Simple SAR(1)_12 model
#' X <- create_sar_matrix(x, p=0, P=1, s=12)
#'
#' # AR(1) + SAR(1)_12 model
#' X <- create_sar_matrix(x, p=1, P=1, s=12)
#'
#' # With multiplicative terms
#' X <- create_sar_matrix(x, p=1, P=1, s=12, multiplicative=TRUE)
#' }
create_sar_matrix <- function(x, p = 0, P = 1, s = 12, multiplicative = FALSE) {
  n <- length(x)

  # Validate inputs
  if (p < 0 || P < 0) {
    stop("AR orders (p, P) must be non-negative")
  }
  if (p == 0 && P == 0) {
    stop("At least one of p or P must be positive")
  }
  if (s <= 1) {
    stop("Seasonal period s must be greater than 1")
  }

  # Calculate maximum lag
  max_nonseasonal_lag <- p
  max_seasonal_lag <- P * s
  max_lag <- max(max_nonseasonal_lag, max_seasonal_lag)

  # Check data sufficiency
  if (n <= max_lag) {
    stop("Insufficient data: need at least ", max_lag + 1, " observations")
  }

  # Number of rows in design matrix
  nr <- n - max_lag

  # Number of columns depends on multiplicative option
  if (multiplicative) {
    ncol <- p + P + (p * P)  # Non-seasonal + Seasonal + Cross-terms
  } else {
    ncol <- p + P  # Additive model
  }

  # Initialize design matrix
  X <- matrix(0, nrow = nr, ncol = ncol)
  col_idx <- 1

  # Add non-seasonal AR lags (if p > 0)
  if (p > 0) {
    for (i in seq_len(p)) {
      X[, col_idx] <- x[(max_lag - i + 1):(n - i)]
      col_idx <- col_idx + 1
    }
  }

  # Add seasonal AR lags (if P > 0)
  if (P > 0) {
    for (j in seq_len(P)) {
      seasonal_lag <- j * s
      X[, col_idx] <- x[(max_lag - seasonal_lag + 1):(n - seasonal_lag)]
      col_idx <- col_idx + 1
    }
  }

  # Add multiplicative cross-terms (if requested)
  if (multiplicative && p > 0 && P > 0) {
    for (i in seq_len(p)) {
      for (j in seq_len(P)) {
        cross_lag <- i + j * s
        X[, col_idx] <- x[(max_lag - cross_lag + 1):(n - cross_lag)]
        col_idx <- col_idx + 1
      }
    }
  }

  # Set column names
  col_names <- character(ncol)
  name_idx <- 1

  if (p > 0) {
    for (i in seq_len(p)) {
      col_names[name_idx] <- paste0("ar", i)
      name_idx <- name_idx + 1
    }
  }

  if (P > 0) {
    for (j in seq_len(P)) {
      col_names[name_idx] <- paste0("sar", j)
      name_idx <- name_idx + 1
    }
  }

  if (multiplicative && p > 0 && P > 0) {
    for (i in seq_len(p)) {
      for (j in seq_len(P)) {
        col_names[name_idx] <- paste0("ar", i, "_sar", j)
        name_idx <- name_idx + 1
      }
    }
  }

  colnames(X) <- col_names
  return(X)
}


# =============================================================================
# 2. ГОЛОВНА ФУНКЦІЯ SAR_PMM2
# =============================================================================

#' Fit Seasonal AR model using PMM2 method
#'
#' @param x Numeric vector of time series data
#' @param order Vector of length 2: c(p, P) where p is non-seasonal AR order
#'   and P is seasonal AR order
#' @param season List with seasonal specification: list(period = s)
#' @param method Estimation method: "pmm2" (default), "ols", "css"
#' @param include.mean Logical, include intercept term (default TRUE)
#' @param multiplicative Logical, use multiplicative form with cross-terms (default FALSE)
#' @param max_iter Maximum iterations for PMM2 algorithm (default 50)
#' @param tol Convergence tolerance (default 1e-6)
#' @param regularize Logical, use regularization (default TRUE)
#' @param reg_lambda Regularization parameter (default 1e-8)
#' @param verbose Logical, print progress (default FALSE)
#'
#' @return S4 object of class SARPMM2 with fitted model
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate synthetic seasonal data
#' n <- 120
#' s <- 12
#' y <- arima.sim(n = n, list(ar = 0.7, sar = 0.5, seasonal = list(period = s)))
#'
#' # Fit SAR model
#' fit <- sar_pmm2(y, order = c(1, 1), season = list(period = 12))
#' summary(fit)
#'
#' # Simple seasonal model (no non-seasonal component)
#' fit_pure_sar <- sar_pmm2(y, order = c(0, 1), season = list(period = 12))
#' }
sar_pmm2 <- function(x,
                     order = c(0, 1),
                     season = list(period = 12),
                     method = "pmm2",
                     include.mean = TRUE,
                     multiplicative = FALSE,
                     max_iter = 50,
                     tol = 1e-6,
                     regularize = TRUE,
                     reg_lambda = 1e-8,
                     verbose = FALSE) {

  # Store original call
  cl <- match.call()

  # Validate and prepare data
  x <- as.numeric(x)
  if (any(is.na(x)) || any(is.infinite(x))) {
    stop("Time series contains NA or infinite values")
  }

  # Parse order specification
  if (length(order) != 2) {
    stop("'order' must be a vector of length 2: c(p, P)")
  }
  p <- as.integer(order[1])  # Non-seasonal AR order
  P <- as.integer(order[2])  # Seasonal AR order

  # Parse seasonal specification
  if (!is.list(season) || is.null(season$period)) {
    stop("'season' must be a list with 'period' element")
  }
  s <- as.integer(season$period)

  if (verbose) {
    cat("Fitting SAR(", p, ",", P, ")_", s, " model\n", sep = "")
    cat("Method:", method, "\n")
    cat("Multiplicative:", multiplicative, "\n")
  }

  # Center data (handle mean)
  orig_x <- x
  if (include.mean) {
    x_mean <- mean(x)
    x_centered <- x - x_mean
  } else {
    x_mean <- 0
    x_centered <- x
  }

  # Create design matrix
  X <- create_sar_matrix(x_centered, p, P, s, multiplicative)
  max_lag <- max(p, P * s)
  y <- x_centered[(max_lag + 1):length(x_centered)]

  if (verbose) {
    cat("Design matrix dimensions:", nrow(X), "x", ncol(X), "\n")
    cat("Effective sample size:", length(y), "\n")
  }

  # Check for sufficient data
  if (nrow(X) < ncol(X) + 5) {
    warning("Very small sample size relative to number of parameters")
  }

  # Step 1: Initial estimation (OLS/CSS)
  if (method %in% c("pmm2", "css")) {
    # Use OLS for initial estimates
    XtX <- t(X) %*% X
    Xty <- t(X) %*% y

    # Add small regularization if needed
    if (det(XtX) < 1e-10) {
      if (verbose) cat("Adding regularization to X'X for initial fit\n")
      diag(XtX) <- diag(XtX) + 1e-6
    }

    b_init <- solve(XtX, Xty)

    if (verbose) {
      cat("Initial coefficients (OLS):\n")
      print(b_init)
    }
  } else if (method == "ols") {
    # Return OLS estimates directly
    b_init <- solve(t(X) %*% X, t(X) %*% y)
    residuals <- y - X %*% b_init
    moments <- compute_moments(residuals)

    # Return SARPMM2 object with OLS estimates
    return(new("SARPMM2",
               coefficients = as.numeric(b_init),
               residuals = as.numeric(residuals),
               m2 = as.numeric(moments$m2),
               m3 = as.numeric(moments$m3),
               m4 = as.numeric(moments$m4),
               convergence = TRUE,
               iterations = 0L,
               call = cl,
               model_type = "sar",
               intercept = if (include.mean) as.numeric(x_mean) else 0,
               original_series = as.numeric(orig_x),
               order = list(ar = p, sar = P, period = s)))
  }

  # Step 2: Compute moments from initial residuals
  residuals_init <- y - X %*% b_init
  moments <- compute_moments(residuals_init)

  if (verbose) {
    cat("\nMoments from initial residuals:\n")
    cat("  m2 (variance):", moments$m2, "\n")
    cat("  m3 (skewness indicator):", moments$m3, "\n")
    cat("  m4 (kurtosis indicator):", moments$m4, "\n")
    cat("  c3 (skewness coef):", moments$c3, "\n")
    cat("  c4 (excess kurtosis):", moments$c4, "\n")
    cat("  g (variance reduction):", moments$g, "\n")
  }

  # Step 3: Apply PMM2 algorithm
  if (method == "pmm2") {
    if (verbose) cat("\nApplying PMM2 refinement...\n")

    # Use existing pmm2_algorithm function (no changes needed!)
    pmm2_result <- pmm2_algorithm(
      b_init = b_init,
      X = X,
      y = y,
      m2 = moments$m2,
      m3 = moments$m3,
      m4 = moments$m4,
      max_iter = max_iter,
      tol = tol,
      regularize = regularize,
      reg_lambda = reg_lambda,
      verbose = verbose
    )

    b_final <- pmm2_result$coefficients
    converged <- pmm2_result$convergence
    iterations <- pmm2_result$iterations

    # Compute final residuals
    residuals_final <- y - X %*% b_final
    moments_final <- compute_moments(residuals_final)

    if (verbose) {
      cat("\nPMM2 converged:", converged, "\n")
      cat("Iterations:", iterations, "\n")
      cat("\nFinal coefficients:\n")
      print(b_final)
      cat("\nFinal moments:\n")
      cat("  m2:", moments_final$m2, "\n")
      cat("  c3:", moments_final$c3, "\n")
      cat("  g:", moments_final$g, "\n")
    }

    # Create SARPMM2 object
    result <- new("SARPMM2",
                  coefficients = as.numeric(b_final),
                  residuals = as.numeric(residuals_final),
                  m2 = as.numeric(moments_final$m2),
                  m3 = as.numeric(moments_final$m3),
                  m4 = as.numeric(moments_final$m4),
                  convergence = converged,
                  iterations = as.integer(iterations),
                  call = cl,
                  model_type = "sar",
                  intercept = if (include.mean) as.numeric(x_mean) else 0,
                  original_series = as.numeric(orig_x),
                  order = list(ar = p, sar = P, period = s))

    return(result)

  } else {
    # CSS method (return initial estimates)
    result <- new("SARPMM2",
                  coefficients = as.numeric(b_init),
                  residuals = as.numeric(residuals_init),
                  m2 = as.numeric(moments$m2),
                  m3 = as.numeric(moments$m3),
                  m4 = as.numeric(moments$m4),
                  convergence = TRUE,
                  iterations = 0L,
                  call = cl,
                  model_type = "sar",
                  intercept = if (include.mean) as.numeric(x_mean) else 0,
                  original_series = as.numeric(orig_x),
                  order = list(ar = p, sar = P, period = s))

    return(result)
  }
}


# =============================================================================
# 3. S4 КЛАС ДЛЯ SAR МОДЕЛЕЙ
# =============================================================================

#' S4 class for SAR model results with PMM2
#'
#' @slot coefficients Numeric vector of estimated parameters
#' @slot residuals Numeric vector of residuals/innovations
#' @slot m2 Second central moment
#' @slot m3 Third central moment
#' @slot m4 Fourth central moment
#' @slot convergence Logical, whether algorithm converged
#' @slot iterations Integer, number of iterations
#' @slot call Original function call
#' @slot model_type Character, model type identifier
#' @slot intercept Numeric, intercept/mean term
#' @slot original_series Numeric vector, original time series
#' @slot order List with model order: list(ar, sar, period)
#'
#' @exportClass SARPMM2
setClass("SARPMM2",
         slots = c(
           coefficients = "numeric",
           residuals = "numeric",
           m2 = "numeric",
           m3 = "numeric",
           m4 = "numeric",
           convergence = "logical",
           iterations = "numeric",
           call = "call",
           model_type = "character",
           intercept = "numeric",
           original_series = "numeric",
           order = "list"
         ),
         contains = "TS2fit")  # Inherits from existing TS2fit class


# =============================================================================
# 4. МЕТОДИ ДЛЯ SARPMM2 КЛАСУ
# =============================================================================

#' Extract coefficients from SARPMM2 object
#' @param object SARPMM2 object
#' @export
setMethod("coef", "SARPMM2",
          function(object) {
            coefs <- object@coefficients
            p <- object@order$ar
            P <- object@order$sar

            # Create proper names
            names_vec <- character(length(coefs))
            idx <- 1

            if (p > 0) {
              for (i in seq_len(p)) {
                names_vec[idx] <- paste0("ar", i)
                idx <- idx + 1
              }
            }

            if (P > 0) {
              for (j in seq_len(P)) {
                names_vec[idx] <- paste0("sar", j)
                idx <- idx + 1
              }
            }

            # If there are more coefficients (multiplicative terms)
            if (idx <= length(coefs)) {
              for (k in idx:length(coefs)) {
                names_vec[k] <- paste0("coef", k)
              }
            }

            names(coefs) <- names_vec
            return(coefs)
          })


#' Summary method for SARPMM2 objects
#' @param object SARPMM2 object
#' @export
setMethod("summary", "SARPMM2",
          function(object) {
            cat("\nSeasonal AR Model fitted with PMM2\n")
            cat("===================================\n\n")

            cat("Call:\n")
            print(object@call)
            cat("\n")

            cat("Model: SAR(", object@order$ar, ",", object@order$sar, ")_",
                object@order$period, "\n", sep = "")
            cat("Observations:", length(object@original_series), "\n")
            cat("Effective sample size:", length(object@residuals), "\n\n")

            cat("Coefficients:\n")
            coefs <- coef(object)
            print(coefs)

            if (object@intercept != 0) {
              cat("\nIntercept:", object@intercept, "\n")
            }

            cat("\nResidual moments:\n")
            cat("  m2 (variance):", object@m2, "\n")
            cat("  m3 (skewness):", object@m3, "\n")
            cat("  m4 (kurtosis):", object@m4, "\n")

            # Calculate distribution characteristics
            c3 <- object@m3 / (object@m2^(3/2))
            c4 <- object@m4 / (object@m2^2) - 3
            g <- 1 - c3^2 / (2 + c4)

            cat("\nDistribution characteristics:\n")
            cat("  Skewness coefficient (c3):", c3, "\n")
            cat("  Excess kurtosis (c4):", c4, "\n")
            cat("  Variance reduction factor (g):", g, "\n")

            if (g < 1) {
              reduction_pct <- (1 - g) * 100
              cat("  => Expected", round(reduction_pct, 1),
                  "% variance reduction vs OLS\n")
            }

            cat("\nAlgorithm:\n")
            cat("  Converged:", object@convergence, "\n")
            cat("  Iterations:", object@iterations, "\n")

            # Residual diagnostics
            cat("\nResidual statistics:\n")
            cat("  Min:", min(object@residuals), "\n")
            cat("  Q1:", quantile(object@residuals, 0.25), "\n")
            cat("  Median:", median(object@residuals), "\n")
            cat("  Q3:", quantile(object@residuals, 0.75), "\n")
            cat("  Max:", max(object@residuals), "\n")

            invisible(object)
          })


# =============================================================================
# 5. ФУНКЦІЯ ПОРІВНЯННЯ МЕТОДІВ
# =============================================================================

#' Compare SAR model estimation methods
#'
#' @param x Time series data
#' @param order Model order c(p, P)
#' @param period Seasonal period
#'
#' @return Data frame with comparison results
#' @export
compare_sar_methods <- function(x, order = c(1, 1), period = 12) {
  cat("Comparing SAR estimation methods\n")
  cat("================================\n\n")

  p <- order[1]
  P <- order[2]

  # Fit with different methods
  fit_ols <- sar_pmm2(x, order = order, season = list(period = period),
                      method = "ols", verbose = FALSE)

  fit_pmm2 <- sar_pmm2(x, order = order, season = list(period = period),
                       method = "pmm2", verbose = FALSE)

  # Also try standard arima for comparison
  fit_arima <- tryCatch({
    arima(x, order = c(p, 0, 0),
          seasonal = list(order = c(P, 0, 0), period = period),
          method = "CSS")
  }, error = function(e) NULL)

  # Extract coefficients
  coef_ols <- coef(fit_ols)
  coef_pmm2 <- coef(fit_pmm2)
  coef_arima <- if (!is.null(fit_arima)) coef(fit_arima) else rep(NA, length(coef_ols))

  # Compute residual standard errors
  sigma_ols <- sqrt(fit_ols@m2)
  sigma_pmm2 <- sqrt(fit_pmm2@m2)
  sigma_arima <- if (!is.null(fit_arima)) sqrt(fit_arima$sigma2) else NA

  # Create comparison table
  result <- data.frame(
    Parameter = names(coef_ols),
    OLS = coef_ols,
    PMM2 = coef_pmm2,
    ARIMA_CSS = coef_arima[1:length(coef_ols)]
  )

  cat("Coefficient estimates:\n")
  print(result, digits = 4)

  cat("\nResidual standard errors:\n")
  cat("  OLS:", sigma_ols, "\n")
  cat("  PMM2:", sigma_pmm2, "\n")
  cat("  ARIMA (CSS):", sigma_arima, "\n")

  # Distribution characteristics
  cat("\nResidual distribution:\n")
  cat("  OLS - Skewness:", fit_ols@m3 / (fit_ols@m2^(3/2)), "\n")
  cat("  PMM2 - Skewness:", fit_pmm2@m3 / (fit_pmm2@m2^(3/2)), "\n")

  invisible(result)
}


# =============================================================================
# 6. ПРИКЛАД ВИКОРИСТАННЯ
# =============================================================================

#' Example: Fit SAR model to synthetic data
#'
#' @export
example_sar_pmm2 <- function() {
  cat("SAR-PMM2 Example with Synthetic Data\n")
  cat("====================================\n\n")

  # Generate synthetic data with seasonal pattern
  set.seed(42)
  n <- 120  # 10 years of monthly data
  s <- 12   # Annual seasonality

  # True parameters
  true_ar <- 0.7     # Non-seasonal AR(1)
  true_sar <- 0.5    # Seasonal AR(1)
  true_mean <- 5

  # Generate seasonal AR process with asymmetric errors
  y <- numeric(n)
  y[1:s] <- rnorm(s, mean = true_mean, sd = 1)

  for (t in (s+1):n) {
    # AR(1) + SAR(1) structure
    ar_component <- true_ar * (y[t-1] - true_mean)
    sar_component <- true_sar * (y[t-s] - true_mean)

    # Asymmetric innovations (gamma distribution)
    innovation <- rgamma(1, shape = 2, scale = 0.5) - 1  # Mean ≈ 0

    y[t] <- true_mean + ar_component + sar_component + innovation
  }

  cat("Generated time series:\n")
  cat("  Length:", n, "\n")
  cat("  Mean:", mean(y), "\n")
  cat("  SD:", sd(y), "\n\n")

  # Plot the series
  plot(y, type = "l", main = "Synthetic SAR(1,1)_12 Series",
       xlab = "Time", ylab = "Value")
  abline(h = mean(y), col = "red", lty = 2)

  # Fit models
  cat("\nFitting SAR(1,1)_12 model with different methods...\n\n")

  fit_ols <- sar_pmm2(y, order = c(1, 1), season = list(period = 12),
                      method = "ols")
  fit_pmm2 <- sar_pmm2(y, order = c(1, 1), season = list(period = 12),
                       method = "pmm2")

  # Compare results
  cat("\n\n=== OLS Results ===\n")
  summary(fit_ols)

  cat("\n\n=== PMM2 Results ===\n")
  summary(fit_pmm2)

  # Compare with true values
  cat("\n\n=== Comparison with True Values ===\n")
  true_coefs <- c(ar1 = true_ar, sar1 = true_sar)
  estimated_ols <- coef(fit_ols)
  estimated_pmm2 <- coef(fit_pmm2)

  comparison <- data.frame(
    Parameter = names(true_coefs),
    True = true_coefs,
    OLS = estimated_ols,
    PMM2 = estimated_pmm2,
    OLS_Error = abs(estimated_ols - true_coefs),
    PMM2_Error = abs(estimated_pmm2 - true_coefs)
  )

  print(comparison, digits = 4)

  cat("\n")
  if (mean(comparison$PMM2_Error) < mean(comparison$OLS_Error)) {
    cat("=> PMM2 has lower average estimation error!\n")
  } else {
    cat("=> OLS has lower average estimation error in this case.\n")
  }

  invisible(list(data = y, fit_ols = fit_ols, fit_pmm2 = fit_pmm2))
}


# =============================================================================
# 7. ЗАПУСК ПРИКЛАДУ
# =============================================================================

# Uncomment to run the example:
# example_sar_pmm2()

cat("\nSAR-PMM2 prototype loaded successfully!\n")
cat("Run example_sar_pmm2() to see a demonstration.\n")
