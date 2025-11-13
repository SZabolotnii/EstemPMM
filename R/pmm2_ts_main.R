# pmm2_ts_main.R - Unified module for time series models with PMM2

#' Fit a time series model using the PMM2 method
#'
#' @param x Numeric vector of time series data
#' @param order Model order specification:
#'        - For AR models: a single integer (AR order)
#'        - For MA models: a single integer (MA order)
#'        - For ARMA models: vector c(p, q) (AR and MA orders)
#'        - For ARIMA models: vector c(p, d, q) (AR, differencing, and MA orders)
#' @param model_type String specifying the model type: "ar", "ma", "arma", or "arima"
#' @param method String: estimation method, one of "pmm2" (default), "css", "ml", "yw", "ols"
#' @param max_iter Integer: maximum number of iterations for the algorithm
#' @param tol Numeric: tolerance for convergence
#' @param include.mean Logical: whether to include a mean (intercept) term
#' @param initial List or vector of initial parameter estimates (optional)
#' @param na.action Function for handling missing values, default is na.fail
#' @param regularize Logical, add small values to diagonal for numerical stability
#' @param reg_lambda Regularization parameter (if regularize=TRUE)
#' @param verbose Logical: whether to print progress information
#'
#' @details
#' The PMM2 algorithm works as follows:
#'
#' 1. Fits an initial model using a standard method (OLS, Yule-Walker, CSS or ML)
#' 2. Computes central moments (m2, m3, m4) from initial residuals/innovations
#' 3. Uses these moments with a specialized solver (pmm2_algorithm) to find
#'    robust parameter estimates
#'
#' @return An S4 object \code{TS2fit} of the corresponding subclass
#' @export
ts_pmm2 <- function(x, order,
                    model_type = c("ar", "ma", "arma", "arima"),
                    method      = "pmm2",
                    max_iter    = 50,
                    tol         = 1e-6,
                    include.mean= TRUE,
                    initial     = NULL,
                    na.action   = na.fail,
                    regularize  = TRUE,
                    reg_lambda  = 1e-8,
                    verbose     = FALSE) {

  model_type <- match.arg(model_type)
  cl <- match.call()

  if (!is.null(na.action)) {
    x <- na.action(x)
  }

  # 1) Validate input data
  model_params <- validate_ts_parameters(x, order, model_type, include.mean)

  if (model_params$model_type == "ma") {
    q <- model_params$ma_order
    css_fit <- ma_css_fit(model_params$original_x, q, include.mean, verbose)

    if (method == "css") {
      moments <- compute_moments(css_fit$residuals)
      return(new("MAPMM2",
                 coefficients    = as.numeric(css_fit$coefficients),
                 residuals       = as.numeric(css_fit$residuals),
                 m2              = as.numeric(moments$m2),
                 m3              = as.numeric(moments$m3),
                 m4              = as.numeric(moments$m4),
                 convergence     = css_fit$convergence,
                 iterations      = as.numeric(css_fit$iterations),
                 call            = cl,
                 model_type      = "ma",
                 intercept       = if (include.mean) as.numeric(css_fit$intercept) else 0,
                 original_series = as.numeric(model_params$original_x),
                 order           = list(ar = 0L, ma = q, d = 0L)))
    }

    pmm2_fit <- ma_pmm2_fit(model_params$original_x, q, css_fit,
                             max_iter = max_iter, tol = tol,
                             verbose = verbose)
    moments <- compute_moments(pmm2_fit$innovations)

    return(new("MAPMM2",
               coefficients    = as.numeric(pmm2_fit$coefficients),
               residuals       = as.numeric(pmm2_fit$innovations),
               m2              = as.numeric(moments$m2),
               m3              = as.numeric(moments$m3),
               m4              = as.numeric(moments$m4),
               convergence     = pmm2_fit$convergence,
               iterations      = as.numeric(pmm2_fit$iterations),
               call            = cl,
               model_type      = "ma",
               intercept       = if (include.mean) as.numeric(pmm2_fit$intercept) else 0,
               original_series = as.numeric(model_params$original_x),
               order           = list(ar = 0L, ma = q, d = 0L)))
  }

  if (model_params$model_type == "arima") {
    p <- model_params$ar_order
    d <- model_params$d
    q <- model_params$ma_order

    css_fit <- tryCatch(
      stats::arima(model_params$original_x,
                   order = c(p, d, q),
                   method = "CSS-ML",
                   include.mean = include.mean),
      error = function(e) NULL
    )

    if (is.null(css_fit)) {
      stop("Failed to estimate ARIMA model using classical method")
    }

    coef_names <- names(css_fit$coef)
    ar_css <- if (p > 0) as.numeric(css_fit$coef[paste0("ar", seq_len(p))]) else numeric(0)
    ma_css <- if (q > 0) as.numeric(css_fit$coef[paste0("ma", seq_len(q))]) else numeric(0)

    intercept_css <- 0
    if (include.mean && d == 0) {
      intercept_name <- setdiff(coef_names,
                                c(paste0("ar", seq_len(p)), paste0("ma", seq_len(q))))
      if (length(intercept_name) > 0) {
        intercept_css <- as.numeric(css_fit$coef[intercept_name[1]])
      }
    }

    residuals_css <- as.numeric(css_fit$residuals)
    residuals_css[is.na(residuals_css)] <- 0

    x_diff <- if (d > 0) diff(model_params$original_x, differences = d)
              else model_params$original_x
    res_diff <- tail(residuals_css, length(x_diff))

    include_intercept_diff <- include.mean && d == 0

    if (method == "css") {
      res_clean <- res_diff[is.finite(res_diff)]
      moments <- compute_moments(res_clean)
      return(new("ARIMAPMM2",
                 coefficients    = as.numeric(c(ar_css, ma_css)),
                 residuals       = as.numeric(residuals_css),
                 m2              = as.numeric(moments$m2),
                 m3              = as.numeric(moments$m3),
                 m4              = as.numeric(moments$m4),
                 convergence     = TRUE,
                 iterations      = 1L,
                 call            = cl,
                 model_type      = "arima",
                 intercept       = if (include_intercept_diff) intercept_css else 0,
                 original_series = as.numeric(model_params$original_x),
                 order           = list(ar = p, ma = q, d = d)))
    }

    design <- arma_build_design(x_diff, res_diff,
                                 p = p, q = q,
                                 intercept = intercept_css,
                                 include_intercept = include_intercept_diff)

    moments <- compute_moments(res_diff[is.finite(res_diff)])
    b_init <- c(if (include_intercept_diff) 0 else NULL, ar_css, ma_css)

    algo_res <- pmm2_algorithm(
      b_init = b_init,
      X = design$X,
      y = design$y,
      m2 = moments$m2,
      m3 = moments$m3,
      m4 = moments$m4,
      max_iter = max_iter,
      tol = tol,
      regularize = regularize,
      reg_lambda = reg_lambda,
      verbose = verbose
    )

    if (include_intercept_diff) {
      intercept_hat <- algo_res$b[1]
      ar_hat <- if (p > 0) algo_res$b[1 + seq_len(p)] else numeric(0)
      ma_hat <- if (q > 0) algo_res$b[1 + p + seq_len(q)] else numeric(0)
    } else {
      intercept_hat <- 0
      ar_hat <- if (p > 0) algo_res$b[seq_len(p)] else numeric(0)
      ma_hat <- if (q > 0) algo_res$b[p + seq_len(q)] else numeric(0)
    }

    x_mean <- if (include_intercept_diff) intercept_hat else 0
    x_centered <- x_diff - x_mean

    model_info <- list(
      ar_order = p,
      ma_order = q,
      d = d,
      model_type = "arima",
      include.mean = include.mean,
      innovations = res_diff,
      x = x_centered,
      x_mean = x_mean,
      original_x = model_params$original_x,
      verbose = verbose
    )

    final_coef <- c(ar_hat, ma_hat)
    final_res <- compute_ts_residuals(final_coef, model_info)
    res_clean <- final_res[is.finite(final_res)]
    moments_final <- compute_moments(res_clean)

    return(new("ARIMAPMM2",
               coefficients    = as.numeric(final_coef),
               residuals       = as.numeric(final_res),
               m2              = as.numeric(moments_final$m2),
               m3              = as.numeric(moments_final$m3),
               m4              = as.numeric(moments_final$m4),
               convergence     = algo_res$convergence,
               iterations      = as.numeric(algo_res$iterations),
               call            = cl,
               model_type      = "arima",
               intercept       = if (include_intercept_diff) as.numeric(intercept_hat) else 0,
               original_series = as.numeric(model_params$original_x),
               order           = list(ar = p, ma = q, d = d)))
  }

  # 2) Obtain initial estimates
  init <- get_initial_estimates(model_params, initial, method, verbose)
  b_init      <- init$b_init
  x_mean      <- init$x_mean
  innovations <- init$innovations
  x_centered  <- init$x_centered
  orig_x      <- init$orig_x
  m2          <- init$m2
  m3          <- init$m3
  m4          <- init$m4

  # 3) Create design matrix
  dm <- create_ts_design_matrix(
    x = x_centered,
    model_info = list(
      ar_order = model_params$ar_order,
      ma_order = model_params$ma_order,
      d = model_params$d,
      model_type = model_params$model_type,
      include.mean = model_params$include.mean
    ),
    innovations = innovations
  )

  # 4) If method == "pmm2", use the PMM2 algorithm
  if (method == "pmm2") {
    if (verbose) cat("Starting PMM2 optimization...\n")

    model_info <- list(
      ar_order = model_params$ar_order,
      ma_order = model_params$ma_order,
      d = model_params$d,
      model_type = model_params$model_type,
      include.mean = model_params$include.mean,
      innovations = innovations,
      x = x_centered,
      x_mean = x_mean,
      original_x = orig_x,
      verbose = verbose
    )

    result <- pmm2_algorithm(
      b_init = b_init,
      X = dm$X,
      y = dm$y,
      m2 = m2,
      m3 = m3,
      m4 = m4,
      max_iter = max_iter,
      tol = tol,
      regularize = regularize,
      reg_lambda = reg_lambda,
      verbose = verbose
    )

    final_coef <- result$b
    converged <- result$convergence
    iterations <- result$iterations

    # Compute final residuals
    final_res <- compute_ts_residuals(final_coef, model_info)
  } else {
    # For other methods, simply use initial estimates
    final_coef <- b_init
    converged <- TRUE
    iterations <- 0
    final_res <- innovations
  }

  # 5) Create the appropriate class object
  if (model_type == "ar") {
    result_class <- "ARPMM2"
  } else if (model_type == "ma") {
    result_class <- "MAPMM2"
  } else if (model_type == "arma") {
    result_class <- "ARMAPMM2"
  } else if (model_type == "arima") {
    result_class <- "ARIMAPMM2"
  } else {
    result_class <- "TS2fit"  # Base class by default
  }

  # Create and return the object of the appropriate class
  new(result_class,
      coefficients    = as.numeric(final_coef),
      residuals       = as.numeric(final_res),
      m2              = as.numeric(m2),
      m3              = as.numeric(m3),
      m4              = as.numeric(m4),
      convergence     = converged,
      iterations      = as.numeric(iterations),
      call            = cl,
      model_type      = model_type,
      intercept       = as.numeric(x_mean),
      original_series = as.numeric(orig_x),
      order           = list(ar = model_params$ar_order,
                             ma = model_params$ma_order,
                             d = model_params$d))
}

#' Fit an AR model using PMM2 (wrapper)
#'
#' @inheritParams ts_pmm2
#' @return An S4 object of class \code{ARPMM2} containing fitted autoregressive
#'   coefficients, residuals, central moment estimates (m2-m4), model order,
#'   intercept, original series, and convergence diagnostics.
#' @export
ar_pmm2 <- function(x, order = 1, method = "pmm2", max_iter = 50, tol = 1e-6,
                    include.mean = TRUE, initial = NULL, na.action = na.fail,
                    regularize = TRUE, reg_lambda = 1e-8, verbose = FALSE) {
  ts_pmm2(x, order = order, model_type = "ar", method = method,
          max_iter = max_iter, tol = tol,
          include.mean = include.mean, initial = initial,
          na.action = na.action, regularize = regularize,
          reg_lambda = reg_lambda, verbose = verbose)
}

#' Fit an MA model using PMM2 (wrapper)
#'
#' @inheritParams ts_pmm2
#' @return An S4 object of class \code{MAPMM2} containing moving-average
#'   coefficients, residual innovations, central moments, model order,
#'   intercept, original series, and convergence diagnostics.
#' @export
ma_pmm2 <- function(x, order = 1, method = "pmm2", max_iter = 50, tol = 1e-6,
                    include.mean = TRUE, initial = NULL, na.action = na.fail,
                    regularize = TRUE, reg_lambda = 1e-8, verbose = FALSE) {
  ts_pmm2(x, order = order, model_type = "ma", method = method,
          max_iter = max_iter, tol = tol,
          include.mean = include.mean, initial = initial,
          na.action = na.action, regularize = regularize,
          reg_lambda = reg_lambda, verbose = verbose)
}

#' Fit an ARMA model using PMM2 (wrapper)
#'
#' @inheritParams ts_pmm2
#' @return An S4 object of class \code{ARMAPMM2} containing fitted AR and MA
#'   coefficients, residuals, central moments, model specification, intercept,
#'   original series, and convergence diagnostics.
#' @export
arma_pmm2 <- function(x, order = c(1, 1), method = "pmm2", max_iter = 50, tol = 1e-6,
                      include.mean = TRUE, initial = NULL, na.action = na.fail,
                      regularize = TRUE, reg_lambda = 1e-8, verbose = FALSE) {
  ts_pmm2(x, order = order, model_type = "arma", method = method,
          max_iter = max_iter, tol = tol,
          include.mean = include.mean, initial = initial,
          na.action = na.action, regularize = regularize,
          reg_lambda = reg_lambda, verbose = verbose)
}

#' Fit an ARIMA model using PMM2 (wrapper)
#'
#' @inheritParams ts_pmm2
#' @return An S4 object of class \code{ARIMAPMM2} containing fitted AR and MA
#'   coefficients, residual series, central moments, differencing order,
#'   intercept, original series, and convergence diagnostics.
#' @export
arima_pmm2 <- function(x, order = c(1, 1, 1), method = "pmm2", max_iter = 50, tol = 1e-6,
                       include.mean = TRUE, initial = NULL, na.action = na.fail,
                       regularize = TRUE, reg_lambda = 1e-8, verbose = FALSE) {
  ts_pmm2(x, order = order, model_type = "arima", method = method,
          max_iter = max_iter, tol = tol,
          include.mean = include.mean, initial = initial,
          na.action = na.action, regularize = regularize,
          reg_lambda = reg_lambda, verbose = verbose)
}

#' Fit a Seasonal MA model using PMM2
#'
#' Fits a Seasonal Moving Average (SMA) model using the Polynomial Maximization
#' Method (PMM2). This is particularly effective when the innovations have
#' asymmetric or non-Gaussian distributions.
#'
#' @param x Numeric vector (time series data)
#' @param order Seasonal MA order (Q)
#' @param season List with seasonal specification: list(period = s)
#' @param method Estimation method: "pmm2" or "css"
#' @param max_iter Maximum iterations for PMM2 algorithm
#' @param tol Convergence tolerance
#' @param include.mean Include intercept in the model
#' @param na.action Function to handle missing values
#' @param regularize Add regularization to Jacobian matrix
#' @param reg_lambda Regularization parameter
#' @param verbose Print diagnostic information
#'
#' @return An S4 object of class \code{SMAPMM2} containing:
#'   \itemize{
#'     \item coefficients: Seasonal MA coefficients (Θ₁, Θ₂, ..., Θ_Q)
#'     \item innovations: Estimated innovations (residuals)
#'     \item m2, m3, m4: Central moments of innovations
#'     \item convergence: Convergence status
#'     \item iterations: Number of iterations performed
#'     \item intercept: Model intercept
#'     \item original_series: Original time series
#'     \item order: Model order list(Q, s)
#'   }
#'
#' @details
#' The SMA(Q)_s model has the form:
#'   y_t = μ + ε_t + Θ₁·ε_{t-s} + Θ₂·ε_{t-2s} + ... + Θ_Q·ε_{t-Qs}
#'
#' Where:
#'   - Q is the seasonal MA order
#'   - s is the seasonal period (12 for monthly, 4 for quarterly)
#'   - ε_t are innovations (errors)
#'
#' The PMM2 method provides more efficient parameter estimates than ML/CSS when
#' the innovation distribution is asymmetric (non-Gaussian). The expected
#' variance reduction is given by g = 1 - c₃²/(2 + c₄), where c₃ and c₄
#' are the skewness and excess kurtosis coefficients.
#'
#' @examples
#' \dontrun{
#' # Generate synthetic seasonal data
#' set.seed(123)
#' n <- 120
#' s <- 12
#' theta <- 0.6
#'
#' # Gamma innovations (asymmetric)
#' innov <- rgamma(n, shape = 2, scale = 1) - 2
#' y <- numeric(n)
#' for (t in 1:n) {
#'   ma_term <- if (t > s) theta * innov[t-s] else 0
#'   y[t] <- innov[t] + ma_term
#' }
#'
#' # Fit SMA(1)_12 model with PMM2
#' fit <- sma_pmm2(y, order = 1, season = list(period = 12))
#' summary(fit)
#'
#' # Compare with CSS
#' fit_css <- sma_pmm2(y, order = 1, season = list(period = 12), method = "css")
#' }
#'
#' @seealso \code{\link{ma_pmm2}}, \code{\link{sar_pmm2}}, \code{\link{arima_pmm2}}
#'
#' @export
sma_pmm2 <- function(x,
                     order = 1,
                     season = list(period = 12),
                     method = "pmm2",
                     max_iter = 50,
                     tol = 1e-6,
                     include.mean = TRUE,
                     na.action = na.fail,
                     regularize = TRUE,
                     reg_lambda = 1e-8,
                     verbose = FALSE) {

  # Validate inputs
  if (!is.numeric(x)) {
    stop("x must be a numeric vector")
  }

  Q <- as.integer(order[1])
  if (Q < 1) {
    stop("Seasonal MA order must be at least 1")
  }

  if (is.null(season$period) || season$period < 2) {
    stop("season$period must be specified and >= 2")
  }

  s <- as.integer(season$period)

  # Handle NA values
  na.action <- match.fun(na.action)
  x <- na.action(x)

  # Store call and original series
  cl <- match.call()
  orig_x <- x

  # Check data sufficiency
  min_obs <- Q * s + 10
  if (length(x) < min_obs) {
    warning(sprintf("Very small sample size. Need at least %d observations, have %d",
                    min_obs, length(x)))
  }

  if (verbose) {
    cat("Fitting SMA(", Q, ")_", s, " model\n", sep = "")
    cat("Method:", method, "\n")
    cat("Sample size:", length(x), "\n")
  }

  # Center data if needed
  x_mean <- if (include.mean) mean(x) else 0
  x_centered <- x - x_mean

  # Step 1: Initial CSS estimation
  css_fit <- sma_css_fit(x, Q, s, include_mean = include.mean, verbose = verbose)

  if (verbose) {
    cat("Initial CSS coefficients:\n")
    print(css_fit$coefficients)
  }

  # Step 2: Apply PMM2 refinement or return CSS
  if (method == "pmm2") {
    if (verbose) cat("\nApplying PMM2 refinement...\n")

    pmm2_fit <- sma_pmm2_fit(x, Q, s, css_fit,
                             max_iter = max_iter,
                             tol = tol,
                             verbose = verbose)

    final_coef <- pmm2_fit$coefficients
    final_innovations <- pmm2_fit$innovations
    converged <- pmm2_fit$convergence
    iterations <- pmm2_fit$iterations

  } else {
    # Return CSS estimates
    final_coef <- css_fit$coefficients
    final_innovations <- css_fit$residuals
    converged <- css_fit$convergence
    iterations <- css_fit$iterations
  }

  # Compute final moments
  final_moments <- compute_moments(final_innovations)

  if (verbose) {
    cat("\nFinal coefficients:\n")
    print(final_coef)
    cat("\nConvergence:", converged, "\n")
    cat("Iterations:", iterations, "\n")
  }

  # Create SMAPMM2 object
  result <- new("SMAPMM2",
                coefficients = as.numeric(final_coef),
                innovations = as.numeric(final_innovations),
                m2 = as.numeric(final_moments$m2),
                m3 = as.numeric(final_moments$m3),
                m4 = as.numeric(final_moments$m4),
                convergence = converged,
                iterations = as.integer(iterations),
                call = cl,
                model_type = "sma",
                intercept = as.numeric(x_mean),
                original_series = as.numeric(orig_x),
                order = list(Q = Q, s = s))

  return(result)
}

# --- MA utilities ---------------------------------------------------------

ma_css_fit <- function(x, q, include_mean = TRUE, verbose = FALSE) {
  fit <- tryCatch(
    stats::arima(x, order = c(0, 0, q), method = "CSS-ML",
                 include.mean = include_mean),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    if (verbose) {
      cat("Failed to estimate MA model via stats::arima: returning zero coefficients\n")
    }
    coef_css <- rep(0, q)
    intercept <- if (include_mean) mean(x) else 0
    residuals <- x - intercept
    residuals[is.na(residuals)] <- 0
    return(list(
      coefficients = coef_css,
      intercept = intercept,
      residuals = residuals,
      convergence = FALSE,
      iterations = 0L
    ))
  }

  names_coef <- names(fit$coef)
  coef_css <- numeric(q)
  for (j in seq_len(q)) {
    coef_name <- paste0("ma", j)
    if (coef_name %in% names_coef) {
      coef_css[j] <- fit$coef[coef_name]
    } else if (length(fit$coef) >= j) {
      coef_css[j] <- fit$coef[j]
    } else {
      coef_css[j] <- 0
    }
  }

  intercept <- 0
  if (include_mean) {
    intercept_name <- setdiff(names_coef, paste0("ma", seq_len(q)))
    if (length(intercept_name) > 0) {
      intercept <- as.numeric(fit$coef[intercept_name[1]])
    }
  }

  residuals <- as.numeric(fit$residuals)
  residuals[is.na(residuals)] <- 0

  list(
    coefficients = coef_css,
    intercept = intercept,
    residuals = residuals,
    convergence = TRUE,
    iterations = 1L
  )
}

ma_pmm2_fit <- function(x, q, css_fit, max_iter = 50, tol = 1e-6, verbose = FALSE) {
  design <- ma_build_design(css_fit$intercept, css_fit$residuals, x, q)
  moments <- compute_moments(css_fit$residuals)

  b_init <- c(0, css_fit$coefficients)
  solve_res <- ma_solve_pmm2(b_init, design$X, design$y,
                             moments$m2, moments$m3, moments$m4,
                             max_iter = max_iter, tol = tol,
                             verbose = verbose)

  theta <- solve_res$coefficients
  innovations <- ma_compute_innovations(x - css_fit$intercept, theta, q)

  list(
    coefficients = theta,
    intercept = css_fit$intercept,
    innovations = innovations,
    convergence = solve_res$convergence,
    iterations = solve_res$iterations
  )
}

ma_build_design <- function(intercept, residuals, x, q) {
  idx <- seq.int(q + 1L, length(x))
  X <- matrix(1, nrow = length(idx), ncol = q + 1L)
  for (j in seq_len(q)) {
    X[, j + 1L] <- residuals[idx - j]
  }
  y <- x[idx] - intercept
  list(X = X, y = y)
}

ma_solve_pmm2 <- function(b_init, X, Y, m2, m3, m4,
                          max_iter = 50, tol = 1e-6,
                          verbose = FALSE) {
  b <- as.numeric(b_init)
  iterations <- 0L
  converged <- FALSE
  for (iter in seq_len(max_iter)) {
    iterations <- iter
    S <- as.vector(X %*% b)
    Z1 <- m3 * S^2 + (m4 - m2^2 - 2 * m3 * Y) * S +
      (m3 * Y^2 - (m4 - m2^2) * Y - m2 * m3)
    Z <- as.numeric(t(X) %*% Z1)
    JZ11 <- 2 * m3 * S + (m4 - m2^2 - 2 * m3 * Y)
    J <- t(X) %*% (X * JZ11)
    step <- tryCatch(solve(J, Z), error = function(e) NULL)
    if (is.null(step)) {
      if (verbose) cat("System is singular at iteration", iter, "\n")
      break
    }
    b_new <- b - step
    if (sqrt(sum((b_new - b)^2)) < tol) {
      b <- b_new
      converged <- TRUE
      break
    }
    b <- b_new
  }
  list(coefficients = b[-1], convergence = converged, iterations = iterations)
}

ma_compute_innovations <- function(x, theta, q) {
  n <- length(x)
  innovations <- numeric(n)
  history <- rep(0, q)
  for (t in seq_len(n)) {
    ma_component <- if (q > 0) sum(theta * history) else 0
    innovations[t] <- x[t] - ma_component
    if (q > 0) {
      history <- c(innovations[t], history)[seq_len(q)]
    }
  }
  innovations
}

# --- SMA (Seasonal Moving Average) utilities --------------------------------

#' Build design matrix for Seasonal MA model
#'
#' @param intercept Model intercept
#' @param residuals Initial residuals from CSS fit
#' @param x Time series data
#' @param Q Seasonal MA order
#' @param s Seasonal period
#'
#' @return List with design matrix X and dependent variable y
#' @keywords internal
sma_build_design <- function(intercept, residuals, x, Q, s) {
  max_lag <- Q * s
  idx <- seq.int(max_lag + 1L, length(x))

  X <- matrix(1, nrow = length(idx), ncol = Q + 1L)

  # Add seasonal MA lags: ε_{t-s}, ε_{t-2s}, ..., ε_{t-Qs}
  for (j in seq_len(Q)) {
    lag_seasonal <- j * s  # KEY DIFFERENCE from MA: j*s instead of j
    X[, j + 1L] <- residuals[idx - lag_seasonal]
  }

  y <- x[idx] - intercept
  list(X = X, y = y)
}

#' Compute innovations for Seasonal MA model
#'
#' @param x Centered time series
#' @param theta Seasonal MA coefficients
#' @param Q Seasonal MA order
#' @param s Seasonal period
#'
#' @return Vector of innovations
#' @keywords internal
sma_compute_innovations <- function(x, theta, Q, s) {
  n <- length(x)
  innovations <- numeric(n)

  for (t in seq_len(n)) {
    ma_component <- 0

    # Add seasonal MA terms: Θ₁·ε_{t-s} + Θ₂·ε_{t-2s} + ...
    for (j in seq_len(Q)) {
      lag_seasonal <- j * s
      if (t > lag_seasonal) {
        ma_component <- ma_component + theta[j] * innovations[t - lag_seasonal]
      }
    }

    innovations[t] <- x[t] - ma_component
  }

  innovations
}

#' Fit Seasonal MA model using PMM2
#'
#' @param x Time series data
#' @param Q Seasonal MA order
#' @param s Seasonal period
#' @param css_fit Initial CSS fit
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @param verbose Print progress
#'
#' @return List with coefficients, intercept, innovations, convergence info
#' @keywords internal
sma_pmm2_fit <- function(x, Q, s, css_fit, max_iter = 50, tol = 1e-6, verbose = FALSE) {
  # Build SMA design matrix
  design <- sma_build_design(css_fit$intercept, css_fit$residuals, x, Q, s)

  # Compute moments
  moments <- compute_moments(css_fit$residuals)

  # Initial parameters: [intercept, Θ₁, Θ₂, ..., Θ_Q]
  b_init <- c(0, css_fit$coefficients)

  # Solve PMM2 (reuse existing MA solver!)
  solve_res <- ma_solve_pmm2(b_init, design$X, design$y,
                             moments$m2, moments$m3, moments$m4,
                             max_iter = max_iter, tol = tol,
                             verbose = verbose)

  theta <- solve_res$coefficients

  # Compute final innovations with SMA structure
  innovations <- sma_compute_innovations(x - css_fit$intercept, theta, Q, s)

  list(
    coefficients = theta,
    intercept = css_fit$intercept,
    innovations = innovations,
    convergence = solve_res$convergence,
    iterations = solve_res$iterations
  )
}

#' CSS fit for Seasonal MA model
#'
#' @param x Time series data
#' @param Q Seasonal MA order
#' @param s Seasonal period
#' @param include_mean Include intercept
#' @param verbose Print diagnostics
#'
#' @return List with coefficients, intercept, residuals, convergence info
#' @keywords internal
sma_css_fit <- function(x, Q, s, include_mean = TRUE, verbose = FALSE) {
  fit <- tryCatch(
    stats::arima(x,
                 order = c(0, 0, 0),
                 seasonal = list(order = c(0, 0, Q), period = s),
                 method = "CSS-ML",
                 include.mean = include_mean),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    if (verbose) {
      cat("Failed to estimate SMA model via stats::arima: returning zero coefficients\n")
    }
    coef_css <- rep(0, Q)
    intercept <- if (include_mean) mean(x) else 0
    residuals <- x - intercept
    residuals[is.na(residuals)] <- 0
    return(list(
      coefficients = coef_css,
      intercept = intercept,
      residuals = residuals,
      convergence = FALSE,
      iterations = 0L
    ))
  }

  # Extract seasonal MA coefficients
  names_coef <- names(fit$coef)
  coef_css <- numeric(Q)
  for (j in seq_len(Q)) {
    coef_name <- paste0("sma", j)
    if (coef_name %in% names_coef) {
      coef_css[j] <- fit$coef[coef_name]
    } else if (length(fit$coef) >= j) {
      coef_css[j] <- fit$coef[j]
    } else {
      coef_css[j] <- 0
    }
  }

  intercept <- 0
  if (include_mean) {
    intercept_name <- setdiff(names_coef, paste0("sma", seq_len(Q)))
    if (length(intercept_name) > 0) {
      intercept <- as.numeric(fit$coef[intercept_name[1]])
    }
  }

  residuals <- as.numeric(fit$residuals)
  residuals[is.na(residuals)] <- 0

  list(
    coefficients = coef_css,
    intercept = intercept,
    residuals = residuals,
    convergence = fit$code == 0,
    iterations = 1L
  )
}

arma_build_design <- function(x, residuals, p, q, intercept = 0, include_intercept = FALSE) {
  n <- length(x)
  max_lag <- max(p, q)
  if (n <= max_lag) {
    stop("Insufficient data to build ARMA design matrix")
  }

  idx <- seq.int(max_lag + 1L, n)
  columns <- list()
  if (include_intercept) {
    columns <- c(columns, list(rep(1, length(idx))))
  }
  if (p > 0) {
    for (j in seq_len(p)) {
      columns <- c(columns, list(x[idx - j]))
    }
  }
  if (q > 0) {
    for (j in seq_len(q)) {
      columns <- c(columns, list(residuals[idx - j]))
    }
  }

  X <- if (length(columns) > 0) do.call(cbind, columns) else matrix(0, length(idx), 0)
  y <- x[idx] - intercept
  list(X = X, y = y)
}


#' Fit Seasonal AR model using PMM2 method
#'
#' Fits a Seasonal Autoregressive (SAR) model using the Polynomial Maximization
#' Method (PMM2). The model can include both non-seasonal and seasonal AR components.
#'
#' @param x Numeric vector of time series data
#' @param order Vector of length 2: c(p, P) where:
#'   \itemize{
#'     \item p: Non-seasonal AR order
#'     \item P: Seasonal AR order
#'   }
#' @param season List with seasonal specification: list(period = s)
#'   where s is the seasonal period (e.g., 12 for monthly data with annual seasonality)
#' @param method Estimation method: "pmm2" (default), "ols", "css"
#' @param include.mean Logical, include intercept/mean term (default TRUE)
#' @param multiplicative Logical, use multiplicative form with cross-terms (default FALSE)
#' @param max_iter Maximum iterations for PMM2 algorithm (default 50)
#' @param tol Convergence tolerance (default 1e-6)
#' @param regularize Logical, use regularization for numerical stability (default TRUE)
#' @param reg_lambda Regularization parameter (default 1e-8)
#' @param verbose Logical, print progress information (default FALSE)
#'
#' @return S4 object of class SARPMM2 containing:
#'   \itemize{
#'     \item coefficients: Estimated AR and SAR parameters
#'     \item residuals: Model residuals/innovations
#'     \item m2, m3, m4: Central moments of residuals
#'     \item convergence: Convergence status
#'     \item iterations: Number of iterations performed
#'   }
#'
#' @details
#' The SAR model has the form:
#'   y_t = φ₁·y_{t-1} + ... + φ_p·y_{t-p} + Φ₁·y_{t-s} + ... + Φ_P·y_{t-Ps} + μ + ε_t
#'
#' Where:
#'   - p is the non-seasonal AR order
#'   - P is the seasonal AR order
#'   - s is the seasonal period
#'   - ε_t are innovations (errors)
#'
#' The PMM2 method provides more efficient parameter estimates than OLS when
#' the innovation distribution is asymmetric (non-Gaussian). The expected
#' variance reduction is given by g = 1 - c₃²/(2 + c₄), where c₃ and c₄
#' are the skewness and excess kurtosis coefficients.
#'
#' @examples
#' \dontrun{
#' # Generate synthetic seasonal data
#' n <- 120
#' y <- arima.sim(n = n, list(ar = 0.7, seasonal = list(sar = 0.5, period = 12)))
#'
#' # Fit SAR(1,1)_12 model with PMM2
#' fit <- sar_pmm2(y, order = c(1, 1), season = list(period = 12))
#' summary(fit)
#'
#' # Simple seasonal model (no non-seasonal component)
#' fit_pure_sar <- sar_pmm2(y, order = c(0, 1), season = list(period = 12))
#'
#' # Compare with OLS
#' fit_ols <- sar_pmm2(y, order = c(1, 1), season = list(period = 12), method = "ols")
#' }
#'
#' @seealso \code{\link{ar_pmm2}}, \code{\link{ts_pmm2}}, \code{\link{compare_sar_methods}}
#'
#' @export
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

  if (p < 0 || P < 0) {
    stop("AR orders must be non-negative")
  }
  if (p == 0 && P == 0) {
    stop("At least one of p or P must be positive")
  }

  # Parse seasonal specification
  if (!is.list(season) || is.null(season$period)) {
    stop("'season' must be a list with 'period' element")
  }
  s <- as.integer(season$period)

  if (s <= 1) {
    stop("Seasonal period must be greater than 1")
  }

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
  X <- create_sar_matrix(x_centered, p = p, P = P, s = s, multiplicative = multiplicative)
  max_lag <- max(p, P * s)

  if (multiplicative && p > 0 && P > 0) {
    max_lag <- max(max_lag, p + P * s)
  }

  y <- x_centered[(max_lag + 1):length(x_centered)]

  if (verbose) {
    cat("Design matrix dimensions:", nrow(X), "x", ncol(X), "\n")
    cat("Effective sample size:", length(y), "\n")
  }

  # Check for sufficient data
  if (nrow(X) < ncol(X) + 5) {
    warning("Very small sample size relative to number of parameters. ",
            "Need at least ", ncol(X) + 5, " observations after lags, have ", nrow(X))
  }

  # Step 1: Initial estimation (OLS/CSS)
  if (method %in% c("pmm2", "css", "ols")) {
    # Use OLS for initial estimates
    XtX <- t(X) %*% X
    Xty <- t(X) %*% y

    # Add small regularization if needed
    if (det(XtX) < 1e-10) {
      if (verbose) cat("Adding regularization to X'X for initial fit\n")
      diag(XtX) <- diag(XtX) + 1e-6
    }

    b_init <- as.numeric(solve(XtX, Xty))

    if (verbose) {
      cat("Initial coefficients (OLS):\n")
      print(b_init)
    }
  }

  # Step 2: Compute moments from initial residuals
  residuals_init <- as.numeric(y - X %*% b_init)
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

  # Step 3: Apply PMM2 algorithm or return initial estimates
  if (method == "pmm2") {
    if (verbose) cat("\nApplying PMM2 refinement...\n")

    # Use existing pmm2_algorithm function
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

    b_final <- pmm2_result$b  # Already converted to numeric by pmm2_algorithm
    converged <- pmm2_result$convergence
    iterations <- pmm2_result$iterations

    # Compute final residuals
    residuals_final <- as.numeric(y - X %*% b_final)
    moments_final <- compute_moments(residuals_final)

    if (verbose) {
      cat("\nPMM2 converged:", converged, "\n")
      cat("Iterations:", iterations, "\n")
      cat("\nFinal coefficients:\n")
      print(b_final)
    }

  } else {
    # Return initial estimates for OLS/CSS methods
    b_final <- b_init  # Already numeric from line 712
    residuals_final <- residuals_init
    moments_final <- moments
    converged <- TRUE
    iterations <- 0L
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
}

#' Fit a Seasonal ARMA model using PMM2 method
#'
#' Fits a Seasonal Autoregressive Moving Average (SARMA) model that combines
#' both seasonal AR and seasonal MA components, optionally with non-seasonal
#' AR and MA terms as well.
#'
#' @param x Numeric vector of time series data
#' @param order Vector of length 4: c(p, P, q, Q) where:
#'   \itemize{
#'     \item p: Non-seasonal AR order
#'     \item P: Seasonal AR order
#'     \item q: Non-seasonal MA order
#'     \item Q: Seasonal MA order
#'   }
#' @param season List with seasonal specification: list(period = s)
#'   where s is the seasonal period (e.g., 12 for monthly data)
#' @param method Estimation method: "pmm2" (default), "css-ml", "css"
#' @param include.mean Logical, include intercept/mean term (default TRUE)
#' @param max_iter Maximum iterations for PMM2 algorithm (default 50)
#' @param tol Convergence tolerance (default 1e-6)
#' @param regularize Logical, use regularization for numerical stability (default TRUE)
#' @param reg_lambda Regularization parameter (default 1e-8)
#' @param verbose Logical, print progress information (default FALSE)
#'
#' @return S4 object of class SARMAPMM2 containing:
#'   \itemize{
#'     \item coefficients: Estimated parameters
#'     \item residuals: Model residuals/innovations
#'     \item m2, m3, m4: Central moments of residuals
#'     \item convergence: Convergence status
#'     \item iterations: Number of iterations performed
#'   }
#'
#' @details
#' The SARMA model combines four components:
#'   - AR(p): φ₁·y_{t-1} + ... + φ_p·y_{t-p}
#'   - SAR(P)_s: Φ₁·y_{t-s} + ... + Φ_P·y_{t-Ps}
#'   - MA(q): θ₁·ε_{t-1} + ... + θ_q·ε_{t-q}
#'   - SMA(Q)_s: Θ₁·ε_{t-s} + ... + Θ_Q·ε_{t-Qs}
#'
#' The PMM2 method provides more efficient parameter estimates than ML/CSS when
#' the innovation distribution is asymmetric (non-Gaussian).
#'
#' @examples
#' \dontrun{
#' # Generate synthetic seasonal data with SARMA structure
#' set.seed(123)
#' n <- 200
#' y <- arima.sim(n = n, list(ar = 0.5, ma = 0.3,
#'                             seasonal = list(sar = 0.6, sma = 0.4, period = 12)))
#'
#' # Fit SARMA(1,1,1,1)_12 model with PMM2
#' fit <- sarma_pmm2(y, order = c(1, 1, 1, 1), season = list(period = 12))
#' summary(fit)
#'
#' # Pure seasonal model (no non-seasonal components)
#' fit_pure <- sarma_pmm2(y, order = c(0, 1, 0, 1), season = list(period = 12))
#' }
#'
#' @seealso \code{\link{sar_pmm2}}, \code{\link{sma_pmm2}}, \code{\link{sarima_pmm2}}
#'
#' @export
sarma_pmm2 <- function(x,
                       order = c(0, 1, 0, 1),
                       season = list(period = 12),
                       method = "pmm2",
                       include.mean = TRUE,
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
  if (length(order) != 4) {
    stop("'order' must be a vector of length 4: c(p, P, q, Q)")
  }
  p <- as.integer(order[1])  # Non-seasonal AR order
  P <- as.integer(order[2])  # Seasonal AR order
  q <- as.integer(order[3])  # Non-seasonal MA order
  Q <- as.integer(order[4])  # Seasonal MA order

  if (p < 0 || P < 0 || q < 0 || Q < 0) {
    stop("All orders must be non-negative")
  }
  if (p == 0 && P == 0 && q == 0 && Q == 0) {
    stop("At least one order must be positive")
  }

  # Parse seasonal specification
  if (!is.list(season) || is.null(season$period)) {
    stop("'season' must be a list with 'period' element")
  }
  s <- as.integer(season$period)

  if (s <= 1) {
    stop("Seasonal period must be greater than 1")
  }

  if (verbose) {
    cat("Fitting SARMA(", p, ",", q, ")×(", P, ",", Q, ")_", s, " model\n", sep = "")
    cat("Method:", method, "\n")
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

  # Step 1: Get initial estimates using stats::arima
  arima_fit <- tryCatch(
    stats::arima(x,
                 order = c(p, 0, q),
                 seasonal = list(order = c(P, 0, Q), period = s),
                 method = "CSS-ML",
                 include.mean = include.mean),
    error = function(e) {
      if (verbose) cat("Warning: Initial ARIMA fit failed:", e$message, "\n")
      NULL
    }
  )

  if (is.null(arima_fit)) {
    stop("Failed to obtain initial estimates using stats::arima")
  }

  # Extract initial coefficients
  coef_names <- names(arima_fit$coef)
  ar_init <- if (p > 0) as.numeric(arima_fit$coef[paste0("ar", 1:p)]) else numeric(0)
  sar_init <- if (P > 0) as.numeric(arima_fit$coef[paste0("sar", 1:P)]) else numeric(0)
  ma_init <- if (q > 0) as.numeric(arima_fit$coef[paste0("ma", 1:q)]) else numeric(0)
  sma_init <- if (Q > 0) as.numeric(arima_fit$coef[paste0("sma", 1:Q)]) else numeric(0)

  # Handle any NAs in initial estimates
  ar_init[is.na(ar_init)] <- 0
  sar_init[is.na(sar_init)] <- 0
  ma_init[is.na(ma_init)] <- 0
  sma_init[is.na(sma_init)] <- 0

  b_init <- c(ar_init, sar_init, ma_init, sma_init)

  # Get initial residuals
  residuals_init <- as.numeric(arima_fit$residuals)
  residuals_init[is.na(residuals_init)] <- 0

  if (verbose) {
    cat("Initial coefficients:\n")
    cat("  AR:", ar_init, "\n")
    cat("  SAR:", sar_init, "\n")
    cat("  MA:", ma_init, "\n")
    cat("  SMA:", sma_init, "\n")
  }

  # Compute moments from initial residuals
  moments <- compute_moments(residuals_init)

  if (verbose) {
    cat("\nMoments from initial residuals:\n")
    cat("  m2:", moments$m2, "\n")
    cat("  m3:", moments$m3, "\n")
    cat("  m4:", moments$m4, "\n")
  }

  # If method is not pmm2, return initial estimates
  if (method != "pmm2") {
    return(new("SARMAPMM2",
               coefficients = as.numeric(b_init),
               residuals = as.numeric(residuals_init),
               m2 = as.numeric(moments$m2),
               m3 = as.numeric(moments$m3),
               m4 = as.numeric(moments$m4),
               convergence = TRUE,
               iterations = 0L,
               call = cl,
               model_type = "sarma",
               intercept = if (include.mean) as.numeric(x_mean) else 0,
               original_series = as.numeric(orig_x),
               order = list(ar = p, sar = P, ma = q, sma = Q, period = s)))
  }

  # Step 2: Build design matrix for PMM2
  X <- create_sarma_matrix(x_centered, residuals_init,
                           p = p, P = P, q = q, Q = Q, s = s)

  max_ar_lag <- max(p, P * s)
  max_ma_lag <- max(q, Q * s)
  max_lag <- max(max_ar_lag, max_ma_lag)

  y <- x_centered[(max_lag + 1):length(x_centered)]

  if (verbose) {
    cat("\nDesign matrix dimensions:", nrow(X), "x", ncol(X), "\n")
    cat("Effective sample size:", length(y), "\n")
  }

  # Check for sufficient data
  if (nrow(X) < ncol(X) + 5) {
    warning("Very small sample size relative to number of parameters.")
  }

  # Step 3: Apply PMM2 algorithm
  if (verbose) cat("\nApplying PMM2 refinement...\n")

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

  b_final <- pmm2_result$b
  converged <- pmm2_result$convergence
  iterations <- pmm2_result$iterations

  # Compute final residuals using stats::arima with fixed parameters
  final_coef_names <- c(
    if (p > 0) paste0("ar", 1:p) else NULL,
    if (P > 0) paste0("sar", 1:P) else NULL,
    if (q > 0) paste0("ma", 1:q) else NULL,
    if (Q > 0) paste0("sma", 1:Q) else NULL
  )

  if (include.mean) {
    final_coef_names <- c(final_coef_names, "intercept")
    final_coef_values <- c(b_final, x_mean)
  } else {
    final_coef_values <- b_final
  }

  names(final_coef_values) <- final_coef_names

  final_fit <- tryCatch(
    stats::arima(x,
                 order = c(p, 0, q),
                 seasonal = list(order = c(P, 0, Q), period = s),
                 fixed = final_coef_values,
                 include.mean = FALSE),  # Mean already included in fixed
    error = function(e) {
      if (verbose) cat("Warning: Final fit failed, using PMM2 residuals\n")
      list(residuals = pmm2_result$residuals)
    }
  )

  residuals_final <- as.numeric(final_fit$residuals)
  if (length(residuals_final) < length(x)) {
    residuals_final <- c(rep(NA, length(x) - length(residuals_final)), residuals_final)
  }

  moments_final <- compute_moments(residuals_final[!is.na(residuals_final)])

  if (verbose) {
    cat("\nPMM2 converged:", converged, "\n")
    cat("Iterations:", iterations, "\n")
    cat("\nFinal coefficients:\n")
    print(b_final)
  }

  # Create SARMAPMM2 object
  result <- new("SARMAPMM2",
                coefficients = as.numeric(b_final),
                residuals = as.numeric(residuals_final),
                m2 = as.numeric(moments_final$m2),
                m3 = as.numeric(moments_final$m3),
                m4 = as.numeric(moments_final$m4),
                convergence = converged,
                iterations = as.integer(iterations),
                call = cl,
                model_type = "sarma",
                intercept = if (include.mean) as.numeric(x_mean) else 0,
                original_series = as.numeric(orig_x),
                order = list(ar = p, sar = P, ma = q, sma = Q, period = s))

  return(result)
}

#' Fit a Seasonal ARIMA model using PMM2 method
#'
#' Fits a Seasonal Autoregressive Integrated Moving Average (SARIMA) model
#' with both non-seasonal and seasonal differencing operators.
#'
#' @param x Numeric vector of time series data
#' @param order Vector of length 4: c(p, P, q, Q) where:
#'   \itemize{
#'     \item p: Non-seasonal AR order
#'     \item P: Seasonal AR order
#'     \item q: Non-seasonal MA order
#'     \item Q: Seasonal MA order
#'   }
#' @param seasonal List with seasonal specification: list(order = c(d, D), period = s)
#'   \itemize{
#'     \item d: Non-seasonal differencing order
#'     \item D: Seasonal differencing order
#'     \item period: Seasonal period (s)
#'   }
#' @param method Estimation method: "pmm2" (default), "css-ml", "css"
#' @param include.mean Logical, include drift term (default TRUE for d+D=0, FALSE otherwise)
#' @param max_iter Maximum iterations for PMM2 algorithm (default 50)
#' @param tol Convergence tolerance (default 1e-6)
#' @param regularize Logical, use regularization for numerical stability (default TRUE)
#' @param reg_lambda Regularization parameter (default 1e-8)
#' @param verbose Logical, print progress information (default FALSE)
#'
#' @return S4 object of class SARIMAPMM2 containing fitted model results
#'
#' @details
#' The SARIMA(p,d,q)×(P,D,Q)_s model is:
#'   (1-φ₁B-...-φ_pB^p)(1-Φ₁B^s-...-Φ_PB^{Ps})(1-B)^d(1-B^s)^D y_t =
#'   (1+θ₁B+...+θ_qB^q)(1+Θ₁B^s+...+Θ_QB^{Qs}) ε_t
#'
#' Where B is the backshift operator. The model combines:
#'   - Non-seasonal differencing: (1-B)^d
#'   - Seasonal differencing: (1-B^s)^D
#'   - Non-seasonal ARMA components
#'   - Seasonal ARMA components
#'
#' @examples
#' \dontrun{
#' # Generate synthetic data
#' set.seed(123)
#' n <- 200
#' y <- arima.sim(n = n, list(order = c(1, 1, 1),
#'                             seasonal = list(order = c(1, 1, 1), period = 12)))
#'
#' # Fit SARIMA(1,1,1)×(1,1,1)_12 model
#' fit <- sarima_pmm2(y, order = c(1, 1, 1, 1),
#'                    seasonal = list(order = c(1, 1), period = 12))
#' summary(fit)
#' }
#'
#' @seealso \code{\link{sarma_pmm2}}, \code{\link{arima_pmm2}}
#'
#' @export
sarima_pmm2 <- function(x,
                        order = c(0, 1, 0, 1),
                        seasonal = list(order = c(1, 1), period = 12),
                        method = "pmm2",
                        include.mean = NULL,
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
  if (length(order) != 4) {
    stop("'order' must be a vector of length 4: c(p, P, q, Q)")
  }
  p <- as.integer(order[1])  # Non-seasonal AR order
  P <- as.integer(order[2])  # Seasonal AR order
  q <- as.integer(order[3])  # Non-seasonal MA order
  Q <- as.integer(order[4])  # Seasonal MA order

  # Parse seasonal specification
  if (!is.list(seasonal) || is.null(seasonal$order) || is.null(seasonal$period)) {
    stop("'seasonal' must be a list with 'order' and 'period' elements")
  }
  if (length(seasonal$order) != 2) {
    stop("'seasonal$order' must be a vector of length 2: c(d, D)")
  }

  d <- as.integer(seasonal$order[1])  # Non-seasonal differencing
  D <- as.integer(seasonal$order[2])  # Seasonal differencing
  s <- as.integer(seasonal$period)    # Seasonal period

  if (p < 0 || P < 0 || q < 0 || Q < 0 || d < 0 || D < 0) {
    stop("All orders must be non-negative")
  }
  if (s <= 1) {
    stop("Seasonal period must be greater than 1")
  }

  # Determine if mean should be included
  if (is.null(include.mean)) {
    include.mean <- (d + D == 0)
  }

  if (verbose) {
    cat("Fitting SARIMA(", p, ",", d, ",", q, ")×(", P, ",", D, ",", Q, ")_", s, " model\n", sep = "")
    cat("Method:", method, "\n")
    cat("Include mean:", include.mean, "\n")
  }

  # Store original series
  orig_x <- x

  # Step 1: Get initial estimates using stats::arima
  arima_fit <- tryCatch(
    stats::arima(x,
                 order = c(p, d, q),
                 seasonal = list(order = c(P, D, Q), period = s),
                 method = "CSS-ML",
                 include.mean = include.mean),
    error = function(e) {
      if (verbose) cat("Warning: Initial ARIMA fit failed:", e$message, "\n")
      NULL
    }
  )

  if (is.null(arima_fit)) {
    stop("Failed to obtain initial estimates using stats::arima")
  }

  # Extract initial coefficients
  coef_names <- names(arima_fit$coef)
  ar_init <- if (p > 0) as.numeric(arima_fit$coef[paste0("ar", 1:p)]) else numeric(0)
  sar_init <- if (P > 0) as.numeric(arima_fit$coef[paste0("sar", 1:P)]) else numeric(0)
  ma_init <- if (q > 0) as.numeric(arima_fit$coef[paste0("ma", 1:q)]) else numeric(0)
  sma_init <- if (Q > 0) as.numeric(arima_fit$coef[paste0("sma", 1:Q)]) else numeric(0)

  # Handle any NAs
  ar_init[is.na(ar_init)] <- 0
  sar_init[is.na(sar_init)] <- 0
  ma_init[is.na(ma_init)] <- 0
  sma_init[is.na(sma_init)] <- 0

  # Extract mean/drift if present
  x_mean <- 0
  if (include.mean && "intercept" %in% coef_names) {
    x_mean <- as.numeric(arima_fit$coef["intercept"])
  } else if (include.mean && "drift" %in% coef_names) {
    x_mean <- as.numeric(arima_fit$coef["drift"])
  }

  b_init <- c(ar_init, sar_init, ma_init, sma_init)

  # Get initial residuals
  residuals_init <- as.numeric(arima_fit$residuals)
  residuals_init[is.na(residuals_init)] <- 0

  if (verbose) {
    cat("Initial coefficients:\n")
    cat("  AR:", ar_init, "\n")
    cat("  SAR:", sar_init, "\n")
    cat("  MA:", ma_init, "\n")
    cat("  SMA:", sma_init, "\n")
    if (include.mean) cat("  Mean/Drift:", x_mean, "\n")
  }

  # Compute moments from initial residuals
  moments <- compute_moments(residuals_init[residuals_init != 0])

  if (verbose) {
    cat("\nMoments from initial residuals:\n")
    cat("  m2:", moments$m2, "\n")
    cat("  m3:", moments$m3, "\n")
    cat("  m4:", moments$m4, "\n")
  }

  # If method is not pmm2, return initial estimates
  if (method != "pmm2") {
    return(new("SARIMAPMM2",
               coefficients = as.numeric(b_init),
               residuals = as.numeric(residuals_init),
               m2 = as.numeric(moments$m2),
               m3 = as.numeric(moments$m3),
               m4 = as.numeric(moments$m4),
               convergence = TRUE,
               iterations = 0L,
               call = cl,
               model_type = "sarima",
               intercept = if (include.mean) as.numeric(x_mean) else 0,
               original_series = as.numeric(orig_x),
               order = list(ar = p, sar = P, ma = q, sma = Q,
                           d = d, D = D, period = s)))
  }

  # Apply differencing to get stationary series
  x_diff <- x
  if (d > 0) {
    x_diff <- diff(x_diff, differences = d)
  }
  if (D > 0) {
    x_diff <- diff(x_diff, lag = s, differences = D)
  }

  # Center if needed
  if (include.mean && d + D == 0) {
    x_diff <- x_diff - mean(x_diff)
  }

  # Align residuals with differenced series
  residuals_diff <- tail(residuals_init, length(x_diff))

  # Step 2: Build design matrix for PMM2
  X <- create_sarma_matrix(x_diff, residuals_diff,
                           p = p, P = P, q = q, Q = Q, s = s)

  max_ar_lag <- max(p, P * s)
  max_ma_lag <- max(q, Q * s)
  max_lag <- max(max_ar_lag, max_ma_lag)

  y <- x_diff[(max_lag + 1):length(x_diff)]

  if (verbose) {
    cat("\nDesign matrix dimensions:", nrow(X), "x", ncol(X), "\n")
    cat("Effective sample size:", length(y), "\n")
  }

  # Step 3: Apply PMM2 algorithm
  if (verbose) cat("\nApplying PMM2 refinement...\n")

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

  b_final <- pmm2_result$b
  converged <- pmm2_result$convergence
  iterations <- pmm2_result$iterations

  # Compute final residuals
  final_coef_names <- c(
    if (p > 0) paste0("ar", 1:p) else NULL,
    if (P > 0) paste0("sar", 1:P) else NULL,
    if (q > 0) paste0("ma", 1:q) else NULL,
    if (Q > 0) paste0("sma", 1:Q) else NULL
  )

  if (include.mean) {
    final_coef_names <- c(final_coef_names, if (d + D == 0) "intercept" else "drift")
    final_coef_values <- c(b_final, x_mean)
  } else {
    final_coef_values <- b_final
  }

  names(final_coef_values) <- final_coef_names

  final_fit <- tryCatch(
    stats::arima(x,
                 order = c(p, d, q),
                 seasonal = list(order = c(P, D, Q), period = s),
                 fixed = final_coef_values,
                 include.mean = FALSE),
    error = function(e) {
      if (verbose) cat("Warning: Final fit failed\n")
      list(residuals = c(rep(NA, length(x) - length(pmm2_result$residuals)),
                        pmm2_result$residuals))
    }
  )

  residuals_final <- as.numeric(final_fit$residuals)
  moments_final <- compute_moments(residuals_final[!is.na(residuals_final)])

  if (verbose) {
    cat("\nPMM2 converged:", converged, "\n")
    cat("Iterations:", iterations, "\n")
  }

  # Create SARIMAPMM2 object
  result <- new("SARIMAPMM2",
                coefficients = as.numeric(b_final),
                residuals = as.numeric(residuals_final),
                m2 = as.numeric(moments_final$m2),
                m3 = as.numeric(moments_final$m3),
                m4 = as.numeric(moments_final$m4),
                convergence = converged,
                iterations = as.integer(iterations),
                call = cl,
                model_type = "sarima",
                intercept = if (include.mean) as.numeric(x_mean) else 0,
                original_series = as.numeric(orig_x),
                order = list(ar = p, sar = P, ma = q, sma = Q,
                            d = d, D = D, period = s))

  return(result)
}
