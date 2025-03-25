# pmm_ts.R

#' An S4 class to store time series PMM2 fit results (base class)
#'
#' @slot coefficients numeric vector of fitted parameters
#' @slot residuals numeric vector of final residuals
#' @slot m2 numeric second central moment from the initial residuals
#' @slot m3 numeric third central moment from the initial residuals
#' @slot m4 numeric fourth central moment from the initial residuals
#' @slot convergence logical or integer code indicating if the algorithm converged
#' @slot iterations integer number of iterations performed
#' @slot call the original function call
#' @slot model_type character string indicating the type of model
#' @slot intercept numeric intercept value
#' @slot original_series numeric original time series data
#' @slot order list of order parameters
#'
#' @exportClass TS2fit
setClass("TS2fit",
         slots = c(coefficients = "numeric",
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
                   order = "list"))

#' Validates the time series model parameters and prepares model information
#'
#' @param x Time series data
#' @param order Model order specification
#' @param model_type Type of model (ar, ma, arma, or arima)
#' @param include.mean Whether to include mean/intercept
#'
#' @return List of validated parameters and model information
#' @keywords internal
validate_ts_parameters <- function(x, order, model_type, include.mean) {
  # Check input validity
  if(missing(x)) {
    stop("'x' must be provided")
  }

  # Convert input to numeric vector
  x <- as.numeric(x)

  if(!is.numeric(x)) {
    stop("'x' must be a numeric vector")
  }

  if(missing(order)) {
    stop("'order' must be provided")
  }

  # Parse order parameter based on model_type
  if(model_type == "ar") {
    if(!is.numeric(order) || length(order) != 1)
      stop("For AR models, 'order' must be a single integer")
    ar_order <- as.integer(order)
    ma_order <- 0
    d <- 0
    if(ar_order <= 0)
      stop("AR order must be positive")
  } else if(model_type == "ma") {
    if(!is.numeric(order) || length(order) != 1)
      stop("For MA models, 'order' must be a single integer")
    ar_order <- 0
    ma_order <- as.integer(order)
    d <- 0
    if(ma_order <= 0)
      stop("MA order must be positive")
  } else if(model_type == "arma") {
    if(!is.numeric(order) || length(order) != 2)
      stop("For ARMA models, 'order' must be a vector of length 2 (AR order, MA order)")
    ar_order <- as.integer(order[1])
    ma_order <- as.integer(order[2])
    d <- 0
    if(ar_order < 0 || ma_order < 0)
      stop("AR and MA orders must be non-negative")
    if(ar_order == 0 && ma_order == 0)
      stop("At least one of AR or MA order must be positive")
  } else if(model_type == "arima") {
    if(!is.numeric(order) || length(order) != 3)
      stop("For ARIMA models, 'order' must be a vector of length 3 (AR order, differencing, MA order)")
    ar_order <- as.integer(order[1])
    d <- as.integer(order[2])
    ma_order <- as.integer(order[3])
    if(ar_order < 0 || ma_order < 0 || d < 0)
      stop("AR, differencing, and MA orders must be non-negative")
    if(ar_order == 0 && ma_order == 0 && d == 0)
      stop("At least one of AR, differencing, or MA order must be positive")
  } else {
    stop("Unknown model type: ", model_type)
  }

  # Store original series
  orig_x <- as.numeric(x)

  # Apply differencing for ARIMA models
  if(model_type == "arima" && d > 0) {
    x_diff <- as.numeric(diff(x, differences = d))
  } else {
    x_diff <- x
  }

  # Compute series mean if requested
  if(include.mean) {
    x_mean <- mean(x_diff)
    x_centered <- x_diff - x_mean
  } else {
    x_mean <- 0
    x_centered <- x_diff
  }

  # Order list
  order_list <- list(
    ar = as.integer(ar_order),
    ma = as.integer(ma_order),
    d = as.integer(d)
  )

  # Return validated and parsed parameters
  return(list(
    original_x = orig_x,
    x_centered = x_centered,
    x_mean = x_mean,
    ar_order = ar_order,
    ma_order = ma_order,
    d = d,
    model_type = model_type,
    order = order_list
  ))
}

#' Get initial parameter estimates for time series models
#'
#' @param model_params Validated model parameters from validate_ts_parameters
#' @param initial Optional user-provided initial estimates
#' @param method Estimation method
#' @param verbose Print verbose output
#'
#' @return List containing initial parameter estimates and innovations
#' @keywords internal
get_initial_estimates <- function(model_params, initial = NULL, method = "pmm2", verbose = FALSE) {
  x_centered <- model_params$x_centered
  x_mean <- model_params$x_mean
  ar_order <- model_params$ar_order
  ma_order <- model_params$ma_order
  d <- model_params$d
  model_type <- model_params$model_type

  if(model_type == "ar") {
    # Create lagged design matrix for AR models
    X <- create_ar_matrix(x_centered, ar_order)
    y <- x_centered[(ar_order + 1):length(x_centered)]

    # Get initial parameter estimates
    if(is.null(initial)) {
      if(method == "yw") {
        # Use Yule-Walker method
        b_init <- get_yw_estimates(x_centered, ar_order)
      } else {
        # Use OLS
        fit_ols <- lm.fit(x = X, y = y)
        b_init <- fit_ols$coefficients
      }
    } else {
      # Use provided initial estimates
      if(length(initial) != ar_order) {
        stop("Length of 'initial' must equal AR order")
      }
      b_init <- initial
    }

    # Calculate initial residuals
    res_init <- as.numeric(y - X %*% b_init)
    ar_part <- b_init
    ma_part <- numeric(0)
    innovations <- res_init

  } else if(model_type %in% c("ma", "arma", "arima")) {
    # For MA, ARMA, and ARIMA models, use arima() to get initial estimates
    arima_order <- c(ar_order, ifelse(model_type == "arima", d, 0), ma_order)

    if(is.null(initial)) {
      # Use arima() for initial estimates
      init_fit <- tryCatch({
        stats::arima(
          ifelse(model_type == "arima", model_params$original_x, x_centered),
          order = arima_order,
          method = ifelse(method == "pmm2", "CSS", method),
          include.mean = (x_mean != 0) & (model_type != "arima"))
      }, error = function(e) {
        stop("Error in initial parameter estimation: ", conditionMessage(e))
      })

      ar_init <- if(ar_order > 0) init_fit$coef[paste0("ar", 1:ar_order)] else numeric(0)
      ma_init <- if(ma_order > 0) init_fit$coef[paste0("ma", 1:ma_order)] else numeric(0)

      # Adjust mean if necessary
      if(x_mean != 0 && "intercept" %in% names(init_fit$coef)) {
        x_mean <- init_fit$coef["intercept"]
      }

      # Get innovations
      innovations <- as.numeric(residuals(init_fit))

    } else {
      # Handle provided initial estimates
      if(is.list(initial)) {
        if(!all(c("ar", "ma") %in% names(initial)) &&
           (ar_order > 0 || ma_order > 0)) {
          stop("'initial' must be a list with components 'ar' and 'ma'")
        }
        ar_init <- if(ar_order > 0) initial$ar else numeric(0)
        ma_init <- if(ma_order > 0) initial$ma else numeric(0)
      } else {
        # Assume it's a vector with AR coefficients followed by MA coefficients
        if(length(initial) != ar_order + ma_order) {
          stop("Length of 'initial' must equal sum of AR and MA orders")
        }
        ar_init <- if(ar_order > 0) initial[1:ar_order] else numeric(0)
        ma_init <- if(ma_order > 0) initial[(ar_order+1):(ar_order+ma_order)] else numeric(0)
      }

      # Validate lengths
      if(ar_order > 0 && length(ar_init) != ar_order) {
        stop("Length of AR initial estimates must match AR order")
      }
      if(ma_order > 0 && length(ma_init) != ma_order) {
        stop("Length of MA initial estimates must match MA order")
      }

      # Get innovations using arima with fixed parameters
      fixed_params <- c(ar_init, ma_init)

      init_fit <- tryCatch({
        stats::arima(
          ifelse(model_type == "arima", model_params$original_x, x_centered),
          order = arima_order,
          include.mean = (x_mean != 0) & (model_type != "arima"),
          fixed = fixed_params)
      }, error = function(e) {
        stop("Error in initial model fitting with provided parameters: ", conditionMessage(e))
      })

      innovations <- as.numeric(residuals(init_fit))
    }

    ar_part <- ar_init
    ma_part <- ma_init
  }

  # Combine parameters for PMM algorithm
  b_init <- c(ar_part, ma_part)

  # Handle NAs in initial coefficients
  if(any(is.na(b_init))) {
    stop("Initial parameter estimation resulted in NA coefficients.")
  }

  # Print verbose output if requested
  if(verbose) {
    cat("Initial parameter estimates:\n")
    if(ar_order > 0) cat("  AR coefficients:", ar_part, "\n")
    if(ma_order > 0) cat("  MA coefficients:", ma_part, "\n")
    if(x_mean != 0) cat("  Intercept:", x_mean, "\n")
  }

  return(list(
    b_init = b_init,
    x_mean = x_mean,
    innovations = innovations
  ))
}

#' Calculate fitted values for AR models
#'
#' @param object A TS2fit object with model_type="ar"
#' @return Vector of fitted values
#' @keywords internal
get_ar_fitted <- function(object) {
  if(object@model_type != "ar") {
    stop("This function is only for AR models")
  }

  # Extract parameters
  x <- object@original_series
  ar_order <- object@order$ar
  ar_coef <- object@coefficients[1:ar_order]
  intercept <- object@intercept

  # Create lagged design matrix
  n <- length(x)
  if(intercept != 0) {
    x_centered <- x - intercept
  } else {
    x_centered <- x
  }

  X <- create_ar_matrix(x_centered, ar_order)

  # Calculate fitted values
  fitted <- as.vector(X %*% ar_coef) + intercept

  return(fitted)
}

#' Generate MA model predictions
#'
#' @param innovations Vector of innovations
#' @param ma_coef Vector of MA coefficients
#' @param n.ahead Number of steps to predict
#'
#' @return Vector of predictions
#' @keywords internal
ma_predictions <- function(innovations, ma_coef, n.ahead) {
  n <- length(innovations)
  q <- length(ma_coef)

  # Initialize predictions
  pred <- numeric(n.ahead)

  # Generate predictions
  for(i in 1:n.ahead) {
    for(j in 1:min(i, q)) {
      if(n - i + j > 0) {
        pred[i] <- pred[i] + ma_coef[j] * innovations[n - i + j]
      }
    }
  }

  return(pred)
}

#' Fit a time series model using the Polynomial Maximization Method (order 2)
#'
#' This function fits an AR, MA, ARMA, or ARIMA model to a time series using the
#' Polynomial Maximization Method (PMM) of order 2, which is robust against
#' non-Gaussian errors.
#'
#' @param x numeric vector of time series data
#' @param order specification of the model order:
#'        - For AR models: a single integer (the AR order)
#'        - For MA models: a single integer (the MA order)
#'        - For ARMA models: a vector c(p, q) (AR and MA orders)
#'        - For ARIMA models: a vector c(p, d, q) (AR, differencing, and MA orders)
#' @param model_type string specifying the model type: "ar", "ma", "arma", or "arima"
#' @param method string: estimation method, one of "pmm2" (default), "css", "ml", "yw", "ols"
#' @param max_iter integer: maximum number of iterations for the algorithm
#' @param tol numeric: tolerance for convergence
#' @param include.mean logical: whether to include a mean (intercept) term
#' @param initial list or vector of initial parameter estimates (optional)
#' @param na.action function to handle missing values, default is na.fail
#' @param regularize logical, add small value to diagonal for numerical stability
#' @param reg_lambda regularization parameter (if regularize=TRUE)
#' @param verbose logical: whether to print progress information
#'
#' @details
#' The PMM2 algorithm works as follows:
#'
#' 1. Fits an initial model using a standard method (OLS, Yule-Walker, CSS, or ML)
#' 2. Computes central moments (m2, m3, m4) from initial residuals/innovations
#' 3. Iteratively improves parameter estimates using PMM approach
#'
#' PMM2 is particularly useful when error terms are non-Gaussian and asymmetric.
#'
#' @return An S4 \code{TS2fit} object
#' @export
#'
#' @examples
#' \dontrun{
#' # AR model example
#' n <- 200
#' ar_coef <- c(0.7, -0.3)
#' x <- arima.sim(model = list(ar = ar_coef), n = n,
#'                rand.gen = function(n) rt(n, df=3))
#' fit_ar <- ts_pmm2(x, order = 2, model_type = "ar")
#'
#' # ARIMA model example
#' y <- arima.sim(model = list(ar = 0.7, ma = 0.4), n = n,
#'                rand.gen = function(n) rgamma(n, shape=2, scale=1) - 2)
#' z <- cumsum(y)  # Create non-stationary series
#' fit_arima <- ts_pmm2(z, order = c(1, 1, 1), model_type = "arima")
#' }
ts_pmm2 <- function(x, order, model_type = c("ar", "ma", "arma", "arima"),
                    method = "pmm2", max_iter = 50, tol = 1e-6,
                    include.mean = TRUE, initial = NULL, na.action = na.fail,
                    regularize = TRUE, reg_lambda = 1e-8, verbose = FALSE) {
  # Capture the call
  call <- match.call()

  # Match model_type argument
  model_type <- match.arg(model_type)

  # Handle missing values
  if(!is.null(na.action)) {
    x <- na.action(x)
  }

  # Validate and preprocess parameters
  model_params <- validate_ts_parameters(x, order, model_type, include.mean)

  # Get initial estimates
  initial_estimates <- get_initial_estimates(model_params, initial, method, verbose)
  b_init <- initial_estimates$b_init
  x_mean <- initial_estimates$x_mean
  innovations <- initial_estimates$innovations

  # Calculate central moments from residuals/innovations
  m2 <- mean(innovations^2)
  m3 <- mean(innovations^3)
  m4 <- mean(innovations^4)

  if(verbose) {
    cat("Initial moments from residuals/innovations:\n")
    cat("  m2 =", m2, "\n")
    cat("  m3 =", m3, "\n")
    cat("  m4 =", m4, "\n")
  }

  # Check for potential issues
  if(m2 <= 0) {
    warning("Second central moment (m2) is non-positive. Results may be unreliable.")
  }

  if(m4 <= m2^2) {
    warning("Fourth central moment (m4) is less than m2^2. This violates a basic inequality for probability distributions.")
  }

  # Call the PMM2 fitting function if requested
  if(method == "pmm2") {
    if(verbose) cat("Starting PMM2 iterations...\n")

    # Prepare model info for .ts_pmm2_fit
    model_info <- list(
      x = model_params$x_centered,
      model_type = model_type,
      ar_order = model_params$ar_order,
      ma_order = model_params$ma_order
    )

    out <- .ts_pmm2_fit(
      b_init = b_init,
      innovations_init = innovations,
      model_info = model_info,
      m2 = m2,
      m3 = m3,
      m4 = m4,
      max_iter = max_iter,
      tol = tol,
      regularize = regularize,
      reg_lambda = reg_lambda,
      verbose = verbose
    )

    b_est <- out$b
    conv <- out$convergence
    iter <- out$iterations
    final_innovations <- out$innovations

    # Extract AR and MA coefficients
    if(model_params$ar_order > 0) {
      ar_est <- b_est[1:model_params$ar_order]
    } else {
      ar_est <- numeric(0)
    }

    if(model_params$ma_order > 0) {
      ma_est <- b_est[(model_params$ar_order+1):(model_params$ar_order+model_params$ma_order)]
    } else {
      ma_est <- numeric(0)
    }

    if(verbose) {
      cat("PMM2 algorithm finished.\n")
      cat("  Converged:", conv, "\n")
      cat("  Iterations:", iter, "\n")
      if(model_params$ar_order > 0) cat("  Final AR parameters:", ar_est, "\n")
      if(model_params$ma_order > 0) cat("  Final MA parameters:", ma_est, "\n")
    }
  } else {
    # For non-PMM2 methods, just use the initial estimates
    ar_est <- if(model_params$ar_order > 0) b_init[1:model_params$ar_order] else numeric(0)
    ma_est <- if(model_params$ma_order > 0) b_init[(model_params$ar_order+1):length(b_init)] else numeric(0)
    conv <- TRUE
    iter <- 0
    final_innovations <- innovations
  }

  # Ensure all values are proper types
  final_innovations <- as.numeric(final_innovations)

  # Return S4 object with results
  ans <- new("TS2fit",
             coefficients = as.numeric(c(ar_est, ma_est)),
             residuals = as.numeric(final_innovations),
             m2 = as.numeric(m2),
             m3 = as.numeric(m3),
             m4 = as.numeric(m4),
             convergence = as.logical(conv),
             iterations = as.numeric(iter),
             call = call,
             model_type = as.character(model_type),
             intercept = as.numeric(x_mean),
             original_series = as.numeric(model_params$original_x),
             order = model_params$order)

  return(ans)
}

#' Fit an AR model using PMM2
#'
#' @inheritParams ts_pmm2
#' @export
ar_pmm2 <- function(x, order = 1, method = "pmm2", max_iter = 50, tol = 1e-6,
                    include.mean = TRUE, initial = NULL, na.action = na.fail,
                    verbose = FALSE) {
  ts_pmm2(x, order = order, model_type = "ar", method = method, max_iter = max_iter,
          tol = tol, include.mean = include.mean, initial = initial,
          na.action = na.action, verbose = verbose)
}

#' Fit a MA model using PMM2
#'
#' @inheritParams ts_pmm2
#' @export
ma_pmm2 <- function(x, order = 1, method = "pmm2", max_iter = 50, tol = 1e-6,
                    include.mean = TRUE, initial = NULL, na.action = na.fail,
                    verbose = FALSE) {
  ts_pmm2(x, order = order, model_type = "ma", method = method, max_iter = max_iter,
          tol = tol, include.mean = include.mean, initial = initial,
          na.action = na.action, verbose = verbose)
}

#' Fit an ARMA model using PMM2
#'
#' @inheritParams ts_pmm2
#' @export
arma_pmm2 <- function(x, order = c(1, 1), method = "pmm2", max_iter = 50, tol = 1e-6,
                      include.mean = TRUE, initial = NULL, na.action = na.fail,
                      verbose = FALSE) {
  ts_pmm2(x, order = order, model_type = "arma", method = method, max_iter = max_iter,
          tol = tol, include.mean = include.mean, initial = initial,
          na.action = na.action, verbose = verbose)
}

#' Fit an ARIMA model using PMM2
#'
#' @inheritParams ts_pmm2
#' @export
arima_pmm2 <- function(x, order = c(1, 1, 1), method = "pmm2", max_iter = 50, tol = 1e-6,
                       include.mean = TRUE, initial = NULL, na.action = na.fail,
                       verbose = FALSE) {
  ts_pmm2(x, order = order, model_type = "arima", method = method, max_iter = max_iter,
          tol = tol, include.mean = include.mean, initial = initial,
          na.action = na.action, verbose = verbose)
}
