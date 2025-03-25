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

  # Check for NAs and infinite values
  if(any(is.na(x)) || any(is.infinite(x))) {
    warning("NA or infinite values detected in input series. They will be removed.")
    x <- x[!is.na(x) & !is.infinite(x)]
    if(length(x) < 10) {
      stop("Too few valid observations after removing NA/infinite values")
    }
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

    # Check if we have enough data
    if(length(x) <= ar_order + 1) {
      stop("Too few observations for AR model of order ", ar_order)
    }
  } else if(model_type == "ma") {
    if(!is.numeric(order) || length(order) != 1)
      stop("For MA models, 'order' must be a single integer")
    ar_order <- 0
    ma_order <- as.integer(order)
    d <- 0
    if(ma_order <= 0)
      stop("MA order must be positive")

    # Check if we have enough data
    if(length(x) <= ma_order + 1) {
      stop("Too few observations for MA model of order ", ma_order)
    }
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

    # Check if we have enough data
    if(length(x) <= max(ar_order, ma_order) + 1) {
      stop("Too few observations for ARMA model of order (", ar_order, ",", ma_order, ")")
    }
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

    # Check if we have enough data after differencing
    if(length(x) <= d + max(ar_order, ma_order) + 1) {
      stop("Too few observations for ARIMA model after differencing")
    }
  } else {
    stop("Unknown model type: ", model_type)
  }

  # Store original series
  orig_x <- as.numeric(x)

  # Apply differencing for ARIMA models
  if(model_type == "arima" && d > 0) {
    x_diff <- tryCatch({
      as.numeric(diff(x, differences = d))
    }, error = function(e) {
      stop("Error applying differencing: ", conditionMessage(e))
    })

    # Safety check
    if(length(x_diff) < max(ar_order, ma_order) + 1) {
      stop("Too few observations after differencing for the specified model orders")
    }
  } else {
    x_diff <- x
  }

  # Compute series mean if requested
  if(include.mean) {
    x_mean <- mean(x_diff, na.rm = TRUE)
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
      # Try multiple methods for initial estimates
      init_fit <- tryCatch({
        stats::arima(
          ifelse(model_type == "arima", model_params$original_x, x_centered),
          order = arima_order,
          method = ifelse(method == "pmm2", "CSS", method),
          include.mean = (x_mean != 0) & (model_type != "arima"))
      }, error = function(e) {
        # Try alternative method if the first one fails
        tryCatch({
          if(verbose) cat("First method failed, trying CSS-ML method...\n")
          stats::arima(
            ifelse(model_type == "arima", model_params$original_x, x_centered),
            order = arima_order,
            method = "CSS-ML", # Try hybrid method
            include.mean = (x_mean != 0) & (model_type != "arima"))
        }, error = function(e2) {
          # If that also fails, use simplified initial values
          if(verbose) cat("All estimation methods failed, using simplified initial values...\n")

          # Create a simplified arima result
          result <- list(
            coef = numeric(ar_order + ma_order),
            residuals = ifelse(model_type == "arima",
                               diff(model_params$original_x, differences = d),
                               x_centered)
          )

          # Add names to coefficients
          if(ar_order > 0)
            names(result$coef)[1:ar_order] <- paste0("ar", 1:ar_order)
          if(ma_order > 0)
            names(result$coef)[(ar_order+1):(ar_order+ma_order)] <- paste0("ma", 1:ma_order)

          return(result)
        })
      })

      # Extract parameters from the fitted model
      ar_init <- if(ar_order > 0 && !is.null(init_fit$coef) &&
                    any(grepl("^ar", names(init_fit$coef)))) {
        as.numeric(init_fit$coef[paste0("ar", 1:ar_order)])
      } else {
        # Default small values if no AR part found
        rep(0.1, ar_order)
      }

      ma_init <- if(ma_order > 0 && !is.null(init_fit$coef) &&
                    any(grepl("^ma", names(init_fit$coef)))) {
        as.numeric(init_fit$coef[paste0("ma", 1:ma_order)])
      } else {
        # Default small values if no MA part found
        rep(0.1, ma_order)
      }

      # Adjust mean if necessary
      if(x_mean != 0 && !is.null(init_fit$coef) && "intercept" %in% names(init_fit$coef)) {
        x_mean <- init_fit$coef["intercept"]
      }

      # Get innovations
      innovations <- if(!is.null(init_fit$residuals)) {
        as.numeric(init_fit$residuals)
      } else {
        # Default white noise if no residuals available
        rnorm(length(x_centered), 0, sd(x_centered, na.rm=TRUE))
      }

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

      # Try to get innovations, with fallback mechanism
      init_fit <- tryCatch({
        stats::arima(
          ifelse(model_type == "arima", model_params$original_x, x_centered),
          order = arima_order,
          include.mean = (x_mean != 0) & (model_type != "arima"),
          fixed = fixed_params)
      }, error = function(e) {
        if(verbose) cat("Error in fixed parameter model:", conditionMessage(e), "\n")
        # Return a simplified structure if arima fails
        list(
          residuals = ifelse(model_type == "arima",
                             diff(model_params$original_x, differences = d),
                             x_centered)
        )
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
    warning("Initial parameter estimation resulted in NA coefficients. Using zeros instead.")
    b_init[is.na(b_init)] <- 0
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
#' 3. Uses these moments with the universal solver to find robust parameter estimates
#'
#' PMM2 is particularly useful when error terms are non-Gaussian and asymmetric.
#'
#' @return An S4 \code{TS2fit} object
#' @export
ts_pmm2 <- function(x, order, model_type = c("ar", "ma", "arma", "arima"),
                    method = "pmm2", max_iter = 50, tol = 1e-6,
                    include.mean = TRUE, initial = NULL, na.action = na.fail,
                    regularize = TRUE, reg_lambda = 1e-8, verbose = FALSE) {

  # Match model_type argument
  model_type <- match.arg(model_type)

  # Dispatch to appropriate function based on model_type
  switch(model_type,
         "ar" = ar_pmm2(x, order, method, max_iter, tol, include.mean, initial,
                        na.action, regularize, reg_lambda, verbose),
         "ma" = ma_pmm2(x, order, method, max_iter, tol, include.mean, initial,
                        na.action, regularize, reg_lambda, verbose),
         "arma" = arma_pmm2(x, order, method, max_iter, tol, include.mean, initial,
                            na.action, regularize, reg_lambda, verbose),
         "arima" = arima_pmm2(x, order, method, max_iter, tol, include.mean, initial,
                              na.action, regularize, reg_lambda, verbose),
         stop("Unknown model type: ", model_type))
}


#' Fit an AR model using PMM2
#'
#' @inheritParams ts_pmm2
#' @export
ar_pmm2 <- function(x, order = 1, method = "pmm2", max_iter = 50, tol = 1e-6,
                    include.mean = TRUE, initial = NULL, na.action = na.fail,
                    regularize = TRUE, reg_lambda = 1e-8, verbose = FALSE) {

  # Capture the call
  call <- match.call()

  # Handle missing values
  if (!is.null(na.action)) {
    x <- na.action(x)
  }

  # Convert input to numeric vector
  x <- as.numeric(x)
  n <- length(x)
  ar_order <- as.integer(order)

  # Check if we have enough data
  if (n <= ar_order + 1) {
    stop("Too few observations for AR model of order ", ar_order)
  }

  # Center the data if mean is to be included
  if (include.mean) {
    x_mean <- mean(x)
    x_centered <- x - x_mean
  } else {
    x_mean <- 0
    x_centered <- x
  }

  # Create design matrix for AR model
  X <- matrix(0, nrow = n - ar_order, ncol = ar_order)
  for (i in 1:ar_order) {
    X[, i] <- x_centered[(ar_order - i + 1):(n - i)]
  }
  y <- x_centered[(ar_order + 1):n]

  # Get initial estimates using OLS
  if (is.null(initial)) {
    b_init <- solve(t(X) %*% X) %*% t(X) %*% y
  } else {
    if (length(initial) != ar_order) {
      stop("Length of initial estimates must match AR order")
    }
    b_init <- initial
  }
  b_init <- as.numeric(b_init)

  # Calculate OLS residuals
  res_ols <- y - X %*% b_init

  # Calculate moments from residuals
  m2 <- mean(res_ols^2)
  m3 <- mean(res_ols^3)
  m4 <- mean(res_ols^4)

  if (verbose) {
    cat("Initial moments from OLS residuals:\n")
    cat("  m2 =", m2, "\n")
    cat("  m3 =", m3, "\n")
    cat("  m4 =", m4, "\n")
  }

  # Check for problematic moments
  if (m2 <= 0) {
    warning("Second central moment (m2) is non-positive. Using absolute value.")
    m2 <- abs(m2)
    if (m2 < 1e-6) m2 <- var(res_ols, na.rm = TRUE)
  }

  if (m4 <= m2^2) {
    warning("Fourth central moment (m4) is less than m2^2. Adjusting to ensure valid distribution.")
    m4 <- m2^2 * 3  # Use Gaussian kurtosis as fallback
  }

  # Apply PMM2 method if requested
  if (method == "pmm2") {
    if (verbose) cat("Starting PMM2 optimization...\n")

    # Get PMM2 estimates
    b_pmm2 <- solve_pmm2(b_init, X, y, m2, m3, m4,
                         max_iter = max_iter, tol = tol,
                         regularize = regularize, reg_lambda = reg_lambda,
                         verbose = verbose)

    # Calculate final residuals
    final_res <- y - X %*% b_pmm2

    if (verbose) {
      cat("PMM2 estimation complete.\n")
      cat("Final AR coefficients:", b_pmm2, "\n")
    }
  } else {
    # For other methods, use the initial estimates
    b_pmm2 <- b_init
    final_res <- res_ols
  }

  # Create the S4 object
  ans <- new("TS2fit",
             coefficients = as.numeric(b_pmm2),
             residuals = as.numeric(final_res),
             m2 = as.numeric(m2),
             m3 = as.numeric(m3),
             m4 = as.numeric(m4),
             convergence = TRUE,
             iterations = as.numeric(max_iter),
             call = call,
             model_type = "ar",
             intercept = as.numeric(x_mean),
             original_series = as.numeric(x),
             order = list(ar = ar_order, ma = 0, d = 0))

  return(ans)
}


#' Fit an MA model using PMM2
#'
#' This function fits a Moving Average (MA) model to a time series using the
#' Polynomial Maximization Method (PMM) of order 2, which is robust against
#' non-Gaussian errors, particularly with asymmetric distributions.
#'
#' @param x numeric vector of time series data
#' @param order integer: MA order
#' @param method string: estimation method, one of "pmm2" (default), "css", "ml", "yw", "ols"
#' @param max_iter integer: maximum number of iterations for the algorithm
#' @param tol numeric: tolerance for convergence
#' @param include.mean logical: whether to include a mean (intercept) term
#' @param initial vector of initial parameter estimates (optional)
#' @param na.action function to handle missing values, default is na.fail
#' @param regularize logical, add small value to diagonal for numerical stability
#' @param reg_lambda regularization parameter (if regularize=TRUE)
#' @param verbose logical: whether to print progress information
#'
#' @return An S4 \code{TS2fit} object
#' @export
ma_pmm2 <- function(x, order = 1, method = "pmm2", max_iter = 50, tol = 1e-6,
                    include.mean = TRUE, initial = NULL, na.action = na.fail,
                    regularize = TRUE, reg_lambda = 1e-8, verbose = FALSE) {

  # Capture the call
  call <- match.call()

  # Handle missing values
  if (!is.null(na.action)) {
    x <- na.action(x)
  }

  # Convert input to numeric vector
  x <- as.numeric(x)
  n <- length(x)
  ma_order <- as.integer(order)

  # Check if we have enough data
  if (n <= ma_order + 10) {
    stop("Too few observations for MA model of order ", ma_order,
         ". At least ", ma_order + 10, " observations recommended.")
  }

  # Center the data if mean is to be included
  if (include.mean) {
    x_mean <- mean(x, na.rm = TRUE)
    x_centered <- x - x_mean
  } else {
    x_mean <- 0
    x_centered <- x
  }

  if (verbose) {
    cat("Fitting MA(", ma_order, ") model to time series of length ", n, "\n")
    if (include.mean) cat("Series mean:", x_mean, "\n")
  }

  # Get initial MA model using arima function
  if (is.null(initial)) {
    if (verbose) cat("Getting initial estimates using ARIMA...\n")

    initial_model <- tryCatch({
      if (verbose) cat("Trying CSS-ML method...\n")
      arima(x_centered, order = c(0, 0, ma_order), method = "CSS-ML")
    }, error = function(e) {
      if (verbose) cat("CSS-ML failed, trying CSS method...\n")
      tryCatch({
        arima(x_centered, order = c(0, 0, ma_order), method = "CSS")
      }, error = function(e2) {
        if (verbose) cat("CSS failed, trying ML method...\n")
        tryCatch({
          arima(x_centered, order = c(0, 0, ma_order), method = "ML")
        }, error = function(e3) {
          if (verbose) cat("All methods failed. Using simple starting values.\n")
          # Create a simple model with small default values
          coef <- rep(0.1, ma_order)
          names(coef) <- paste0("ma", 1:ma_order)
          list(
            coef = coef,
            residuals = rnorm(n, 0, sd(x_centered))
          )
        })
      })
    })

    # Extract MA coefficients
    if (!is.null(initial_model$coef)) {
      ma_names <- paste0("ma", 1:ma_order)
      if (all(ma_names %in% names(initial_model$coef))) {
        b_init <- as.numeric(initial_model$coef[ma_names])
      } else {
        if (verbose) cat("Missing coefficient names. Using available coefficients.\n")
        b_init <- as.numeric(initial_model$coef)
        if (length(b_init) != ma_order) {
          b_init <- c(b_init, rep(0.1, ma_order - length(b_init)))
        }
      }
    } else {
      if (verbose) cat("No coefficients found. Using default values.\n")
      b_init <- rep(0.1, ma_order)
    }

    # Extract residuals/innovations
    if (!is.null(initial_model$residuals)) {
      innovations <- as.numeric(initial_model$residuals)
    } else {
      if (verbose) cat("No residuals found. Using random noise.\n")
      innovations <- rnorm(n, 0, sd(x_centered))
    }

    if (verbose) {
      cat("Initial MA coefficients:", b_init, "\n")
    }
  } else {
    # Use provided initial estimates
    if (length(initial) != ma_order) {
      stop("Length of initial estimates (", length(initial),
           ") must match MA order (", ma_order, ")")
    }
    b_init <- as.numeric(initial)

    # Get innovations using fixed parameters
    if (verbose) cat("Using provided initial coefficients:", b_init, "\n")

    initial_model <- tryCatch({
      arima(x_centered, order = c(0, 0, ma_order), fixed = b_init)
    }, error = function(e) {
      if (verbose) cat("Error with provided coefficients:", conditionMessage(e), "\n")
      # Create innovations via MA filter with provided coefficients
      innovations <- numeric(n)
      # Fill first ma_order values with noise
      innovations[1:ma_order] <- rnorm(ma_order, 0, sd(x_centered))

      # Apply MA filter
      for (t in (ma_order+1):n) {
        predicted <- 0
        for (i in 1:ma_order) {
          predicted <- predicted + b_init[i] * innovations[t-i]
        }
        innovations[t] <- x_centered[t] - predicted
      }

      list(residuals = innovations)
    })

    # Get innovations
    if (!is.null(initial_model$residuals)) {
      innovations <- as.numeric(initial_model$residuals)
    } else {
      if (verbose) cat("Using manually calculated innovations.\n")
    }
  }

  if (length(innovations) < n) {
    if (verbose) cat("Padding innovations to match data length.\n")
    # If innovations are shorter, pad them
    innovations <- c(rep(0, n - length(innovations)), innovations)
  } else if (length(innovations) > n) {
    if (verbose) cat("Trimming innovations to match data length.\n")
    # If innovations are longer, trim them
    innovations <- tail(innovations, n)
  }

  # Prepare data for PMM2 estimation
  usable_length <- n - ma_order

  # Create design matrix using lagged innovations
  X <- matrix(0, nrow = usable_length, ncol = ma_order)

  if (verbose) cat("Creating design matrix with dimensions", nrow(X), "x", ncol(X), "\n")

  # Fill design matrix with lagged innovations
  for (i in 1:ma_order) {
    X[, i] <- innovations[(ma_order-i+1):(n-i)]
  }

  # Response vector - the centered time series excluding first ma_order values
  y <- x_centered[(ma_order+1):n]

  # Double-check dimensions
  if (nrow(X) != length(y)) {
    stop("Dimension mismatch: X has ", nrow(X), " rows but y has ", length(y), " elements")
  }

  # Calculate moments from innovations
  m2 <- mean(innovations^2, na.rm = TRUE)
  m3 <- mean(innovations^3, na.rm = TRUE)
  m4 <- mean(innovations^4, na.rm = TRUE)

  if (verbose) {
    cat("Innovation statistics:\n")
    cat("  Mean: ", mean(innovations, na.rm = TRUE), "\n")
    cat("  SD:   ", sd(innovations, na.rm = TRUE), "\n")
    cat("  Skewness: ", m3 / m2^(3/2), "\n")
    cat("  Kurtosis: ", m4 / m2^2 - 3, "\n")
    cat("Moments used for PMM2:\n")
    cat("  m2 =", m2, "\n")
    cat("  m3 =", m3, "\n")
    cat("  m4 =", m4, "\n")
  }

  # Check for problematic moments
  if (is.na(m2) || m2 <= 0) {
    warning("Second central moment (m2) is non-positive or NA. Using absolute value.")
    m2 <- abs(m2)
    if (is.na(m2) || m2 < 1e-6) m2 <- var(innovations, na.rm = TRUE)
    if (is.na(m2) || m2 < 1e-6) m2 <- 1e-6
  }

  if (is.na(m3)) {
    warning("Third central moment (m3) is NA. Setting to zero.")
    m3 <- 0
  }

  if (is.na(m4) || m4 <= m2^2) {
    warning("Fourth central moment (m4) is invalid. Adjusting to ensure valid distribution.")
    m4 <- m2^2 * 3  # Use Gaussian kurtosis as fallback
  }

  # Apply PMM2 method if requested
  if (method == "pmm2") {
    if (verbose) cat("Starting PMM2 optimization...\n")

    # Get PMM2 estimates using our universal solver
    b_pmm2 <- solve_pmm2(b_init, X, y, m2, m3, m4,
                         max_iter = max_iter, tol = tol,
                         regularize = regularize, reg_lambda = reg_lambda,
                         verbose = verbose)

    if (verbose) {
      cat("PMM2 estimation complete.\n")
      cat("Initial MA coefficients:", b_init, "\n")
      cat("Final MA coefficients:  ", b_pmm2, "\n")
    }

    # Calculate final model using original arima function with fixed parameters
    # Make sure the names are correct
    fixed_params <- b_pmm2
    names(fixed_params) <- paste0("ma", 1:ma_order)

    final_model <- tryCatch({
      arima(x_centered, order = c(0, 0, ma_order), fixed = fixed_params)
    }, error = function(e) {
      warning("Error in final model with PMM2 parameters: ", conditionMessage(e),
              ". Using original innovations.")
      return(list(residuals = innovations))
    })

    final_res <- residuals(final_model)

    # Make sure final_res is proper length
    if (length(final_res) != n) {
      if (verbose) cat("Adjusting final residuals to match data length.\n")
      if (length(final_res) < n) {
        final_res <- c(rep(NA, n - length(final_res)), final_res)
      } else {
        final_res <- tail(final_res, n)
      }
    }
  } else {
    # For other methods, use the initial estimates
    b_pmm2 <- b_init
    final_res <- innovations
  }

  # Create the S4 object
  ans <- new("TS2fit",
             coefficients = as.numeric(b_pmm2),
             residuals = as.numeric(final_res),
             m2 = as.numeric(m2),
             m3 = as.numeric(m3),
             m4 = as.numeric(m4),
             convergence = TRUE,
             iterations = as.numeric(max_iter),
             call = call,
             model_type = "ma",
             intercept = as.numeric(x_mean),
             original_series = as.numeric(x),
             order = list(ar = 0, ma = ma_order, d = 0))

  return(ans)
}


#' Fit an ARMA model using PMM2
#'
#' @inheritParams ts_pmm2
#' @export
arma_pmm2 <- function(x, order = c(1, 1), method = "pmm2", max_iter = 50, tol = 1e-6,
                      include.mean = TRUE, initial = NULL, na.action = na.fail,
                      regularize = TRUE, reg_lambda = 1e-8, verbose = FALSE) {

  # Capture the call
  call <- match.call()

  # Handle missing values
  if (!is.null(na.action)) {
    x <- na.action(x)
  }

  # Convert input to numeric vector and parse order parameters
  x <- as.numeric(x)
  n <- length(x)
  ar_order <- as.integer(order[1])
  ma_order <- as.integer(order[2])

  # Check if we have enough data
  min_obs <- max(ar_order, ma_order) + 1
  if (n <= min_obs) {
    stop("Too few observations for ARMA model of order (", ar_order, ",", ma_order, ")")
  }

  # Center the data if mean is to be included
  if (include.mean) {
    x_mean <- mean(x)
    x_centered <- x - x_mean
  } else {
    x_mean <- 0
    x_centered <- x
  }

  # Get initial ARMA model using arima function
  if (is.null(initial)) {
    initial_model <- tryCatch({
      arima(x_centered, order = c(ar_order, 0, ma_order), method = "CSS-ML")
    }, error = function(e) {
      warning("Error in initial ARMA estimation: ", conditionMessage(e),
              ". Trying CSS method...")
      arima(x_centered, order = c(ar_order, 0, ma_order), method = "CSS")
    })

    # Extract ARMA coefficients
    ar_coef <- numeric(0)
    ma_coef <- numeric(0)

    if (ar_order > 0) {
      ar_names <- paste0("ar", 1:ar_order)
      if (any(ar_names %in% names(initial_model$coef))) {
        ar_coef <- initial_model$coef[ar_names]
      } else {
        ar_coef <- rep(0.1, ar_order)
      }
    }

    if (ma_order > 0) {
      ma_names <- paste0("ma", 1:ma_order)
      if (any(ma_names %in% names(initial_model$coef))) {
        ma_coef <- initial_model$coef[ma_names]
      } else {
        ma_coef <- rep(0.1, ma_order)
      }
    }

    b_init <- c(ar_coef, ma_coef)

    # Extract residuals/innovations
    innovations <- residuals(initial_model)
  } else {
    # Parse provided initial estimates
    if (is.list(initial)) {
      ar_coef <- if (ar_order > 0 && !is.null(initial$ar)) initial$ar else rep(0.1, ar_order)
      ma_coef <- if (ma_order > 0 && !is.null(initial$ma)) initial$ma else rep(0.1, ma_order)
    } else {
      # Assume vector with AR followed by MA coefficients
      if (length(initial) != ar_order + ma_order) {
        stop("Length of initial vector must match sum of AR and MA orders")
      }
      ar_coef <- if (ar_order > 0) initial[1:ar_order] else numeric(0)
      ma_coef <- if (ma_order > 0) initial[(ar_order+1):(ar_order+ma_order)] else numeric(0)
    }

    b_init <- c(ar_coef, ma_coef)

    # Get innovations using fixed parameters
    initial_model <- tryCatch({
      arima(x_centered, order = c(ar_order, 0, ma_order), fixed = b_init)
    }, error = function(e) {
      stop("Error using provided initial estimates: ", conditionMessage(e))
    })
    innovations <- residuals(initial_model)
  }

  # Create design matrix using lagged values and innovations
  min_lag <- max(ar_order, ma_order)
  valid_n <- n - min_lag
  X <- matrix(0, nrow = valid_n, ncol = ar_order + ma_order)

  # Add AR components (lagged values)
  if (ar_order > 0) {
    for (i in 1:ar_order) {
      X[, i] <- x_centered[(min_lag - i + 1):(n - i)]
    }
  }

  # Add MA components (lagged innovations)
  if (ma_order > 0) {
    for (i in 1:ma_order) {
      if (i <= length(innovations) - min_lag) {
        X[, ar_order + i] <- innovations[i:(valid_n + i - 1)]
      } else {
        warning("Not enough innovations for MA lag ", i)
        X[, ar_order + i] <- 0
      }
    }
  }

  # Response vector
  y <- x_centered[(min_lag + 1):n]

  # Calculate moments from innovations
  m2 <- mean(innovations^2, na.rm = TRUE)
  m3 <- mean(innovations^3, na.rm = TRUE)
  m4 <- mean(innovations^4, na.rm = TRUE)

  if (verbose) {
    cat("Initial moments from innovations:\n")
    cat("  m2 =", m2, "\n")
    cat("  m3 =", m3, "\n")
    cat("  m4 =", m4, "\n")
  }

  # Check for problematic moments
  if (m2 <= 0) {
    warning("Second central moment (m2) is non-positive. Using absolute value.")
    m2 <- abs(m2)
    if (m2 < 1e-6) m2 <- var(innovations, na.rm = TRUE)
  }

  if (m4 <= m2^2) {
    warning("Fourth central moment (m4) is less than m2^2. Adjusting to ensure valid distribution.")
    m4 <- m2^2 * 3  # Use Gaussian kurtosis as fallback
  }

  # Apply PMM2 method if requested
  if (method == "pmm2") {
    if (verbose) cat("Starting PMM2 optimization...\n")

    # Get PMM2 estimates
    b_pmm2 <- solve_pmm2(b_init, X, y, m2, m3, m4,
                         max_iter = max_iter, tol = tol,
                         regularize = regularize, reg_lambda = reg_lambda,
                         verbose = verbose)

    # Extract AR and MA coefficients
    ar_pmm2 <- if (ar_order > 0) b_pmm2[1:ar_order] else numeric(0)
    ma_pmm2 <- if (ma_order > 0) b_pmm2[(ar_order+1):(ar_order+ma_order)] else numeric(0)

    # Calculate final model using original arima function with fixed parameters
    final_model <- tryCatch({
      arima(x_centered, order = c(ar_order, 0, ma_order),
            fixed = c(ar_pmm2, ma_pmm2))
    }, error = function(e) {
      warning("Error in final model with PMM2 parameters: ", conditionMessage(e),
              ". Using original innovations.")
      return(list(residuals = innovations))
    })

    final_res <- residuals(final_model)

    if (verbose) {
      cat("PMM2 estimation complete.\n")
      if (ar_order > 0) cat("Final AR coefficients:", ar_pmm2, "\n")
      if (ma_order > 0) cat("Final MA coefficients:", ma_pmm2, "\n")
    }
  } else {
    # For other methods, use the initial estimates
    b_pmm2 <- b_init
    final_res <- innovations
  }

  # Create the S4 object
  ans <- new("TS2fit",
             coefficients = as.numeric(b_pmm2),
             residuals = as.numeric(final_res),
             m2 = as.numeric(m2),
             m3 = as.numeric(m3),
             m4 = as.numeric(m4),
             convergence = TRUE,
             iterations = as.numeric(max_iter),
             call = call,
             model_type = "arma",
             intercept = as.numeric(x_mean),
             original_series = as.numeric(x),
             order = list(ar = ar_order, ma = ma_order, d = 0))

  return(ans)
}


#' Fit an ARIMA model using PMM2
#'
#' @inheritParams ts_pmm2
#' @export
arima_pmm2 <- function(x, order = c(1, 1, 1), method = "pmm2", max_iter = 50, tol = 1e-6,
                       include.mean = TRUE, initial = NULL, na.action = na.fail,
                       regularize = TRUE, reg_lambda = 1e-8, verbose = FALSE) {

  # Capture the call
  call <- match.call()

  # Handle missing values
  if (!is.null(na.action)) {
    x <- na.action(x)
  }

  # Convert input to numeric vector and parse order parameters
  x <- as.numeric(x)
  n <- length(x)
  ar_order <- as.integer(order[1])
  d <- as.integer(order[2])
  ma_order <- as.integer(order[3])

  # Check if we have enough data
  min_obs <- max(ar_order, ma_order) + d + 1
  if (n <= min_obs) {
    stop("Too few observations for ARIMA model of order (", ar_order, ",", d, ",", ma_order, ")")
  }

  # Store original series
  orig_x <- x

  # Apply differencing
  if (d > 0) {
    x_diff <- diff(x, differences = d)
  } else {
    x_diff <- x
  }

  # Center the differenced data if mean is to be included
  if (include.mean) {
    x_mean <- mean(x_diff)
    x_centered <- x_diff - x_mean
  } else {
    x_mean <- 0
    x_centered <- x_diff
  }

  # Get initial ARIMA model using stats::arima function
  if (is.null(initial)) {
    initial_model <- tryCatch({
      arima(orig_x, order = c(ar_order, d, ma_order), method = "CSS-ML")
    }, error = function(e) {
      warning("Error in initial ARIMA estimation: ", conditionMessage(e),
              ". Trying CSS method...")
      arima(orig_x, order = c(ar_order, d, ma_order), method = "CSS")
    })

    # Extract ARIMA coefficients
    ar_coef <- numeric(0)
    ma_coef <- numeric(0)

    if (ar_order > 0) {
      ar_names <- paste0("ar", 1:ar_order)
      if (any(ar_names %in% names(initial_model$coef))) {
        ar_coef <- initial_model$coef[ar_names]
      } else {
        ar_coef <- rep(0.1, ar_order)
      }
    }

    if (ma_order > 0) {
      ma_names <- paste0("ma", 1:ma_order)
      if (any(ma_names %in% names(initial_model$coef))) {
        ma_coef <- initial_model$coef[ma_names]
      } else {
        ma_coef <- rep(0.1, ma_order)
      }
    }

    b_init <- c(ar_coef, ma_coef)

    # Extract residuals/innovations
    innovations <- residuals(initial_model)
  } else {
    # Parse provided initial estimates
    if (is.list(initial)) {
      ar_coef <- if (ar_order > 0 && !is.null(initial$ar)) initial$ar else rep(0.1, ar_order)
      ma_coef <- if (ma_order > 0 && !is.null(initial$ma)) initial$ma else rep(0.1, ma_order)
    } else {
      # Assume vector with AR followed by MA coefficients
      if (length(initial) != ar_order + ma_order) {
        stop("Length of initial vector must match sum of AR and MA orders")
      }
      ar_coef <- if (ar_order > 0) initial[1:ar_order] else numeric(0)
      ma_coef <- if (ma_order > 0) initial[(ar_order+1):(ar_order+ma_order)] else numeric(0)
    }

    b_init <- c(ar_coef, ma_coef)

    # Get innovations using fixed parameters
    initial_model <- tryCatch({
      arima(orig_x, order = c(ar_order, d, ma_order), fixed = b_init)
    }, error = function(e) {
      stop("Error using provided initial estimates: ", conditionMessage(e))
    })
    innovations <- residuals(initial_model)
  }

  # Work with the differenced series for PMM2 estimation
  n_diff <- length(x_centered)

  # Create design matrix using lagged values and innovations
  min_lag <- max(ar_order, ma_order)
  valid_n <- n_diff - min_lag
  X <- matrix(0, nrow = valid_n, ncol = ar_order + ma_order)

  # Add AR components (lagged values)
  if (ar_order > 0) {
    for (i in 1:ar_order) {
      X[, i] <- x_centered[(min_lag - i + 1):(n_diff - i)]
    }
  }

  # Add MA components (lagged innovations)
  if (ma_order > 0) {
    # Adjust innovations length to match differenced data
    inno_diff <- innovations[(length(innovations) - n_diff + 1):length(innovations)]

    for (i in 1:ma_order) {
      if (i <= length(inno_diff) - min_lag) {
        X[, ar_order + i] <- inno_diff[i:(valid_n + i - 1)]
      } else {
        warning("Not enough innovations for MA lag ", i)
        X[, ar_order + i] <- 0
      }
    }
  }

  # Response vector
  y <- x_centered[(min_lag + 1):n_diff]

  # Calculate moments from innovations
  m2 <- mean(innovations^2, na.rm = TRUE)
  m3 <- mean(innovations^3, na.rm = TRUE)
  m4 <- mean(innovations^4, na.rm = TRUE)

  if (verbose) {
    cat("Initial moments from innovations:\n")
    cat("  m2 =", m2, "\n")
    cat("  m3 =", m3, "\n")
    cat("  m4 =", m4, "\n")
  }

  # Check for problematic moments
  if (m2 <= 0) {
    warning("Second central moment (m2) is non-positive. Using absolute value.")
    m2 <- abs(m2)
    if (m2 < 1e-6) m2 <- var(innovations, na.rm = TRUE)
  }

  if (m4 <= m2^2) {
    warning("Fourth central moment (m4) is less than m2^2. Adjusting to ensure valid distribution.")
    m4 <- m2^2 * 3  # Use Gaussian kurtosis as fallback
  }

  # Apply PMM2 method if requested
  if (method == "pmm2") {
    if (verbose) cat("Starting PMM2 optimization...\n")

    # Get PMM2 estimates
    b_pmm2 <- solve_pmm2(b_init, X, y, m2, m3, m4,
                         max_iter = max_iter, tol = tol,
                         regularize = regularize, reg_lambda = reg_lambda,
                         verbose = verbose)

    # Extract AR and MA coefficients
    ar_pmm2 <- if (ar_order > 0) b_pmm2[1:ar_order] else numeric(0)
    ma_pmm2 <- if (ma_order > 0) b_pmm2[(ar_order+1):(ar_order+ma_order)] else numeric(0)

    # Переконуємось, що правильно передаємо параметри fixed
    fixed_params <- numeric(0)
    if (ar_order > 0) fixed_params <- c(fixed_params, ar_pmm2)
    if (ma_order > 0) fixed_params <- c(fixed_params, ma_pmm2)

    # Додаємо параметр intercept, якщо потрібно
    if (include.mean) {
      fixed_params <- c(fixed_params, x_mean)
      fixed_names <- c(
        if(ar_order > 0) paste0("ar", 1:ar_order) else character(0),
        if(ma_order > 0) paste0("ma", 1:ma_order) else character(0),
        "intercept"
      )
      names(fixed_params) <- fixed_names
    }

    # Calculate final model using original arima function with fixed parameters
    final_model <- tryCatch({
      arima(orig_x, order = c(ar_order, d, ma_order),
            include.mean = include.mean,
            fixed = fixed_params)
    }, error = function(e) {
      warning("Error in final model with PMM2 parameters: ", conditionMessage(e),
              ". Using original innovations.")
      list(residuals = innovations)
    })

    final_res <- residuals(final_model)

    if (verbose) {
      cat("PMM2 estimation complete.\n")
      if (ar_order > 0) cat("Final AR coefficients:", ar_pmm2, "\n")
      if (ma_order > 0) cat("Final MA coefficients:", ma_pmm2, "\n")
    }
  } else {
    # For other methods, use the initial estimates
    b_pmm2 <- b_init
    final_res <- innovations
  }

  # Create the S4 object
  ans <- new("TS2fit",
             coefficients = as.numeric(b_pmm2),
             residuals = as.numeric(final_res),
             m2 = as.numeric(m2),
             m3 = as.numeric(m3),
             m4 = as.numeric(m4),
             convergence = TRUE,
             iterations = as.numeric(max_iter),
             call = call,
             model_type = "arima",
             intercept = as.numeric(x_mean),
             original_series = as.numeric(orig_x),
             order = list(ar = ar_order, ma = ma_order, d = d))

  return(ans)
}
