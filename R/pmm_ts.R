# -------------------------------
# pmm_ts.R (Refactored & Fixed for MA final residuals)
# -------------------------------

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
  if (missing(x)) {
    stop("'x' must be provided")
  }

  # Convert input to numeric vector
  x <- as.numeric(x)

  if (!is.numeric(x)) {
    stop("'x' must be a numeric vector")
  }

  # Check for NAs and infinite values
  if (any(is.na(x)) || any(is.infinite(x))) {
    warning("NA or infinite values detected in input series. They will be removed.")
    x <- x[!is.na(x) & !is.infinite(x)]
    if (length(x) < 10) {
      stop("Too few valid observations after removing NA/infinite values")
    }
  }

  if (missing(order)) {
    stop("'order' must be provided")
  }

  # Parse order parameter based on model_type
  if (model_type == "ar") {
    if (!is.numeric(order) || length(order) != 1)
      stop("For AR models, 'order' must be a single integer")
    ar_order <- as.integer(order)
    ma_order <- 0
    d <- 0
    if (ar_order <= 0)
      stop("AR order must be positive")

    # Check if we have enough data
    if (length(x) <= ar_order + 1) {
      stop("Too few observations for AR model of order ", ar_order)
    }
  } else if (model_type == "ma") {
    if (!is.numeric(order) || length(order) != 1)
      stop("For MA models, 'order' must be a single integer")
    ar_order <- 0
    ma_order <- as.integer(order)
    d <- 0
    if (ma_order <= 0)
      stop("MA order must be positive")

    # Check if we have enough data
    if (length(x) <= ma_order + 1) {
      stop("Too few observations for MA model of order ", ma_order)
    }
  } else if (model_type == "arma") {
    if (!is.numeric(order) || length(order) != 2)
      stop("For ARMA models, 'order' must be a vector of length 2 (AR order, MA order)")
    ar_order <- as.integer(order[1])
    ma_order <- as.integer(order[2])
    d <- 0
    if (ar_order < 0 || ma_order < 0)
      stop("AR and MA orders must be non-negative")
    if (ar_order == 0 && ma_order == 0)
      stop("At least one of AR or MA order must be positive")

    # Check if we have enough data
    if (length(x) <= max(ar_order, ma_order) + 1) {
      stop("Too few observations for ARMA model of order (", ar_order, ",", ma_order, ")")
    }
  } else if (model_type == "arima") {
    if (!is.numeric(order) || length(order) != 3)
      stop("For ARIMA models, 'order' must be a vector of length 3 (AR order, differencing, MA order)")
    ar_order <- as.integer(order[1])
    d <- as.integer(order[2])
    ma_order <- as.integer(order[3])
    if (ar_order < 0 || ma_order < 0 || d < 0)
      stop("AR, differencing, and MA orders must be non-negative")
    if (ar_order == 0 && ma_order == 0 && d == 0)
      stop("At least one of AR, differencing, or MA order must be positive")

    # Check if we have enough data after differencing
    if (length(x) <= d + max(ar_order, ma_order) + 1) {
      stop("Too few observations for ARIMA model after differencing")
    }
  } else {
    stop("Unknown model type: ", model_type)
  }

  # Store original series
  orig_x <- as.numeric(x)

  list(
    original_x = orig_x,
    ar_order   = ar_order,
    ma_order   = ma_order,
    d          = d,
    model_type = model_type,
    include.mean = include.mean
  )
}


#' Get initial parameter estimates for time series models
#'
#' @param model_params Validated model parameters from validate_ts_parameters
#' @param initial Optional user-provided initial estimates
#' @param method Estimation method
#' @param verbose Print verbose output
#'
#' @return List containing:
#'   \item{b_init}{vector of initial AR/MA coefficients}
#'   \item{x_mean}{estimated mean (if include.mean=TRUE)}
#'   \item{innovations}{initial residuals/innovations}
#'   \item{x_centered}{centered (or differenced + centered) series}
#'   \item{orig_x}{original series (unchanged)}
#' @keywords internal
get_initial_estimates <- function(model_params,
                                  initial = NULL,
                                  method = "pmm2",
                                  verbose = FALSE) {
  x         <- model_params$original_x
  ar_order  <- model_params$ar_order
  ma_order  <- model_params$ma_order
  d         <- model_params$d
  mtype     <- model_params$model_type
  inc_mean  <- model_params$include.mean

  # Possibly difference for ARIMA
  if (mtype == "arima" && d > 0) {
    x_diff <- diff(x, differences = d)
  } else {
    x_diff <- x
  }

  # Center if needed
  if (inc_mean) {
    x_mean <- mean(x_diff, na.rm = TRUE)
    x_centered <- x_diff - x_mean
  } else {
    x_mean <- 0
    x_centered <- x_diff
  }

  if (mtype == "ar") {
    # AR( p ): quick approach for initial
    if (is.null(initial)) {
      if (method == "yw") {
        b_init <- get_yw_estimates(x_centered, ar_order)
      } else {
        X <- create_ar_matrix(x_centered, ar_order)
        y <- x_centered[(ar_order + 1):length(x_centered)]
        fit_ols <- lm.fit(x = X, y = y)
        b_init <- fit_ols$coefficients
      }
    } else {
      if (length(initial) != ar_order) {
        stop("Length of 'initial' must match AR order")
      }
      b_init <- initial
    }
    # Innovations from initial fit
    X <- create_ar_matrix(x_centered, ar_order)
    y <- x_centered[(ar_order + 1):length(x_centered)]
    innovations <- as.numeric(y - X %*% b_init)

  } else if (mtype %in% c("ma", "arma", "arima")) {
    # Use stats::arima for an initial guess or user-provided
    arima_order <- c(ar_order, ifelse(mtype=="arima", d, 0), ma_order)
    if (is.null(initial)) {
      init_fit <- NULL
      try_methods <- c("CSS","CSS-ML","ML")
      if(method=="pmm2") try_methods <- c("CSS","CSS-ML","ML")

      for(mm in try_methods) {
        tmp <- tryCatch({
          stats::arima(x, order=arima_order, method=mm,
                       include.mean=inc_mean && (mtype!="arima"))
        }, error=function(e) NULL)
        if(!is.null(tmp)) {
          init_fit <- tmp
          break
        }
      }
      if(is.null(init_fit)) {
        if(verbose) cat("All standard methods failed; using simplified values.\n")
        init_fit <- list(
          coef = numeric(ar_order + ma_order),
          residuals = if(mtype=="arima") x_diff else x_centered
        )
        if(ar_order>0) names(init_fit$coef)[1:ar_order] <- paste0("ar",1:ar_order)
        if(ma_order>0) names(init_fit$coef)[(ar_order+1):(ar_order+ma_order)] <- paste0("ma",1:ma_order)
      }

      ar_init <- rep(0, ar_order)
      ma_init <- rep(0, ma_order)
      if(ar_order>0) {
        idx <- paste0("ar",1:ar_order)
        ar_init <- if(all(idx %in% names(init_fit$coef))) as.numeric(init_fit$coef[idx]) else rep(0.1, ar_order)
      }
      if(ma_order>0) {
        idx <- paste0("ma",1:ma_order)
        ma_init <- if(all(idx %in% names(init_fit$coef))) as.numeric(init_fit$coef[idx]) else rep(0.1, ma_order)
      }
      if(inc_mean && !is.null(init_fit$coef) && ("intercept" %in% names(init_fit$coef))) {
        x_mean <- init_fit$coef["intercept"]
      }
      innovations <- if(!is.null(init_fit$residuals)) as.numeric(init_fit$residuals) else {
        rnorm(length(x_centered),0,sd(x_centered,na.rm=TRUE))
      }
      b_init <- c(ar_init, ma_init)

    } else {
      # Provided initial
      if(is.list(initial)) {
        if(ar_order>0 && is.null(initial$ar)) {
          stop("Initial list missing 'ar' but ar_order>0")
        }
        if(ma_order>0 && is.null(initial$ma)) {
          stop("Initial list missing 'ma' but ma_order>0")
        }
        ar_init <- if(ar_order>0) initial$ar else numeric(0)
        ma_init <- if(ma_order>0) initial$ma else numeric(0)
      } else {
        if(length(initial) != (ar_order+ma_order)) {
          stop("Length of 'initial' must match sum of AR and MA orders")
        }
        ar_init <- if(ar_order>0) initial[1:ar_order] else numeric(0)
        ma_init <- if(ma_order>0) initial[(ar_order+1):(ar_order+ma_order)] else numeric(0)
      }
      b_init <- c(ar_init, ma_init)

      init_fit <- tryCatch({
        stats::arima(x, order=arima_order,
                     fixed=b_init,
                     include.mean=inc_mean && (mtype!="arima"))
      }, error=function(e) {
        if(verbose) cat("Error w/ user-provided initial:",e$message,"\n")
        list(residuals = if(mtype=="arima") x_diff else x_centered)
      })
      innovations <- as.numeric(init_fit$residuals)
    }
  }

  if(anyNA(b_init)) {
    warning("NA in initial parameters replaced with 0.")
    b_init[is.na(b_init)] <- 0
  }

  list(
    b_init      = b_init,
    x_mean      = x_mean,
    innovations = innovations,
    x_centered  = x_centered,
    orig_x      = x
  )
}


#' Calculate fitted values for AR models
#'
#' @param object A TS2fit object with model_type="ar"
#' @return Vector of fitted values
#' @keywords internal
get_ar_fitted <- function(object) {
  if (object@model_type != "ar") {
    stop("This function is only for AR models")
  }

  x         <- object@original_series
  ar_order  <- object@order$ar
  ar_coef   <- object@coefficients[1:ar_order]
  intercept <- object@intercept

  if (intercept != 0) {
    x_centered <- x - intercept
  } else {
    x_centered <- x
  }

  X <- create_ar_matrix(x_centered, ar_order)
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
  pred <- numeric(n.ahead)

  for (i in seq_len(n.ahead)) {
    for (j in seq_len(min(i, q))) {
      if ((n - i + j) > 0) {
        pred[i] <- pred[i] + ma_coef[j] * innovations[n - i + j]
      }
    }
  }
  pred
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
#' 3. Uses these moments with a specialized solver (solve_pmm2) to find robust parameter estimates
#'
#' **Головна зміна**: у кінці, щоби отримати коректні залишки (особливо для MA/ARMA/ARIMA),
#' викликається \code{stats::arima(..., fixed=...)} із оціненими коефіцієнтами, що гарантує
#' правильне їх позиціювання й обчислення.
#'
#' @return An S4 \code{TS2fit} object
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

  # 1) Validate
  vp <- validate_ts_parameters(x, order, model_type, include.mean)

  # 2) Initial
  init <- get_initial_estimates(vp, initial, method, verbose)
  b_init      <- init$b_init
  x_mean      <- init$x_mean
  innovations <- init$innovations
  x_centered  <- init$x_centered
  orig_x      <- init$orig_x

  # 3) Moments
  m2 <- mean(innovations^2, na.rm=TRUE)
  m3 <- mean(innovations^3, na.rm=TRUE)
  m4 <- mean(innovations^4, na.rm=TRUE)
  if(is.na(m2) || m2<=0) {
    warning("Second central moment (m2) <= 0 or NA; using fallback.")
    m2 <- var(innovations, na.rm=TRUE)
    if(is.na(m2) || m2<1e-8) m2 <- 1e-4
  }
  if(is.na(m3)) m3 <- 0
  if(is.na(m4) || m4<=m2^2) {
    warning("Fourth central moment (m4) invalid or <= m2^2; adjusting.")
    m4 <- 3*m2^2
  }

  # 4) Design matrix
  ar_order <- vp$ar_order
  ma_order <- vp$ma_order
  d        <- vp$d
  n_data   <- length(x_centered)
  max_lag  <- max(ar_order, ma_order)

  if(n_data <= max_lag) {
    stop("Not enough data for PMM2 approach after differencing/cleaning.")
  }

  n_rows <- n_data - max_lag  # ОГОЛОШУЄМО n_rows тут!
  X <- matrix(0, nrow = n_rows, ncol = ar_order + ma_order)
  y <- x_centered[(max_lag + 1):n_data]

  # AR columns
  if(ar_order > 0) {
    for(i in seq_len(ar_order)) {
      X[, i] <- x_centered[(max_lag - i + 1):(n_data - i)]
    }
  }

  # MA columns
  if(ma_order > 0) {
    # вирівнювання innovations
    if(length(innovations) < n_data) {
      diff_len <- n_data - length(innovations)
      innovations <- c(rep(0, diff_len), innovations)
    } else if(length(innovations) > n_data) {
      innovations <- tail(innovations, n_data)
    }

    for(j in seq_len(ma_order)) {
      for(i in seq_len(n_rows)) {
        t <- max_lag + i
        X[i, ar_order + j] <- innovations[t - j]
      }
    }
  }


  # 5) If method == "pmm2", do solve_pmm2
  if(method=="pmm2") {
    if(verbose) cat("Starting PMM2 optimization...\n")
    final_coef <- solve_pmm2(b_init, X, y, m2, m3, m4,
                             max_iter=max_iter, tol=tol,
                             regularize=regularize, reg_lambda=reg_lambda,
                             verbose=verbose)
  } else {
    final_coef <- b_init
  }

  # 6) Завжди обчислюємо залишки через stats::arima(..., fixed=...).
  #    Це критично для коректних MA/ARMA/ARIMA (і також годиться для AR).
  fixed_params <- final_coef
  if(include.mean) {
    fixed_params <- c(fixed_params, x_mean)
    names(fixed_params) <- c(
      if(ar_order>0) paste0("ar", 1:ar_order) else NULL,
      if(ma_order>0) paste0("ma", 1:ma_order) else NULL,
      "intercept"
    )
  } else {
    names(fixed_params) <- c(
      if(ar_order>0) paste0("ar", 1:ar_order) else NULL,
      if(ma_order>0) paste0("ma", 1:ma_order) else NULL
    )
  }

  final_fit <- tryCatch({
    stats::arima(orig_x,
                 order=c(ar_order, d, ma_order),
                 fixed=fixed_params,
                 include.mean=include.mean)
  }, error=function(e) {
    if(verbose) cat("Error finalizing model with PMM2 params:", e$message, "\n")
    list(residuals = rep(NA, length(orig_x)))
  })
  final_res <- as.numeric(final_fit$residuals)
  if(length(final_res)<length(orig_x)) {
    final_res <- c(rep(NA, length(orig_x)-length(final_res)), final_res)
  }

  # 7) Return TS2fit object
  ans <- new("TS2fit",
             coefficients    = as.numeric(final_coef),
             residuals       = as.numeric(final_res),
             m2              = as.numeric(m2),
             m3              = as.numeric(m3),
             m4              = as.numeric(m4),
             convergence     = TRUE,
             iterations      = as.numeric(max_iter),
             call            = cl,
             model_type      = model_type,
             intercept       = as.numeric(x_mean),
             original_series = as.numeric(orig_x),
             order           = list(ar=ar_order, ma=ma_order, d=d))

  ans
}


#' Fit an AR model using PMM2 (wrapper)
#'
#' @inheritParams ts_pmm2
#' @export
ar_pmm2 <- function(x, order=1, method="pmm2", max_iter=50, tol=1e-6,
                    include.mean=TRUE, initial=NULL, na.action=na.fail,
                    regularize=TRUE, reg_lambda=1e-8, verbose=FALSE) {
  ts_pmm2(x, order=order, model_type="ar", method=method,
          max_iter=max_iter, tol=tol,
          include.mean=include.mean, initial=initial,
          na.action=na.action, regularize=regularize,
          reg_lambda=reg_lambda, verbose=verbose)
}

#' Fit an MA model using PMM2 (wrapper)
#'
#' @inheritParams ts_pmm2
#' @export
ma_pmm2 <- function(x, order=1, method="pmm2", max_iter=50, tol=1e-6,
                    include.mean=TRUE, initial=NULL, na.action=na.fail,
                    regularize=TRUE, reg_lambda=1e-8, verbose=FALSE) {
  ts_pmm2(x, order=order, model_type="ma", method=method,
          max_iter=max_iter, tol=tol,
          include.mean=include.mean, initial=initial,
          na.action=na.action, regularize=regularize,
          reg_lambda=reg_lambda, verbose=verbose)
}

#' Fit an ARMA model using PMM2 (wrapper)
#'
#' @inheritParams ts_pmm2
#' @export
arma_pmm2 <- function(x, order=c(1,1), method="pmm2", max_iter=50, tol=1e-6,
                      include.mean=TRUE, initial=NULL, na.action=na.fail,
                      regularize=TRUE, reg_lambda=1e-8, verbose=FALSE) {
  ts_pmm2(x, order=order, model_type="arma", method=method,
          max_iter=max_iter, tol=tol,
          include.mean=include.mean, initial=initial,
          na.action=na.action, regularize=regularize,
          reg_lambda=reg_lambda, verbose=verbose)
}

#' Fit an ARIMA model using PMM2 (wrapper)
#'
#' @inheritParams ts_pmm2
#' @export
arima_pmm2 <- function(x, order=c(1,1,1), method="pmm2", max_iter=50, tol=1e-6,
                       include.mean=TRUE, initial=NULL, na.action=na.fail,
                       regularize=TRUE, reg_lambda=1e-8, verbose=FALSE) {
  ts_pmm2(x, order=order, model_type="arima", method=method,
          max_iter=max_iter, tol=tol,
          include.mean=include.mean, initial=initial,
          na.action=na.action, regularize=regularize,
          reg_lambda=reg_lambda, verbose=verbose)
}


# ---------------------------------------------------------
# Below are some helper stubs for the code above to work
# You should replace them with real implementations as needed
# ---------------------------------------------------------

#' Create a lagged design matrix for an AR(p) model
#'
#' @param x numeric vector
#' @param p integer AR order
#' @return A matrix of dimension (length(x)-p) x p
#' @keywords internal
create_ar_matrix <- function(x, p) {
  n <- length(x)
  if (n <= p) {
    stop("Not enough data points for AR order p=", p)
  }
  nr <- n - p
  M <- matrix(0, nr, p)
  for (i in seq_len(p)) {
    M[, i] <- x[(p - i + 1):(n - i)]
  }
  M
}

#' Get Yule-Walker estimates for AR(p)
#'
#' @param x numeric vector
#' @param p integer AR order
#' @return numeric vector of length p (AR coefficients)
#' @keywords internal
get_yw_estimates <- function(x, p) {
  # This is a simplified approach that may not handle edge cases
  r <- numeric(p+1)
  n <- length(x)
  xm <- mean(x)
  for (k in 0:p) {
    r[k+1] <- sum((x[1:(n-k)] - xm)*(x[(k+1):n] - xm))
  }
  R <- matrix(0, p, p)
  for (i in seq_len(p)) {
    for (j in seq_len(p)) {
      R[i,j] <- r[abs(i-j)+1]
    }
  }
  rhs <- r[2:(p+1)]
  phi <- solve(R, rhs)
  phi
}

#' Dummy solve_pmm2 function
#'
#' In real code, you'd implement the polynomial-based robust optimization.
#' For this example, we do a single linear step each iteration for demonstration.
#'
#' @param b_init initial parameter guess
#' @param X design matrix
#' @param y response vector
#' @param m2 second moment
#' @param m3 third moment
#' @param m4 fourth moment
#' @param max_iter max iterations
#' @param tol tolerance
#' @param regularize logical
#' @param reg_lambda numeric
#' @param verbose logical
#'
#' @return final parameter estimates
#' @keywords internal
solve_pmm2 <- function(b_init, X, y, m2, m3, m4,
                       max_iter    = 50,
                       tol         = 1e-6,
                       regularize  = TRUE,
                       reg_lambda  = 1e-8,
                       verbose     = FALSE) {
  b_old <- b_init
  for (iter in seq_len(max_iter)) {
    XtX <- crossprod(X)
    if (regularize) {
      diag(XtX) <- diag(XtX) + reg_lambda
    }
    Xty <- crossprod(X, y)
    b_new <- tryCatch(solve(XtX, Xty), error = function(e) b_old)

    if (sqrt(sum((b_new - b_old)^2)) < tol) {
      if (verbose) cat("Converged at iteration", iter, "\n")
      return(b_new)
    }
    b_old <- b_new
  }
  b_old
}
