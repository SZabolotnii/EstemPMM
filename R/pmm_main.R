# pmm_main.R

#' pmm2: Master function for PMM2 (S=2)
#'
#' Fits a linear model using the Polynomial Maximization Method (order 2)
#' which is robust against non-Gaussian errors.
#'
#' @param formula R formula for the model
#' @param data data.frame containing the variables in the formula
#' @param max_iter integer: maximum number of iterations for the algorithm
#' @param tol numeric: tolerance for convergence
#' @param regularize logical: add small value to diagonal for numerical stability
#' @param reg_lambda numeric: regularization parameter (if regularize=TRUE)
#' @param na.action function to handle missing values, default is na.fail
#' @param weights optional vector of weights (not yet implemented)
#' @param verbose logical: whether to print progress information
#'
#' @details
#' The PMM2 algorithm works as follows:
#'
#' 1. Fits ordinary least squares (OLS) regression to get initial estimates
#' 2. Computes central moments (m2, m3, m4) from OLS residuals
#' 3. Iteratively improves parameter estimates using a gradient-based approach
#'
#' PMM2 is particularly useful when error terms are non-Gaussian.
#'
#' @return An S4 \code{PMM2fit} object
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate sample data with t-distributed errors
#' n <- 100
#' x <- rnorm(n)
#' y <- 2 + 3*x + rt(n, df=3)
#' dat <- data.frame(y=y, x=x)
#'
#' # Fit model using PMM2
#' fit <- lm_pmm2(y ~ x, data=dat)
#'
#' # Summary and inference
#' summary(fit, formula=y~x, data=dat)
#' }
lm_pmm2 <- function(formula, data,
                    max_iter=50, tol=1e-6,
                    regularize=TRUE, reg_lambda=1e-8,
                    na.action=na.fail, weights=NULL,
                    verbose=FALSE)
{
  # Capture the call
  call <- match.call()

  # Check input validity
  if(missing(formula) || missing(data)) {
    stop("Both 'formula' and 'data' must be provided")
  }

  if(!is.data.frame(data)) {
    stop("'data' must be a data.frame")
  }

  if(max_iter <= 0) {
    stop("'max_iter' must be positive")
  }

  if(tol <= 0) {
    stop("'tol' must be positive")
  }

  # Handle missing values
  if(!is.null(na.action)) {
    data <- na.action(data)
  }

  # Check for weights
  if(!is.null(weights)) {
    warning("Weights are not yet implemented in PMM2. Ignoring weights.")
  }

  # 1) OLS
  if(verbose) cat("Fitting initial OLS model...\n")

  mf <- model.frame(formula, data)
  X  <- model.matrix(formula, mf)
  y  <- model.response(mf)

  # Check for rank deficiency
  qr_X <- qr(X)
  if(qr_X$rank < ncol(X)) {
    warning("Design matrix is rank deficient, some coefficients may be inestimable")
  }

  fit_ols <- lm.fit(x=X, y=y)
  b_ols   <- fit_ols$coefficients

  # Handle NAs in OLS coefficients
  if(any(is.na(b_ols))) {
    stop("OLS fitting resulted in NA coefficients. Check for multicollinearity.")
  }

  # 2) OLS residuals => m2, m3, m4
  res_ols <- y - (X %*% b_ols)
  m2 <- mean(res_ols^2)
  m3 <- mean(res_ols^3)
  m4 <- mean(res_ols^4)

  if(verbose) {
    cat("Initial moments from OLS residuals:\n")
    cat("  m2 =", m2, "\n")
    cat("  m3 =", m3, "\n")
    cat("  m4 =", m4, "\n")
  }

  # Check for potential issues with moments
  if(m2 <= 0) {
    warning("Second central moment (m2) is non-positive. Results may be unreliable.")
  }

  if(m4 <= m2^2) {
    warning("Fourth central moment (m4) is less than m2^2. This violates a basic inequality for probability distributions.")
  }

  # 3) Call the main fitting function
  if(verbose) cat("Starting PMM2 iterations...\n")

  out <- .pmm2_fit(b_ols, X, y, m2, m3, m4, max_iter, tol, regularize, reg_lambda, verbose)

  # Extract results
  b_est   <- out$b
  conv    <- out$convergence
  iter    <- out$iterations

  if(verbose) {
    cat("PMM2 algorithm finished.\n")
    cat("  Converged:", conv, "\n")
    cat("  Iterations:", iter, "\n")
  }

  # Final residuals
  final_res <- as.numeric(y - X %*% b_est)

  # Return the S4 object with all results
  ans <- new("PMM2fit",
             coefficients = b_est,
             residuals = final_res,
             m2 = m2,
             m3 = m3,
             m4 = m4,
             convergence = conv,
             iterations = iter,
             call = call)

  return(ans)
}

#' Extract coefficients from a PMM2fit object
#'
#' @param object A PMM2fit object
#' @param ... Additional arguments (not used)
#'
#' @return A vector of coefficients
#' @export
setMethod("coef", "PMM2fit",
          function(object, ...) {
            object@coefficients
          })

#' Extract residuals from a PMM2fit object
#'
#' @param object A PMM2fit object
#' @param ... Additional arguments (not used)
#'
#' @return A vector of residuals
#' @export
setMethod("residuals", "PMM2fit",
          function(object, ...) {
            object@residuals
          })

#' Extract fitted values from a PMM2fit object
#'
#' @param object A PMM2fit object
#' @param ... Additional arguments (not used)
#'
#' @return A vector of fitted values
#' @export
setMethod("fitted", "PMM2fit",
          function(object, ...) {
            fitted_values(object)
          })

#' Calculate AIC for a PMM2fit object
#'
#' @param object A PMM2fit object
#' @param ... Additional arguments (not used)
#' @param k The penalty per parameter to be used; default is 2
#'
#' @return The AIC value
#' @export
setMethod("AIC", "PMM2fit",
          function(object, ..., k = 2) {
            res <- object@residuals
            n <- length(res)
            p <- length(object@coefficients)

            # Log-likelihood approximation
            ll <- -n/2 * log(sum(res^2)/n) - n/2 * (1 + log(2*pi))

            # AIC
            -2 * ll + k * p
          })
