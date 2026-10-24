# pmm3_main.R - Main module for PMM3 linear models

#' PMM3: Fit linear model using Polynomial Maximization Method (S=3)
#'
#' Fits a linear model using PMM3, which is designed for symmetric
#' platykurtic error distributions. Uses a cubic stochastic polynomial
#' with Newton-Raphson solver.
#'
#' @param formula R formula for the model
#' @param data data.frame containing variables in the formula
#' @param max_iter integer: maximum number of NR iterations (default 100)
#' @param tol numeric: convergence tolerance (default 1e-6)
#' @param adaptive logical: re-estimate kappa each iteration (default FALSE)
#' @param step_max numeric: maximum NR step, measured as the root-mean-square
#'   change of the fitted values in units of the OLS residual standard
#'   deviation (default 5.0)
#' @param na.action function for handling missing values (default na.fail)
#' @param verbose logical: print progress information (default FALSE)
#' @param symmetric logical: \code{TRUE} (default) fits the symmetric PMM3
#'   estimator described below; \code{FALSE} fits the general degree-3
#'   estimator, which also uses the odd moments and is valid for skewed
#'   errors.
#'
#' @details
#' The PMM3 algorithm works as follows:
#'
#' 1. Fits OLS regression to obtain initial estimates
#' 2. Computes central moments (m2, m4, m6) from OLS residuals
#' 3. Checks symmetry: warns if |gamma3| > 0.3 (PMM2 may be more appropriate)
#' 4. Computes moment ratio kappa = (m6 - 3*m4*m2) / (m4 - 3*m2^2)
#' 5. Iterates Newton-Raphson with score Z = X'(eps*(kappa - eps^2))
#'
#' PMM3 achieves lower variance than OLS when errors are symmetric and
#' platykurtic (gamma4 < 0), with theoretical efficiency
#' g3 = 1 - gamma4^2 / (6 + 9*gamma4 + gamma6).
#'
#' The symmetric estimator assumes \eqn{\gamma_3 = \gamma_5 = 0}. On skewed
#' errors it remains consistent for the slopes but loses most of its gain.
#' With chi-squared(3), gamma, exponential and lognormal errors its OLS-to-PMM3
#' MSE ratio is 1.1 to 1.5, against 1.7 to 2.4 for PMM2 and 2.1 to 3.8 for the
#' general estimator (\code{dev_scripts/check_pmm3_general.R}). With
#' \code{symmetric = FALSE} the polynomial is
#' \eqn{\psi(e) = \sum_{k=1}^{3} h_k (e^k - \alpha_k)}, where \eqn{\alpha_k}
#' are the central moments of the OLS residuals and \eqn{h = F^{-1} b} solves
#' the full Kunchenko normal system with \eqn{F_{pq} = \alpha_{p+q} - \alpha_p
#' \alpha_q} and \eqn{b_k = k \alpha_{k-1}}. Its variance-reduction factor
#' is \eqn{g_3 = 1 / (\alpha_2\, b^\top F^{-1} b)}. This factor never
#' exceeds the PMM2 factor \eqn{g_2}, and it equals the symmetric
#' \eqn{g_3} when the odd moments vanish. The Newton iteration uses
#' backtracking on the score norm and falls back to OLS with a warning when
#' it fails. \code{adaptive = TRUE} re-estimates the moments and \eqn{h}
#' at every iteration.
#'
#' @return S4 object of class \code{PMM3fit}
#' @export
#'
#' @examples
#' set.seed(123)
#' n <- 100
#' x <- rnorm(n)
#' y <- 2 + 3 * x + runif(n, -sqrt(3), sqrt(3))
#' dat <- data.frame(y = y, x = x)
#'
#' fit <- lm_pmm3(y ~ x, data = dat)
#' summary(fit)
#'
#' # Skewed errors: the general degree-3 estimator
#' dat$y2 <- 2 + 3 * x + rgamma(n, 2, 1) - 2
#' fit_g <- lm_pmm3(y2 ~ x, data = dat, symmetric = FALSE)
#' coef(fit_g)
lm_pmm3 <- function(formula, data,
                     max_iter = 100, tol = 1e-6,
                     adaptive = FALSE, step_max = 5.0,
                     na.action = na.fail,
                     verbose = FALSE,
                     symmetric = TRUE) {
  # Capture call
  call <- match.call()

  # Validate input
  if (missing(formula) || missing(data)) {
    stop("Both 'formula' and 'data' must be provided")
  }

  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame")
  }

  if (max_iter <= 0) {
    stop("'max_iter' must be positive")
  }

  if (tol <= 0) {
    stop("'tol' must be positive")
  }

  # Handle missing values
  if (!is.null(na.action)) {
    data <- na.action(data)
  }

  # 1) OLS
  if (verbose) cat("Fitting initial OLS model...\n")

  mf <- model.frame(formula, data)
  X  <- model.matrix(formula, mf)
  y  <- model.response(mf)

  # Check for rank deficiency
  qr_X <- qr(X)
  if (qr_X$rank < ncol(X)) {
    warning("Design matrix has rank deficiency, some coefficients may not be estimable")
  }

  fit_ols <- lm.fit(x = X, y = y)
  b_ols   <- as.numeric(fit_ols$coefficients)

  if (any(is.na(b_ols))) {
    stop("OLS fit resulted in NA coefficients. Check for multicollinearity.")
  }

  # 2) OLS residuals -> moments
  res_ols <- as.numeric(y - X %*% b_ols)
  moments <- compute_moments_pmm3(res_ols)

  if (verbose) {
    cat("Initial moments from OLS residuals:\n")
    cat("  m2 =", moments$m2, "\n")
    cat("  m4 =", moments$m4, "\n")
    cat("  m6 =", moments$m6, "\n")
    cat("  gamma4 =", moments$gamma4, "\n")
    cat("  kappa =", moments$kappa, "\n")
  }

  if (moments$m2 <= 0) {
    warning("Second central moment (m2) is not positive. Results may be unreliable.")
  }

  if (!isTRUE(symmetric)) {
    return(.lm_pmm3_general(X, y, b_ols, res_ols, adaptive = adaptive,
                            tol = tol, max_iter = max_iter, call = call,
                            mf = mf, data = data))
  }

  # Symmetry warning
  if (abs(moments$gamma3) > 0.3) {
    warning(sprintf(
      "Residuals appear asymmetric (|gamma3| = %.3f > 0.3). Consider lm_pmm2() or lm_pmm3(symmetric = FALSE) instead.",
      abs(moments$gamma3)
    ))
  }

  # 3) Newton-Raphson
  if (verbose) cat("Starting PMM3 iterations...\n")

  nr <- .pmm3_nr_solver(
    B_ols    = b_ols,
    X        = X,
    Y        = y,
    kappa    = moments$kappa,
    adaptive = adaptive,
    tol      = tol,
    max_iter = max_iter,
    step_max = step_max
  )

  b_est <- nr$B
  names(b_est) <- colnames(X)
  final_res <- as.numeric(y - X %*% b_est)

  if (verbose) {
    cat("PMM3 algorithm completed.\n")
    cat("  Converged:", nr$converged, "\n")
    cat("  Iterations:", nr$iter, "\n")
  }

  # Compute variance factor from initial moments
  vf <- pmm3_variance_factor(moments$m2, moments$m4, moments$m6)

  # Return S4 object
  ans <- new("PMM3fit",
             coefficients  = b_est,
             residuals     = final_res,
             m2            = moments$m2,
             m4            = moments$m4,
             m6            = moments$m6,
             gamma4        = vf$gamma4,
             gamma6        = vf$gamma6,
             g_coefficient = if (is.na(vf$g3)) 1.0 else vf$g3,
             kappa         = if (is.na(moments$kappa)) NA_real_ else moments$kappa,
             convergence   = nr$converged,
             iterations    = as.numeric(nr$iter),
             call          = call)
  attr(ans, "model_matrix") <- X
  attr(ans, "model_frame")  <- mf
  attr(ans, "response")     <- as.numeric(y)
  attr(ans, "data")         <- data
  # Near-Gaussian residuals leave kappa undefined and the OLS values are
  # returned without a PMM3 step, so g3 does not apply to them either
  attr(ans, "fallback")     <- isTRUE(nr$fallback) || isTRUE(nr$near_gaussian)

  return(ans)
}

# General (asymmetric) degree-3 branch of lm_pmm3(symmetric = FALSE).
.lm_pmm3_general <- function(X, y, b_ols, res_ols, adaptive, tol, max_iter,
                             call, mf, data) {
  r <- res_ols - mean(res_ols)
  alpha <- vapply(1:6, function(k) mean(r^k), numeric(1))
  sol <- .pmm3_general_solver(b_ols, X, y, alpha, adaptive = adaptive,
                              tol = tol, max_iter = max_iter)
  b_est <- sol$B
  names(b_est) <- colnames(X)

  m2 <- alpha[2]
  gamma3 <- alpha[3] / m2^1.5
  gamma4 <- alpha[4] / m2^2 - 3
  gamma6 <- alpha[6] / m2^3 - 15 * alpha[4] / m2^2 - 10 * gamma3^2 + 30
  g3 <- .pmm_poly_coef(alpha, 3)$g

  ans <- new("PMM3fit",
             coefficients  = b_est,
             residuals     = as.numeric(y - X %*% b_est),
             m2            = m2,
             m4            = alpha[4],
             m6            = alpha[6],
             gamma4        = gamma4,
             gamma6        = gamma6,
             g_coefficient = if (is.finite(g3)) g3 else 1.0,
             kappa         = NA_real_,
             convergence   = sol$converged,
             iterations    = as.numeric(sol$iter),
             call          = call,
             symmetric     = FALSE,
             poly_coef     = sol$h)
  attr(ans, "model_matrix") <- X
  attr(ans, "model_frame")  <- mf
  attr(ans, "response")     <- as.numeric(y)
  attr(ans, "data")         <- data
  attr(ans, "moments")      <- alpha
  attr(ans, "jacobian_min_eigen") <- sol$jacobian_min_eigen
  attr(ans, "fallback")     <- isTRUE(sol$fallback)
  ans
}

# =============================================================================
# S4 Methods for PMM3fit
# =============================================================================

#' Extract coefficients from PMM3fit object
#'
#' @param object PMM3fit object
#' @param ... Additional arguments (not used)
#'
#' @return Vector of coefficients
#' @export
setMethod("coef", "PMM3fit",
          function(object, ...) {
            object@coefficients
          })

#' Extract residuals from PMM3fit object
#'
#' @param object PMM3fit object
#' @param ... Additional arguments (not used)
#'
#' @return Vector of residuals
#' @export
setMethod("residuals", "PMM3fit",
          function(object, ...) {
            object@residuals
          })

#' Extract fitted values from PMM3fit object
#'
#' @param object PMM3fit object
#' @param ... Additional arguments (not used)
#'
#' @return Vector of fitted values
#' @export
setMethod("fitted", "PMM3fit",
          function(object, ...) {
            stored_X <- attr(object, "model_matrix")
            if (!is.null(stored_X)) {
              return(as.vector(stored_X %*% object@coefficients))
            }
            stored_response <- attr(object, "response")
            if (!is.null(stored_response)) {
              return(as.vector(stored_response - object@residuals))
            }
            stop("Cannot compute fitted values: no stored model matrix or response")
          })

#' Predict method for PMM3fit objects
#'
#' @param object PMM3fit object
#' @param newdata Data frame with predictor variables
#' @param ... Additional arguments (not used)
#'
#' @return Numeric vector of predicted values
#' @export
setMethod("predict", "PMM3fit",
          function(object, newdata = NULL, ...) {
            if (is.null(newdata)) {
              return(fitted(object))
            }
            formula <- eval(object@call$formula)
            rhs <- formula[[3]]
            design_formula <- as.formula(paste("~", deparse(rhs)))
            X <- model.matrix(design_formula, newdata)

            if (is.null(names(object@coefficients)) ||
                all(names(object@coefficients) == "")) {
              if (length(colnames(X)) == length(object@coefficients)) {
                names(object@coefficients) <- colnames(X)
              } else {
                stop("Number of coefficients does not match design matrix columns.")
              }
            }

            if (!identical(names(object@coefficients), colnames(X))) {
              if (all(colnames(X) %in% names(object@coefficients))) {
                object@coefficients <- object@coefficients[colnames(X)]
              } else {
                stop("Design matrix columns do not match coefficient names.")
              }
            }

            as.vector(X %*% object@coefficients)
          })

#' Summary method for PMM3fit objects
#'
#' @param object PMM3fit object
#' @param ... Additional arguments (not used)
#'
#' @return Prints summary to console; returns object (invisibly)
#' @export
setMethod("summary", "PMM3fit",
          function(object, ...) {
            general <- methods::.hasSlot(object, "symmetric") && !isTRUE(object@symmetric)
            if (general) {
              cat("PMM3 estimation results (S = 3, general polynomial)\n")
            } else {
              cat("PMM3 estimation results (S = 3, symmetric errors)\n")
            }
            if (!is.null(object@call)) {
              cat("Call:\n")
              print(object@call)
              cat("\n")
            }

            cat("Coefficients:\n")
            print(object@coefficients)

            cat("\nCentral moments of initial residuals:\n")
            cat("  m2 =", object@m2, "\n")
            cat("  m4 =", object@m4, "\n")
            cat("  m6 =", object@m6, "\n\n")

            cat("Theoretical characteristics of PMM3 (S = 3):\n")
            cat("  gamma4 =", object@gamma4, "\n")
            cat("  gamma6 =", object@gamma6, "\n")
            cat("  g3     =", object@g_coefficient,
                " (expected ratio Var[PMM3]/Var[OLS])\n")
            if (!is.na(object@kappa)) {
              cat("  kappa  =", object@kappa, "\n")
            }
            if (general) {
              cat("  h      =", format(object@poly_coef, digits = 4), "\n")
            }
            cat("\n")

            cat("Algorithm information:\n")
            cat("  Convergence status:", object@convergence, "\n")
            cat("  Iterations:", object@iterations, "\n\n")

            invisible(object)
          })

# AIC and BIC for PMM3fit are obtained via the default dispatch through
# logLik(), which returns a "logLik" object with df and nobs attributes.
# See setMethod("logLik", "PMM3fit", ...) below.

#' Extract log-likelihood from PMM3fit object
#'
#' Returns a Gaussian approximate log-likelihood. Defined as an S3 method
#' so that \code{stats::AIC} and \code{stats::BIC} (which dispatch
#' \code{logLik} via \code{UseMethod}) work without explicit S4 methods.
#'
#' @param object PMM3fit object
#' @param ... Additional arguments (not used)
#'
#' @return Object of class \code{logLik}
#' @method logLik PMM3fit
#' @export
logLik.PMM3fit <- function(object, ...) {
  res <- object@residuals
  n <- length(res)
  p <- length(object@coefficients)
  ll <- -n/2 * log(sum(res^2)/n) - n/2 * (1 + log(2*pi))
  attr(ll, "df") <- p
  attr(ll, "nobs") <- n
  class(ll) <- "logLik"
  ll
}

#' Number of observations in PMM3fit object
#'
#' @param object PMM3fit object
#' @param ... Additional arguments (not used)
#'
#' @return Integer number of observations
#' @export
setMethod("nobs", "PMM3fit",
          function(object, ...) {
            length(object@residuals)
          })

#' Variance-covariance matrix for PMM3fit objects
#'
#' Returns the asymptotic covariance matrix
#' \eqn{V_{\mathrm{PMM3}} = g_3\,\sigma^2 (X^\top X)^{-1}} where
#' \eqn{g_3 = 1 - \gamma_4^{2} / (6 + 9\gamma_4 + \gamma_6)}
#' is the PMM3 variance reduction factor (Kunchenko's polynomial
#' maximization method, s = 3 specialised to symmetric platykurtic
#' errors; see \code{\link{pmm3_variance_matrices}}). Under Gaussian
#' errors \eqn{g_3 = 1} and the result equals the OLS covariance.
#' For \code{lm_pmm3(symmetric = FALSE)}, \eqn{g_3} is the general factor
#' \eqn{1 / (\alpha_2\, b^\top F^{-1} b)}, and the matrix gains the term
#' \eqn{(1 - g_3)\sigma^2 n^{-1} a a^\top}, \eqn{a = (X^\top X)^{-1} X^\top 1},
#' for the moments estimated from the OLS residuals (see
#' \code{\link{vcov,PMM2fit-method}}); with an intercept only the intercept
#' variance changes. The symmetric estimator needs no such term. When the
#' iteration failed and
#' the OLS estimate was returned, the OLS covariance is returned.
#'
#' @param object PMM3fit object
#' @param ... Additional arguments (not used)
#'
#' @return Numeric matrix of the same dimension as the number of coefficients
#' @export
setMethod("vcov", "PMM3fit",
          function(object, ...) {
            X <- attr(object, "model_matrix")
            if (is.null(X))
              stop("model_matrix not found; refit the model with lm_pmm3()")
            if (isTRUE(attr(object, "fallback"))) {
              V <- object@m2 * solve(crossprod(X))
            } else if (methods::.hasSlot(object, "symmetric") && !isTRUE(object@symmetric)) {
              V <- object@g_coefficient * object@m2 * solve(crossprod(X)) +
                .pmm_location_term(X, object@m2, object@g_coefficient)
            } else {
              V <- pmm3_variance_matrices(X, object@m2, object@m4, object@m6)$pmm3
            }
            nms <- names(object@coefficients)
            if (!is.null(nms) && length(nms) == nrow(V)) {
              rownames(V) <- nms
              colnames(V) <- nms
            }
            V
          })

#' Confidence intervals for PMM3fit coefficients
#'
#' Computes normal-approximation confidence intervals using the
#' asymptotic covariance matrix from \code{\link{vcov,PMM3fit-method}}.
#'
#' @param object PMM3fit object
#' @param parm character or integer vector of parameter names/indices to include
#' @param level confidence level (default 0.95)
#' @param ... Additional arguments (not used)
#'
#' @return Matrix with lower and upper confidence limits
#' @export
setMethod("confint", "PMM3fit",
          function(object, parm, level = 0.95, ...) {
            cf <- object@coefficients
            se <- sqrt(diag(vcov(object)))
            a  <- (1 - level) / 2
            fac <- stats::qnorm(c(a, 1 - a))
            ci  <- cf + outer(se, fac)
            colnames(ci) <- paste0(format(100 * c(a, 1 - a), trim = TRUE), " %")
            if (!missing(parm)) ci <- ci[parm, , drop = FALSE]
            ci
          })

#' Plot diagnostic plots for PMM3fit object
#'
#' @param x PMM3fit object
#' @param y Not used (compatibility with generic)
#' @param which Set of plots to display (values 1-4)
#' @param ... Additional arguments passed to plotting functions
#'
#' @return Invisibly returns the input object
#' @export
setMethod("plot", signature(x = "PMM3fit", y = "missing"),
          function(x, y, which = 1:4, ...) {
            res <- as.numeric(x@residuals)
            fitted_vals <- tryCatch({
              fitted(x)
            }, error = function(e) {
              stored_X <- attr(x, "model_matrix")
              if (!is.null(stored_X)) {
                as.vector(stored_X %*% x@coefficients)
              } else {
                seq_along(res)
              }
            })

            which <- intersect(unique(which), 1:4)
            if (length(which) == 0) which <- 1:4

            old_par <- graphics::par(no.readonly = TRUE)
            on.exit(graphics::par(old_par))
            graphics::par(mfrow = c(2, 2))

            for (idx in which) {
              switch(idx,
                     {
                       graphics::plot(fitted_vals, res,
                                      main = "Residuals vs Fitted",
                                      xlab = "Fitted values",
                                      ylab = "Residuals", ...)
                       graphics::abline(h = 0, col = "red", lty = 2)
                     },
                     {
                       stats::qqnorm(res, main = "Normal Q-Q", ...)
                       stats::qqline(res, col = "red", lty = 2)
                     },
                     {
                       graphics::plot(seq_along(res), res, type = "l",
                                      main = "Residuals over Index",
                                      xlab = "Observation",
                                      ylab = "Residual", ...)
                       graphics::abline(h = 0, col = "red", lty = 2)
                     },
                     {
                       graphics::hist(res,
                                      main = "Residual Histogram",
                                      xlab = "Residuals",
                                      breaks = "FD", ...)
                     })
            }

            invisible(x)
          })
