# wpmm_main.R - Weak-moment (windowed) PMM2 and PMM3 for linear models
#
# Classical PMM2/PMM3 coefficients are built from raw residual moments, which a
# single heavy-tail outlier destabilises (raw m6 need not even exist). The weak
# variants replace raw averaging E[g(r)] by a value-domain windowed functional
#   E_w[g(r)] = sum_i w(r_i) g(r_i) / sum_i w(r_i),
#   w(r) = exp(-0.5 ((r - median(r)) / (sigma_mult * MAD(r)))^2),
# and solve the factored estimating equation sum_i W_i r_i x_i = 0 by IRLS:
#   WPMM2: W_i = w_i (1 + a2 rc_i),  a2 = -c3^w / (2 (c2^w)^2 + c4^w)
#   WPMM3: W_i = w_i (kappa^w - rc_i^2) / |kappa^w|,
#          kappa^w = (m6^w - 3 m4^w m2^w) / (m4^w - 3 (m2^w)^2)
# with rc_i = r_i - E_w[r] and negative weights clipped to zero. The window and
# the weak moments are re-evaluated at every step, so the estimator is a
# redescending M-estimator. Ported from gmdhpmm (EstemPMM-GMDH, R/weak_pmm.R),
# where it is the inner estimator of the GMDH tournament; the IRLS steps and
# their order are kept so that both packages give the same estimates.

#' @include pmm_base_classes.R
NULL

#' S4 class for weak-moment PMM regression fits
#'
#' Concrete class returned by \code{\link{lm_wpmm2}} and
#' \code{\link{lm_wpmm3}}. Inherits the common slots from
#' \code{\linkS4class{PMMfit}}.
#'
#' @slot degree        numeric: 2 for WPMM2, 3 for WPMM3
#' @slot sigma_mult    numeric: window width as a multiple of the robust
#'   residual scale
#' @slot weights       numeric: final IRLS weights \eqn{W_i}
#' @slot window_center numeric: median of the final residuals
#' @slot window_scale  numeric: window standard deviation
#'   (\code{sigma_mult} times the MAD of the final residuals)
#' @slot weak_moments  numeric: windowed central moments of the final
#'   residuals (\code{c2, c3, c4} cumulants for WPMM2;
#'   \code{m2, m4, m6} moments for WPMM3)
#' @slot correction    numeric: the polynomial coefficient, \eqn{a_2} for
#'   WPMM2 or \eqn{\kappa^w} for WPMM3; \code{NA} when the windowed law is
#'   mesokurtic and WPMM3 reduces to the Welsch estimator
#'
#' @exportClass WPMMfit
setClass("WPMMfit",
         contains = "PMMfit",
         slots = c(degree        = "numeric",
                   sigma_mult    = "numeric",
                   weights       = "numeric",
                   window_center = "numeric",
                   window_scale  = "numeric",
                   weak_moments  = "numeric",
                   correction    = "numeric"))

# IRLS weights W_i at residuals r, together with d(W_i r_i)/dr_i for the
# sandwich bread. The derivative holds the window centre, scale and weak
# moments fixed but keeps the window slope w'(r), which does not vanish by
# symmetry. Returns NULL when the window carries no mass.
.wpmm_weights <- function(r, degree, sigma_mult) {
  mu <- stats::median(r)
  s <- stats::mad(r, constant = 1.4826)
  if (!is.finite(s) || s <= 1e-12) s <- stats::sd(r)
  if (!is.finite(s) || s <= 1e-12) s <- 1
  sw <- sigma_mult * s
  w <- exp(-0.5 * ((r - mu) / sw)^2)
  sw_sum <- sum(w)
  if (!is.finite(sw_sum) || sw_sum <= 0) return(NULL)
  dw <- -w * (r - mu) / sw^2
  rc <- r - sum(w * r) / sw_sum

  if (degree == 2) {
    c2 <- sum(w * rc^2) / sw_sum
    c3 <- sum(w * rc^3) / sw_sum
    c4 <- sum(w * rc^4) / sw_sum - 3 * c2^2
    denom <- 2 * c2^2 + c4
    a2 <- if (!is.finite(c2) || c2 <= 1e-12 || !is.finite(denom) || abs(denom) < 1e-12)
            0 else -c3 / denom
    a2 <- max(min(a2, 0.5 / sqrt(c2 + 1e-12)), -0.5 / sqrt(c2 + 1e-12))
    W <- w * (1 + a2 * rc)
    dpsi <- dw * (1 + a2 * rc) * r + w * a2 * r + W
    moments <- c(c2 = c2, c3 = c3, c4 = c4)
    correction <- a2
  } else {
    m2 <- sum(w * rc^2) / sw_sum
    m4 <- sum(w * rc^4) / sw_sum
    m6 <- sum(w * rc^6) / sw_sum
    dk <- m4 - 3 * m2^2
    moments <- c(m2 = m2, m4 = m4, m6 = m6)
    if (!is.finite(dk) || abs(dk) < 1e-3 * max(m2^2, 1e-12)) {
      W <- w
      dpsi <- dw * r + w
      correction <- NA_real_
    } else {
      kappa <- (m6 - 3 * m4 * m2) / dk
      ak <- max(abs(kappa), 1e-8)
      W <- w * (kappa - rc^2) / ak
      dpsi <- (dw * (kappa - rc^2) - 2 * w * rc) * r / ak + W
      correction <- kappa
    }
  }
  clipped <- !is.finite(W) | W < 0
  W[clipped] <- 0
  dpsi[clipped] <- 0
  if (degree == 3 && sum(W) <= 0) {
    W <- w
    dpsi <- dw * r + w
  }
  list(W = W, dpsi = dpsi, center = mu, scale = sw,
       moments = moments, correction = correction)
}

# Shared fitter behind lm_wpmm2() and lm_wpmm3().
.lm_wpmm <- function(formula, data, degree, sigma_mult, max_iter, tol,
                     na.action, call) {
  if (missing(formula) || missing(data)) {
    stop("Both 'formula' and 'data' must be provided")
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame")
  }
  if (!is.numeric(sigma_mult) || length(sigma_mult) != 1L || !(sigma_mult > 0)) {
    stop("'sigma_mult' must be a single positive number")
  }
  if (max_iter <= 0) {
    stop("'max_iter' must be positive")
  }
  if (tol <= 0) {
    stop("'tol' must be positive")
  }
  if (!is.null(na.action)) {
    data <- na.action(data)
  }

  mf <- stats::model.frame(formula, data)
  X  <- stats::model.matrix(formula, mf)
  y  <- as.numeric(stats::model.response(mf))

  theta <- stats::lm.fit(x = X, y = y)$coefficients
  if (any(is.na(theta))) {
    stop("OLS fit resulted in NA coefficients. Check for multicollinearity.")
  }
  theta <- as.numeric(theta)

  converged <- FALSE
  failure <- NULL
  iter <- 0L
  for (iter in seq_len(max_iter)) {
    r <- as.numeric(y - X %*% theta)
    ww <- .wpmm_weights(r, degree, sigma_mult)
    if (is.null(ww)) {
      failure <- "the window carries no mass"
      break
    }
    theta_new <- stats::lm.wfit(x = X, y = y, w = ww$W)$coefficients
    if (any(!is.finite(theta_new))) {
      failure <- "the weighted fit is not estimable"
      break
    }
    theta_new <- as.numeric(theta_new)
    step <- sqrt(sum((theta_new - theta)^2)) / max(1, sqrt(sum(theta^2)))
    theta <- theta_new
    if (step < tol) {
      converged <- TRUE
      break
    }
  }
  label <- if (degree == 2) "WPMM2" else "WPMM3"
  if (!is.null(failure)) {
    warning(label, " iteration stopped because ", failure,
            "; returning the last finite estimate.")
  } else if (!converged) {
    warning(label, " algorithm did not converge after ", max_iter,
            " iterations.")
  }

  names(theta) <- colnames(X)
  res <- as.numeric(y - X %*% theta)
  ww <- .wpmm_weights(res, degree, sigma_mult)
  if (is.null(ww)) {
    ww <- list(W = rep(NA_real_, length(res)), center = NA_real_,
               scale = NA_real_, moments = NA_real_, correction = NA_real_)
  }

  ans <- new("WPMMfit",
             coefficients  = theta,
             residuals     = res,
             convergence   = converged,
             iterations    = as.numeric(iter),
             call          = call,
             degree        = degree,
             sigma_mult    = sigma_mult,
             weights       = ww$W,
             window_center = ww$center,
             window_scale  = ww$scale,
             weak_moments  = ww$moments,
             correction    = ww$correction)
  attr(ans, "model_matrix") <- X
  attr(ans, "model_frame")  <- mf
  attr(ans, "response")     <- y
  attr(ans, "data")         <- data
  ans
}

#' Weak-moment PMM2 and PMM3 for linear models
#'
#' Fit a linear model by the polynomial maximization method with windowed
#' (weak) moments in place of raw residual moments. \code{lm_wpmm2} is the
#' degree-2 estimator for skewed errors, \code{lm_wpmm3} the degree-3
#' estimator for symmetric heavy-tailed errors. Both are meant for residuals
#' whose raw cumulants are unstable: heavy tails, outliers, or moments that do
#' not exist.
#'
#' @param formula R formula for the model.
#' @param data data.frame containing the variables in \code{formula}.
#' @param sigma_mult window width as a multiple of the robust residual scale
#'   (default 2.5). A wider window moves the estimator towards classical
#'   PMM2/PMM3 and OLS; a narrower one makes it more robust.
#' @param max_iter maximum number of IRLS iterations (default 50).
#' @param tol convergence tolerance on the relative change of the
#'   coefficient vector (default 1e-7).
#' @param na.action function for handling missing values (default
#'   \code{na.fail}).
#'
#' @details
#' Starting from OLS, each iteration
#' \enumerate{
#'   \item centres a Gaussian window \eqn{w_i = \exp(-\frac12 ((r_i - \mathrm{med}(r)) / (\sigma\, s))^2)}
#'     on the residual median, with \eqn{s} the MAD of the residuals and
#'     \eqn{\sigma} = \code{sigma_mult};
#'   \item computes windowed central moments of the residuals around the
#'     windowed mean, \eqn{r^c_i = r_i - \sum w_i r_i / \sum w_i};
#'   \item solves the weighted least-squares problem with weights
#'     \eqn{W_i = w_i (1 + a_2 r^c_i)}, \eqn{a_2 = -c_3^w / (2 (c_2^w)^2 + c_4^w)}
#'     for WPMM2, or \eqn{W_i = w_i (\kappa^w - (r^c_i)^2) / |\kappa^w|},
#'     \eqn{\kappa^w = (m_6^w - 3 m_4^w m_2^w) / (m_4^w - 3 (m_2^w)^2)} for
#'     WPMM3. Negative weights are set to zero, so both estimators are
#'     redescending.
#' }
#' The estimates agree with the inner estimators \code{fit_kg2_weak_pmm2}
#' and \code{fit_kg2_weak_pmm3} of the gmdhpmm package (EstemPMM-GMDH), which
#' add a ridge penalty of 1e-8.
#'
#' Two limitations follow from the construction.
#' \itemize{
#'   \item Under skewed errors the intercept estimates the location at which
#'     the windowed score vanishes. That location depends on
#'     \code{sigma_mult}, does not shrink as \eqn{n} grows, and is not the
#'     error mean. The slopes are unaffected. \code{lm_wpmm2} does not remove
#'     this offset.
#'   \item \eqn{\kappa^w} is the coefficient of the windowed law. It ignores
#'     the derivative of the window, so it is not the variance-minimising
#'     coefficient, and the gap grows as the window narrows.
#' }
#' \code{lm_wpmm3} targets symmetric leptokurtic errors. For symmetric
#' platykurtic errors use \code{\link{lm_pmm3}}. When the windowed law is
#' mesokurtic (\eqn{m_4^w \approx 3 (m_2^w)^2}), \code{lm_wpmm3} reduces to
#' the Welsch M-estimator.
#'
#' Standard errors (\code{vcov}) come from the M-estimator sandwich
#' \eqn{A^{-1} B A^{-1}} with \eqn{A = \sum \psi'(r_i) x_i x_i^\top},
#' \eqn{B = \sum \psi(r_i)^2 x_i x_i^\top} and \eqn{\psi(r) = W(r)\,r}.
#' The window centre, scale and weak moments are treated as fixed. For
#' errors independent of the regressors this leaves the slope standard
#' errors asymptotically correct. It can understate the intercept standard
#' error.
#'
#' @return S4 object of class \code{\linkS4class{WPMMfit}}.
#'
#' @seealso \code{\link{lm_pmm2}}, \code{\link{lm_pmm3}},
#'   \code{\link{pmm_lm}}
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' x <- rnorm(n)
#' dat <- data.frame(x = x, y = 1 + 2 * x + rt(n, df = 3))
#' fit3 <- lm_wpmm3(y ~ x, data = dat)
#' summary(fit3)
#'
#' dat$y2 <- 1 + 2 * x + (rexp(n) - 1) + 3 * rbinom(n, 1, 0.05)
#' fit2 <- lm_wpmm2(y2 ~ x, data = dat)
#' coef(fit2)
#' @export
lm_wpmm2 <- function(formula, data, sigma_mult = 2.5,
                     max_iter = 50, tol = 1e-7, na.action = na.fail) {
  .lm_wpmm(formula, data, degree = 2, sigma_mult = sigma_mult,
           max_iter = max_iter, tol = tol, na.action = na.action,
           call = match.call())
}

#' @rdname lm_wpmm2
#' @export
lm_wpmm3 <- function(formula, data, sigma_mult = 2.5,
                     max_iter = 50, tol = 1e-7, na.action = na.fail) {
  .lm_wpmm(formula, data, degree = 3, sigma_mult = sigma_mult,
           max_iter = max_iter, tol = tol, na.action = na.action,
           call = match.call())
}

# =============================================================================
# S4 Methods for WPMMfit
# =============================================================================

#' Methods for WPMMfit objects
#'
#' Accessors, prediction and inference for fits returned by
#' \code{\link{lm_wpmm2}} and \code{\link{lm_wpmm3}}.
#'
#' @param object WPMMfit object.
#' @param newdata optional data frame with the predictor variables.
#' @param parm character or integer vector of parameters to include.
#' @param level confidence level (default 0.95).
#' @param ... additional arguments (not used).
#'
#' @return \code{coef}, \code{residuals}, \code{fitted} and \code{predict}
#'   return numeric vectors; \code{vcov} the sandwich covariance matrix
#'   described in \code{\link{lm_wpmm2}}; \code{confint} a matrix of
#'   normal-approximation limits; \code{nobs} an integer; \code{summary}
#'   and \code{show} print and return the object invisibly.
#'
#' @name WPMMfit-methods
NULL

#' @rdname WPMMfit-methods
#' @export
setMethod("coef", "WPMMfit",
          function(object, ...) {
            object@coefficients
          })

#' @rdname WPMMfit-methods
#' @export
setMethod("residuals", "WPMMfit",
          function(object, ...) {
            object@residuals
          })

#' @rdname WPMMfit-methods
#' @export
setMethod("fitted", "WPMMfit",
          function(object, ...) {
            as.vector(attr(object, "model_matrix") %*% object@coefficients)
          })

#' @rdname WPMMfit-methods
#' @export
setMethod("predict", "WPMMfit",
          function(object, newdata = NULL, ...) {
            if (is.null(newdata)) {
              return(fitted(object))
            }
            tt <- stats::delete.response(stats::terms(attr(object, "model_frame")))
            X <- stats::model.matrix(tt, stats::model.frame(tt, newdata))
            as.vector(X %*% object@coefficients[colnames(X)])
          })

#' @rdname WPMMfit-methods
#' @export
setMethod("nobs", "WPMMfit",
          function(object, ...) {
            length(object@residuals)
          })

#' @rdname WPMMfit-methods
#' @export
setMethod("vcov", "WPMMfit",
          function(object, ...) {
            X <- attr(object, "model_matrix")
            if (is.null(X))
              stop("model_matrix not found; refit the model with lm_wpmm2() or lm_wpmm3()")
            r <- object@residuals
            ww <- .wpmm_weights(r, object@degree, object@sigma_mult)
            if (is.null(ww))
              stop("the window carries no mass at the final residuals")
            A <- crossprod(X, X * ww$dpsi)
            B <- crossprod(X * (ww$W * r))
            A_inv <- tryCatch(solve(A), error = function(e) NULL)
            if (is.null(A_inv)) {
              warning("Sandwich bread is singular; returning NA covariance.")
              V <- matrix(NA_real_, ncol(X), ncol(X))
            } else {
              V <- A_inv %*% B %*% t(A_inv)
            }
            dimnames(V) <- list(names(object@coefficients),
                                names(object@coefficients))
            V
          })

#' @rdname WPMMfit-methods
#' @export
setMethod("confint", "WPMMfit",
          function(object, parm, level = 0.95, ...) {
            cf <- object@coefficients
            se <- sqrt(diag(vcov(object)))
            a  <- (1 - level) / 2
            ci <- cf + outer(se, stats::qnorm(c(a, 1 - a)))
            colnames(ci) <- paste0(format(100 * c(a, 1 - a), trim = TRUE), " %")
            if (!missing(parm)) ci <- ci[parm, , drop = FALSE]
            ci
          })

#' @rdname WPMMfit-methods
#' @export
setMethod("summary", "WPMMfit",
          function(object, ...) {
            label <- if (object@degree == 2) "WPMM2" else "WPMM3"
            cat(label, " estimation results (weak-moment PMM, S = ",
                object@degree, ")\n", sep = "")
            cat("Call:\n")
            print(object@call)
            cat("\n")

            se <- sqrt(diag(vcov(object)))
            z  <- object@coefficients / se
            tab <- cbind(Estimate     = object@coefficients,
                         `Std. Error` = se,
                         `z value`    = z,
                         `Pr(>|z|)`   = 2 * stats::pnorm(-abs(z)))
            cat("Coefficients (sandwich standard errors):\n")
            stats::printCoefmat(tab, P.values = TRUE, has.Pvalue = TRUE)

            cat("\nWindow: centre =", format(object@window_center, digits = 4),
                ", sd =", format(object@window_scale, digits = 4),
                paste0("(sigma_mult = ", object@sigma_mult, ")"), "\n")
            cat("Windowed moments of final residuals:\n")
            print(object@weak_moments, digits = 4)
            coef_name <- if (object@degree == 2) "a2" else "kappa^w"
            if (is.na(object@correction)) {
              cat(coef_name, "= NA (mesokurtic window: Welsch estimator)\n")
            } else {
              cat(coef_name, "=", format(object@correction, digits = 4), "\n")
            }
            cat("Zero-weight (rejected) observations:",
                sum(object@weights == 0, na.rm = TRUE), "of",
                length(object@weights), "\n\n")

            cat("Algorithm information:\n")
            cat("  Convergence status:", object@convergence, "\n")
            cat("  Iterations:", object@iterations, "\n\n")
            invisible(object)
          })

#' @rdname WPMMfit-methods
#' @export
setMethod("show", "WPMMfit",
          function(object) {
            digits <- max(3L, getOption("digits") - 3L)
            cat("\nCall:\n",
                paste(deparse(object@call), sep = "\n", collapse = "\n"),
                "\n\n", sep = "")
            cat("Coefficients:\n")
            print.default(format(object@coefficients, digits = digits),
                          print.gap = 2L, quote = FALSE)
            conv_str <- if (isTRUE(object@convergence)) "converged" else "did NOT converge"
            cat("\nWPMM", object@degree, " algorithm: ", conv_str,
                " in ", object@iterations, " iteration",
                if (object@iterations == 1L) "" else "s", ".\n", sep = "")
            invisible(object)
          })
