# pmm_unified_api.R - Unified, method-dispatching user-facing entry
# points for the package. These thin wrappers consolidate the older
# family of separate `lm_pmm2()` / `lm_pmm3()` / `ar_pmm2()` / ... /
# `arima_pmm3()` functions into a single dispatcher per model class,
# selecting the polynomial order via a `method` argument.
#
# Introduced in EstemPMM 0.4.0. As of 0.5.0 the old names are still
# exported and emit no deprecation message; they are documented as the
# lower-level fitters these wrappers call. The 0.4.0 NEWS entry that
# scheduled warnings for 0.4.1 and removal for 0.5.0 was not acted on,
# and the JSS manuscript's examples call them directly.

#' Fit a linear model with the polynomial maximization method
#'
#' Unified entry point for PMM2 and PMM3 linear regression. Dispatches
#' to \code{\link{lm_pmm2}} or \code{\link{lm_pmm3}} based on the
#' \code{method} argument. With \code{method = "auto"}, the
#' \code{\link{pmm_dispatch}} rule is applied to the residuals of an
#' initial OLS fit to select the recommended method.
#'
#' @param formula R formula for the model.
#' @param data data.frame containing the variables in \code{formula}.
#' @param method One of \code{"auto"}, \code{"pmm2"}, \code{"pmm3"},
#'   \code{"wpmm2"} or \code{"wpmm3"}.
#'   \code{"auto"} consults \code{\link{pmm_dispatch}} and falls back to
#'   PMM2 when the dispatcher recommends OLS, so the return value is
#'   always a \code{\linkS4class{PMMfit}} object. In that case a message
#'   says so. The dispatcher's decision is stored in the attribute
#'   \code{"dispatch"} of the fit. \code{"auto"} never selects the
#'   weak-moment estimators.
#' @param ... Additional arguments forwarded to the underlying fitter
#'   (\code{\link{lm_pmm2}}, \code{\link{lm_pmm3}}, \code{\link{lm_wpmm2}}
#'   or \code{\link{lm_wpmm3}}). With \code{method = "auto"}, named
#'   arguments the selected fitter does not take (for example
#'   \code{reg_lambda} when PMM3 is selected) are dropped with a message.
#'
#' @return An S4 object of class \code{\linkS4class{PMM2fit}},
#'   \code{\linkS4class{PMM3fit}} or \code{\linkS4class{WPMMfit}}, all
#'   inheriting from \code{\linkS4class{PMMfit}}.
#'
#' @seealso \code{\link{pmm_dispatch}}, \code{\link{pmm_arima}},
#'   \code{\link{lm_pmm2}}, \code{\link{lm_pmm3}}.
#'
#' @examples
#' set.seed(123)
#' n <- 80
#' x <- rnorm(n)
#' y <- 2 + 3 * x + rgamma(n, 2, 1) - 2
#' dat <- data.frame(y = y, x = x)
#' fit <- pmm_lm(y ~ x, data = dat, method = "pmm2")
#' fit
#'
#' @export
pmm_lm <- function(formula, data,
                   method = c("auto", "pmm2", "pmm3", "wpmm2", "wpmm3"),
                   ...) {
  method <- match.arg(method)
  outer_call <- match.call()
  dec <- NULL
  if (method == "auto") {
    ols    <- stats::lm(formula, data = data)
    dec    <- pmm_dispatch(stats::residuals(ols), verbose = FALSE)
    chosen <- tolower(dec$method)
    method <- if (chosen == "pmm3") "pmm3" else "pmm2"
    if (chosen == "ols") {
      message("pmm_dispatch() recommends OLS (", dec$reasoning,
              "); fitting PMM2 as in earlier versions. The decision is in ",
              "attr(fit, \"dispatch\").")
    }
  }
  fitter <- switch(method, pmm2 = lm_pmm2, pmm3 = lm_pmm3,
                   wpmm2 = lm_wpmm2, wpmm3 = lm_wpmm3)
  args <- list(...)
  if (!is.null(dec) && length(args) > 0) {
    # With "auto" the caller cannot know which fitter runs, so an argument
    # of the other fitter (reg_lambda for lm_pmm2, step_max for lm_pmm3)
    # is dropped with a message instead of failing as unused.
    unused <- setdiff(names(args), names(formals(fitter)))
    unused <- unused[nzchar(unused)]
    if (length(unused) > 0) {
      message("pmm_lm(method = \"auto\") selected ", method, "; ignoring ",
              "argument(s) it does not take: ",
              paste(unused, collapse = ", "), ".")
      args <- args[setdiff(names(args), unused)]
    }
  }
  fit <- do.call(fitter, c(list(formula = formula, data = data), args))
  # Overwrite the inner call so that predict()/show() see the formula the
  # user actually wrote, not the symbol "formula" forwarded by this wrapper.
  fit@call <- outer_call
  if (!is.null(dec)) attr(fit, "dispatch") <- dec
  fit
}

#' Fit an AR model with the polynomial maximization method
#'
#' Unified entry point for PMM2 and PMM3 autoregressive estimation.
#'
#' @param x univariate time series or numeric vector.
#' @param order non-negative integer giving the autoregressive order.
#' @param method One of \code{"pmm2"} or \code{"pmm3"}.
#' @param ... Additional arguments forwarded to \code{\link{ar_pmm2}}
#'   or \code{\link{ar_pmm3}}.
#'
#' @return An S4 object of class \code{\linkS4class{ARPMM2}} (for
#'   \code{method = "pmm2"}) or \code{\linkS4class{ARPMM3}}.
#'
#' @export
pmm_ar <- function(x, order, method = c("pmm2", "pmm3"), ...) {
  method <- match.arg(method)
  outer_call <- match.call()
  fit <- switch(method,
    pmm2 = ar_pmm2(x, order = order, ...),
    pmm3 = ar_pmm3(x, order = order, ...)
  )
  fit@call <- outer_call
  fit
}

#' Fit an MA model with the polynomial maximization method
#'
#' Unified entry point for PMM2 and PMM3 moving-average estimation.
#'
#' @param x univariate time series or numeric vector.
#' @param order non-negative integer giving the moving-average order.
#' @param method One of \code{"pmm2"} or \code{"pmm3"}.
#' @param ... Additional arguments forwarded to \code{\link{ma_pmm2}}
#'   or \code{\link{ma_pmm3}}.
#'
#' @return An S4 object of class \code{\linkS4class{MAPMM2}} or
#'   \code{\linkS4class{MAPMM3}}.
#'
#' @export
pmm_ma <- function(x, order, method = c("pmm2", "pmm3"), ...) {
  method <- match.arg(method)
  outer_call <- match.call()
  fit <- switch(method,
    pmm2 = ma_pmm2(x, order = order, ...),
    pmm3 = ma_pmm3(x, order = order, ...)
  )
  fit@call <- outer_call
  fit
}

#' Fit an ARMA model with the polynomial maximization method
#'
#' Unified entry point for PMM2 and PMM3 ARMA estimation.
#'
#' @param x univariate time series or numeric vector.
#' @param order integer vector \code{c(p, q)} giving the AR and MA orders.
#' @param method One of \code{"pmm2"} or \code{"pmm3"}.
#' @param ... Additional arguments forwarded to \code{\link{arma_pmm2}}
#'   or \code{\link{arma_pmm3}}.
#'
#' @return An S4 object of class \code{\linkS4class{ARMAPMM2}} or
#'   \code{\linkS4class{ARMAPMM3}}.
#'
#' @export
pmm_arma <- function(x, order, method = c("pmm2", "pmm3"), ...) {
  method <- match.arg(method)
  outer_call <- match.call()
  fit <- switch(method,
    pmm2 = arma_pmm2(x, order = order, ...),
    pmm3 = arma_pmm3(x, order = order, ...)
  )
  fit@call <- outer_call
  fit
}

#' Fit an ARIMA model with the polynomial maximization method
#'
#' Unified entry point for PMM2 and PMM3 ARIMA estimation. This is the
#' recommended primary API for non-Gaussian ARIMA modelling and the
#' interface featured in the JSS submission accompanying the package.
#'
#' @param x univariate time series or numeric vector.
#' @param order integer vector \code{c(p, d, q)} giving the AR order,
#'   degree of differencing, and MA order.
#' @param method One of \code{"pmm2"} or \code{"pmm3"}.
#' @param ... Additional arguments forwarded to \code{\link{arima_pmm2}}
#'   or \code{\link{arima_pmm3}}.
#'
#' @return An S4 object of class \code{\linkS4class{ARIMAPMM2}} or
#'   \code{\linkS4class{ARIMAPMM3}}.
#'
#' @examples
#' set.seed(2026)
#' x <- arima.sim(list(ar = 0.6), n = 200,
#'                rand.gen = function(k) rgamma(k, 2, 1) - 2)
#' pmm_arima(x, order = c(1, 0, 0), method = "pmm2")
#'
#' @export
pmm_arima <- function(x, order, method = c("pmm2", "pmm3"), ...) {
  method <- match.arg(method)
  outer_call <- match.call()
  fit <- switch(method,
    pmm2 = arima_pmm2(x, order = order, ...),
    pmm3 = arima_pmm3(x, order = order, ...)
  )
  fit@call <- outer_call
  fit
}

#' Fit a seasonal ARIMA model with the polynomial maximization method
#'
#' Unified entry point for SARIMA-PMM2 estimation. PMM3 seasonal
#' estimation is not yet implemented; \code{method = "pmm3"} therefore
#' raises an error. All arguments other than \code{method} are
#' forwarded verbatim to \code{\link{sarima_pmm2}} -- see that function
#' for the seasonal order specification.
#'
#' @param x univariate time series or numeric vector.
#' @param ... Arguments forwarded to \code{\link{sarima_pmm2}} (notably
#'   \code{order} and \code{seasonal}).
#' @param method Currently only \code{"pmm2"} is supported.
#'
#' @return An S4 object of class \code{\linkS4class{SARIMAPMM2}}.
#'
#' @export
pmm_sarima <- function(x, ..., method = c("pmm2", "pmm3")) {
  method <- match.arg(method)
  if (method == "pmm3") {
    stop("Seasonal PMM3 estimation is not implemented in this release.")
  }
  outer_call <- match.call()
  fit <- sarima_pmm2(x, ...)
  fit@call <- outer_call
  fit
}
