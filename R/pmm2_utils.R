# pmm2_utils.R

#' Calculate kurtosis from data
#'
#' @param x numeric vector
#' @param excess logical, whether to return excess kurtosis (kurtosis - 3)
#'
#' @return Kurtosis value
#' @export
pmm_kurtosis <- function(x, excess = TRUE) {
  # Remove NAs
  x <- x[!is.na(x)]
  n <- length(x)

  if(n < 4) {
    warning("At least 4 non-missing values are needed to compute kurtosis")
    return(NA)
  }

  # Center the data
  x_centered <- x - mean(x)

  # Calculate moments
  m2 <- mean(x_centered^2)
  m4 <- mean(x_centered^4)

  # Calculate kurtosis
  kurt <- m4 / (m2^2)

  # Return excess kurtosis if requested
  if(excess) {
    kurt <- kurt - 3
  }

  return(kurt)
}

#' Calculate skewness from data
#'
#' @param x numeric vector
#'
#' @return Skewness value
#' @export
pmm_skewness <- function(x) {
  # Remove NAs
  x <- x[!is.na(x)]
  n <- length(x)

  if(n < 3) {
    warning("At least 3 non-missing values are needed to compute skewness")
    return(NA)
  }

  # Center the data
  x_centered <- x - mean(x)

  # Calculate moments
  m2 <- mean(x_centered^2)
  m3 <- mean(x_centered^3)

  # Calculate skewness
  skew <- m3 / (m2^(3/2))

  return(skew)
}


#' Calculate moments and cumulants of error distribution
#'
#' @param errors numeric vector of errors
#' @return list with moments, cumulants and theoretical variance reduction coefficient
#' @export
compute_moments <- function(errors) {
  m2 <- mean(errors^2)
  m3 <- mean(errors^3)
  m4 <- mean(errors^4)

  c3 <- m3 / m2^(3/2)  # Skewness coefficient
  c4 <- m4 / m2^2 - 3  # Excess kurtosis coefficient

  # Theoretical variance reduction coefficient
  g <- 1 - c3^2 / (2 + c4)

  return(list(m2 = m2, m3 = m3, m4 = m4,
              c3 = c3, c4 = c4,
              g = g))
}

#' Calculate theoretical skewness, kurtosis coefficients and variance reduction factor
#'
#' @param m2,m3,m4 central moments of second, third and fourth orders
#'
#' @return List with fields `c3`, `c4` and `g`
#' @export
pmm2_variance_factor <- function(m2, m3, m4) {
  if (is.na(m2) || m2 <= 0) {
    return(list(c3 = NA_real_, c4 = NA_real_, g = NA_real_))
  }
  c3 <- m3 / m2^(3/2)
  c4 <- m4 / m2^2 - 3
  denom <- 2 + c4
  g <- if (is.na(denom) || denom == 0) NA_real_ else 1 - (c3^2) / denom
  list(c3 = c3, c4 = c4, g = g)
}

#' Calculate theoretical variance matrices for OLS and PMM2
#'
#' @param X Design matrix with column of ones
#' @param m2,m3,m4 central moments of OLS residuals
#'
#' @return List with fields `ols`, `pmm2`, `c3`, `c4`, `g`
#' @export
pmm2_variance_matrices <- function(X, m2, m3, m4) {
  X <- as.matrix(X)
  XtX <- crossprod(X)
  V1 <- tryCatch({
    m2 * solve(XtX)
  }, error = function(e) {
    stop("Failed to invert matrix X'X: ", conditionMessage(e), call. = FALSE)
  })

  vf <- pmm2_variance_factor(m2, m3, m4)
  if (is.na(vf$g)) {
    V2 <- matrix(NA_real_, nrow = nrow(V1), ncol = ncol(V1))
  } else {
    V2 <- vf$g * V1
  }

  list(ols = V1, pmm2 = V2, c3 = vf$c3, c4 = vf$c4, g = vf$g)
}

# Covariance term that the plug-in moments add to a PMM regression fit.
# The polynomial is built from the central moments of the OLS residuals, and
# the error in alpha_2-hat (which is O(n^-1/2)) shifts the fit along
# a = (X'X)^-1 X'1, i.e. the intercept. Expanding the estimating equation
# around the true moments gives
#   beta-hat - beta = (X'X)^-1 X'u + a (eps-bar - u-bar),  u = psi(eps) / E psi',
# with Var(u) = g m2, Cov(u, eps) = g m2 and hence
#   Var = g m2 (X'X)^-1 + (1 - g) m2 / n * a a'.
# The same holds for every degree S. With an intercept a is its unit vector,
# so the slopes keep g m2 (X'X)^-1 and the intercept is close to OLS.
.pmm_location_term <- function(X, m2, g) {
  X <- as.matrix(X)
  a <- solve(crossprod(X), colSums(X))
  (1 - g) * m2 / nrow(X) * tcrossprod(a)
}
