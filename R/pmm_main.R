# pmm_main.R

#' pmm2: Master function for PMM2 (S=2)
#'
#' @param formula R formula
#' @param data data.frame
#' @param max_iter integer: iteration limit
#' @param tol numeric: tolerance for convergence
#' @param regularize logical: add small value to diagonal for numerical stability
#'
#' @details
#' 1) Fits OLS => b_ols
#' 2) computes m2,m3,m4 from OLS residuals
#' 3) calls internal iteration approach
#'
#' @return An S4 \code{PMM2fit} object
#' @export
lm_pmm2 <- function(formula, data,
                 max_iter=50, tol=1e-6,
                 regularize=TRUE)
{

  # 1) OLS
  mf <- model.frame(formula, data)
  X  <- model.matrix(formula, mf)
  y  <- model.response(mf)
  fit_ols <- lm.fit(x=X, y=y)
  b_ols   <- fit_ols$coefficients

  # 2) OLS residuals => m2,m3,m4
  res_ols <- y - (X %*% b_ols)
  m2 <- mean(res_ols^2)
  m3 <- mean(res_ols^3)
  m4 <- mean(res_ols^4)

  # 3) Викликаємо єдиний уніфікований метод
  out <- .pmm2_fit(b_ols, X, y, m2, m3, m4, max_iter, tol, regularize)

  # out => list(b, conv)
  b_est   <- out$b
  conv    <- out$convergence

  final_res <- as.numeric(y - X %*% b_est)

  # wrap into S4
  ans <- new("PMM2fit",
             coefficients=b_est,
             residuals=final_res,
             m2=m2, m3=m3, m4=m4,
             convergence=conv)
  ans
}
