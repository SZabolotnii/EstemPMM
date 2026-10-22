# pmm3_solver.R - Internal Newton-Raphson solver for PMM3

#' PMM3 Newton-Raphson solver (vectorised)
#'
#' Solves the PMM3 score equations using Newton-Raphson iteration.
#' The score is \eqn{Z = X'(\varepsilon(\kappa - \varepsilon^2))} and the
#' Jacobian is \eqn{J = X' \mathrm{diag}(3\varepsilon^2 - \kappa) X}.
#'
#' @param B_ols     OLS starting values (length P)
#' @param X         Design matrix (n x P)
#' @param Y         Response vector (length n)
#' @param kappa     Moment ratio (scalar)
#' @param adaptive  Re-estimate kappa each iteration? (default FALSE)
#' @param tol       Convergence tolerance on ||delta||_2 (default 1e-6)
#' @param max_iter  Maximum NR iterations (default 100)
#' @param step_max  Maximum step size (default 5.0)
#'
#' @return List: B (estimates), converged, iter, kappa, near_gaussian
#' @keywords internal
.pmm3_nr_solver <- function(B_ols, X, Y, kappa,
                             adaptive = FALSE,
                             tol      = 1e-6,
                             max_iter = 100,
                             step_max = 5.0) {

  P <- length(B_ols)
  B <- B_ols

  # Near-Gaussian: kappa undefined -> return OLS
  if (is.na(kappa)) {
    return(list(B = B_ols, converged = TRUE, iter = 0L,
                kappa = NA_real_, near_gaussian = TRUE))
  }

  converged <- FALSE
  iter      <- 0L

  for (i in seq_len(max_iter)) {
    iter <- as.integer(i)
    eps  <- as.numeric(Y - X %*% B)

    # Adaptive mode: re-estimate kappa from current residuals
    if (adaptive) {
      eps_c <- eps - mean(eps)
      m2_c  <- mean(eps_c^2)
      m4_c  <- mean(eps_c^4)
      m6_c  <- mean(eps_c^6)
      d     <- m4_c - 3 * m2_c^2
      if (abs(d) > .Machine$double.eps * 1e6) {
        kappa <- (m6_c - 3 * m4_c * m2_c) / d
      }
    }

    # Score: Z = X'(eps * (kappa - eps^2))
    Z1 <- eps * (kappa - eps^2)
    Z  <- as.numeric(crossprod(X, Z1))

    # Jacobian: J = X' diag(3*eps^2 - kappa) X
    w  <- 3 * eps^2 - kappa
    J  <- crossprod(X * w, X)

    # Regularise if near-singular
    if (rcond(J) < 1e-12) J <- J + diag(1e-8, P)

    # Solve
    delta <- tryCatch(as.numeric(solve(J, Z)),
                      error = function(e) rep(NA_real_, P))
    if (anyNA(delta)) break

    # Step-size limiting
    step_norm <- sqrt(sum(delta^2))
    if (step_norm > step_max) delta <- delta * step_max / step_norm

    B <- B - delta

    if (step_norm < tol) {
      converged <- TRUE
      break
    }
  }

  # Divergence guard: revert to OLS if too far
  if (sqrt(sum((B - B_ols)^2)) > 10 || anyNA(B)) {
    B         <- B_ols
    converged <- FALSE
  }

  list(B = B, converged = converged, iter = iter,
       kappa = kappa, near_gaussian = FALSE)
}


#' Coefficients of the general degree-S PMM polynomial
#'
#' Solves the Kunchenko normal system \eqn{F h = b} with
#' \eqn{F_{pq} = \alpha_{p+q} - \alpha_p \alpha_q} and
#' \eqn{b_k = k \alpha_{k-1}} (\eqn{\alpha_0 = 1}), where \eqn{\alpha_k} are
#' the central moments of the errors, and returns \eqn{h} together with the
#' variance-reduction factor \eqn{g_S = 1 / (\alpha_2\, b^\top F^{-1} b)}.
#'
#' @param alpha central moments \eqn{\alpha_1, \dots, \alpha_{2S}}
#'   (\eqn{\alpha_1 = 0}).
#' @param S polynomial degree.
#'
#' @return List: h (length S, or NULL if F is singular) and g.
#' @keywords internal
.pmm_poly_coef <- function(alpha, S) {
  a <- function(k) if (k == 0) 1 else alpha[k]
  i <- seq_len(S)
  F <- outer(i, i, Vectorize(function(p, q) a(p + q) - a(p) * a(q)))
  b <- i * vapply(i - 1, a, numeric(1))
  h <- tryCatch(as.numeric(solve(F, b)), error = function(e) NULL)
  if (is.null(h) || any(!is.finite(h))) return(list(h = NULL, g = NA_real_))
  list(h = h, g = 1 / (alpha[2] * sum(b * h)))
}

#' Newton solver for the general (asymmetric) degree-3 PMM estimator
#'
#' Solves \eqn{\sum_i x_i \psi(\varepsilon_i) = 0} with
#' \eqn{\psi(e) = \sum_{k=1}^{3} h_k (e^k - \alpha_k)}, safeguarded like
#' \code{pmm2_algorithm()}: backtracking on the score norm, and a return to
#' the starting values when the line search stalls away from a root or the
#' iterate is not finite.
#'
#' @param B_ols    starting values (OLS).
#' @param X        design matrix.
#' @param Y        response vector.
#' @param alpha    central moments \eqn{\alpha_1..\alpha_6} of the OLS
#'   residuals.
#' @param adaptive re-estimate the moments and \eqn{h} from the current
#'   residuals at every iteration.
#' @param tol      convergence tolerance on the Newton step norm.
#' @param max_iter maximum number of iterations.
#'
#' @return List: B, converged, iter, h, alpha, jacobian_min_eigen.
#' @keywords internal
.pmm3_general_solver <- function(B_ols, X, Y, alpha, adaptive = FALSE,
                                 tol = 1e-6, max_iter = 100) {
  central <- function(e) { e <- e - mean(e); vapply(1:6, function(k) mean(e^k), numeric(1)) }
  psi  <- function(e, h, al) h[1] * (e - al[1]) + h[2] * (e^2 - al[2]) + h[3] * (e^3 - al[3])
  dpsi <- function(e, h) h[1] + 2 * h[2] * e + 3 * h[3] * e^2
  score <- function(B, h, al) as.numeric(crossprod(X, psi(as.numeric(Y - X %*% B), h, al)))

  h <- .pmm_poly_coef(alpha, 3)$h
  if (is.null(h)) {
    warning("PMM3 normal system is singular; returning the initial estimate.")
    return(list(B = B_ols, converged = FALSE, iter = 0L, h = rep(NA_real_, 3),
                alpha = alpha, jacobian_min_eigen = NA_real_))
  }

  B <- B_ols
  converged <- FALSE
  failed <- NULL
  iter <- 0L
  norm_init <- NA_real_
  for (i in seq_len(max_iter)) {
    iter <- i
    e <- as.numeric(Y - X %*% B)
    if (adaptive) {
      alpha_i <- central(e)
      h_i <- .pmm_poly_coef(alpha_i, 3)$h
      if (!is.null(h_i)) { alpha <- alpha_i; h <- h_i }
    }
    G <- as.numeric(crossprod(X, psi(e, h, alpha)))
    J <- crossprod(X * dpsi(e, h), X)
    delta <- tryCatch(as.numeric(solve(J, G)), error = function(err) NULL)
    if (is.null(delta) || any(!is.finite(delta))) {
      failed <- "singular Jacobian"
      break
    }
    norm_cur <- sqrt(sum(G^2))
    if (i == 1) norm_init <- norm_cur
    step_norm <- sqrt(sum(delta^2))
    t_step <- 1
    B_new <- B + delta
    if (step_norm >= tol) {
      repeat {
        norm_new <- sqrt(sum(score(B_new, h, alpha)^2))
        if (all(is.finite(B_new)) && is.finite(norm_new) &&
            norm_new <= (1 - 1e-4 * t_step) * norm_cur) break
        t_step <- t_step / 2
        if (t_step < 2^-30) break
        B_new <- B + t_step * delta
      }
    }
    if (t_step < 2^-30) {
      if (norm_cur <= 1e-10 * norm_init) converged <- TRUE else failed <- "line search stalled"
      break
    }
    B <- B_new
    if (step_norm < tol) {
      converged <- TRUE
      break
    }
  }
  if (is.null(failed) && any(!is.finite(B))) failed <- "non-finite estimate"
  if (!is.null(failed)) {
    warning("PMM3 iteration failed (", failed, "); returning the initial estimate.")
    return(list(B = B_ols, converged = FALSE, iter = iter, h = h, alpha = alpha,
                jacobian_min_eigen = NA_real_))
  }
  if (!converged) {
    warning("PMM3 algorithm did not converge after ", max_iter, " iterations.")
  }
  e <- as.numeric(Y - X %*% B)
  J <- crossprod(X * dpsi(e, h), X) / length(e)
  min_eig <- min(eigen((J + t(J)) / 2, symmetric = TRUE, only.values = TRUE)$values)
  if (converged && min_eig <= 0) {
    warning("PMM3 root has a Jacobian that is not positive definite; ",
            "it may not be the consistent root.")
  }
  list(B = B, converged = converged, iter = iter, h = h, alpha = alpha,
       jacobian_min_eigen = min_eig)
}


#' PMM3 Newton-Raphson solver for nonlinear models
#'
#' Uses a residual function and numerical Jacobian (via \code{numDeriv})
#' instead of a fixed design matrix. This is necessary for ARIMA and other
#' models where the linearized design matrix does not capture full dynamics.
#'
#' The PMM3 estimating equation is:
#' \deqn{\sum_t \varepsilon_t(\kappa - \varepsilon_t^2) \cdot
#'       \frac{\partial \varepsilon_t}{\partial \theta} = 0}
#'
#' @param theta_init Initial parameter estimates (from MLE/CSS)
#' @param fn_residuals Function(theta) returning residual vector
#' @param kappa Moment ratio (scalar)
#' @param adaptive Re-estimate kappa each iteration? (default FALSE)
#' @param tol Convergence tolerance (default 1e-6)
#' @param max_iter Maximum NR iterations (default 100)
#' @param step_max Maximum step size (default 5.0)
#'
#' @return List: B (estimates), converged, iter, kappa, near_gaussian
#' @keywords internal
.pmm3_nr_solver_nonlinear <- function(theta_init, fn_residuals, kappa,
                                       adaptive = FALSE,
                                       tol      = 1e-6,
                                       max_iter = 100,
                                       step_max = 5.0) {

  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    stop("Package 'numDeriv' is required for ARIMA PMM3. Please install it.",
         call. = FALSE)
  }

  P <- length(theta_init)
  B <- theta_init

  # Near-Gaussian: kappa undefined -> return initial
  if (is.na(kappa)) {
    return(list(B = theta_init, converged = TRUE, iter = 0L,
                kappa = NA_real_, near_gaussian = TRUE))
  }

  converged <- FALSE
  iter      <- 0L

  # Determine valid residual indices from initial residuals
  eps_full <- fn_residuals(B)
  valid_idx <- which(is.finite(eps_full))
  n_valid   <- length(valid_idx)

  # Wrapper that always returns fixed-length vector (only valid indices)
  fn_res_valid <- function(th) {
    r <- fn_residuals(th)
    out <- r[valid_idx]
    out[!is.finite(out)] <- 0  # safety: replace any new NAs with 0
    out
  }

  for (i in seq_len(max_iter)) {
    iter <- as.integer(i)

    # Current residuals (fixed length)
    eps <- fn_res_valid(B)

    # Adaptive mode: re-estimate kappa
    if (adaptive) {
      eps_c <- eps - mean(eps)
      m2_c  <- mean(eps_c^2)
      m4_c  <- mean(eps_c^4)
      m6_c  <- mean(eps_c^6)
      d     <- m4_c - 3 * m2_c^2
      if (abs(d) > .Machine$double.eps * 1e6) {
        kappa <- (m6_c - 3 * m4_c * m2_c) / d
      }
    }

    # Numerical Jacobian: J_res[i,j] = d(eps_i)/d(theta_j)
    J_res <- numDeriv::jacobian(fn_res_valid, B, method = "Richardson")

    # Score: gradient = J_res' * eps*(kappa - eps^2)
    Z1       <- eps * (kappa - eps^2)
    gradient <- as.numeric(crossprod(J_res, Z1))

    # Hessian approximation: H = J_res' * diag(3*eps^2 - kappa) * J_res
    w <- 3 * eps^2 - kappa
    H <- crossprod(J_res * w, J_res)

    # Regularise
    if (rcond(H) < 1e-12) H <- H + diag(1e-8, P)

    # Solve for delta
    delta <- tryCatch(as.numeric(solve(H, gradient)),
                      error = function(e) rep(NA_real_, P))
    if (anyNA(delta)) break

    # Step-size limiting
    step_norm <- sqrt(sum(delta^2))
    if (step_norm > step_max) delta <- delta * step_max / step_norm

    # NB: + not - because J_res = d(eps)/d(theta), and the design matrix
    # analog D = -J_res. The NR step is theta - H^{-1}*D'*Z1,
    # which equals theta + H^{-1}*J_res'*Z1 = theta + delta.
    B <- B + delta

    if (step_norm < tol) {
      converged <- TRUE
      break
    }
  }

  # Divergence guard
  if (sqrt(sum((B - theta_init)^2)) > 10 || anyNA(B)) {
    B         <- theta_init
    converged <- FALSE
  }

  list(B = B, converged = converged, iter = iter,
       kappa = kappa, near_gaussian = FALSE)
}
