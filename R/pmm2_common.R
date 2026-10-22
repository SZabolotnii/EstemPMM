# pmm2_common.R - Common utilities for all PMM2 models


#' Universal PMM2 algorithm for all model types
#'
#' @param b_init Initial parameter estimates
#' @param X Design matrix
#' @param y Response vector
#' @param m2,m3,m4 Central moments
#' @param max_iter Maximum number of iterations
#' @param tol Tolerance for convergence
#' @param regularize Whether to add regularization
#' @param reg_lambda Regularization parameter
#' @param verbose Whether to print progress information
#' @param poly_terms Pre-computed polynomial coefficients (list with elements \code{A}, \code{B}, \code{C});
#'   allows passing custom values for special scenarios, otherwise they are computed from moments
#'
#' @return List with estimation results
#' @keywords internal
pmm2_algorithm <- function(b_init, X, y, m2, m3, m4,
                           max_iter = 50, tol = 1e-6,
                           regularize = TRUE, reg_lambda = 1e-8,
                           verbose = FALSE,
                           poly_terms = NULL) {
  # Current parameter estimates
  b_cur <- b_init
  converged <- FALSE
  iterations <- 0

  # Calculate PMM2 polynomial coefficients
  # Derivation corresponds to equation (10) from the paper:
  # A = c3 * sigma^3, B = (c4 + 3) * sigma^4 - sigma^4 - 2 c3 * sigma^3 * y,
  # C = c3 * y^2 * sigma^3 - ((c4 + 3) * sigma^4 - sigma^4) * y - sigma^2 * c3 * sigma^3.
  # Substituting sigma^2 = m2, m3 = c3 sigma^3, m4 = (c4 + 3) sigma^4 gives the formulas below.
  if (is.null(poly_terms)) {
    A <- m3
    B <- m4 - m2^2 - 2 * m3 * y
    C <- m3 * y^2 - (m4 - m2^2) * y - m2 * m3
  } else {
    A <- poly_terms$A
    B <- poly_terms$B
    C <- poly_terms$C
    if (length(B) != nrow(X) || length(C) != nrow(X)) {
      stop("poly_terms$B and poly_terms$C must have length equal to the number of rows in X")
    }
  }

  # Track convergence history if verbose
  if (verbose) {
    conv_history <- numeric(max_iter)
  }

  p <- length(b_cur)

  # Estimating function Z(b) = sum_t Z1_t x_t, Z1 = A*y_pred^2 + B*y_pred + C
  score <- function(b) {
    y_pred <- as.vector(X %*% b)
    Z1 <- A*(y_pred^2) + B*y_pred + C
    Z <- numeric(p)
    for (r in 1:p) {
      Z[r] <- sum(Z1 * X[, r])
    }
    Z
  }

  # Jacobian dZ/db = sum_t (2*A*y_pred + B) x_t x_t'
  jacobian <- function(b) {
    y_pred <- as.vector(X %*% b)
    JZ11 <- 2 * A * y_pred + B
    JZs <- matrix(0, p, p)
    for (i in 1:p) {
      for (j in 1:p) {
        JZs[i, j] <- sum(JZ11 * X[, i] * X[, j])
      }
    }
    JZs
  }

  stalled <- FALSE

  # Main iteration loop
  for (iter in seq_len(max_iter)) {
    iterations <- iter

    Z <- score(b_cur)
    JZs <- jacobian(b_cur)

    # Add regularization if needed
    if (regularize) {
      diag(JZs) <- diag(JZs) + reg_lambda
    }

    # Solve system JZs * delta = Z
    delta <- tryCatch({
      solve(JZs, Z)
    }, error = function(e) {
      if (verbose) {
        cat("Error solving linear system:", conditionMessage(e), "\n")
        cat("Adding stronger regularization\n")
      }
      diag(JZs) <- diag(JZs) + 1e-4
      solve(JZs, Z)
    })

    # Update parameters. The full Newton step is taken whenever it reduces
    # ||Z||, so a well-behaved fit follows exactly the 0.5.0 trajectory;
    # otherwise the step is halved (backtracking on the norm of the
    # estimating function, as in the recursive MA solver).
    norm_cur <- sqrt(sum(Z^2))
    if (iter == 1) norm_init <- norm_cur
    t_step <- 1
    b_new <- b_cur - delta
    # A step already below tol passes the 0.5.0 convergence test; near the
    # root ||Z|| sits at rounding level and need not decrease, so it is
    # taken as is.
    if (sqrt(sum(delta^2)) >= tol) {
      repeat {
        norm_new <- sqrt(sum(score(b_new)^2))
        if (all(is.finite(b_new)) && is.finite(norm_new) &&
            norm_new <= (1 - 1e-4 * t_step) * norm_cur) break
        t_step <- t_step / 2
        if (t_step < 2^-30) break
        b_new <- b_cur - t_step * delta
      }
    }
    if (t_step < 2^-30) {
      if (norm_cur <= 1e-10 * norm_init) {
        # ||Z|| is already at rounding level: this is the root
        converged <- TRUE
      } else {
        # No step along the Newton direction reduces ||Z||: stop here
        stalled <- TRUE
        if (verbose) cat("Line search failed at iteration", iter, "\n")
      }
      break
    }
    diff_par <- sqrt(sum((b_new - b_cur)^2))

    # Store convergence history if verbose
    if (verbose) {
      conv_history[iter] <- diff_par
      if (iter %% 5 == 0 || iter == 1) {
        cat("Iteration", iter, ": Parameter change =",
            formatC(diff_par, digits = 8), "\n")
      }
    }

    b_cur <- b_new

    # Check convergence
    if (diff_par < tol) {
      converged <- TRUE
      if (verbose) cat("Converged after", iter, "iterations\n")
      break
    }
  }

  # A stalled line search or a non-finite iterate has no estimate to offer:
  # return the initial (consistent) estimate, as the PMM3 solver does.
  if (stalled || !all(is.finite(b_cur))) {
    warning("PMM2 iteration failed (",
            if (stalled) "line search stalled" else "non-finite estimate",
            "); returning the initial estimate.", call. = FALSE)
    b_cur <- b_init
    converged <- FALSE
  } else if (!converged) {
    warning("PMM2 algorithm did not converge after ", iterations,
            " iterations.", call. = FALSE)
  }

  # Root check. At the consistent root the Jacobian
  # sum_t ((m4 - m2^2) - 2*m3*eps_t) x_t x_t' is positive definite in
  # expectation; the estimating equation can have other roots where it is
  # not. In small samples the consistent root can fail this test too, so it
  # is reported, not acted on.
  min_eig <- tryCatch(
    min(eigen(jacobian(b_cur), symmetric = TRUE, only.values = TRUE)$values),
    error = function(e) NA_real_)
  root_pd <- is.finite(min_eig) && min_eig > 0
  if (converged && !root_pd) {
    warning("PMM2 root has a Jacobian that is not positive definite; it may ",
            "not be the consistent root. Compare with the OLS fit.",
            call. = FALSE)
  }

  # Calculate final residuals
  final_res <- as.numeric(y - X %*% b_cur)

  # Plot convergence history if verbose
  if (verbose && iterations > 1) {
    if (requireNamespace("graphics", quietly = TRUE)) {
      graphics::plot(1:iterations, conv_history[1:iterations], type = "b",
                     xlab = "Iteration", ylab = "Parameter change",
                     main = "Convergence history")
      graphics::abline(h = tol, col = "red", lty = 2)
    }
  }

  # Return results
  list(
    b = as.numeric(b_cur),
    convergence = converged,
    iterations = iterations,
    residuals = final_res,
    jacobian_min_eigen = min_eig
  )
}
