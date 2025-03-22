# pmm_utils.R

#' PMM2 fitting algorithm - unified implementation
#'
#' Core iterative algorithm for PMM2 parameter estimation
#'
#' @param b_init initial parameter estimates (typically from OLS)
#' @param X design matrix
#' @param y response vector
#' @param m2,m3,m4 central moments
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
#' @param regularize logical, add small value to diagonal for numerical stability
#' @param reg_lambda regularization parameter (if regularize=TRUE)
#' @param verbose logical, whether to print progress information
#'
#' @return A list with components:
#'   \item{b}{estimated parameters}
#'   \item{convergence}{logical convergence status}
#'   \item{iterations}{number of iterations performed}
#'
#' @keywords internal
.pmm2_fit <- function(b_init, X, y,
                      m2, m3, m4,
                      max_iter=50, tol=1e-6,
                      regularize=TRUE, reg_lambda=1e-8,
                      verbose=FALSE) {
  n <- nrow(X)
  p <- ncol(X)

  # Compute A, B, C for the entire dataset at once
  A <- m3
  B <- (m4 - m2^2) - 2*m3*y
  C <- m3*(y^2) - ((m4 - m2^2)*y) - m2*m3

  # Current parameter estimates
  b_cur <- b_init
  converged <- FALSE
  iterations <- 0

  # Track convergence history if verbose
  if(verbose) {
    conv_history <- numeric(max_iter)
  }

  for (iter in seq_len(max_iter)) {
    iterations <- iter

    # Compute Yx = X %*% b_cur (predictor)
    # More efficient matrix multiplication
    Yx <- X %*% b_cur

    # Compute Z1 = A*Yx^2 + B*Yx + C
    Z1 <- A*(Yx^2) + B*Yx + C

    # Form Z - vector of equations
    # More efficient matrix operations
    Z <- c(sum(Z1), crossprod(X[,-1, drop=FALSE], Z1))

    # Form JZs - Jacobian matrix
    JZ11 <- 2*A*Yx + B

    # Vectorized computation
    JZs <- matrix(NA, p, p)
    JZs[1,1] <- sum(JZ11)

    for (ii in 2:p) {
      tmp <- JZ11*X[,ii]
      JZs[1,ii] <- sum(tmp)
      JZs[ii,1] <- JZs[1,ii]
    }

    for (ii in 2:p) {
      for (jj in 2:p) {
        tmp <- JZ11 * X[,ii] * X[,jj]
        JZs[ii,jj] <- sum(tmp)
      }
    }

    # Apply regularization if needed to avoid singularity
    if (regularize) {
      diag(JZs) <- diag(JZs) + reg_lambda
    }

    # Solve system JZs * delta = Z
    step <- tryCatch({
      solve(JZs, Z)
    }, error = function(e) {
      if(verbose) {
        cat("Error solving linear system in iteration", iter, ":", conditionMessage(e), "\n")
      }
      # Fallback to pseudoinverse or more stable solver
      if(requireNamespace("MASS", quietly = TRUE)) {
        MASS::ginv(JZs) %*% Z
      } else {
        warning("Failed to solve linear system. Consider installing 'MASS' package.")
        rep(NA, length(Z))
      }
    })

    # Check for numerical problems
    if (any(is.na(step)) || any(is.infinite(step))) {
      warning("Numerical problems encountered in iteration ", iter,
              ". Algorithm stopped.")
      break
    }

    # Update parameters
    b_new <- b_cur - step
    diff_par <- sqrt(sum((b_new - b_cur)^2))

    # Store convergence history if verbose
    if(verbose) {
      conv_history[iter] <- diff_par
      if(iter %% 5 == 0 || iter == 1) {
        cat("Iteration", iter, ": Parameter change =",
            formatC(diff_par, digits=8), "\n")
      }
    }

    b_cur <- b_new

    # Check convergence
    if (diff_par < tol) {
      converged <- TRUE
      if(verbose) cat("Converged after", iter, "iterations\n")
      break
    }
  }

  # Warning if max iterations reached without convergence
  if(!converged && verbose) {
    cat("Warning: Algorithm did not converge after", max_iter, "iterations\n")
  }

  # Plot convergence history if verbose
  if(verbose && iterations > 1) {
    if(requireNamespace("graphics", quietly = TRUE)) {
      graphics::plot(1:iterations, conv_history[1:iterations], type="b",
                     xlab="Iteration", ylab="Parameter change",
                     main="Convergence history")
      graphics::abline(h=tol, col="red", lty=2)
    }
  }

  # Return results
  list(
    b = b_cur,
    convergence = converged,
    iterations = iterations
  )
}

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

#' Compare PMM2 with OLS
#'
#' @param formula Model formula
#' @param data Data frame
#' @param pmm2_args List of arguments to pass to lm_pmm2()
#'
#' @return A list with OLS and PMM2 fit objects
#' @export
compare_with_ols <- function(formula, data, pmm2_args = list()) {
  # Fit OLS model
  fit_ols <- lm(formula, data)

  # Fit PMM2 model with default or specified arguments
  args <- c(list(formula = formula, data = data), pmm2_args)
  fit_pmm2 <- do.call(lm_pmm2, args)

  # Extract and compare coefficients
  coef_ols <- coef(fit_ols)
  coef_pmm2 <- coef(fit_pmm2)

  # Compute residual statistics
  res_ols <- residuals(fit_ols)
  res_pmm2 <- residuals(fit_pmm2)

  res_stats <- data.frame(
    Method = c("OLS", "PMM2"),
    RSS = c(sum(res_ols^2), sum(res_pmm2^2)),
    MAE = c(mean(abs(res_ols)), mean(abs(res_pmm2))),
    Skewness = c(pmm_skewness(res_ols), pmm_skewness(res_pmm2)),
    Kurtosis = c(pmm_kurtosis(res_ols), pmm_kurtosis(res_pmm2))
  )

  # Create comparison table of coefficients
  coef_table <- data.frame(
    Coefficient = names(coef_ols),
    OLS = coef_ols,
    PMM2 = coef_pmm2[names(coef_ols)],
    Diff_Percent = 100 * (coef_pmm2[names(coef_ols)] - coef_ols) / abs(coef_ols)
  )

  return(list(
    ols = fit_ols,
    pmm2 = fit_pmm2,
    coefficients = coef_table,
    residual_stats = res_stats
  ))
}
