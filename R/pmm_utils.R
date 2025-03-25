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
    Yx <- as.vector(X %*% b_cur)  # Ensure Yx is a vector

    # Compute Z1 = A*Yx^2 + B*Yx + C
    Z1 <- A*(Yx^2) + B*Yx + C

    # Form Z - vector of equations
    # Improved handling for multi-variable cases
    if (p > 1) {
      # Виконуємо crossprod для всіх предикторів, окрім перетину
      Z_rest <- crossprod(X[, -1, drop=FALSE], Z1)

      # Перевіряємо, чи Z_rest є матрицею (більше ніж 1 предиктор)
      if (is.matrix(Z_rest)) {
        Z_rest <- as.vector(Z_rest)
      }

      # Об'єднуємо результати
      Z <- c(sum(Z1), Z_rest)
    } else {
      # Якщо лише перетин, Z - це просто сума
      Z <- sum(Z1)
    }

    # Form JZs - Jacobian matrix
    JZ11 <- 2*A*Yx + B

    # Vectorized computation of Jacobian matrix
    JZs <- matrix(NA, p, p)
    JZs[1,1] <- sum(JZ11)

    # Перший рядок і перший стовпець Якобіана
    for (ii in 2:p) {
      tmp <- JZ11 * X[,ii]
      JZs[1,ii] <- sum(tmp)
      JZs[ii,1] <- JZs[1,ii]  # Якобіан симетричний
    }

    # Решта елементів Якобіана
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

  # Return results - ensure b is a numeric vector, not a matrix
  list(
    b = as.numeric(b_cur),
    convergence = converged,
    iterations = iterations
  )
}

#' PMM2 fitting algorithm for time series models
#'
#' Implementation of PMM2 algorithm for time series models with various approaches
#' depending on the model structure (AR, MA, ARMA, ARIMA)
#'
#' @param b_init initial parameter estimates (AR followed by MA)
#' @param innovations_init initial innovations (errors)
#' @param model_info list with model structure information
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
#'   \item{innovations}{final innovations (errors)}
#'
#' @keywords internal
.ts_pmm2_fit <- function(b_init, innovations_init, model_info,
                         m2, m3, m4,
                         max_iter = 50, tol = 1e-6,
                         regularize = TRUE, reg_lambda = 1e-8,
                         verbose = FALSE) {

  # Extract model information
  x <- model_info$x
  model_type <- model_info$model_type
  ar_order <- model_info$ar_order
  ma_order <- model_info$ma_order

  n <- length(x)
  total_params <- ar_order + ma_order

  # Current parameter estimates
  b_cur <- b_init
  converged <- FALSE
  iterations <- 0

  # Current innovations
  innovations <- innovations_init

  # Track convergence history if verbose
  if(verbose) {
    conv_history <- numeric(max_iter)
  }

  for(iter in seq_len(max_iter)) {
    iterations <- iter

    # Update innovations based on model type and current parameters
    if(model_type == "ar") {
      # For AR models, calculate predicted values
      X <- create_ar_matrix(x, ar_order)
      y <- x[(ar_order + 1):n]
      y_pred <- as.vector(X %*% b_cur)

      # Form matrices for the PMM algorithm
      A <- m3
      B <- (m4 - m2^2) - 2*m3*y
      C <- m3*(y^2) - ((m4 - m2^2)*y) - m2*m3

      # Calculate Z1 = A*yhat^2 + B*yhat + C
      Z1 <- A*(y_pred^2) + B*y_pred + C

      # Form Z vector for each parameter
      Z <- numeric(ar_order)
      for(i in 1:ar_order) {
        Z[i] <- sum(Z1 * X[, i])
      }

      # Form Jacobian matrix
      JZ11 <- 2*A*y_pred + B
      JZs <- matrix(0, ar_order, ar_order)

      for(i in 1:ar_order) {
        for(j in 1:ar_order) {
          JZs[i, j] <- sum(JZ11 * X[, i] * X[, j])
        }
      }

      # Update innovations
      innovations <- y - y_pred

    } else if(model_type %in% c("ma", "arma", "arima")) {
      # Extract AR and MA parts from current parameters
      ar_part <- if(ar_order > 0) b_cur[1:ar_order] else numeric(0)
      ma_part <- if(ma_order > 0) b_cur[(ar_order+1):total_params] else numeric(0)

      # For MA-based models, we need to update innovations using ARMA framework
      if(model_type == "ma") {
        # For pure MA, update innovations directly
        innovations <- update_ma_innovations(x, ma_part)
      } else {
        # For ARMA and ARIMA, use arima with fixed parameters
        arima_order <- c(ar_order, 0, ma_order) # d=0 as we're working with (differenced) x

        # Handle potential errors in model fitting
        curr_model <- tryCatch({
          stats::arima(x, order = arima_order,
                       include.mean = FALSE,
                       fixed = c(ar_part, ma_part))
        }, error = function(e) {
          if(verbose) {
            cat("Error in model fitting:", conditionMessage(e), "\n")
          }
          return(NULL)
        })

        if(is.null(curr_model)) {
          warning("Failed to update model in iteration ", iter, ". Algorithm stopped.")
          break
        }

        innovations <- residuals(curr_model)
      }

      # Prepare design matrices for both AR and MA parts
      J <- matrix(0, nrow = n, ncol = total_params)

      # For AR part, use lagged values of the series
      if(ar_order > 0) {
        for(i in 1:ar_order) {
          J[(i+1):n, i] <- x[1:(n-i)]
        }
      }

      # For MA part, use lagged innovations
      if(ma_order > 0) {
        for(i in 1:ma_order) {
          J[(i+1):n, ar_order+i] <- innovations[1:(n-i)]
        }
      }

      # Remove initial rows that can't be used for all lags
      start_idx <- max(ar_order, ma_order) + 1
      if(start_idx > 1) {
        J_valid <- J[start_idx:n, , drop = FALSE]
        x_valid <- x[start_idx:n]
        innovations_valid <- innovations[start_idx:n]
      } else {
        J_valid <- J
        x_valid <- x
        innovations_valid <- innovations
      }

      # Compute PMM2 coefficients
      A <- m3
      B <- (m4 - m2^2) - 2*m3*x_valid
      C <- m3*(x_valid^2) - ((m4 - m2^2)*x_valid) - m2*m3

      # Form Z vector for each parameter
      Z <- numeric(total_params)
      for(j in 1:total_params) {
        Z1 <- A*(innovations_valid^2) + B*innovations_valid + C
        Z[j] <- sum(Z1 * J_valid[, j])
      }

      # Form Jacobian matrix
      JZ11 <- 2*A*innovations_valid + B
      JZs <- matrix(0, total_params, total_params)

      for(i in 1:total_params) {
        for(j in 1:total_params) {
          JZs[i, j] <- sum(JZ11 * J_valid[, i] * J_valid[, j])
        }
      }
    }

    # Apply regularization if needed
    if(regularize) {
      diag(JZs) <- diag(JZs) + reg_lambda
    }

    # Solve system JZs * delta = Z
    step <- tryCatch({
      solve(JZs, Z)
    }, error = function(e) {
      if(verbose) {
        cat("Error solving linear system in iteration", iter, ":", conditionMessage(e), "\n")
      }
      # Fallback to pseudoinverse
      if(requireNamespace("MASS", quietly = TRUE)) {
        MASS::ginv(JZs) %*% Z
      } else {
        warning("Failed to solve linear system. Consider installing 'MASS' package.")
        rep(NA, length(Z))
      }
    })

    # Check for numerical problems
    if(any(is.na(step)) || any(is.infinite(step))) {
      warning("Numerical problems encountered in iteration ", iter, ". Algorithm stopped.")
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
            formatC(diff_par, digits = 8), "\n")
      }
    }

    b_cur <- b_new

    # Check convergence
    if(diff_par < tol) {
      converged <- TRUE
      if(verbose) cat("Converged after", iter, "iterations\n")
      break
    }
  }

  # Warning if max iterations reached without convergence
  if(!converged && verbose) {
    cat("Warning: Algorithm did not converge after", max_iter, "iterations\n")
  }

  # Final innovation update based on model type
  if(model_type == "ar") {
    X <- create_ar_matrix(x, ar_order)
    y <- x[(ar_order + 1):n]
    final_innovations <- y - X %*% b_cur
  } else {
    # Extract AR and MA parts
    ar_part <- if(ar_order > 0) b_cur[1:ar_order] else numeric(0)
    ma_part <- if(ma_order > 0) b_cur[(ar_order+1):total_params] else numeric(0)

    # Arima order for the final model
    arima_order <- c(ar_order, 0, ma_order)

    # Fit final model to get innovations
    final_model <- tryCatch({
      stats::arima(x, order = arima_order,
                   include.mean = FALSE,
                   fixed = c(ar_part, ma_part))
    }, error = function(e) {
      if(verbose) {
        cat("Error in final model fitting:", conditionMessage(e), "\n")
      }
      NULL
    })

    if(!is.null(final_model)) {
      final_innovations <- residuals(final_model)
    } else {
      final_innovations <- innovations
      warning("Failed to compute final innovations. Using last iteration values.")
    }
  }

  # Plot convergence history if verbose
  if(verbose && iterations > 1) {
    if(requireNamespace("graphics", quietly = TRUE)) {
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
    innovations = final_innovations
  )
}

#' Create design matrix for AR models
#'
#' @param x centered time series
#' @param p AR order
#' @return Design matrix with lagged values
#' @keywords internal
create_ar_matrix <- function(x, p) {
  n <- length(x)
  X <- matrix(NA, nrow = n - p, ncol = p)
  for(i in 1:p) {
    X[, i] <- x[(p - i + 1):(n - i)]
  }
  return(X)
}

#' Get Yule-Walker estimates for AR model
#'
#' @param x centered time series
#' @param p AR order
#' @return Vector of AR coefficients
#' @keywords internal
get_yw_estimates <- function(x, p) {
  acf_values <- stats::acf(x, lag.max = p, plot = FALSE)$acf[-1]
  if(p == 1) {
    return(acf_values)
  } else {
    # Solve Yule-Walker equations
    R <- stats::toeplitz(c(1, acf_values[1:(p-1)]))
    return(solve(R, acf_values))
  }
}

#' Update innovations for MA model
#'
#' @param x centered time series
#' @param ma_coef vector of MA coefficients
#' @return vector of innovations
#' @keywords internal
update_ma_innovations <- function(x, ma_coef) {
  n <- length(x)
  q <- length(ma_coef)

  # Initialize innovations as zeros
  innovations <- numeric(n)

  # Iteratively compute innovations
  for(t in 1:n) {
    # Compute expected value based on previous innovations
    expected <- 0
    for(j in 1:q) {
      if(t - j > 0) {
        expected <- expected + ma_coef[j] * innovations[t - j]
      }
    }

    # Compute current innovation
    innovations[t] <- x[t] - expected
  }

  return(innovations)
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

#' Calculate moments and cumulants of error distribution
#'
#' @param errors numeric vector of errors
#' @return list with moments, cumulants and theoretical variance reduction coefficient
#' @export
compute_moments <- function(errors) {
  m2 <- mean(errors^2)
  m3 <- mean(errors^3)
  m4 <- mean(errors^4)

  c3 <- m3 / m2^(3/2)  # Коефіцієнт асиметрії
  c4 <- m4 / m2^2 - 3  # Коефіцієнт ексцесу

  # Теоретичний коефіцієнт зменшення дисперсії
  g <- 1 - c3^2 / (2 + c4)

  return(list(m2 = m2, m3 = m3, m4 = m4,
              c3 = c3, c4 = c4,
              g = g))
}
