#!/usr/bin/env Rscript
# ==============================================================================
# EstemPMM-стиль PMM2 для MA/SMA параметрів
# ==============================================================================
#
# Базується на аналізі EstemPMM коду:
# https://github.com/SZabolotnii/EstemPMM/blob/main/R/pmm2_ts_main.R
#
# Ключовий підхід:
#   1. CSS fit → початкові оцінки + залишки
#   2. Використати залишки як ФІКСОВАНІ регресори
#   3. Застосувати лінійний PMM2 з формулою градієнту
#
# Формули PMM2 (з EstemPMM):
#   Z1 = m3*S^2 + (m4 - m2^2 - 2*m3*Y)*S + (m3*Y^2 - (m4-m2^2)*Y - m2*m3)
#   Z = t(X) %*% Z1
#   JZ11 = 2*m3*S + (m4 - m2^2 - 2*m3*Y)
#   J = t(X) %*% (X * JZ11)
#   b_new = b - solve(J, Z)
#
# ==============================================================================


#' Обчислення MA інновацій (forward recursion)
#'
#' @param x Time series (centered)
#' @param theta MA coefficients
#' @param q MA order
#' @return Vector of innovations
ma_compute_innovations <- function(x, theta, q) {
  n <- length(x)
  innovations <- numeric(n)
  history <- rep(0, q)

  for (t in seq_len(n)) {
    ma_component <- if (q > 0) sum(theta * history) else 0
    innovations[t] <- x[t] - ma_component

    if (q > 0) {
      history <- c(innovations[t], history)[seq_len(q)]
    }
  }

  innovations
}


#' Обчислення SMA інновацій (forward recursion)
#'
#' @param x Time series (centered)
#' @param Theta SMA coefficients
#' @param Q SMA order
#' @param s Seasonal period
#' @return Vector of innovations
sma_compute_innovations <- function(x, Theta, Q, s) {
  n <- length(x)
  innovations <- numeric(n)

  for (t in seq_len(n)) {
    sma_component <- 0
    if (Q > 0) {
      for (J in 1:Q) {
        lag <- J * s
        if (t - lag >= 1) {
          sma_component <- sma_component + Theta[J] * innovations[t - lag]
        }
      }
    }
    innovations[t] <- x[t] - sma_component
  }

  innovations
}


#' Побудова дизайн матриці для MA
#'
#' @param intercept Intercept (mean)
#' @param residuals Initial residuals from CSS
#' @param x Time series
#' @param q MA order
#' @return List with X (design matrix) and y (response)
ma_build_design <- function(intercept, residuals, x, q) {
  idx <- seq.int(q + 1L, length(x))
  X <- matrix(1, nrow = length(idx), ncol = q + 1L)

  for (j in seq_len(q)) {
    X[, j + 1L] <- residuals[idx - j]
  }

  y <- x[idx] - intercept

  list(X = X, y = y)
}


#' Побудова дизайн матриці для SMA
#'
#' @param intercept Intercept (mean)
#' @param residuals Initial residuals from CSS
#' @param x Time series
#' @param Q SMA order
#' @param s Seasonal period
#' @return List with X (design matrix) and y (response)
sma_build_design <- function(intercept, residuals, x, Q, s) {
  max_lag <- Q * s
  idx <- seq.int(max_lag + 1L, length(x))
  X <- matrix(1, nrow = length(idx), ncol = Q + 1L)

  for (J in seq_len(Q)) {
    lag_seasonal <- J * s
    X[, J + 1L] <- residuals[idx - lag_seasonal]
  }

  y <- x[idx] - intercept

  list(X = X, y = y)
}


#' PMM2 solver (EstemPMM формула)
#'
#' @param b_init Initial coefficients [intercept, theta1, ..., thetaq]
#' @param X Design matrix
#' @param Y Response vector
#' @param m2 Second central moment
#' @param m3 Third central moment
#' @param m4 Fourth central moment
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @param verbose Print diagnostics
#'
#' @return List with coefficients, convergence, iterations
ma_solve_pmm2 <- function(b_init, X, Y, m2, m3, m4,
                          max_iter = 50, tol = 1e-6,
                          verbose = FALSE) {
  b <- as.numeric(b_init)
  iterations <- 0L
  converged <- FALSE

  for (iter in seq_len(max_iter)) {
    iterations <- iter

    # S = X %*% b (predicted values)
    S <- as.vector(X %*% b)

    # Градієнт Z1 (формула з EstemPMM)
    Z1 <- m3 * S^2 +
          (m4 - m2^2 - 2 * m3 * Y) * S +
          (m3 * Y^2 - (m4 - m2^2) * Y - m2 * m3)

    # Z = t(X) %*% Z1
    Z <- as.numeric(t(X) %*% Z1)

    # Jacobian JZ11
    JZ11 <- 2 * m3 * S + (m4 - m2^2 - 2 * m3 * Y)

    # J = t(X) %*% (X * JZ11)
    J <- t(X) %*% (X * JZ11)

    # Newton step
    step <- tryCatch(solve(J, Z), error = function(e) NULL)

    if (is.null(step)) {
      if (verbose) cat("  System is singular at iteration", iter, "\n")
      break
    }

    # Update
    b_new <- b - step

    # Check convergence
    if (sqrt(sum((b_new - b)^2)) < tol) {
      b <- b_new
      converged <- TRUE
      if (verbose) cat("  Converged at iteration", iter, "\n")
      break
    }

    b <- b_new

    if (verbose && (iter %% 10 == 0 || iter <= 3)) {
      cat("  Iter", iter, ": ||Δb|| =", sqrt(sum(step^2)), "\n")
    }
  }

  list(
    coefficients = b[-1],  # Видалити intercept
    intercept = b[1],
    convergence = converged,
    iterations = iterations
  )
}


#' EstemPMM-стиль оцінювач для MA
#'
#' @param x Time series
#' @param q MA order
#' @param include.mean Include intercept
#' @param max_iter Maximum PMM2 iterations
#' @param verbose Print diagnostics
#'
#' @return List with ma_coef, mean, innovations, convergence, method
#'
#' @export
estpmm_style_ma <- function(x,
                            q = 1,
                            include.mean = TRUE,
                            max_iter = 50,
                            verbose = FALSE) {

  n <- length(x)

  if (verbose) {
    cat("=== EstemPMM-Style MA Estimator ===\n")
    cat("Model: MA(", q, ")\n")
    cat("Sample size:", n, "\n\n")
  }

  # ===========================================================================
  # STEP 1: CSS fit для початкових оцінок
  # ===========================================================================
  if (verbose) cat("Step 1: CSS initialization...\n")

  css_fit <- tryCatch({
    stats::arima(x, order = c(0, 0, q),
                 method = "CSS",
                 include.mean = include.mean)
  }, error = function(e) {
    if (verbose) cat("  CSS failed, trying CSS-ML\n")
    stats::arima(x, order = c(0, 0, q),
                 method = "CSS-ML",
                 include.mean = include.mean)
  })

  # Витягти початкові параметри
  coefs_css <- coef(css_fit)

  theta_init <- numeric(q)
  for (j in 1:q) {
    coef_name <- paste0("ma", j)
    if (coef_name %in% names(coefs_css)) {
      theta_init[j] <- coefs_css[coef_name]
    }
  }

  intercept_init <- if (include.mean && "intercept" %in% names(coefs_css)) {
    coefs_css["intercept"]
  } else {
    0
  }

  # Обчислити залишки з CSS
  residuals_css <- as.numeric(residuals(css_fit))

  if (verbose) {
    cat("  CSS estimates: θ =", paste(round(theta_init, 4), collapse = ", "), "\n")
    cat("  Intercept:", round(intercept_init, 4), "\n\n")
  }

  # ===========================================================================
  # STEP 2: Побудувати дизайн матрицю з ФІКСОВАНИХ залишків
  # ===========================================================================
  if (verbose) cat("Step 2: Build design matrix from fixed CSS residuals...\n")

  design <- ma_build_design(intercept_init, residuals_css, x, q)
  X <- design$X
  y <- design$y

  if (verbose) {
    cat("  Design matrix: ", nrow(X), "×", ncol(X), "\n")
    cat("  Effective sample size:", nrow(X), "\n\n")
  }

  # ===========================================================================
  # STEP 3: Обчислити моменти з CSS залишків
  # ===========================================================================
  if (verbose) cat("Step 3: Compute moments from CSS residuals...\n")

  # Використати ефективні залишки (після видалення початкових q спостережень)
  eff_residuals <- residuals_css[(q+1):n]

  m2 <- mean(eff_residuals^2)
  m3 <- mean(eff_residuals^3)
  m4 <- mean(eff_residuals^4)

  gamma3 <- m3 / (m2^(3/2))
  gamma4 <- m4 / (m2^2) - 3

  if (verbose) {
    cat("  m2:", round(m2, 4), "\n")
    cat("  m3:", round(m3, 4), "\n")
    cat("  m4:", round(m4, 4), "\n")
    cat("  γ₃:", round(gamma3, 4), "\n")
    cat("  γ₄:", round(gamma4, 4), "\n\n")
  }

  # ===========================================================================
  # STEP 4: PMM2 оптимізація
  # ===========================================================================
  if (verbose) cat("Step 4: PMM2 optimization...\n")

  b_init <- c(intercept_init, theta_init)

  pmm2_result <- ma_solve_pmm2(b_init, X, y, m2, m3, m4,
                               max_iter = max_iter,
                               tol = 1e-6,
                               verbose = verbose)

  # Фінальні параметри
  theta_final <- pmm2_result$coefficients
  intercept_final <- pmm2_result$intercept

  if (verbose) {
    cat("\nFinal estimates:\n")
    cat("  θ:", paste(round(theta_final, 4), collapse = ", "), "\n")
    cat("  Intercept:", round(intercept_final, 4), "\n")
  }

  # ===========================================================================
  # STEP 5: Обчислити фінальні інновації
  # ===========================================================================
  x_centered <- x - intercept_final
  innovations <- ma_compute_innovations(x_centered, theta_final, q)

  list(
    ma_coef = theta_final,
    mean = intercept_final,
    innovations = innovations,
    convergence = pmm2_result$convergence,
    iterations = pmm2_result$iterations,
    css_estimates = theta_init,
    method = "EstemPMM-style PMM2"
  )
}


#' EstemPMM-стиль оцінювач для SMA
#'
#' @param x Time series
#' @param Q SMA order
#' @param s Seasonal period
#' @param include.mean Include intercept
#' @param max_iter Maximum PMM2 iterations
#' @param verbose Print diagnostics
#'
#' @return List with sma_coef, mean, innovations, convergence, method
#'
#' @export
estpmm_style_sma <- function(x,
                             Q = 1,
                             s = 4,
                             include.mean = TRUE,
                             max_iter = 50,
                             verbose = FALSE) {

  n <- length(x)

  if (verbose) {
    cat("=== EstemPMM-Style SMA Estimator ===\n")
    cat("Model: SMA(", Q, ")_", s, "\n")
    cat("Sample size:", n, "\n\n")
  }

  # ===========================================================================
  # STEP 1: CSS fit
  # ===========================================================================
  if (verbose) cat("Step 1: CSS initialization...\n")

  css_fit <- tryCatch({
    stats::arima(x, order = c(0, 0, 0),
                 seasonal = list(order = c(0, 0, Q), period = s),
                 method = "CSS",
                 include.mean = include.mean)
  }, error = function(e) {
    if (verbose) cat("  CSS failed, trying CSS-ML\n")
    stats::arima(x, order = c(0, 0, 0),
                 seasonal = list(order = c(0, 0, Q), period = s),
                 method = "CSS-ML",
                 include.mean = include.mean)
  })

  coefs_css <- coef(css_fit)

  Theta_init <- numeric(Q)
  for (J in 1:Q) {
    coef_name <- paste0("sma", J)
    if (coef_name %in% names(coefs_css)) {
      Theta_init[J] <- coefs_css[coef_name]
    }
  }

  intercept_init <- if (include.mean && "intercept" %in% names(coefs_css)) {
    coefs_css["intercept"]
  } else {
    0
  }

  residuals_css <- as.numeric(residuals(css_fit))

  if (verbose) {
    cat("  CSS estimates: Θ =", paste(round(Theta_init, 4), collapse = ", "), "\n")
    cat("  Intercept:", round(intercept_init, 4), "\n\n")
  }

  # ===========================================================================
  # STEP 2: Дизайн матриця
  # ===========================================================================
  if (verbose) cat("Step 2: Build design matrix...\n")

  design <- sma_build_design(intercept_init, residuals_css, x, Q, s)
  X <- design$X
  y <- design$y

  if (verbose) {
    cat("  Design matrix:", nrow(X), "×", ncol(X), "\n\n")
  }

  # ===========================================================================
  # STEP 3: Моменти
  # ===========================================================================
  if (verbose) cat("Step 3: Compute moments...\n")

  max_lag <- Q * s
  eff_residuals <- residuals_css[(max_lag+1):n]

  m2 <- mean(eff_residuals^2)
  m3 <- mean(eff_residuals^3)
  m4 <- mean(eff_residuals^4)

  gamma3 <- m3 / (m2^(3/2))
  gamma4 <- m4 / (m2^2) - 3

  if (verbose) {
    cat("  γ₃:", round(gamma3, 4), ", γ₄:", round(gamma4, 4), "\n\n")
  }

  # ===========================================================================
  # STEP 4: PMM2
  # ===========================================================================
  if (verbose) cat("Step 4: PMM2 optimization...\n")

  b_init <- c(intercept_init, Theta_init)

  pmm2_result <- ma_solve_pmm2(b_init, X, y, m2, m3, m4,
                               max_iter = max_iter,
                               tol = 1e-6,
                               verbose = verbose)

  Theta_final <- pmm2_result$coefficients
  intercept_final <- pmm2_result$intercept

  if (verbose) {
    cat("\nFinal estimates:\n")
    cat("  Θ:", paste(round(Theta_final, 4), collapse = ", "), "\n")
  }

  # ===========================================================================
  # STEP 5: Інновації
  # ===========================================================================
  x_centered <- x - intercept_final
  innovations <- sma_compute_innovations(x_centered, Theta_final, Q, s)

  list(
    sma_coef = Theta_final,
    mean = intercept_final,
    innovations = innovations,
    convergence = pmm2_result$convergence,
    iterations = pmm2_result$iterations,
    css_estimates = Theta_init,
    method = "EstemPMM-style PMM2"
  )
}
