#!/usr/bin/env Rscript
# ==============================================================================
# Нелінійний PMM2 для MA/SMA параметрів
# ==============================================================================
#
# Реалізація справжнього нелінійного PMM2 для MA моделей через оптимізацію
# цільової функції, як описано в:
# "Application of the Polynomial Maximization Method for Estimating
#  Nonlinear Regression Parameters with Non-gaussian Asymmetric Errors"
#
# Ключова ідея:
#   - MA інновації залежать від параметрів: ε_t = ε_t(θ)
#   - Оптимізуємо цільову функцію PMM2: γ₃²/(2 + γ₄)
#   - де γ₃, γ₄ обчислюються з інновацій ε(θ)
#
# ==============================================================================

#' Обчислення MA інновацій (forward recursion)
#'
#' @param x Time series (centered)
#' @param ma_coef MA coefficients (length q)
#' @param sma_coef SMA coefficients (length Q)
#' @param s Seasonal period
#' @param init_method Initialization: "zero" або "mle"
#'
#' @return Vector of innovations
compute_ma_sma_innovations <- function(x, ma_coef = NULL, sma_coef = NULL,
                                       s = 1, init_method = "zero") {
  n <- length(x)
  q <- if (is.null(ma_coef)) 0 else length(ma_coef)
  Q <- if (is.null(sma_coef)) 0 else length(sma_coef)

  innovations <- numeric(n)

  # Ініціалізація
  if (init_method == "zero") {
    # Перші max(q, Q*s) інновацій = 0
    init_len <- max(q, Q * s)
    if (init_len > 0) {
      innovations[1:min(init_len, n)] <- 0
    }
  } else {
    # MLE ініціалізація (можна покращити)
    innovations[1:min(max(q, Q*s), n)] <- 0
  }

  # Forward recursion: ε_t = x_t - Σθ_j*ε_{t-j} - ΣΘ_J*ε_{t-s*J}
  start_idx <- max(q, Q * s) + 1
  if (start_idx <= n) {
    for (t in start_idx:n) {
      eps_t <- x[t]

      # MA частина
      if (q > 0) {
        for (j in 1:q) {
          if (t - j >= 1) {
            eps_t <- eps_t - ma_coef[j] * innovations[t - j]
          }
        }
      }

      # SMA частина
      if (Q > 0) {
        for (J in 1:Q) {
          lag <- J * s
          if (t - lag >= 1) {
            eps_t <- eps_t - sma_coef[J] * innovations[t - lag]
          }
        }
      }

      innovations[t] <- eps_t
    }
  }

  innovations
}


#' Обчислення моментів з інновацій
#'
#' @param innovations Vector of innovations
#' @param remove_init Remove first k observations (default: TRUE)
#' @param k Number of initial observations to remove
#'
#' @return List with m2, m3, m4, gamma3, gamma4
compute_moments <- function(innovations, remove_init = TRUE, k = NULL) {
  if (remove_init) {
    if (is.null(k)) {
      k <- min(10, floor(length(innovations) * 0.05))
    }
    if (k > 0 && k < length(innovations)) {
      innovations <- innovations[(k+1):length(innovations)]
    }
  }

  # Видалити NA
  innovations <- innovations[!is.na(innovations)]

  if (length(innovations) < 10) {
    return(list(m2 = NA, m3 = NA, m4 = NA, gamma3 = NA, gamma4 = NA, valid = FALSE))
  }

  # Центрування (якщо потрібно)
  # innovations <- innovations - mean(innovations)

  m2 <- mean(innovations^2)
  m3 <- mean(innovations^3)
  m4 <- mean(innovations^4)

  if (m2 <= 0) {
    return(list(m2 = m2, m3 = m3, m4 = m4, gamma3 = NA, gamma4 = NA, valid = FALSE))
  }

  gamma3 <- m3 / (m2^(3/2))
  gamma4 <- m4 / (m2^2) - 3

  list(
    m2 = m2,
    m3 = m3,
    m4 = m4,
    gamma3 = gamma3,
    gamma4 = gamma4,
    valid = TRUE
  )
}


#' PMM2 цільова функція для MA/SMA параметрів
#'
#' @param theta Vector of MA/SMA parameters [MA(1), ..., MA(q), SMA(1), ..., SMA(Q)]
#' @param x Time series (centered)
#' @param q MA order
#' @param Q SMA order
#' @param s Seasonal period
#' @param verbose Print diagnostics
#'
#' @return PMM2 objective value: γ₃²/(2 + γ₄)
objective_pmm2_ma <- function(theta, x, q = 0, Q = 0, s = 1, verbose = FALSE) {
  # Розділити параметри
  ma_coef <- if (q > 0) theta[1:q] else NULL
  sma_coef <- if (Q > 0) theta[(q+1):(q+Q)] else NULL

  # Перевірка стаціонарності/invertibility (приблизна)
  # Для MA(1): |θ| < 1, для SMA(1): |Θ| < 1
  if (q == 1 && !is.null(ma_coef)) {
    if (abs(ma_coef[1]) > 0.99) return(-1e10)
  }
  if (Q == 1 && !is.null(sma_coef)) {
    if (abs(sma_coef[1]) > 0.99) return(-1e10)
  }

  # Обчислити інновації
  innovations <- compute_ma_sma_innovations(x, ma_coef, sma_coef, s,
                                            init_method = "zero")

  # Обчислити моменти
  moments <- compute_moments(innovations, remove_init = TRUE,
                             k = max(q, Q * s))

  if (!moments$valid) {
    return(-1e10)
  }

  gamma3 <- moments$gamma3
  gamma4 <- moments$gamma4

  # Цільова функція PMM2: γ₃²/(2 + γ₄)
  # Потрібно уникнути ділення на 0 або від'ємний знаменник
  denominator <- 2 + gamma4

  if (is.na(denominator) || denominator <= 0.1) {
    return(-1e10)
  }

  objective <- (gamma3^2) / denominator

  if (verbose) {
    cat("  θ =", round(theta, 4), "| γ₃ =", round(gamma3, 4),
        ", γ₄ =", round(gamma4, 4), ", obj =", round(objective, 6), "\n")
  }

  # Повертаємо objective (для максимізації)
  objective
}


#' Нелінійний PMM2 оцінювач для MA/SMA
#'
#' @param x Time series
#' @param q MA order
#' @param Q SMA order (default: 0)
#' @param s Seasonal period (default: 1)
#' @param include.mean Include intercept (default: TRUE)
#' @param init_method Initialization: "mle" (default) або "zero"
#' @param optim_method Optimization method: "BFGS", "Nelder-Mead", "L-BFGS-B"
#' @param max_iter Maximum iterations (default: 100)
#' @param verbose Print diagnostics
#'
#' @return List with:
#'   \item{ma_coef}{MA coefficients}
#'   \item{sma_coef}{SMA coefficients}
#'   \item{mean}{Intercept}
#'   \item{innovations}{Final innovations}
#'   \item{convergence}{Convergence status}
#'   \item{objective}{Final objective value}
#'   \item{iterations}{Number of iterations}
#'   \item{init_estimates}{Initial estimates}
#'   \item{method}{"Nonlinear PMM2"}
#'
#' @export
nonlinear_pmm2_ma <- function(x,
                              q = 1,
                              Q = 0,
                              s = 1,
                              include.mean = TRUE,
                              init_method = c("mle", "zero"),
                              optim_method = c("BFGS", "Nelder-Mead", "L-BFGS-B"),
                              max_iter = 100,
                              verbose = FALSE) {

  init_method <- match.arg(init_method)
  optim_method <- match.arg(optim_method)

  n <- length(x)

  if (verbose) {
    cat("=== Nonlinear PMM2 for MA/SMA ===\n")
    cat("Model: MA(", q, "), SMA(", Q, ")_", s, "\n")
    cat("Initialization:", init_method, "\n")
    cat("Optimization:", optim_method, "\n\n")
  }

  # ===========================================================================
  # STEP 1: Центрування та оцінка середнього
  # ===========================================================================
  x_mean <- if (include.mean) mean(x) else 0
  x_centered <- x - x_mean

  # ===========================================================================
  # STEP 2: Початкові оцінки
  # ===========================================================================
  if (init_method == "mle") {
    # MLE ініціалізація
    if (verbose) cat("Step 1: MLE initialization...\n")

    mle_fit <- tryCatch({
      stats::arima(x, order = c(0, 0, q),
                   seasonal = list(order = c(0, 0, Q), period = s),
                   include.mean = include.mean,
                   method = "CSS-ML")
    }, error = function(e) {
      if (verbose) cat("  MLE failed, using zeros\n")
      return(NULL)
    })

    if (!is.null(mle_fit)) {
      coefs_mle <- coef(mle_fit)

      theta_init <- numeric(q + Q)

      # MA coefficients
      if (q > 0) {
        for (j in 1:q) {
          coef_name <- paste0("ma", j)
          if (coef_name %in% names(coefs_mle)) {
            theta_init[j] <- coefs_mle[coef_name]
          }
        }
      }

      # SMA coefficients
      if (Q > 0) {
        for (J in 1:Q) {
          coef_name <- paste0("sma", J)
          if (coef_name %in% names(coefs_mle)) {
            theta_init[q + J] <- coefs_mle[coef_name]
          }
        }
      }

      if (verbose) {
        cat("  MLE initial estimates:", paste(round(theta_init, 4), collapse = ", "), "\n\n")
      }
    } else {
      theta_init <- rep(0.1, q + Q)
    }
  } else {
    # Zero initialization
    theta_init <- rep(0.1, q + Q)
    if (verbose) {
      cat("  Initial estimates:", paste(round(theta_init, 4), collapse = ", "), "\n\n")
    }
  }

  # ===========================================================================
  # STEP 3: Нелінійна оптимізація PMM2
  # ===========================================================================
  if (verbose) cat("Step 2: Nonlinear PMM2 optimization...\n")

  # Bounds для L-BFGS-B
  lower <- rep(-0.99, q + Q)
  upper <- rep(0.99, q + Q)

  optim_result <- tryCatch({
    if (optim_method == "L-BFGS-B") {
      optim(
        par = theta_init,
        fn = objective_pmm2_ma,
        method = "L-BFGS-B",
        lower = lower,
        upper = upper,
        control = list(fnscale = -1, maxit = max_iter),
        x = x_centered,
        q = q,
        Q = Q,
        s = s,
        verbose = verbose
      )
    } else {
      optim(
        par = theta_init,
        fn = objective_pmm2_ma,
        method = optim_method,
        control = list(fnscale = -1, maxit = max_iter),
        x = x_centered,
        q = q,
        Q = Q,
        s = s,
        verbose = verbose
      )
    }
  }, error = function(e) {
    if (verbose) cat("  Optimization failed:", e$message, "\n")
    return(NULL)
  })

  if (is.null(optim_result)) {
    return(list(
      ma_coef = if (q > 0) rep(NA, q) else numeric(0),
      sma_coef = if (Q > 0) rep(NA, Q) else numeric(0),
      mean = x_mean,
      innovations = rep(NA, n),
      convergence = FALSE,
      objective = NA,
      iterations = 0,
      init_estimates = theta_init,
      method = "Nonlinear PMM2 (failed)"
    ))
  }

  # ===========================================================================
  # STEP 4: Витягти результати
  # ===========================================================================
  theta_final <- optim_result$par

  ma_coef <- if (q > 0) theta_final[1:q] else numeric(0)
  sma_coef <- if (Q > 0) theta_final[(q+1):(q+Q)] else numeric(0)

  # Фінальні інновації
  innovations <- compute_ma_sma_innovations(x_centered, ma_coef, sma_coef, s)

  # Convergence
  converged <- (optim_result$convergence == 0)

  if (verbose) {
    cat("\nResults:\n")
    cat("  Convergence:", converged, "(code:", optim_result$convergence, ")\n")
    cat("  Objective value:", round(optim_result$value, 6), "\n")
    if (q > 0) cat("  MA coefficients:", paste(round(ma_coef, 4), collapse = ", "), "\n")
    if (Q > 0) cat("  SMA coefficients:", paste(round(sma_coef, 4), collapse = ", "), "\n")

    # Фінальні моменти
    final_moments <- compute_moments(innovations, remove_init = TRUE, k = max(q, Q*s))
    if (final_moments$valid) {
      cat("  Final γ₃:", round(final_moments$gamma3, 4), "\n")
      cat("  Final γ₄:", round(final_moments$gamma4, 4), "\n")
    }
  }

  list(
    ma_coef = ma_coef,
    sma_coef = sma_coef,
    mean = x_mean,
    innovations = innovations,
    convergence = converged,
    objective = optim_result$value,
    iterations = if (!is.null(optim_result$counts)) optim_result$counts[1] else NA,
    init_estimates = theta_init,
    method = "Nonlinear PMM2"
  )
}
