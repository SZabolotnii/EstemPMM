#!/usr/bin/env Rscript
# =============================================================================
# mc_pmm3_comprehensive.R
# Масштабне Монте-Карло тестування PMM3 для всіх реалізованих варіацій:
#   1. lm_pmm3()   — лінійна регресія
#   2. ar_pmm3()   — AR моделі
#   3. ma_pmm3()   — MA моделі
#   4. arma_pmm3() — ARMA моделі
#   5. arima_pmm3() — ARIMA моделі
#
# Розподіли похибок:
#   - Uniform(-sqrt(3), sqrt(3))  — gamma4 ≈ -1.2, сильно платикуртичний
#   - Beta(2,2) rescaled          — gamma4 ≈ -0.86, помірно платикуртичний
#   - TN (Two-piece Normal)       — gamma4 залежить від lambda
#
# Порівняння: PMM3 vs OLS/MLE, з теоретичним g3
# Розміри вибірки: n = 100, 200, 500
# Реплікацій: M = 500
# =============================================================================

cat("============================================================\n")
cat("  PMM3 Comprehensive Monte Carlo Simulation\n")
cat("============================================================\n\n")

# Завантаження пакету
suppressMessages(library(EstemPMM))

# ── Конфігурація ─────────────────────────────────────────────────────────────

MC_CFG <- list(
  M          = 500L,
  n_grid     = c(100L, 200L, 500L),
  seed       = 20260319L,
  output_dir = "mc_results_pmm3"
)

dir.create(MC_CFG$output_dir, showWarnings = FALSE, recursive = TRUE)

# ── Генератори розподілів ────────────────────────────────────────────────────

gen_uniform <- function(n) runif(n, -sqrt(3), sqrt(3))

gen_beta22 <- function(n) {
  # Beta(2,2) rescaled to mean=0, var=1
  x <- rbeta(n, 2, 2) - 0.5
  x / sd(x)
}

gen_tn <- function(n, lambda = 1.5) {
  U <- rnorm(n, mean = lambda, sd = 1)
  S <- sample(c(-1, 1), n, replace = TRUE)
  eps <- S * abs(U)
  eps - mean(eps)  # center
}

# Теоретичні характеристики розподілів
dist_info <- list(
  uniform = list(
    name = "Uniform",
    gen  = gen_uniform,
    # Theoretical: m2=1, m4=9/5, m6=27/7 => gamma4=-1.2, gamma6=6.857
    gamma4_theory = -1.2,
    g3_theory     = 1 - (-1.2)^2 / (6 + 9*(-1.2) + 6.857)  # ≈ 0.379
  ),
  beta22 = list(
    name = "Beta(2,2)",
    gen  = gen_beta22,
    gamma4_theory = -0.857,
    g3_theory     = NA  # обчислимо емпірично
  ),
  tn_1.5 = list(
    name = "TN(lambda=1.5)",
    gen  = function(n) gen_tn(n, lambda = 1.5),
    gamma4_theory = NA,
    g3_theory     = NA
  )
)

# ── Допоміжні функції ────────────────────────────────────────────────────────

#' Обчислити g3 теоретичний з вибірки розподілу (MC оцінка)
estimate_dist_g3 <- function(gen_fn, n_samples = 100000) {
  x <- gen_fn(n_samples)
  mom <- compute_moments_pmm3(x)
  list(gamma4 = mom$gamma4, gamma6 = mom$gamma6, g3 = mom$g3)
}

#' Безпечний виклик функції
safe_call <- function(expr) {
  tryCatch(suppressWarnings(expr), error = function(e) NULL)
}

# Оновити теоретичні g3 для розподілів де NA
set.seed(MC_CFG$seed)
for (nm in names(dist_info)) {
  if (is.na(dist_info[[nm]]$g3_theory)) {
    info <- estimate_dist_g3(dist_info[[nm]]$gen)
    dist_info[[nm]]$gamma4_theory <- info$gamma4
    dist_info[[nm]]$g3_theory     <- info$g3
  }
}

cat("Distribution characteristics (estimated from 100k samples):\n")
for (nm in names(dist_info)) {
  cat(sprintf("  %-15s  gamma4 = %+.4f  g3 = %.4f\n",
              dist_info[[nm]]$name,
              dist_info[[nm]]$gamma4_theory,
              dist_info[[nm]]$g3_theory))
}
cat("\n")

# =============================================================================
# ЧАСТИНА 1: ЛІНІЙНА РЕГРЕСІЯ (lm_pmm3 vs OLS)
# =============================================================================

run_mc_linear <- function(dist, n, M, seed) {
  set.seed(seed)
  beta_true <- c(1.0, 2.0)
  x_fix <- rnorm(n)

  ols_b1  <- numeric(M)
  pmm3_b1 <- numeric(M)
  pmm3_conv <- logical(M)

  for (i in seq_len(M)) {
    eps <- dist$gen(n)
    y   <- beta_true[1] + beta_true[2] * x_fix + eps
    dat <- data.frame(y = y, x = x_fix)

    # OLS
    ols_b1[i] <- coef(lm(y ~ x, data = dat))[2]

    # PMM3
    fit3 <- safe_call(lm_pmm3(y ~ x, data = dat))
    if (!is.null(fit3)) {
      pmm3_b1[i]   <- coef(fit3)[2]
      pmm3_conv[i] <- fit3@convergence
    } else {
      pmm3_b1[i]   <- NA_real_
      pmm3_conv[i] <- FALSE
    }
  }

  valid <- !is.na(pmm3_b1)
  list(
    var_ols  = var(ols_b1),
    var_pmm3 = var(pmm3_b1[valid]),
    bias_ols  = mean(ols_b1) - beta_true[2],
    bias_pmm3 = mean(pmm3_b1[valid]) - beta_true[2],
    mse_ols  = mean((ols_b1 - beta_true[2])^2),
    mse_pmm3 = mean((pmm3_b1[valid] - beta_true[2])^2),
    conv_rate = mean(pmm3_conv),
    n_valid   = sum(valid)
  )
}

cat("================================================================\n")
cat("  PART 1: LINEAR REGRESSION (lm_pmm3 vs OLS)\n")
cat("================================================================\n\n")

lm_results <- list()
idx <- 0

for (dist_name in names(dist_info)) {
  for (nn in MC_CFG$n_grid) {
    idx <- idx + 1
    cat(sprintf("[LM %d] %-15s n=%3d M=%d ... ",
                idx, dist_info[[dist_name]]$name, nn, MC_CFG$M))

    res <- run_mc_linear(dist_info[[dist_name]], nn, MC_CFG$M,
                          MC_CFG$seed + idx)

    g3_emp <- res$var_pmm3 / res$var_ols
    cat(sprintf("g3_emp=%.4f  g3_th=%.4f  conv=%.1f%%\n",
                g3_emp, dist_info[[dist_name]]$g3_theory, res$conv_rate * 100))

    lm_results[[idx]] <- data.frame(
      model      = "lm_pmm3",
      dist       = dist_info[[dist_name]]$name,
      n          = nn,
      M          = MC_CFG$M,
      gamma4_th  = dist_info[[dist_name]]$gamma4_theory,
      g3_theory  = dist_info[[dist_name]]$g3_theory,
      var_ols    = res$var_ols,
      var_pmm3   = res$var_pmm3,
      g3_emp     = g3_emp,
      are_pmm3   = res$var_ols / res$var_pmm3,
      bias_ols   = res$bias_ols,
      bias_pmm3  = res$bias_pmm3,
      mse_ols    = res$mse_ols,
      mse_pmm3   = res$mse_pmm3,
      mse_ratio  = res$mse_pmm3 / res$mse_ols,
      conv_rate  = res$conv_rate,
      n_valid    = res$n_valid,
      stringsAsFactors = FALSE
    )
  }
}

lm_df <- do.call(rbind, lm_results)

# =============================================================================
# ЧАСТИНА 2: AR МОДЕЛІ (ar_pmm3 vs OLS/MLE)
# =============================================================================

run_mc_ar <- function(dist, n, M, seed, ar_true = 0.7, order = 1) {
  set.seed(seed)

  ols_coef  <- matrix(NA_real_, M, order)
  pmm3_coef <- matrix(NA_real_, M, order)
  pmm3_conv <- logical(M)

  for (i in seq_len(M)) {
    eps <- dist$gen(n + 100)  # burn-in buffer
    ar_list <- if (order == 1) list(ar = ar_true) else list(ar = ar_true)
    x <- tryCatch(
      arima.sim(n = n, model = ar_list, innov = eps[1:n],
                start.innov = eps[(n+1):(n+100)]),
      error = function(e) NULL
    )
    if (is.null(x)) next

    # OLS (via ar_pmm2 internals: YW or lm.fit on lagged matrix)
    x_c <- as.numeric(x) - mean(x)
    X_ar <- EstemPMM:::create_ar_matrix(x_c, order)
    y_ar <- x_c[(order + 1):length(x_c)]
    ols_fit <- lm.fit(X_ar, y_ar)
    ols_coef[i, ] <- ols_fit$coefficients

    # PMM3
    fit3 <- safe_call(ar_pmm3(as.numeric(x), order = order))
    if (!is.null(fit3)) {
      pmm3_coef[i, ] <- fit3@coefficients
      pmm3_conv[i]   <- fit3@convergence
    }
  }

  valid_ols  <- complete.cases(ols_coef)
  valid_pmm3 <- complete.cases(pmm3_coef)

  list(
    var_ols   = apply(ols_coef[valid_ols, , drop = FALSE], 2, var),
    var_pmm3  = apply(pmm3_coef[valid_pmm3, , drop = FALSE], 2, var),
    bias_ols  = colMeans(ols_coef[valid_ols, , drop = FALSE]) - ar_true,
    bias_pmm3 = colMeans(pmm3_coef[valid_pmm3, , drop = FALSE]) - ar_true,
    mse_ols   = colMeans((ols_coef[valid_ols, , drop = FALSE] -
                          matrix(ar_true, sum(valid_ols), order, byrow = TRUE))^2),
    mse_pmm3  = colMeans((pmm3_coef[valid_pmm3, , drop = FALSE] -
                          matrix(ar_true, sum(valid_pmm3), order, byrow = TRUE))^2),
    conv_rate = mean(pmm3_conv),
    n_valid_ols = sum(valid_ols),
    n_valid_pmm3 = sum(valid_pmm3)
  )
}

cat("\n================================================================\n")
cat("  PART 2: AR(1) MODELS (ar_pmm3 vs OLS)\n")
cat("================================================================\n\n")

ar_results <- list()
idx <- 0

for (dist_name in names(dist_info)) {
  for (nn in MC_CFG$n_grid) {
    idx <- idx + 1
    cat(sprintf("[AR %d] %-15s n=%3d M=%d ... ",
                idx, dist_info[[dist_name]]$name, nn, MC_CFG$M))

    res <- run_mc_ar(dist_info[[dist_name]], nn, MC_CFG$M,
                      MC_CFG$seed + 100 + idx, ar_true = 0.7, order = 1)

    g3_emp <- res$var_pmm3[1] / res$var_ols[1]
    cat(sprintf("g3_emp=%.4f  MSE_ratio=%.4f  conv=%.1f%%\n",
                g3_emp, res$mse_pmm3[1] / res$mse_ols[1], res$conv_rate * 100))

    ar_results[[idx]] <- data.frame(
      model      = "ar_pmm3",
      dist       = dist_info[[dist_name]]$name,
      n          = nn,
      M          = MC_CFG$M,
      gamma4_th  = dist_info[[dist_name]]$gamma4_theory,
      g3_theory  = dist_info[[dist_name]]$g3_theory,
      var_ols    = res$var_ols[1],
      var_pmm3   = res$var_pmm3[1],
      g3_emp     = g3_emp,
      are_pmm3   = res$var_ols[1] / res$var_pmm3[1],
      bias_ols   = res$bias_ols[1],
      bias_pmm3  = res$bias_pmm3[1],
      mse_ols    = res$mse_ols[1],
      mse_pmm3   = res$mse_pmm3[1],
      mse_ratio  = res$mse_pmm3[1] / res$mse_ols[1],
      conv_rate  = res$conv_rate,
      n_valid    = res$n_valid_pmm3,
      stringsAsFactors = FALSE
    )
  }
}

ar_df <- do.call(rbind, ar_results)

# =============================================================================
# ЧАСТИНА 3: MA МОДЕЛІ (ma_pmm3 vs MLE)
# =============================================================================

run_mc_ma <- function(dist, n, M, seed, ma_true = 0.6, order = 1) {
  set.seed(seed)

  mle_coef  <- numeric(M)
  pmm3_coef <- numeric(M)
  pmm3_conv <- logical(M)
  mle_valid <- logical(M)
  pmm3_valid <- logical(M)

  for (i in seq_len(M)) {
    eps <- dist$gen(n + 50)
    x <- tryCatch(
      arima.sim(n = n, model = list(ma = ma_true), innov = eps[1:n],
                start.innov = eps[(n+1):(n+50)]),
      error = function(e) NULL
    )
    if (is.null(x)) next

    # MLE via stats::arima
    mle_fit <- tryCatch(
      stats::arima(as.numeric(x), order = c(0, 0, order), method = "CSS-ML"),
      error = function(e) NULL
    )
    if (!is.null(mle_fit)) {
      mle_coef[i]  <- mle_fit$coef[paste0("ma", 1)]
      mle_valid[i] <- TRUE
    }

    # PMM3
    fit3 <- safe_call(ma_pmm3(as.numeric(x), order = order))
    if (!is.null(fit3)) {
      pmm3_coef[i]  <- fit3@coefficients[1]
      pmm3_conv[i]  <- fit3@convergence
      pmm3_valid[i] <- TRUE
    }
  }

  list(
    var_mle   = var(mle_coef[mle_valid]),
    var_pmm3  = var(pmm3_coef[pmm3_valid]),
    bias_mle  = mean(mle_coef[mle_valid]) - ma_true,
    bias_pmm3 = mean(pmm3_coef[pmm3_valid]) - ma_true,
    mse_mle   = mean((mle_coef[mle_valid] - ma_true)^2),
    mse_pmm3  = mean((pmm3_coef[pmm3_valid] - ma_true)^2),
    conv_rate = mean(pmm3_conv),
    n_valid_mle  = sum(mle_valid),
    n_valid_pmm3 = sum(pmm3_valid)
  )
}

cat("\n================================================================\n")
cat("  PART 3: MA(1) MODELS (ma_pmm3 vs MLE)\n")
cat("================================================================\n\n")

ma_results <- list()
idx <- 0

for (dist_name in names(dist_info)) {
  for (nn in MC_CFG$n_grid) {
    idx <- idx + 1
    cat(sprintf("[MA %d] %-15s n=%3d M=%d ... ",
                idx, dist_info[[dist_name]]$name, nn, MC_CFG$M))

    res <- run_mc_ma(dist_info[[dist_name]], nn, MC_CFG$M,
                      MC_CFG$seed + 200 + idx, ma_true = 0.6, order = 1)

    g3_emp <- res$var_pmm3 / res$var_mle
    cat(sprintf("Var_ratio=%.4f  MSE_ratio=%.4f  conv=%.1f%%\n",
                g3_emp, res$mse_pmm3 / res$mse_mle, res$conv_rate * 100))

    ma_results[[idx]] <- data.frame(
      model      = "ma_pmm3",
      dist       = dist_info[[dist_name]]$name,
      n          = nn,
      M          = MC_CFG$M,
      gamma4_th  = dist_info[[dist_name]]$gamma4_theory,
      g3_theory  = dist_info[[dist_name]]$g3_theory,
      var_ols    = res$var_mle,  # baseline is MLE for MA
      var_pmm3   = res$var_pmm3,
      g3_emp     = g3_emp,
      are_pmm3   = res$var_mle / res$var_pmm3,
      bias_ols   = res$bias_mle,
      bias_pmm3  = res$bias_pmm3,
      mse_ols    = res$mse_mle,
      mse_pmm3   = res$mse_pmm3,
      mse_ratio  = res$mse_pmm3 / res$mse_mle,
      conv_rate  = res$conv_rate,
      n_valid    = res$n_valid_pmm3,
      stringsAsFactors = FALSE
    )
  }
}

ma_df <- do.call(rbind, ma_results)

# =============================================================================
# ЧАСТИНА 4: ARMA МОДЕЛІ (arma_pmm3 vs MLE)
# =============================================================================

run_mc_arma <- function(dist, n, M, seed, ar_true = 0.5, ma_true = 0.3) {
  set.seed(seed)

  mle_ar  <- numeric(M)
  mle_ma  <- numeric(M)
  pmm3_ar <- numeric(M)
  pmm3_ma <- numeric(M)
  pmm3_conv <- logical(M)
  mle_valid <- logical(M)
  pmm3_valid <- logical(M)

  for (i in seq_len(M)) {
    eps <- dist$gen(n + 100)
    x <- tryCatch(
      arima.sim(n = n, model = list(ar = ar_true, ma = ma_true),
                innov = eps[1:n], start.innov = eps[(n+1):(n+100)]),
      error = function(e) NULL
    )
    if (is.null(x)) next

    # MLE
    mle_fit <- tryCatch(
      stats::arima(as.numeric(x), order = c(1, 0, 1), method = "CSS-ML"),
      error = function(e) NULL
    )
    if (!is.null(mle_fit)) {
      mle_ar[i]    <- mle_fit$coef["ar1"]
      mle_ma[i]    <- mle_fit$coef["ma1"]
      mle_valid[i] <- TRUE
    }

    # PMM3
    fit3 <- safe_call(arma_pmm3(as.numeric(x), order = c(1, 1)))
    if (!is.null(fit3)) {
      pmm3_ar[i]    <- fit3@coefficients[1]
      pmm3_ma[i]    <- fit3@coefficients[2]
      pmm3_conv[i]  <- fit3@convergence
      pmm3_valid[i] <- TRUE
    }
  }

  list(
    var_mle_ar   = var(mle_ar[mle_valid]),
    var_mle_ma   = var(mle_ma[mle_valid]),
    var_pmm3_ar  = var(pmm3_ar[pmm3_valid]),
    var_pmm3_ma  = var(pmm3_ma[pmm3_valid]),
    mse_mle_ar   = mean((mle_ar[mle_valid] - ar_true)^2),
    mse_mle_ma   = mean((mle_ma[mle_valid] - ma_true)^2),
    mse_pmm3_ar  = mean((pmm3_ar[pmm3_valid] - ar_true)^2),
    mse_pmm3_ma  = mean((pmm3_ma[pmm3_valid] - ma_true)^2),
    conv_rate    = mean(pmm3_conv),
    n_valid_pmm3 = sum(pmm3_valid)
  )
}

cat("\n================================================================\n")
cat("  PART 4: ARMA(1,1) MODELS (arma_pmm3 vs MLE)\n")
cat("================================================================\n\n")

arma_results <- list()
idx <- 0

for (dist_name in names(dist_info)) {
  for (nn in MC_CFG$n_grid) {
    idx <- idx + 1
    cat(sprintf("[ARMA %d] %-15s n=%3d M=%d ... ",
                idx, dist_info[[dist_name]]$name, nn, MC_CFG$M))

    res <- run_mc_arma(dist_info[[dist_name]], nn, MC_CFG$M,
                        MC_CFG$seed + 300 + idx)

    mse_ratio_ar <- res$mse_pmm3_ar / res$mse_mle_ar
    mse_ratio_ma <- res$mse_pmm3_ma / res$mse_mle_ma
    cat(sprintf("MSE_ratio(AR)=%.4f  MSE_ratio(MA)=%.4f  conv=%.1f%%\n",
                mse_ratio_ar, mse_ratio_ma, res$conv_rate * 100))

    arma_results[[idx]] <- data.frame(
      model      = "arma_pmm3",
      dist       = dist_info[[dist_name]]$name,
      n          = nn,
      M          = MC_CFG$M,
      gamma4_th  = dist_info[[dist_name]]$gamma4_theory,
      g3_theory  = dist_info[[dist_name]]$g3_theory,
      var_ols    = (res$var_mle_ar + res$var_mle_ma) / 2,
      var_pmm3   = (res$var_pmm3_ar + res$var_pmm3_ma) / 2,
      g3_emp     = (res$var_pmm3_ar + res$var_pmm3_ma) / (res$var_mle_ar + res$var_mle_ma),
      are_pmm3   = (res$var_mle_ar + res$var_mle_ma) / (res$var_pmm3_ar + res$var_pmm3_ma),
      bias_ols   = NA_real_,
      bias_pmm3  = NA_real_,
      mse_ols    = (res$mse_mle_ar + res$mse_mle_ma) / 2,
      mse_pmm3   = (res$mse_pmm3_ar + res$mse_pmm3_ma) / 2,
      mse_ratio  = (res$mse_pmm3_ar + res$mse_pmm3_ma) / (res$mse_mle_ar + res$mse_mle_ma),
      conv_rate  = res$conv_rate,
      n_valid    = res$n_valid_pmm3,
      stringsAsFactors = FALSE
    )
  }
}

arma_df <- do.call(rbind, arma_results)

# =============================================================================
# ЧАСТИНА 5: ARIMA МОДЕЛІ (arima_pmm3 vs MLE)
# =============================================================================

run_mc_arima <- function(dist, n, M, seed, ar_true = 0.6) {
  set.seed(seed)

  mle_coef  <- numeric(M)
  pmm3_coef <- numeric(M)
  pmm3_conv <- logical(M)
  mle_valid <- logical(M)
  pmm3_valid <- logical(M)

  for (i in seq_len(M)) {
    eps <- dist$gen(n + 100)
    # Generate ARIMA(1,1,0): integrated AR(1)
    ar_series <- tryCatch(
      arima.sim(n = n, model = list(ar = ar_true),
                innov = eps[1:n], start.innov = eps[(n+1):(n+100)]),
      error = function(e) NULL
    )
    if (is.null(ar_series)) next
    x <- cumsum(as.numeric(ar_series))

    # MLE via stats::arima
    mle_fit <- tryCatch(
      stats::arima(x, order = c(1, 1, 0), method = "CSS-ML"),
      error = function(e) NULL
    )
    if (!is.null(mle_fit)) {
      mle_coef[i]  <- mle_fit$coef["ar1"]
      mle_valid[i] <- TRUE
    }

    # PMM3
    fit3 <- safe_call(arima_pmm3(x, order = c(1, 1, 0)))
    if (!is.null(fit3)) {
      pmm3_coef[i]  <- fit3@coefficients[1]
      pmm3_conv[i]  <- fit3@convergence
      pmm3_valid[i] <- TRUE
    }
  }

  list(
    var_mle   = var(mle_coef[mle_valid]),
    var_pmm3  = var(pmm3_coef[pmm3_valid]),
    bias_mle  = mean(mle_coef[mle_valid]) - ar_true,
    bias_pmm3 = mean(pmm3_coef[pmm3_valid]) - ar_true,
    mse_mle   = mean((mle_coef[mle_valid] - ar_true)^2),
    mse_pmm3  = mean((pmm3_coef[pmm3_valid] - ar_true)^2),
    conv_rate = mean(pmm3_conv),
    n_valid_pmm3 = sum(pmm3_valid)
  )
}

cat("\n================================================================\n")
cat("  PART 5: ARIMA(1,1,0) MODELS (arima_pmm3 vs MLE)\n")
cat("================================================================\n\n")

arima_results <- list()
idx <- 0

for (dist_name in names(dist_info)) {
  for (nn in MC_CFG$n_grid) {
    idx <- idx + 1
    cat(sprintf("[ARIMA %d] %-15s n=%3d M=%d ... ",
                idx, dist_info[[dist_name]]$name, nn, MC_CFG$M))

    res <- run_mc_arima(dist_info[[dist_name]], nn, MC_CFG$M,
                         MC_CFG$seed + 400 + idx)

    g3_emp <- res$var_pmm3 / res$var_mle
    cat(sprintf("Var_ratio=%.4f  MSE_ratio=%.4f  conv=%.1f%%\n",
                g3_emp, res$mse_pmm3 / res$mse_mle, res$conv_rate * 100))

    arima_results[[idx]] <- data.frame(
      model      = "arima_pmm3",
      dist       = dist_info[[dist_name]]$name,
      n          = nn,
      M          = MC_CFG$M,
      gamma4_th  = dist_info[[dist_name]]$gamma4_theory,
      g3_theory  = dist_info[[dist_name]]$g3_theory,
      var_ols    = res$var_mle,
      var_pmm3   = res$var_pmm3,
      g3_emp     = g3_emp,
      are_pmm3   = res$var_mle / res$var_pmm3,
      bias_ols   = res$bias_mle,
      bias_pmm3  = res$bias_pmm3,
      mse_ols    = res$mse_mle,
      mse_pmm3   = res$mse_pmm3,
      mse_ratio  = res$mse_pmm3 / res$mse_mle,
      conv_rate  = res$conv_rate,
      n_valid    = res$n_valid_pmm3,
      stringsAsFactors = FALSE
    )
  }
}

arima_df <- do.call(rbind, arima_results)

# =============================================================================
# ЗВЕДЕНА ТАБЛИЦЯ
# =============================================================================

all_df <- rbind(lm_df, ar_df, ma_df, arma_df, arima_df)

cat("\n\n")
cat("================================================================\n")
cat("  COMBINED RESULTS TABLE\n")
cat("================================================================\n\n")

# Форматування для виводу
print_table <- all_df[, c("model", "dist", "n", "gamma4_th", "g3_theory",
                           "g3_emp", "mse_ratio", "conv_rate")]
print_table$gamma4_th <- round(print_table$gamma4_th, 3)
print_table$g3_theory <- round(print_table$g3_theory, 4)
print_table$g3_emp    <- round(print_table$g3_emp, 4)
print_table$mse_ratio <- round(print_table$mse_ratio, 4)
print_table$conv_rate <- round(print_table$conv_rate * 100, 1)

print(print_table, row.names = FALSE)

# =============================================================================
# ЗБЕРЕЖЕННЯ РЕЗУЛЬТАТІВ
# =============================================================================

write.csv(all_df,
          file.path(MC_CFG$output_dir, "mc_pmm3_all_results.csv"),
          row.names = FALSE)

# Зведена таблиця по моделях (усереднення по n)
summary_by_model <- aggregate(
  cbind(g3_emp, mse_ratio, conv_rate) ~ model + dist,
  data = all_df,
  FUN = mean
)
summary_by_model$g3_emp    <- round(summary_by_model$g3_emp, 4)
summary_by_model$mse_ratio <- round(summary_by_model$mse_ratio, 4)
summary_by_model$conv_rate <- round(summary_by_model$conv_rate * 100, 1)

write.csv(summary_by_model,
          file.path(MC_CFG$output_dir, "mc_pmm3_summary_by_model.csv"),
          row.names = FALSE)

cat("\n\n================================================================\n")
cat("  SUMMARY BY MODEL (averaged over sample sizes)\n")
cat("================================================================\n\n")
print(summary_by_model, row.names = FALSE)

# =============================================================================
# ВЕРИФІКАЦІЯ
# =============================================================================

cat("\n\n================================================================\n")
cat("  VERIFICATION CHECKS\n")
cat("================================================================\n\n")

n_checks <- 0
n_pass   <- 0

for (i in seq_len(nrow(all_df))) {
  row <- all_df[i, ]

  # 1. Convergence >= 95%
  n_checks <- n_checks + 1
  ok <- row$conv_rate >= 0.95
  if (ok) n_pass <- n_pass + 1
  if (!ok) {
    cat(sprintf("  [FAIL] conv=%.1f%% < 95%% (%s, %s, n=%d)\n",
                row$conv_rate * 100, row$model, row$dist, row$n))
  }

  # 2. g3 empirical within 0.20 of theory (for n >= 200, linear models)
  if (row$model == "lm_pmm3" && row$n >= 200 &&
      !is.na(row$g3_emp) && !is.na(row$g3_theory)) {
    n_checks <- n_checks + 1
    delta <- abs(row$g3_emp - row$g3_theory)
    ok <- delta < 0.20
    if (ok) n_pass <- n_pass + 1
    if (!ok) {
      cat(sprintf("  [FAIL] |g3_emp - g3_th| = %.3f >= 0.20 (%s, n=%d)\n",
                  delta, row$dist, row$n))
    }
  }

  # 3. MSE ratio < 1.10 (PMM3 should not be much worse than baseline)
  if (!is.na(row$mse_ratio)) {
    n_checks <- n_checks + 1
    ok <- row$mse_ratio < 1.10
    if (ok) n_pass <- n_pass + 1
    if (!ok) {
      cat(sprintf("  [FAIL] mse_ratio=%.4f >= 1.10 (%s, %s, n=%d)\n",
                  row$mse_ratio, row$model, row$dist, row$n))
    }
  }
}

cat(sprintf("\n  %d / %d checks passed\n", n_pass, n_checks))

# Зберегти фінальну інформацію
cat(sprintf("\nResults saved to %s/\n", MC_CFG$output_dir))
cat("  mc_pmm3_all_results.csv\n")
cat("  mc_pmm3_summary_by_model.csv\n")

cat("\n============================================================\n")
cat("  PMM3 Monte Carlo Simulation Complete\n")
cat("============================================================\n")
