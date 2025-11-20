#!/usr/bin/env Rscript
# ==============================================================================
# Гібридний SARIMA оцінювач з EstemPMM-Style PMM2 для MA/SMA
# ==============================================================================
#
# Стратегія:
#   - MA/SMA: MLE (за замовчуванням) або EstemPMM-style PMM2 (опція)
#   - AR/SAR: PMM2 (за замовчуванням) або MLE (опція)
#
# Параметри:
#   ma_method = c("mle", "pmm2")      # MLE або EstemPMM-style PMM2
#   ar_method = c("pmm2", "mle")      # PMM2 або MLE
#
# ==============================================================================

source("R/experimental/06_estpmm_style_ma.R")

#' Гібридний SARIMA оцінювач з EstemPMM-style PMM2
#'
#' @param x Time series
#' @param order Vector c(p, d, q)
#' @param seasonal List with order c(P, D, Q) and period
#' @param include.mean Include intercept
#' @param ma_method Method for MA/SMA: "mle" or "pmm2"
#' @param ar_method Method for AR/SAR: "pmm2" or "mle"
#' @param pmm_method PMM variant for AR: "pmm1" or "pmm2"
#' @param pmm2_max_iter Max iterations for PMM2
#' @param verbose Print diagnostics
#'
#' @return List with coefficients, innovations, method used
#'
#' @export
hybrid_sarima_estimator_v2 <- function(x,
                                       order = c(0, 0, 0),
                                       seasonal = list(order = c(0, 0, 0), period = 1),
                                       include.mean = TRUE,
                                       ma_method = c("mle", "pmm2"),
                                       ar_method = c("pmm2", "mle"),
                                       pmm_method = "pmm2",
                                       pmm2_max_iter = 30,
                                       verbose = FALSE) {

  ma_method <- match.arg(ma_method)
  ar_method <- match.arg(ar_method)

  p <- order[1]
  d <- order[2]
  q <- order[3]

  P <- seasonal$order[1]
  D <- seasonal$order[2]
  Q <- seasonal$order[3]
  s <- seasonal$period

  n <- length(x)

  if (verbose) {
    cat("=== Гібридний SARIMA оцінювач v2 ===\n")
    cat("Модель: SARIMA(", p, ",", d, ",", q, ")(", P, ",", D, ",", Q, ")_", s, "\n", sep = "")
    cat("MA method:", ma_method, "\n")
    cat("AR method:", ar_method, "\n")
    cat("Sample size:", n, "\n\n")
  }

  # Ініціалізація результатів
  ma_coef_final <- if (q > 0) numeric(q) else NULL
  sma_coef_final <- if (Q > 0) numeric(Q) else NULL
  ar_coef_final <- if (p > 0) numeric(p) else NULL
  sar_coef_final <- if (P > 0) numeric(P) else NULL
  mean_final <- 0
  innovations <- NULL
  ma_method_used <- ma_method
  ar_method_used <- ar_method

  # ===========================================================================
  # КРОК 1: Оцінка MA/SMA параметрів
  # ===========================================================================

  if (q == 0 && Q == 0) {
    # Немає MA компонентів - пропустити
    if (verbose) cat("Немає MA/SMA параметрів для оцінки\n\n")

    # Потрібно отримати залишки для AR оцінки
    if (p > 0 || P > 0) {
      # Якщо є AR, отримаємо залишки з MLE
      mle_fit <- stats::arima(x, order = order,
                              seasonal = seasonal,
                              include.mean = include.mean,
                              method = "ML")
      innovations <- as.numeric(residuals(mle_fit))

      if (include.mean) {
        mean_final <- coef(mle_fit)["intercept"]
      }

      # Витягти AR коефіцієнти з MLE (будуть перезаписані PMM2 якщо потрібно)
      if (p > 0) {
        for (i in 1:p) {
          coef_name <- paste0("ar", i)
          if (coef_name %in% names(coef(mle_fit))) {
            ar_coef_final[i] <- coef(mle_fit)[coef_name]
          }
        }
      }

      if (P > 0) {
        for (I in 1:P) {
          coef_name <- paste0("sar", I)
          if (coef_name %in% names(coef(mle_fit))) {
            sar_coef_final[I] <- coef(mle_fit)[coef_name]
          }
        }
      }
    }

  } else if (ma_method == "mle") {
    # MLE для MA/SMA
    if (verbose) cat("КРОК 1: MLE оцінка MA/SMA параметрів\n")

    mle_fit <- stats::arima(x, order = order,
                            seasonal = seasonal,
                            include.mean = include.mean,
                            method = "ML")

    innovations <- as.numeric(residuals(mle_fit))

    if (include.mean) {
      mean_final <- coef(mle_fit)["intercept"]
    }

    # Витягти MA коефіцієнти
    if (q > 0) {
      for (j in 1:q) {
        coef_name <- paste0("ma", j)
        if (coef_name %in% names(coef(mle_fit))) {
          ma_coef_final[j] <- coef(mle_fit)[coef_name]
        }
      }
      if (verbose) cat("  MA (MLE):", paste(round(ma_coef_final, 4), collapse = ", "), "\n")
    }

    # Витягти SMA коефіцієнти
    if (Q > 0) {
      for (J in 1:Q) {
        coef_name <- paste0("sma", J)
        if (coef_name %in% names(coef(mle_fit))) {
          sma_coef_final[J] <- coef(mle_fit)[coef_name]
        }
      }
      if (verbose) cat("  SMA (MLE):", paste(round(sma_coef_final, 4), collapse = ", "), "\n")
    }

    # Витягти AR коефіцієнти з MLE (будуть перезаписані PMM2 якщо потрібно)
    if (p > 0) {
      for (i in 1:p) {
        coef_name <- paste0("ar", i)
        if (coef_name %in% names(coef(mle_fit))) {
          ar_coef_final[i] <- coef(mle_fit)[coef_name]
        }
      }
    }

    if (P > 0) {
      for (I in 1:P) {
        coef_name <- paste0("sar", I)
        if (coef_name %in% names(coef(mle_fit))) {
          sar_coef_final[I] <- coef(mle_fit)[coef_name]
        }
      }
    }

  } else {
    # EstemPMM-style PMM2 для MA/SMA
    if (verbose) cat("КРОК 1: EstemPMM-style PMM2 оцінка MA/SMA параметрів\n")

    # Випадок 1: Тільки MA (без SAR/SMA)
    if (q > 0 && Q == 0 && P == 0) {
      pmm2_fit <- estpmm_style_ma(x, q = q,
                                   include.mean = include.mean,
                                   max_iter = pmm2_max_iter,
                                   verbose = verbose)

      if (pmm2_fit$convergence) {
        ma_coef_final <- pmm2_fit$ma_coef
        mean_final <- pmm2_fit$mean
        innovations <- pmm2_fit$innovations
        ma_method_used <- "pmm2"

        if (verbose) {
          cat("  MA (PMM2):", paste(round(ma_coef_final, 4), collapse = ", "), "\n")
        }
      } else {
        if (verbose) cat("  PMM2 не збігся, використовую MLE\n")
        ma_method_used <- "mle (pmm2 failed)"

        # Fallback до MLE
        mle_fit <- stats::arima(x, order = order,
                                seasonal = seasonal,
                                include.mean = include.mean,
                                method = "ML")

        innovations <- as.numeric(residuals(mle_fit))

        if (q > 0) {
          for (j in 1:q) {
            coef_name <- paste0("ma", j)
            if (coef_name %in% names(coef(mle_fit))) {
              ma_coef_final[j] <- coef(mle_fit)[coef_name]
            }
          }
        }

        if (include.mean) {
          mean_final <- coef(mle_fit)["intercept"]
        }
      }

    } else if (q == 0 && Q > 0 && p == 0) {
      # Випадок 2: Тільки SMA (без AR/MA)
      pmm2_fit <- estpmm_style_sma(x, Q = Q, s = s,
                                    include.mean = include.mean,
                                    max_iter = pmm2_max_iter,
                                    verbose = verbose)

      if (pmm2_fit$convergence) {
        sma_coef_final <- pmm2_fit$sma_coef
        mean_final <- pmm2_fit$mean
        innovations <- pmm2_fit$innovations
        ma_method_used <- "pmm2"

        if (verbose) {
          cat("  SMA (PMM2):", paste(round(sma_coef_final, 4), collapse = ", "), "\n")
        }
      } else {
        if (verbose) cat("  PMM2 не збігся, використовую MLE\n")
        ma_method_used <- "mle (pmm2 failed)"

        # Fallback до MLE
        mle_fit <- stats::arima(x, order = order,
                                seasonal = seasonal,
                                include.mean = include.mean,
                                method = "ML")

        innovations <- as.numeric(residuals(mle_fit))

        if (Q > 0) {
          for (J in 1:Q) {
            coef_name <- paste0("sma", J)
            if (coef_name %in% names(coef(mle_fit))) {
              sma_coef_final[J] <- coef(mle_fit)[coef_name]
            }
          }
        }

        if (include.mean) {
          mean_final <- coef(mle_fit)["intercept"]
        }
      }

    } else {
      # Випадок 3: Змішані моделі або моделі з AR
      # Поки що fallback до MLE
      if (verbose) cat("  Змішані моделі - використовую MLE (PMM2 для змішаних поки не реалізовано)\n")
      ma_method_used <- "mle (mixed model)"

      mle_fit <- stats::arima(x, order = order,
                              seasonal = seasonal,
                              include.mean = include.mean,
                              method = "ML")

      innovations <- as.numeric(residuals(mle_fit))

      if (q > 0) {
        for (j in 1:q) {
          coef_name <- paste0("ma", j)
          if (coef_name %in% names(coef(mle_fit))) {
            ma_coef_final[j] <- coef(mle_fit)[coef_name]
          }
        }
      }

      if (Q > 0) {
        for (J in 1:Q) {
          coef_name <- paste0("sma", J)
          if (coef_name %in% names(coef(mle_fit))) {
            sma_coef_final[J] <- coef(mle_fit)[coef_name]
          }
        }
      }

      if (p > 0) {
        for (i in 1:p) {
          coef_name <- paste0("ar", i)
          if (coef_name %in% names(coef(mle_fit))) {
            ar_coef_final[i] <- coef(mle_fit)[coef_name]
          }
        }
      }

      if (P > 0) {
        for (I in 1:P) {
          coef_name <- paste0("sar", I)
          if (coef_name %in% names(coef(mle_fit))) {
            sar_coef_final[I] <- coef(mle_fit)[coef_name]
          }
        }
      }

      if (include.mean) {
        mean_final <- coef(mle_fit)["intercept"]
      }
    }
  }

  # ===========================================================================
  # КРОК 2: Оцінка AR/SAR параметрів (якщо є і якщо ar_method = "pmm2")
  # ===========================================================================

  if ((p > 0 || P > 0) && ar_method == "pmm2") {
    if (verbose) cat("\nКРОК 2: PMM2 оцінка AR/SAR параметрів\n")

    # Завантажити EstemPMM для AR оцінки
    if (!requireNamespace("EstemPMM", quietly = TRUE)) {
      if (verbose) cat("  EstemPMM не встановлено, використовую MLE для AR\n")
      ar_method_used <- "mle (no EstemPMM)"
    } else {
      # Побудувати дизайн матрицю для AR з фіксованих інновацій
      x_centered <- x - mean_final

      # Використати EstemPMM::sarima_pmm для AR оцінки
      tryCatch({
        # Створити тимчасову модель тільки з AR компонентами
        ar_only_order <- c(p, d, 0)
        ar_only_seasonal <- list(order = c(P, D, 0), period = s)

        ar_pmm2_fit <- EstemPMM::sarima_pmm(
          x,
          order = ar_only_order,
          seasonal = ar_only_seasonal,
          include.mean = include.mean,
          method = pmm_method,
          max_iter = pmm2_max_iter
        )

        # Витягти AR коефіцієнти
        if (!is.null(ar_pmm2_fit) && isTRUE(slot(ar_pmm2_fit, "convergence"))) {
          coeffs <- slot(ar_pmm2_fit, "coefficients")

          if (p > 0 && length(coeffs) >= p) {
            ar_coef_final <- coeffs[1:p]
            if (verbose) cat("  AR (PMM2):", paste(round(ar_coef_final, 4), collapse = ", "), "\n")
          }

          if (P > 0 && length(coeffs) >= (p + P)) {
            sar_coef_final <- coeffs[(p+1):(p+P)]
            if (verbose) cat("  SAR (PMM2):", paste(round(sar_coef_final, 4), collapse = ", "), "\n")
          }

          ar_method_used <- "pmm2"
        } else {
          if (verbose) cat("  PMM2 для AR не збігся, залишаю MLE оцінки\n")
          ar_method_used <- "mle (pmm2 failed)"
        }
      }, error = function(e) {
        if (verbose) cat("  Помилка PMM2 для AR:", e$message, "\n")
        ar_method_used <- "mle (pmm2 error)"
      })
    }
  }

  # ===========================================================================
  # Повернути результат
  # ===========================================================================

  method_desc <- paste0(
    "MA: ", ma_method_used,
    ", AR: ", ar_method_used
  )

  if (verbose) {
    cat("\n=== Фінальні оцінки ===\n")
    if (!is.null(ar_coef_final)) cat("AR:", paste(round(ar_coef_final, 4), collapse = ", "), "\n")
    if (!is.null(sar_coef_final)) cat("SAR:", paste(round(sar_coef_final, 4), collapse = ", "), "\n")
    if (!is.null(ma_coef_final)) cat("MA:", paste(round(ma_coef_final, 4), collapse = ", "), "\n")
    if (!is.null(sma_coef_final)) cat("SMA:", paste(round(sma_coef_final, 4), collapse = ", "), "\n")
    if (include.mean) cat("Mean:", round(mean_final, 4), "\n")
    cat("Method:", method_desc, "\n")
  }

  list(
    ar_coef = ar_coef_final,
    sar_coef = sar_coef_final,
    ma_coef = ma_coef_final,
    sma_coef = sma_coef_final,
    mean = mean_final,
    innovations = innovations,
    ma_method = ma_method_used,
    ar_method = ar_method_used,
    method = method_desc,
    convergence = TRUE
  )
}
