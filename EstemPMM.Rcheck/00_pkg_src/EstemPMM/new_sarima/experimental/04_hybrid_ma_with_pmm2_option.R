# ==============================================================================
# Розширений гібридний MA Support: MLE + опціональний PMM2 для SMA
# ==============================================================================
#
# Базується на 03_hybrid_ma_support.R, але додає можливість:
#   1. Використати MLE залишки як початкові інновації
#   2. Застосувати нелінійний PMM2 для SMA параметрів (експериментально)
#   3. Порівняти MLE vs PMM2 для SMA
#
# ==============================================================================

library(EstemPMM)  # Для sarima_pmm2()

#' Розширений гібридний SARIMA оцінювач з опцією PMM2 для MA/SMA
#'
#' @param x Numeric vector (time series data)
#' @param order Vector c(p, d, q)
#' @param seasonal List with order c(P, D, Q) and period
#' @param include.mean Logical, include intercept
#' @param ma_method Method for MA/SMA: "mle" (default) or "pmm2" (experimental)
#' @param ar_method Method for AR/SAR: "mle" or "pmm2" (default)
#' @param pmm2_max_iter Max iterations for PMM2 optimization
#' @param verbose Logical, print diagnostics
#'
#' @return List with:
#'   \item{ma_coef}{MA coefficients}
#'   \item{sma_coef}{SMA coefficients}
#'   \item{ar_coef}{AR coefficients}
#'   \item{sar_coef}{SAR coefficients}
#'   \item{mean}{Intercept}
#'   \item{innovations}{Innovations}
#'   \item{method}{Method description}
#'   \item{ma_method_used}{Method actually used for MA}
#'   \item{ar_method_used}{Method actually used for AR}
#'   \item{mle_fit}{Initial MLE fit}
#'   \item{pmm2_fit}{PMM2 fit (if used)}
#'
#' @details
#' Алгоритм:
#'   1. MLE fit → початкові оцінки + innovations
#'   2. IF ma_method == "pmm2" AND є сезонні MA:
#'        Застосувати sarima_pmm2() з MLE як початкові значення
#'   3. IF ar_method == "pmm2":
#'        Re-estimate AR/SAR через PMM2
#'   4. Combine results
#'
#' @export
hybrid_sarima_estimator_extended <- function(x,
                                              order = c(0, 0, 0),
                                              seasonal = list(order = c(0, 0, 0), period = 1),
                                              include.mean = TRUE,
                                              ma_method = c("mle", "pmm2"),
                                              ar_method = c("pmm2", "mle"),
                                              pmm2_max_iter = 30,
                                              verbose = FALSE) {
  ma_method <- match.arg(ma_method)
  ar_method <- match.arg(ar_method)

  # Extract components
  p <- order[1]; d <- order[2]; q <- order[3]
  P <- seasonal$order[1]; D <- seasonal$order[2]; Q <- seasonal$order[3]
  s <- seasonal$period

  n <- length(x)

  if (verbose) {
    cat("=== Extended Hybrid SARIMA Estimator ===\n")
    cat("Model: SARIMA(", p, ",", d, ",", q, ")(", P, ",", D, ",", Q, ")_", s, "\n")
    cat("MA method:", ma_method, "\n")
    cat("AR method:", ar_method, "\n\n")
  }

  # ===========================================================================
  # STEP 1: Initial MLE fit
  # ===========================================================================
  if (verbose) cat("Step 1: Initial MLE fit...\n")

  mle_fit <- tryCatch({
    stats::arima(x, order = order, seasonal = seasonal,
                 method = "CSS-ML", include.mean = include.mean)
  }, error = function(e) {
    warning("MLE initialization failed: ", e$message)
    return(NULL)
  })

  if (is.null(mle_fit)) {
    stop("MLE initialization failed. Cannot proceed.")
  }

  # Extract MLE coefficients
  coefs_mle <- coef(mle_fit)

  # Extract MA/SMA from MLE
  ma_coef_mle <- numeric(q)
  if (q > 0) {
    for (j in 1:q) {
      coef_name <- paste0("ma", j)
      if (coef_name %in% names(coefs_mle)) {
        ma_coef_mle[j] <- coefs_mle[coef_name]
      }
    }
  }

  sma_coef_mle <- numeric(Q)
  if (Q > 0) {
    for (j in 1:Q) {
      coef_name <- paste0("sma", j)
      if (coef_name %in% names(coefs_mle)) {
        sma_coef_mle[j] <- coefs_mle[coef_name]
      }
    }
  }

  # Extract AR/SAR from MLE
  ar_coef_mle <- numeric(p)
  if (p > 0) {
    for (j in 1:p) {
      coef_name <- paste0("ar", j)
      if (coef_name %in% names(coefs_mle)) {
        ar_coef_mle[j] <- coefs_mle[coef_name]
      }
    }
  }

  sar_coef_mle <- numeric(P)
  if (P > 0) {
    for (j in 1:P) {
      coef_name <- paste0("sar", j)
      if (coef_name %in% names(coefs_mle)) {
        sar_coef_mle[j] <- coefs_mle[coef_name]
      }
    }
  }

  # Extract mean
  mean_est <- if (include.mean && "intercept" %in% names(coefs_mle)) {
    coefs_mle["intercept"]
  } else {
    0
  }

  # Extract innovations
  innovations <- as.numeric(residuals(mle_fit))

  if (verbose) {
    cat("  MLE estimates:\n")
    if (p > 0) cat("    AR:", paste(round(ar_coef_mle, 4), collapse = ", "), "\n")
    if (q > 0) cat("    MA:", paste(round(ma_coef_mle, 4), collapse = ", "), "\n")
    if (P > 0) cat("    SAR:", paste(round(sar_coef_mle, 4), collapse = ", "), "\n")
    if (Q > 0) cat("    SMA:", paste(round(sma_coef_mle, 4), collapse = ", "), "\n")
    cat("  Innovations: μ =", round(mean(innovations), 4),
        ", σ =", round(sd(innovations), 4), "\n\n")
  }

  # ===========================================================================
  # STEP 2: Apply PMM2 for SMA if requested
  # ===========================================================================
  ma_coef_final <- ma_coef_mle
  sma_coef_final <- sma_coef_mle
  pmm2_fit <- NULL
  ma_method_used <- "mle"

  if (ma_method == "pmm2" && (Q > 0 || q > 0)) {
    if (verbose) cat("Step 2: Applying PMM2 for MA/SMA parameters...\n")

    # Використовуємо sarima_pmm2() з EstemPMM
    # Формат: order = c(p, P, q, Q), seasonal = list(order = c(d, D), period = s)

    # Перевіряємо чи є сезонність
    has_seasonal <- (s > 1) && (P > 0 || Q > 0)

    if (!has_seasonal && s == 1) {
      # Несезонна модель - sarima_pmm2 потребує s > 1
      if (verbose) cat("  Note: PMM2 requires seasonal period > 1, using MLE for non-seasonal MA\n")
      ma_method_used <- "mle (pmm2 unavailable)"
    } else {
      # Застосовуємо PMM2
      pmm2_result <- tryCatch({
        sarima_pmm2(x,
                    order = c(p, P, q, Q),
                    seasonal = list(order = c(d, D), period = s),
                    method = "pmm2",
                    max_iter = pmm2_max_iter,
                    verbose = verbose)
      }, error = function(e) {
        if (verbose) cat("  PMM2 failed:", e$message, "\n")
        return(NULL)
      })

      if (!is.null(pmm2_result)) {
        # Перевірка збіжності (sarima_pmm2 повертає S4 з convergence = TRUE/FALSE)
        converged <- FALSE
        coeffs_pmm2 <- NULL

        # Отримуємо convergence (це logical, не string!)
        conv_slot <- tryCatch(slot(pmm2_result, "convergence"), error = function(e) {
          if (verbose) cat("  Warning: Cannot access convergence slot:", e$message, "\n")
          return(FALSE)
        })

        converged <- isTRUE(conv_slot)

        if (converged) {
          # Отримуємо coefficients (це вектор!)
          coeffs_pmm2 <- tryCatch(slot(pmm2_result, "coefficients"), error = function(e) {
            if (verbose) cat("  Warning: Cannot access coefficients slot:", e$message, "\n")
            return(NULL)
          })

          # Отримуємо order для розуміння структури
          order_slot <- tryCatch(slot(pmm2_result, "order"), error = function(e) NULL)

          if (!is.null(coeffs_pmm2) && length(coeffs_pmm2) > 0 && !is.null(order_slot)) {
            # Розбираємо coefficients vector згідно з order
            # Вектор містить: [AR coeffs, SAR coeffs, MA coeffs, SMA coeffs]
            # (без intercept в цьому векторі)

            idx <- 1
            n_ar <- if (is.null(order_slot$ar)) 0 else order_slot$ar
            n_sar <- if (is.null(order_slot$sar)) 0 else order_slot$sar
            n_ma <- if (is.null(order_slot$ma)) 0 else order_slot$ma
            n_sma <- if (is.null(order_slot$sma)) 0 else order_slot$sma

            # Skip AR coefficients
            if (n_ar > 0) {
              idx <- idx + n_ar
            }

            # Skip SAR coefficients
            if (n_sar > 0) {
              idx <- idx + n_sar
            }

            # Extract MA coefficients from PMM2
            if (n_ma > 0 && idx <= length(coeffs_pmm2)) {
              ma_end <- min(idx + n_ma - 1, length(coeffs_pmm2))
              ma_coef_final <- coeffs_pmm2[idx:ma_end]
              idx <- ma_end + 1
            }

            # Extract SMA coefficients from PMM2
            if (n_sma > 0 && idx <= length(coeffs_pmm2)) {
              sma_end <- min(idx + n_sma - 1, length(coeffs_pmm2))
              sma_coef_final <- coeffs_pmm2[idx:sma_end]
            }

            pmm2_fit <- pmm2_result
            ma_method_used <- "pmm2"

            if (verbose) {
              cat("  PMM2 converged successfully\n")
              if (n_ma > 0) cat("    MA (PMM2):", paste(round(ma_coef_final, 4), collapse = ", "), "\n")
              if (n_sma > 0) cat("    SMA (PMM2):", paste(round(sma_coef_final, 4), collapse = ", "), "\n")
            }
          } else {
            if (verbose) cat("  PMM2 coefficients not available, using MLE estimates\n")
            ma_method_used <- "mle (pmm2 no coefs)"
          }
        } else {
          if (verbose) cat("  PMM2 did not converge, using MLE estimates\n")
          ma_method_used <- "mle (pmm2 failed)"
        }
      } else {
        if (verbose) cat("  Using MLE estimates (PMM2 failed)\n")
        ma_method_used <- "mle (pmm2 failed)"
      }
    }
    cat("\n")
  }

  # ===========================================================================
  # STEP 3: Apply PMM2 for AR/SAR if requested
  # ===========================================================================
  ar_coef_final <- ar_coef_mle
  sar_coef_final <- sar_coef_mle
  ar_method_used <- "mle"

  if (ar_method == "pmm2" && (p > 0 || P > 0)) {
    if (verbose) cat("Step 3: Re-estimating AR/SAR via PMM2...\n")

    # TODO: Реалізувати PMM2 для AR на фіксованих інноваціях
    # Поки що використовуємо MLE
    if (verbose) cat("  Note: AR PMM2 re-estimation not yet implemented, using MLE\n\n")
    ar_method_used <- "mle (ar pmm2 not implemented)"
  }

  # ===========================================================================
  # STEP 4: Determine final method description
  # ===========================================================================
  if (p == 0 && P == 0) {
    # Pure MA model
    method_desc <- paste0(toupper(ma_method_used), "(MA)-only")
  } else if (q == 0 && Q == 0) {
    # Pure AR model
    method_desc <- paste0(toupper(ar_method_used), "(AR)-only")
  } else {
    # Mixed model
    method_desc <- paste0(toupper(ma_method_used), "(MA) + ",
                         toupper(ar_method_used), "(AR)")
  }

  # ===========================================================================
  # STEP 5: Return results
  # ===========================================================================
  list(
    ma_coef = ma_coef_final,
    sma_coef = sma_coef_final,
    ar_coef = ar_coef_final,
    sar_coef = sar_coef_final,
    mean = mean_est,
    innovations = innovations,
    method = method_desc,
    ma_method_used = ma_method_used,
    ar_method_used = ar_method_used,
    mle_fit = mle_fit,
    pmm2_fit = pmm2_fit
  )
}


#' Порівняльний тест: MLE vs PMM2 для SMA
#'
#' @param x Time series data
#' @param order ARIMA order
#' @param seasonal Seasonal specification
#' @param include.mean Include intercept
#' @param verbose Print details
#'
#' @return List with both MLE and PMM2 results
#'
#' @export
compare_mle_pmm2_for_sma <- function(x,
                                      order = c(0, 0, 0),
                                      seasonal = list(order = c(0, 0, 1), period = 4),
                                      include.mean = FALSE,
                                      verbose = TRUE) {
  if (verbose) {
    cat("=== Comparing MLE vs PMM2 for SMA ===\n\n")
  }

  # MLE approach
  if (verbose) cat("--- MLE Approach ---\n")
  result_mle <- hybrid_sarima_estimator_extended(x, order, seasonal,
                                                  include.mean = include.mean,
                                                  ma_method = "mle",
                                                  verbose = verbose)

  # PMM2 approach
  if (verbose) cat("\n--- PMM2 Approach ---\n")
  result_pmm2 <- hybrid_sarima_estimator_extended(x, order, seasonal,
                                                   include.mean = include.mean,
                                                   ma_method = "pmm2",
                                                   verbose = verbose)

  # Comparison
  if (verbose) {
    cat("\n=== Comparison ===\n")

    Q <- seasonal$order[3]
    if (Q > 0) {
      cat("SMA coefficients:\n")
      cat(sprintf("  MLE:  %s\n", paste(round(result_mle$sma_coef, 4), collapse = ", ")))
      cat(sprintf("  PMM2: %s\n", paste(round(result_pmm2$sma_coef, 4), collapse = ", ")))

      diff <- result_pmm2$sma_coef - result_mle$sma_coef
      cat(sprintf("  Diff: %s\n", paste(sprintf("%+.4f", diff), collapse = ", ")))
    }

    q <- order[3]
    if (q > 0) {
      cat("\nMA coefficients:\n")
      cat(sprintf("  MLE:  %s\n", paste(round(result_mle$ma_coef, 4), collapse = ", ")))
      cat(sprintf("  PMM2: %s\n", paste(round(result_pmm2$ma_coef, 4), collapse = ", ")))

      diff <- result_pmm2$ma_coef - result_mle$ma_coef
      cat(sprintf("  Diff: %s\n", paste(sprintf("%+.4f", diff), collapse = ", ")))
    }
  }

  list(
    mle = result_mle,
    pmm2 = result_pmm2
  )
}
