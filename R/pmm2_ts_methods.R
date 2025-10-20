# pmm2_ts_methods.R - Методи для роботи з об'єктами моделей часових рядів

#' Витягнути коефіцієнти з об'єкта TS2fit
#'
#' @param object Об'єкт TS2fit
#' @param ... Додаткові аргументи (не використовуються)
#'
#' @return Іменований вектор коефіцієнтів
#' @export
setMethod("coef", "TS2fit",
          function(object, ...) {
            # Отримати параметри моделі
            model_type <- object@model_type
            ar_order <- object@order$ar
            ma_order <- object@order$ma

            # Витягнути та іменувати AR коефіцієнти
            if(ar_order > 0) {
              ar_coefs <- object@coefficients[1:ar_order]
              names(ar_coefs) <- paste0("ar", 1:ar_order)
            } else {
              ar_coefs <- numeric(0)
            }

            # Витягнути та іменувати MA коефіцієнти
            if(ma_order > 0) {
              ma_coefs <- object@coefficients[(ar_order+1):(ar_order+ma_order)]
              names(ma_coefs) <- paste0("ma", 1:ma_order)
            } else {
              ma_coefs <- numeric(0)
            }

            # Об'єднати коефіцієнти
            result <- c(ar_coefs, ma_coefs)

            # Додати перехоплення, якщо присутнє
            if(object@intercept != 0) {
              result <- c(intercept = object@intercept, result)
            }

            return(result)
          })

#' Витягнути залишки з об'єкта TS2fit
#'
#' @param object Об'єкт TS2fit
#' @param ... Додаткові аргументи (не використовуються)
#'
#' @return Вектор залишків (інновацій)
#' @export
setMethod("residuals", "TS2fit",
          function(object, ...) {
            object@residuals
          })

#' Отримати підігнані значення для AR моделі
#'
#' @param object Об'єкт TS2fit з model_type="ar"
#' @return Вектор підігнаних значень
#' @keywords internal
get_ar_fitted <- function(object) {
  if(object@model_type != "ar") {
    stop("Ця функція лише для AR моделей")
  }

  x <- object@original_series
  ar_order <- object@order$ar
  ar_coef <- object@coefficients[1:ar_order]
  intercept <- object@intercept

  if(intercept != 0) {
    x_centered <- x - intercept
  } else {
    x_centered <- x
  }

  # Створити матрицю дизайну та обчислити підігнані значення
  X <- create_ar_matrix(x_centered, ar_order)
  fitted <- as.vector(X %*% ar_coef) + intercept
  return(fitted)
}

#' Витягнути підігнані значення з об'єкта TS2fit
#'
#' @param object Об'єкт TS2fit
#' @param ... Додаткові аргументи (не використовуються)
#'
#' @return Вектор підігнаних значень
#' @export
setMethod("fitted", "TS2fit",
          function(object, ...) {
            # Отримати тип моделі
            model_type <- object@model_type

            # Обчислити підігнані значення в залежності від типу моделі
            if(model_type == "ar") {
              # Для AR моделей використовувати пряме обчислення
              fitted_values <- get_ar_fitted(object)
            } else {
              # Для інших моделей: підігнані = оригінальні мінус залишки
              orig <- object@original_series
              resid <- object@residuals

              # Вирівняти довжини (часто залишки коротші через початкові значення)
              len_diff <- length(orig) - length(resid)
              if(len_diff > 0 && !all(is.na(resid))) {
                # Знайти перше не-NA значення в resid
                first_valid <- min(which(!is.na(resid)))

                # Побудувати вектор fitted з NA в початку
                fitted_values <- rep(NA, length(orig))
                valid_indices <- first_valid:length(resid)

                # Встановити дійсні значення
                fitted_values[(len_diff + valid_indices)] <-
                  orig[(len_diff + valid_indices)] - resid[valid_indices]
              } else {
                fitted_values <- orig - resid
              }
            }

            return(fitted_values)
          })

#' Побудувати діагностичні графіки для об'єктів TS2fit
#'
#' @param x Об'єкт TS2fit
#' @param y Не використовується (для сумісності методу S4)
#' @param which Цілочисельний вектор, що вказує, які графіки виробляти
#' @param ... додаткові аргументи, передані функціям побудови графіків
#'
#' @return Невидимо повертає x
#'
#' @export
setMethod("plot", signature(x = "TS2fit", y = "missing"),
          function(x, y, which = c(1:4), ...) {
            op <- par(no.readonly = TRUE)
            on.exit(par(op))

            # Отримати параметри моделі
            model_type <- x@model_type
            ar_order <- x@order$ar
            ma_order <- x@order$ma
            d <- x@order$d

            # Макет графіка за замовчуванням
            par(mfrow = c(2, 2))

            # Для моделей ARIMA ми можемо захотіти побудувати оригінальний/диференційований ряд також
            if(model_type == "arima" && length(which) > 4) {
              par(mfrow = c(3, 2))
            }

            # Обчислити підігнані значення та залишки
            residuals <- as.numeric(x@residuals)
            fitted <- fitted(x)

            # Визначити, які графіки відображати
            plot_idx <- 1
            n_plots <- min(length(which), 6) # Максимум 6 графіків

            # Для моделей ARIMA, ми можемо хотіти інші графіки
            if(model_type == "arima") {
              # Графік 1: Оригінальний часовий ряд (тільки ARIMA)
              if(1 %in% which && plot_idx <= n_plots) {
                plot(x@original_series, type = "l",
                     main = "Оригінальний часовий ряд",
                     xlab = "Час",
                     ylab = "Значення",
                     ...)
                plot_idx <- plot_idx + 1
              }

              # Графік 2: Диференційований часовий ряд (тільки ARIMA)
              if(2 %in% which && plot_idx <= n_plots && d > 0) {
                diff_series <- diff(x@original_series, differences = d)
                plot(diff_series, type = "l",
                     main = paste0("Диференційований ряд (d=", d, ")"),
                     xlab = "Час",
                     ylab = "Значення",
                     ...)
                plot_idx <- plot_idx + 1
              }
            }

            # Стандартні графіки для всіх типів моделей
            # Графік: Залишки vs Підігнані
            if(3 %in% which && plot_idx <= n_plots) {
              plot(fitted, residuals,
                   main = "Залишки vs Підігнані",
                   xlab = "Підігнані значення",
                   ylab = "Залишки",
                   ...)
              abline(h = 0, lty = 2)
              lines(lowess(fitted, residuals), col = "red")
              plot_idx <- plot_idx + 1
            }

            # Графік: Нормальний Q-Q графік
            if(4 %in% which && plot_idx <= n_plots) {
              qqnorm(residuals, main = "Нормальний Q-Q графік", ...)
              qqline(residuals)
              plot_idx <- plot_idx + 1
            }

            # Графік: ACF залишків
            if(5 %in% which && plot_idx <= n_plots) {
              acf(residuals, main = "ACF залишків", ...)
              plot_idx <- plot_idx + 1
            }

            # Графік: Гістограма залишків
            if(6 %in% which && plot_idx <= n_plots) {
              hist(residuals,
                   main = "Гістограма залишків",
                   xlab = "Залишки",
                   breaks = "FD",
                   ...)
              plot_idx <- plot_idx + 1
            }

            invisible(x)
          })

#' Метод прогнозування для об'єктів TS2fit
#'
#' @param object Об'єкт TS2fit
#' @param n.ahead Кількість кроків вперед для прогнозування
#' @param ... додаткові аргументи (не використовуються)
#'
#' @return Вектор або список прогнозів, залежно від типу моделі
#'
#' @export
setMethod("predict", "TS2fit",
          function(object, n.ahead = 1, ...) {
            # Отримати параметри моделі
            model_type <- object@model_type
            ar_order <- object@order$ar
            ma_order <- object@order$ma
            d <- object@order$d
            intercept <- object@intercept

            # Витягнути коефіцієнти
            if(ar_order > 0) {
              ar_coef <- object@coefficients[1:ar_order]
            } else {
              ar_coef <- numeric(0)
            }

            if(ma_order > 0) {
              ma_coef <- object@coefficients[(ar_order+1):(ar_order+ma_order)]
            } else {
              ma_coef <- numeric(0)
            }

            # Для AR моделей, реалізувати пряме прогнозування
            if(model_type == "ar") {
              x <- object@original_series
              n <- length(x)
              pred <- numeric(n.ahead)

              # Генерувати прогнози
              for(i in 1:n.ahead) {
                # Використовувати оригінальні дані та попередні прогнози за потреби
                lags <- numeric(ar_order)
                for(j in 1:ar_order) {
                  if(i - j <= 0) {
                    # Використовувати оригінальні дані
                    lags[j] <- x[n - j + i]
                  } else {
                    # Використовувати попередні прогнози
                    lags[j] <- pred[i - j]
                  }
                }

                # Обчислити прогноз
                pred[i] <- sum(ar_coef * lags) + intercept
              }

              return(pred)

            } else if(model_type == "ma") {
              # Для MA моделей, прогнози за межами порядку - це просто середнє
              innovations <- object@residuals

              # Генерувати прогнози MA
              ma_pred <- function(innovations, ma_coef, n.ahead) {
                n <- length(innovations)
                q <- length(ma_coef)
                pred <- numeric(n.ahead)

                for(i in 1:n.ahead) {
                  for(j in 1:min(i, q)) {
                    if((n - i + j) > 0) {
                      pred[i] <- pred[i] + ma_coef[j] * innovations[n - i + j]
                    }
                  }
                }
                return(pred)
              }

              if(n.ahead > ma_order) {
                return(c(ma_pred(innovations, ma_coef, ma_order),
                         rep(intercept, n.ahead - ma_order)))
              } else {
                return(ma_pred(innovations, ma_coef, n.ahead))
              }

            } else {
              # Для ARMA та ARIMA моделей, використовувати прогнози stats::arima
              # які правильно обробляють обидва компоненти

              # Налаштувати модель arima з фіксованими параметрами
              arima_order <- c(ar_order, ifelse(model_type == "arima", d, 0), ma_order)

              # Використовувати функцію прогнозування з пакету stats
              arima_pred <- stats::predict(
                stats::arima(object@original_series,
                             order = arima_order,
                             include.mean = (intercept != 0),
                             fixed = c(ar_coef, ma_coef, if(intercept != 0) intercept else NULL)),
                n.ahead = n.ahead
              )

              return(arima_pred)
            }
          })

#' Порівняти PMM2 з класичними методами оцінювання часових рядів
#'
#' @param x Числовий вектор даних часового ряду
#' @param order Специфікація порядку моделі (див. ts_pmm2 для формату)
#' @param model_type Тип моделі: "ar", "ma", "arma", або "arima"
#' @param include.mean Логічне, чи включати член перехоплення
#' @param pmm2_args Список додаткових аргументів для передачі в ts_pmm2()
#'
#' @return Список з підігнаними моделями та таблицями порівняння
#' @export
compare_ts_methods <- function(x, order, model_type = c("ar", "ma", "arma", "arima"),
                               include.mean = TRUE, pmm2_args = list()) {
  # Вибрати аргумент model_type
  model_type <- match.arg(model_type)

  # Підготувати порівняння моделей на основі model_type
  if(model_type == "ar") {
    # Для AR моделей
    # Підігнати AR модель за допомогою методу Юла-Волкера
    yw_fit <- stats::ar(x, order.max = order, aic = FALSE, method = "yw",
                        demean = include.mean)

    # Підігнати AR модель за допомогою методу OLS
    ols_fit <- stats::ar(x, order.max = order, aic = FALSE, method = "ols",
                         demean = include.mean)

    # Підігнати AR модель за допомогою методу MLE
    mle_fit <- stats::ar(x, order.max = order, aic = FALSE, method = "mle",
                         demean = include.mean)

    # Підігнати AR модель за допомогою PMM2
    pmm2_args <- c(list(x = x, order = order, model_type = "ar",
                        include.mean = include.mean), pmm2_args)
    pmm2_fit <- do.call(ts_pmm2, pmm2_args)

    # Витягнути коефіцієнти
    coef_yw <- yw_fit$ar
    coef_ols <- ols_fit$ar
    coef_mle <- mle_fit$ar
    coef_pmm2 <- pmm2_fit@coefficients

    # Обчислити залишки
    res_yw <- yw_fit$resid[!is.na(yw_fit$resid)]
    res_ols <- ols_fit$resid[!is.na(ols_fit$resid)]
    res_mle <- mle_fit$resid[!is.na(mle_fit$resid)]
    res_pmm2 <- pmm2_fit@residuals

    methods <- c("YW", "OLS", "MLE", "PMM2")

    result_list <- list(
      yw = yw_fit,
      ols = ols_fit,
      mle = mle_fit,
      pmm2 = pmm2_fit
    )

  } else if(model_type %in% c("ma", "arma", "arima")) {
    # Для MA, ARMA та ARIMA моделей

    # Підготувати порядок arima на основі типу моделі
    if(model_type == "ma") {
      arima_order <- c(0, 0, order)
    } else if(model_type == "arma") {
      arima_order <- c(order[1], 0, order[2])
    } else {
      arima_order <- order
    }

    # Підігнати модель за допомогою методу CSS
    css_fit <- arima(x, order = arima_order, method = "CSS", include.mean = include.mean)

    # Підігнати модель за допомогою методу ML
    ml_fit <- arima(x, order = arima_order, method = "ML", include.mean = include.mean)

    # Підігнати модель за допомогою PMM2
    pmm2_args <- c(list(x = x, order = order, model_type = model_type,
                        include.mean = include.mean), pmm2_args)
    pmm2_fit <- do.call(ts_pmm2, pmm2_args)

    # Витягнути назви коефіцієнтів AR та MA на основі типу моделі
    if(model_type == "ma") {
      ar_names <- character(0)
      ma_names <- paste0("ma", 1:order)
    } else if(model_type == "arma") {
      ar_names <- paste0("ar", 1:order[1])
      ma_names <- paste0("ma", 1:order[2])
    } else {
      ar_names <- if(order[1] > 0) paste0("ar", 1:order[1]) else character(0)
      ma_names <- if(order[3] > 0) paste0("ma", 1:order[3]) else character(0)
    }

    coef_names <- c(ar_names, ma_names)

    # Витягнути коефіцієнти
    coef_css <- as.numeric(css_fit$coef[coef_names])
    coef_ml <- as.numeric(ml_fit$coef[coef_names])
    coef_pmm2 <- pmm2_fit@coefficients

    # Обчислити залишки
    res_css <- residuals(css_fit)
    res_ml <- residuals(ml_fit)
    res_pmm2 <- pmm2_fit@residuals

    methods <- c("CSS", "ML", "PMM2")

    result_list <- list(
      css = css_fit,
      ml = ml_fit,
      pmm2 = pmm2_fit
    )
  }

  # Обчислити статистику залишків для всіх методів
  residuals_list <- if(model_type == "ar") {
    list(res_yw, res_ols, res_mle, res_pmm2)
  } else {
    list(res_css, res_ml, res_pmm2)
  }

  compute_res_stats <- function(res) {
    m2 <- mean(res^2, na.rm = TRUE)
    m3 <- mean(res^3, na.rm = TRUE)
    m4 <- mean(res^4, na.rm = TRUE)

    c(RSS = sum(res^2, na.rm = TRUE),
      MAE = mean(abs(res), na.rm = TRUE),
      Skewness = m3 / m2^(3/2),
      Kurtosis = m4 / m2^2)
  }

  res_stats <- data.frame(
    Method = methods,
    do.call(rbind, lapply(residuals_list, compute_res_stats))
  )

  # Створити таблицю порівняння коефіцієнтів
  if(model_type == "ar") {
    coef_names <- paste0("ar", 1:order)
    coef_values <- list(coef_yw, coef_ols, coef_mle, coef_pmm2)
  } else if(model_type == "ma") {
    coef_names <- paste0("ma", 1:order)
    coef_values <- list(coef_css, coef_ml, coef_pmm2)
  } else {
    coef_values <- list(coef_css, coef_ml, coef_pmm2)
  }

  coef_table <- data.frame(
    Coefficient = coef_names,
    do.call(cbind, lapply(seq_along(methods), function(i) {
      result <- coef_values[[i]]
      names(result) <- methods[i]
      return(result)
    }))
  )

  # Повернути результати
  result_list$coefficients <- coef_table
  result_list$residual_stats <- res_stats

  return(result_list)
}

#' Порівняти методи AR
#'
#' @inheritParams compare_ts_methods
#' @export
compare_ar_methods <- function(x, order = 1, include.mean = TRUE, pmm2_args = list()) {
  compare_ts_methods(x, order = order, model_type = "ar",
                     include.mean = include.mean, pmm2_args = pmm2_args)
}

#' Порівняти методи MA
#'
#' @inheritParams compare_ts_methods
#' @export
compare_ma_methods <- function(x, order = 1, include.mean = TRUE, pmm2_args = list()) {
  compare_ts_methods(x, order = order, model_type = "ma",
                     include.mean = include.mean, pmm2_args = pmm2_args)
}

#' Порівняти методи ARMA
#'
#' @inheritParams compare_ts_methods
#' @export
compare_arma_methods <- function(x, order = c(1, 1), include.mean = TRUE, pmm2_args = list()) {
  compare_ts_methods(x, order = order, model_type = "arma",
                     include.mean = include.mean, pmm2_args = pmm2_args)
}

#' Порівняти методи ARIMA
#'
#' @inheritParams compare_ts_methods
#' @export
compare_arima_methods <- function(x, order = c(1, 1, 1), include.mean = TRUE, pmm2_args = list()) {
  compare_ts_methods(x, order = order, model_type = "arima",
                     include.mean = include.mean, pmm2_args = pmm2_args)
}
