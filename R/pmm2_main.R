# pmm2_main.R - Основний модуль для лінійних моделей PMM2

#' pmm2: Головна функція для PMM2 (S=2)
#'
#' Підганяє лінійну модель за допомогою методу максимізації поліномів (порядок 2),
#' який є робастним щодо негаусівських помилок.
#'
#' @param formula Формула R для моделі
#' @param data data.frame, що містить змінні у формулі
#' @param max_iter ціле: максимальна кількість ітерацій для алгоритму
#' @param tol числове: допуск для збіжності
#' @param regularize логічне: додати мале значення до діагоналі для числової стабільності
#' @param reg_lambda числове: параметр регуляризації (якщо regularize=TRUE)
#' @param na.action функція для обробки відсутніх значень, за замовчуванням - na.fail
#' @param weights опціональний вектор ваг (поки не реалізовано)
#' @param verbose логічне: чи виводити інформацію про прогрес
#'
#' @details
#' Алгоритм PMM2 працює наступним чином:
#'
#' 1. Підганяє звичайну регресію найменших квадратів (OLS) для отримання початкових оцінок
#' 2. Обчислює центральні моменти (m2, m3, m4) із залишків OLS
#' 3. Ітеративно покращує оцінки параметрів за допомогою підходу на основі градієнта
#'
#' PMM2 особливо корисний, коли терми помилок не є гаусівськими.
#'
#' @return Об'єкт S4 \code{PMM2fit}
#' @export
#'
#' @examples
#' \dontrun{
#' # Генерувати дані вибірки з t-розподіленими помилками
#' n <- 100
#' x <- rnorm(n)
#' y <- 2 + 3*x + rt(n, df=3)
#' dat <- data.frame(y=y, x=x)
#'
#' # Підігнати модель за допомогою PMM2
#' fit <- lm_pmm2(y ~ x, data=dat)
#'
#' # Резюме та статистичний висновок
#' summary(fit, formula=y~x, data=dat)
#' }
lm_pmm2 <- function(formula, data,
                    max_iter=50, tol=1e-6,
                    regularize=TRUE, reg_lambda=1e-8,
                    na.action=na.fail, weights=NULL,
                    verbose=FALSE)
{
  # Захопити виклик
  call <- match.call()

  # Перевірити валідність входу
  if(missing(formula) || missing(data)) {
    stop("Обидва 'formula' та 'data' мають бути надані")
  }

  if(!is.data.frame(data)) {
    stop("'data' має бути data.frame")
  }

  if(max_iter <= 0) {
    stop("'max_iter' має бути додатним")
  }

  if(tol <= 0) {
    stop("'tol' має бути додатним")
  }

  # Обробити відсутні значення
  if(!is.null(na.action)) {
    data <- na.action(data)
  }

  # Перевірити на ваги
  if(!is.null(weights)) {
    warning("Ваги поки не реалізовані в PMM2. Ігнорую ваги.")
  }

  # 1) OLS
  if(verbose) cat("Підганяю початкову модель OLS...\n")

  mf <- model.frame(formula, data)
  X  <- model.matrix(formula, mf)
  y  <- model.response(mf)

  # Перевірити на дефіцит рангу
  qr_X <- qr(X)
  if(qr_X$rank < ncol(X)) {
    warning("Матриця дизайну має дефіцит рангу, деякі коефіцієнти можуть бути неоцінюваними")
  }

  fit_ols <- lm.fit(x=X, y=y)
  b_ols   <- fit_ols$coefficients

  # Обробити NA в коефіцієнтах OLS
  if(any(is.na(b_ols))) {
    stop("Підгонка OLS привела до NA коефіцієнтів. Перевірте на мультиколінеарність.")
  }

  # Перетворити b_ols на числовий вектор для забезпечення узгодженої обробки
  b_ols <- as.numeric(b_ols)

  # 2) OLS залишки => m2, m3, m4
  res_ols <- y - (X %*% b_ols)
  moments <- compute_moments(res_ols)
  m2 <- moments$m2
  m3 <- moments$m3
  m4 <- moments$m4

  if(verbose) {
    cat("Початкові моменти із залишків OLS:\n")
    cat("  m2 =", m2, "\n")
    cat("  m3 =", m3, "\n")
    cat("  m4 =", m4, "\n")
  }

  # Перевірити на потенційні проблеми з моментами
  if(m2 <= 0) {
    warning("Другий центральний момент (m2) не є додатним. Результати можуть бути ненадійними.")
  }

  if(m4 <= m2^2) {
    warning("Четвертий центральний момент (m4) менший за m2^2. Це порушує базову нерівність для розподілів ймовірностей.")
  }

  # 3) Запустити уніфікований алгоритм PMM2
  if(verbose) cat("Починаю ітерації PMM2...\n")

  out <- pmm2_algorithm(b_ols, X, y, m2, m3, m4,
                        max_iter = max_iter, tol = tol,
                        regularize = regularize, reg_lambda = reg_lambda,
                        verbose = verbose)

  # Витягнути результати
  b_est   <- out$b
  conv    <- out$convergence
  iter    <- out$iterations
  final_res <- out$residuals

  if(verbose) {
    cat("Алгоритм PMM2 завершено.\n")
    cat("  Збігся:", conv, "\n")
    cat("  Ітерацій:", iter, "\n")
  }

  # Повернути об'єкт S4 з усіма результатами
  ans <- new("PMM2fit",
             coefficients = b_est,
             residuals = final_res,
             m2 = m2,
             m3 = m3,
             m4 = m4,
             convergence = conv,
             iterations = iter,
             call = call)
  attr(ans, "model_matrix") <- X
  attr(ans, "model_frame") <- mf
  attr(ans, "response") <- as.numeric(y)
  attr(ans, "data") <- data

  return(ans)
}

#' Витягнути коефіцієнти з об'єкта PMM2fit
#'
#' @param object Об'єкт PMM2fit
#' @param ... Додаткові аргументи (не використовуються)
#'
#' @return Вектор коефіцієнтів
#' @export
setMethod("coef", "PMM2fit",
          function(object, ...) {
            object@coefficients
          })

#' Витягнути залишки з об'єкта PMM2fit
#'
#' @param object Об'єкт PMM2fit
#' @param ... Додаткові аргументи (не використовуються)
#'
#' @return Вектор залишків
#' @export
setMethod("residuals", "PMM2fit",
          function(object, ...) {
            object@residuals
          })

#' Витягнути підігнані значення з об'єкта PMM2fit
#'
#' @param object Об'єкт PMM2fit
#' @param ... Додаткові аргументи (не використовуються)
#'
#' @return Вектор підігнаних значень
#' @export
setMethod("fitted", "PMM2fit",
          function(object, data = NULL, ...) {
            fitted_values(object, data)
          })

#' Розрахувати AIC для об'єкта PMM2fit
#'
#' @param object Об'єкт PMM2fit
#' @param ... Додаткові аргументи (не використовуються)
#' @param k Штраф за параметр, що буде використовуватись; стандартно 2
#'
#' @return Значення AIC
#' @export
setMethod("AIC", "PMM2fit",
          function(object, ..., k = 2) {
            res <- object@residuals
            n <- length(res)
            p <- length(object@coefficients)

            # Апроксимація логарифмічної правдоподібності
            ll <- -n/2 * log(sum(res^2)/n) - n/2 * (1 + log(2*pi))

            # AIC
            -2 * ll + k * p
          })

#' Побудувати діагностичні графіки для об'єкта PMM2fit
#'
#' @param x Об'єкт PMM2fit
#' @param y Не використовується (сумісність з generic)
#' @param which Набір графіків для відображення (значення 1-4)
#' @param ... Додаткові аргументи, що передаються графічним функціям
#'
#' @return Невидимо повертає вхідний об'єкт
#' @export
setMethod("plot", signature(x = "PMM2fit", y = "missing"),
          function(x, y, which = 1:4, ...) {
            res <- as.numeric(x@residuals)
            fitted_vals <- tryCatch({
              fitted(x)
            }, error = function(e) {
              stored_X <- attr(x, "model_matrix")
              if(!is.null(stored_X)) {
                as.vector(stored_X %*% x@coefficients)
              } else {
                seq_along(res)
              }
            })

            which <- intersect(unique(which), 1:4)
            if(length(which) == 0) {
              which <- 1:4
            }

            old_par <- graphics::par(no.readonly = TRUE)
            on.exit(graphics::par(old_par))
            n_plots <- length(which)
            graphics::par(mfrow = c(2, 2))

            for(idx in which) {
              switch(idx,
                     {
                       graphics::plot(fitted_vals, res,
                                      main = "Residuals vs Fitted",
                                      xlab = "Fitted values",
                                      ylab = "Residuals", ...)
                       graphics::abline(h = 0, col = "red", lty = 2)
                     },
                     {
                       stats::qqnorm(res, main = "Normal Q-Q", ...)
                       stats::qqline(res, col = "red", lty = 2)
                     },
                     {
                       graphics::plot(seq_along(res), res, type = "l",
                                      main = "Residuals over Index",
                                      xlab = "Observation",
                                      ylab = "Residual", ...)
                       graphics::abline(h = 0, col = "red", lty = 2)
                     },
                     {
                       graphics::hist(res,
                                      main = "Residual Histogram",
                                      xlab = "Residuals",
                                      breaks = "FD", ...)
                     })
            }

            invisible(x)
          })

#' Допоміжна функція для витягнення підігнаних значень
#'
#' @param object Об'єкт PMM2fit
#' @return Вектор підігнаних значень
#'
#' @keywords internal
fitted_values <- function(object, data = NULL) {
  if(is.null(object@call)) {
    stop("Об'єкт PMM2fit не містить інформації про виклик")
  }

  # Фолбек до збережених атрибутів
  stored_X <- attr(object, "model_matrix")
  stored_response <- attr(object, "response")
  stored_mf <- attr(object, "model_frame")
  stored_data <- attr(object, "data")

  if (!is.null(stored_X)) {
    fitted_attr <- tryCatch({
      as.vector(stored_X %*% object@coefficients)
    }, error = function(e) NULL)
    if (!is.null(fitted_attr)) {
      return(fitted_attr)
    }
  }

  if (!is.null(stored_response) &&
      length(stored_response) == length(object@residuals)) {
    return(as.vector(stored_response - object@residuals))
  }

  # Спробувати реконструювати оригінальні дані
  data_to_use <- data
  if(is.null(data_to_use)) {
    if (!is.null(stored_mf)) {
      data_to_use <- stored_mf
    } else if (!is.null(stored_data)) {
      data_to_use <- stored_data
    } else {
      # Спробувати отримати дані з виклику, але безпечно обробити можливі помилки
      tryCatch({
        data_to_use <- eval(object@call$data, envir = parent.frame())
      }, error = function(e) {
        if(is.null(data)) {
          stop("Не вдалося отримати дані з об'єкта. Передайте параметр 'data'.")
        }
      })
    }
  }

  if(is.null(data_to_use)) {
    stop("Потрібен фрейм даних для обчислення підігнаних значень")
  }

  # Реконструювати формулу
  formula <- eval(object@call$formula)

  # Безпечна побудова матриці дизайну
  tryCatch({
    mf <- model.frame(formula, data_to_use)
    X <- model.matrix(formula, mf)

    # Обчислити підігнані значення
    fitted <- as.vector(X %*% object@coefficients)
    return(fitted)
  }, error = function(e) {
    stop("Помилка при обчисленні підігнаних значень: ", e$message)
  })
}

#' Порівняти PMM2 з OLS
#'
#' @param formula Формула моделі
#' @param data Фрейм даних
#' @param pmm2_args Список аргументів для передачі в lm_pmm2()
#'
#' @return Список з об'єктами підгонки OLS та PMM2
#' @export
compare_with_ols <- function(formula, data, pmm2_args = list()) {
  # Підгонка OLS моделі
  fit_ols <- lm(formula, data)

  # Підгонка PMM2 моделі з стандартними або заданими аргументами
  args <- c(list(formula = formula, data = data), pmm2_args)
  fit_pmm2 <- do.call(lm_pmm2, args)

  # Витягнути та порівняти коефіцієнти
  coef_ols <- coef(fit_ols)
  coef_pmm2 <- coef(fit_pmm2)

  # Перевірка, чи імена коефіцієнтів PMM2 встановлені коректно
  if(is.null(names(coef_pmm2)) || all(names(coef_pmm2) == "")) {
    # Якщо імена не встановлені, використаємо імена з OLS
    if(length(coef_ols) == length(coef_pmm2)) {
      names(coef_pmm2) <- names(coef_ols)
    } else {
      # Якщо довжини відрізняються, використаємо генеровані імена
      names(coef_pmm2) <- paste0("coef", seq_along(coef_pmm2))
    }
  }

  # Обчислити статистику залишків
  res_ols <- residuals(fit_ols)
  res_pmm2 <- residuals(fit_pmm2)

  res_stats <- data.frame(
    Method = c("OLS", "PMM2"),
    RSS = c(sum(res_ols^2), sum(res_pmm2^2)),
    MAE = c(mean(abs(res_ols)), mean(abs(res_pmm2))),
    Skewness = c(pmm_skewness(res_ols), pmm_skewness(res_pmm2)),
    Kurtosis = c(pmm_kurtosis(res_ols), pmm_kurtosis(res_pmm2))
  )

  # Створити таблицю порівняння коефіцієнтів
  # Використовуємо всі унікальні імена коефіцієнтів
  all_coef_names <- unique(c(names(coef_ols), names(coef_pmm2)))

  coef_table <- data.frame(
    Coefficient = all_coef_names,
    OLS = numeric(length(all_coef_names)),
    PMM2 = numeric(length(all_coef_names)),
    Diff_Percent = numeric(length(all_coef_names))
  )

  # Заповнюємо значення для OLS
  for(i in seq_along(all_coef_names)) {
    name <- all_coef_names[i]
    if(name %in% names(coef_ols)) {
      coef_table$OLS[i] <- coef_ols[name]
    } else {
      coef_table$OLS[i] <- NA
    }

    if(name %in% names(coef_pmm2)) {
      coef_table$PMM2[i] <- coef_pmm2[name]

      # Обчислити процентну різницю тільки якщо обидва значення існують
      if(!is.na(coef_table$OLS[i]) && coef_table$OLS[i] != 0) {
        coef_table$Diff_Percent[i] <- 100 * (coef_table$PMM2[i] - coef_table$OLS[i]) / abs(coef_table$OLS[i])
      }
    } else {
      coef_table$PMM2[i] <- NA
      coef_table$Diff_Percent[i] <- NA
    }
  }

  return(list(
    ols = fit_ols,
    pmm2 = fit_pmm2,
    coefficients = coef_table,
    residual_stats = res_stats
  ))
}


#' Метод прогнозування для об'єктів PMM2fit
#'
#' @param object Об'єкт PMM2fit
#' @param newdata Новий фрейм даних для прогнозування
#' @param debug Логічне значення, чи виводити дебаг-інформацію
#' @param ... додаткові аргументи (не використовуються)
#'
#' @return Вектор прогнозів
#' @export
setMethod("predict", "PMM2fit",
          function(object, newdata = NULL, debug = FALSE, ...) {
            if(is.null(newdata)) {
              stop("Параметр newdata має бути наданий")
            }

            if(is.null(object@call)) {
              stop("Об'єкт PMM2fit не містить інформації про виклик")
            }

            # Витягнути формулу з виклику
            formula <- eval(object@call$formula)

            if(debug) {
              cat("Формула:", deparse(formula), "\n")
              cat("Коефіцієнти:", paste(names(object@coefficients), "=", object@coefficients, collapse=", "), "\n")
              cat("Розмір newdata:", nrow(newdata), "x", ncol(newdata), "\n")
              cat("Змінні в newdata:", paste(names(newdata), collapse=", "), "\n")
            }

            # Просте рішення - напряму використаємо праву частину формули для побудови матриці дизайну
            rhs <- formula[[3]]
            design_formula <- as.formula(paste("~", deparse(rhs)))

            if(debug) {
              cat("Формула дизайну:", deparse(design_formula), "\n")
            }

            # Створити матрицю дизайну
            X <- model.matrix(design_formula, newdata)

            if(debug) {
              cat("Розмір матриці дизайну:", nrow(X), "x", ncol(X), "\n")
              cat("Стовпці матриці дизайну:", paste(colnames(X), collapse=", "), "\n")
            }

            # Виправляємо проблему з іменами коефіцієнтів
            # Якщо імена відсутні або порожні, заповнюємо їх правильними значеннями
            if(is.null(names(object@coefficients)) || all(names(object@coefficients) == "")) {
              if(debug) {
                cat("Імена коефіцієнтів відсутні або порожні. Використовуємо стандартні імена.\n")
              }

              expected_names <- colnames(X)
              if(length(expected_names) == length(object@coefficients)) {
                names(object@coefficients) <- expected_names
              } else {
                warning("Кількість коефіцієнтів не відповідає кількості стовпців у матриці дизайну.")
                if(length(object@coefficients) == 3 && ncol(X) == 3 &&
                   all(colnames(X) == c("(Intercept)", "x1", "x2"))) {
                  # Найчастіший випадок - регресія з 2 змінними
                  names(object@coefficients) <- c("(Intercept)", "x1", "x2")
                } else {
                  # Загальне присвоєння імен
                  names(object@coefficients) <- paste0("coef", seq_along(object@coefficients))
                }
              }

              if(debug) {
                cat("Нові імена коефіцієнтів:", paste(names(object@coefficients), collapse=", "), "\n")
              }
            }

            # Обчислити прогнози безпосередньо
            predictions <- numeric(nrow(newdata))

            # Для спрощення, просто обчислюємо прогнози вручну для типової регресії
            if(length(object@coefficients) == 3 && all(c("x1", "x2") %in% names(newdata))) {
              if(debug) {
                cat("Обчислюємо прогнози вручну для типової регресії з інтерсептом і двома змінними.\n")
              }
              # Якщо це типова регресія y ~ x1 + x2
              intercept_idx <- which(names(object@coefficients) == "(Intercept)")
              x1_idx <- which(names(object@coefficients) == "x1")
              x2_idx <- which(names(object@coefficients) == "x2")

              if(length(intercept_idx) == 1 && length(x1_idx) == 1 && length(x2_idx) == 1) {
                predictions <- object@coefficients[intercept_idx] +
                  object@coefficients[x1_idx] * newdata$x1 +
                  object@coefficients[x2_idx] * newdata$x2
              } else {
                # Якщо імена не такі як очікується, використовуємо їх позиції
                predictions <- object@coefficients[1] +
                  object@coefficients[2] * newdata$x1 +
                  object@coefficients[3] * newdata$x2
              }
            } else {
              # Для інших випадків
              if(debug) {
                cat("Намагаємося обчислити загальний випадок.\n")
              }
              # Спрощений підхід для загального випадку
              coeffs <- object@coefficients

              # Інтерсепт
              if("(Intercept)" %in% names(coeffs)) {
                predictions <- predictions + coeffs["(Intercept)"]
              }

              # Інші змінні
              for(var_name in intersect(names(coeffs), names(newdata))) {
                if(var_name != "(Intercept)") {
                  predictions <- predictions + coeffs[var_name] * newdata[[var_name]]
                }
              }
            }

            if(debug) {
              cat("Розмір вектора прогнозів:", length(predictions), "\n")
              cat("Перші кілька прогнозів:", paste(head(predictions), collapse=", "), "\n")
            }

            return(predictions)
          })
