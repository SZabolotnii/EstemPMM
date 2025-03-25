# pmm_ts_examples.R

#' Приклади застосування методу PMM2 для аналізу часових рядів
#'
#' Цей скрипт демонструє використання функцій з пакету EstemPMM
#' для моделювання часових рядів з негаусівськими розподілами
#'
#' Модифікована версія: для всіх моделей (AR, MA, ARMA, ARIMA) інновації
#' генеруються як:
#'   inn <- (rgamma(n = n, shape = shape_gamma) - shape_gamma) / sqrt(shape_gamma)
#' що створює центровану, асиметрично розподілену випадкову величину.

# Завантаження пакетів
library(EstemPMM)
library(moments)
library(stats)

#' Демонстрація PMM2 для аналізу моделей AR, MA, ARMA та ARIMA
#' з загальною моделлю інновацій на основі гамма-розподілу
run_ts_examples <- function() {
  # Встановлення seed для відтворюваності
  set.seed(123)

  # Визначення параметру гамма-розподілу для інновацій
  shape_gamma <- 2

  cat("\n============ PMM2 ДЕМОНСТРАЦІЯ ДЛЯ ЧАСОВИХ РЯДІВ ============\n\n")

  # 1. AR(2) модель з гамма-розподіленими інноваціями
  cat("1. AR(2) модель з гамма-розподіленими інноваціями\n")
  ar_coef <- c(0.7, -0.3)
  n <- 100

  tryCatch({
    ar_series <- as.numeric(arima.sim(model = list(ar = ar_coef), n = n,
                                      rand.gen = function(n) (rgamma(n, shape = shape_gamma) - shape_gamma) / sqrt(shape_gamma)))

    # Порівняння методів оцінки AR моделі
    ar_comparison <- compare_ar_methods(ar_series, order = 2)

    # Вивід результатів
    cat("\nПорівняння оцінених коефіцієнтів AR(2) моделі:\n")
    print(ar_comparison$coefficients)

    cat("\nСтатистика залишків:\n")
    print(ar_comparison$residual_stats)

    # Створення графіка
    par(mfrow = c(2,2))
    plot(ar_comparison$pmm2)

  }, error = function(e) {
    cat("\nПомилка в аналізі AR моделі:", conditionMessage(e), "\n")
  })

  # 2. MA(2) модель з гамма-розподіленими інноваціями
  cat("\n\n2. MA(2) модель з гамма-розподіленими інноваціями\n")

  tryCatch({
    ma_coef <- c(0.6, -0.2)
    # Генерування інновацій та серії MA(2)
    ma_innov <- (rgamma(n, shape = shape_gamma) - shape_gamma) / sqrt(shape_gamma)
    ma_series <- as.numeric(arima.sim(model = list(ma = ma_coef), n = n,
                                      innov = ma_innov))

    # Оцінювання MA(2) класичними методами
    css_fit <- tryCatch({
      arima(ma_series, order = c(0,0,2), method = "CSS")
    }, error = function(e) {
      cat("\nПомилка CSS методу для MA(2):", conditionMessage(e), "\n")
      NULL
    })

    ml_fit <- tryCatch({
      arima(ma_series, order = c(0,0,2), method = "ML")
    }, error = function(e) {
      cat("\nПомилка ML методу для MA(2):", conditionMessage(e), "\n")
      NULL
    })

    # Якщо хоча б один класичний метод працює, виведемо результати
    if(!is.null(css_fit) || !is.null(ml_fit)) {
      # Підготовка таблиці коефіцієнтів
      coef_table <- data.frame(
        Coefficient = c("ma1", "ma2"),
        True = c(ma_coef[1], ma_coef[2])
      )

      if(!is.null(css_fit)) {
        coef_table$CSS <- css_fit$coef[c("ma1", "ma2")]
      }

      if(!is.null(ml_fit)) {
        coef_table$ML <- ml_fit$coef[c("ma1", "ma2")]
      }

      # Вивід результатів
      cat("\nКласичні методи оцінки MA(2) моделі:\n")
      print(coef_table)

      # Статистика залишків для класичних методів
      if(!is.null(css_fit) && !is.null(ml_fit)) {
        res_css <- residuals(css_fit)
        res_ml <- residuals(ml_fit)

        res_table <- data.frame(
          Method = c("CSS", "ML"),
          RSS = c(sum(res_css^2), sum(res_ml^2)),
          MAE = c(mean(abs(res_css)), mean(abs(res_ml))),
          Skewness = c(skewness(res_css), skewness(res_ml)),
          Kurtosis = c(kurtosis(res_css)-3, kurtosis(res_ml)-3)
        )

        cat("\nСтатистика залишків для класичних методів:\n")
        print(res_table)
      }
    }

    # Тепер спробуємо PMM2 метод
    cat("\nСпроба використання PMM2 для MA(2) моделі...\n")
    pmm2_fit <- tryCatch({
      ts_pmm2(ma_series, order = 2, model_type = "ma", method = "pmm2")
    }, error = function(e) {
      cat("PMM2 помилка:", conditionMessage(e), "\n")
      NULL
    })

    if(!is.null(pmm2_fit)) {
      cat("\nPMM2 оцінки для MA(2):\n")
      print(pmm2_fit@coefficients)

      # Графік для PMM2
      par(mfrow = c(2,2))
      plot(pmm2_fit)
    } else {
      cat("Неможливо оцінити MA(2) модель за допомогою PMM2\n")
    }

  }, error = function(e) {
    cat("\nЗагальна помилка в аналізі MA моделі:", conditionMessage(e), "\n")
  })

  # 3. ARMA(1,1) модель з гамма-розподіленими інноваціями
  cat("\n\n3. ARMA(1,1) модель з гамма-розподіленими інноваціями\n")

  tryCatch({
    # Визначення параметрів ARMA(1,1)
    arma_ar <- 0.7
    arma_ma <- 0.4

    # Генерування інновацій та ARMA серії
    arma_innov <- (rgamma(n, shape = shape_gamma) - shape_gamma) / sqrt(shape_gamma)
    arma_series <- as.numeric(arima.sim(model = list(ar = arma_ar, ma = arma_ma),
                                        n = n, innov = arma_innov))

    # Оцінювання ARMA(1,1) класичними методами
    css_fit <- tryCatch({
      arima(arma_series, order = c(1,0,1), method = "CSS")
    }, error = function(e) {
      cat("\nПомилка CSS методу для ARMA(1,1):", conditionMessage(e), "\n")
      NULL
    })

    ml_fit <- tryCatch({
      arima(arma_series, order = c(1,0,1), method = "ML")
    }, error = function(e) {
      cat("\nПомилка ML методу для ARMA(1,1):", conditionMessage(e), "\n")
      NULL
    })

    # Якщо хоча б один класичний метод працює, виведемо результати
    if(!is.null(css_fit) || !is.null(ml_fit)) {
      # Підготовка таблиці коефіцієнтів
      coef_table <- data.frame(
        Coefficient = c("ar1", "ma1"),
        True = c(arma_ar, arma_ma)
      )

      if(!is.null(css_fit)) {
        coef_table$CSS <- css_fit$coef[c("ar1", "ma1")]
      }

      if(!is.null(ml_fit)) {
        coef_table$ML <- ml_fit$coef[c("ar1", "ma1")]
      }

      # Вивід результатів
      cat("\nКласичні методи оцінки ARMA(1,1) моделі:\n")
      print(coef_table)

      # Статистика залишків для класичних методів
      if(!is.null(css_fit) && !is.null(ml_fit)) {
        res_css <- residuals(css_fit)
        res_ml <- residuals(ml_fit)

        res_table <- data.frame(
          Method = c("CSS", "ML"),
          RSS = c(sum(res_css^2), sum(res_ml^2)),
          MAE = c(mean(abs(res_css)), mean(abs(res_ml))),
          Skewness = c(skewness(res_css), skewness(res_ml)),
          Kurtosis = c(kurtosis(res_css)-3, kurtosis(res_ml)-3)
        )

        cat("\nСтатистика залишків для класичних методів:\n")
        print(res_table)
      }
    }

    # Тепер спробуємо PMM2 метод
    cat("\nСпроба використання PMM2 для ARMA(1,1) моделі...\n")
    pmm2_fit <- tryCatch({
      ts_pmm2(arma_series, order = c(1,1), model_type = "arma", method = "pmm2")
    }, error = function(e) {
      cat("PMM2 помилка:", conditionMessage(e), "\n")
      NULL
    })

    if(!is.null(pmm2_fit)) {
      cat("\nPMM2 оцінки для ARMA(1,1):\n")
      print(pmm2_fit@coefficients)

      # Графік для PMM2
      par(mfrow = c(2,2))
      plot(pmm2_fit)
    } else {
      cat("Неможливо оцінити ARMA(1,1) модель за допомогою PMM2\n")
    }

  }, error = function(e) {
    cat("\nЗагальна помилка в аналізі ARMA моделі:", conditionMessage(e), "\n")
  })

  # 4. ARIMA(1,1,1) модель з гамма-розподіленими інноваціями
  cat("\n\n4. ARIMA(1,1,1) модель з гамма-розподіленими інноваціями\n")

  tryCatch({
    # Параметри ARIMA моделі
    arima_ar <- 0.7
    arima_ma <- 0.4
    n_arima <- 400  # Більша вибірка для ARIMA

    # Генерування інновацій
    arima_innov <- (rgamma(n_arima, shape = shape_gamma) - shape_gamma) / sqrt(shape_gamma)

    # Генерування ARMA серії
    arma_base <- as.numeric(arima.sim(model = list(ar = arima_ar, ma = arima_ma),
                                      n = n_arima, innov = arima_innov))

    # Створення нестаціонарного ряду (інтеграція)
    arima_series <- cumsum(arma_base)

    # Оцінювання ARIMA(1,1,1) класичними методами
    css_fit <- tryCatch({
      arima(arima_series, order = c(1,1,1), method = "CSS")
    }, error = function(e) {
      cat("\nПомилка CSS методу для ARIMA(1,1,1):", conditionMessage(e), "\n")
      NULL
    })

    ml_fit <- tryCatch({
      arima(arima_series, order = c(1,1,1), method = "ML")
    }, error = function(e) {
      cat("\nПомилка ML методу для ARIMA(1,1,1):", conditionMessage(e), "\n")
      NULL
    })

    # Якщо хоча б один класичний метод працює, виведемо результати
    if(!is.null(css_fit) || !is.null(ml_fit)) {
      # Підготовка таблиці коефіцієнтів
      coef_table <- data.frame(
        Coefficient = c("ar1", "ma1"),
        True = c(arima_ar, arima_ma)
      )

      if(!is.null(css_fit)) {
        coef_table$CSS <- css_fit$coef[c("ar1", "ma1")]
      }

      if(!is.null(ml_fit)) {
        coef_table$ML <- ml_fit$coef[c("ar1", "ma1")]
      }

      # Вивід результатів
      cat("\nКласичні методи оцінки ARIMA(1,1,1) моделі:\n")
      print(coef_table)

      # Статистика залишків для класичних методів
      if(!is.null(css_fit) && !is.null(ml_fit)) {
        res_css <- residuals(css_fit)
        res_ml <- residuals(ml_fit)

        res_table <- data.frame(
          Method = c("CSS", "ML"),
          RSS = c(sum(res_css^2), sum(res_ml^2)),
          MAE = c(mean(abs(res_css)), mean(abs(res_ml))),
          Skewness = c(skewness(res_css), skewness(res_ml)),
          Kurtosis = c(kurtosis(res_css)-3, kurtosis(res_ml)-3)
        )

        cat("\nСтатистика залишків для класичних методів:\n")
        print(res_table)
      }
    }

    # Тепер спробуємо PMM2 метод
    cat("\nСпроба використання PMM2 для ARIMA(1,1,1) моделі...\n")
    pmm2_fit <- tryCatch({
      ts_pmm2(arima_series, order = c(1,1,1), model_type = "arima", method = "pmm2")
    }, error = function(e) {
      cat("PMM2 помилка:", conditionMessage(e), "\n")
      NULL
    })

    if(!is.null(pmm2_fit)) {
      cat("\nPMM2 оцінки для ARIMA(1,1,1):\n")
      print(pmm2_fit@coefficients)

      # Графік для PMM2
      par(mfrow = c(3,2))
      plot(pmm2_fit, which = 1:6)
    } else {
      cat("Неможливо оцінити ARIMA(1,1,1) модель за допомогою PMM2\n")
    }

  }, error = function(e) {
    cat("\nЗагальна помилка в аналізі ARIMA моделі:", conditionMessage(e), "\n")
  })

  cat("\n============ КІНЕЦЬ ДЕМОНСТРАЦІЇ ============\n")
}

# Запуск демонстрації
run_ts_examples()
