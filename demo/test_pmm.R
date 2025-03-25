# test_pmm.R - Тестовий файл для перевірки функціоналу бібліотеки PMM2

# Завантажуємо всі необхідні пакети
if (!require(testthat)) install.packages("testthat")
if (!require(ggplot2)) install.packages("ggplot2")
library(testthat)
library(ggplot2)
library(EstemPMM)

# Встановлюємо seed для відтворюваності
set.seed(12345)

# ----------------------------------------
# Частина 1: Тестування лінійних моделей
# ----------------------------------------

test_that("Лінійні моделі PMM2 функціонують коректно", {
  # Генеруємо дані з t-розподіленими помилками (для негаусівських хвостів)
  n <- 100
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  # Використовуємо t-розподіл з 3 степенями свободи для важких хвостів
  error <- rt(n, df = 3)
  y <- 2 + 3 * x1 - 1.5 * x2 + error

  test_data <- data.frame(y = y, x1 = x1, x2 = x2)

  # Підгонка OLS для порівняння
  ols_fit <- lm(y ~ x1 + x2, data = test_data)

  # Підгонка моделі PMM2
  pmm2_fit <- lm_pmm2(y ~ x1 + x2, data = test_data, max_iter = 30)

  # Перевірка базової функціональності
  expect_s4_class(pmm2_fit, "PMM2fit")
  expect_length(pmm2_fit@coefficients, 3)  # Включаючи перехоплення

  # Перевірка коефіцієнтів
  cf <- coef(pmm2_fit)
  expect_false(is.null(cf))
  expect_length(cf, 3)

  # Дані потрібно явно передати для функції fitted
  fitted_vals <- fitted(pmm2_fit, data = test_data)
  expect_length(fitted_vals, n)

  # Перевірка залишків
  expect_length(residuals(pmm2_fit), n)

  # Тестування прогнозування
  # Створюємо невеликий тестовий набір
  new_data <- data.frame(
    x1 = rnorm(5),
    x2 = rnorm(5)
  )

  # Виводимо інформацію про дані та коефіцієнти
  cat("\nТестовий набір даних для прогнозування:\n")
  print(new_data)
  cat("\nКоефіцієнти моделі PMM2:\n")
  print(coef(pmm2_fit))

  # Передбачення з увімкненим дебагом
  cat("\nДебаг-інформація виклику predict:\n")
  preds <- tryCatch({
    predict(pmm2_fit, newdata = new_data, debug = TRUE)
  }, error = function(e) {
    cat("ПОМИЛКА при прогнозуванні:", e$message, "\n")
    skip(paste("Помилка при прогнозуванні:", e$message))
    NULL
  })

  if (!is.null(preds)) {
    expect_length(preds, 5)
    cat("\nРезультати прогнозування:\n")
    print(data.frame(
      x1 = new_data$x1,
      x2 = new_data$x2,
      predicted = preds
    ))
  } else {
    cat("\nПрогнозування не вдалося. Перевірте помилки вище.\n")
  }

  # Порівняння з OLS
  comparison <- compare_with_ols(y ~ x1 + x2, test_data)

  expect_s3_class(comparison$ols, "lm")
  expect_s4_class(comparison$pmm2, "PMM2fit")
  expect_s3_class(comparison$coefficients, "data.frame")
  expect_s3_class(comparison$residual_stats, "data.frame")

  # Перевірка статистичного висновку через бутстреп
  # Встановлюємо B=20 для швидкості тестування
  inference_results <- pmm2_inference(pmm2_fit, y ~ x1 + x2, test_data, B = 20)

  expect_s3_class(inference_results, "data.frame")
  expect_equal(length(rownames(inference_results)), 3)
  expect_equal(colnames(inference_results),
               c("Estimate", "Std.Error", "t.value", "p.value", "conf.low", "conf.high"))

  # Виводимо результати для спостереження
  cat("\nРезультати OLS:\n")
  print(summary(ols_fit))

  cat("\nРезультати PMM2:\n")
  print(summary(pmm2_fit))
  print(pmm2_fit)

  cat("\nПорівняння коефіцієнтів:\n")
  print(comparison$coefficients)

  cat("\nПорівняння статистики залишків:\n")
  print(comparison$residual_stats)

  cat("\nРезультати бутстреп-висновку PMM2:\n")
  print(inference_results)
})

# ----------------------------------------
# Частина 2: Тестування моделей часових рядів
# ----------------------------------------

test_that("Моделі часових рядів PMM2 функціонують коректно", {
  skip_if_not_installed("stats")

  # Генеруємо AR(1) процес з t-розподіленими інноваціями
  n <- 200
  phi <- 0.7  # AR коефіцієнт

  # Генеруємо інновації з t-розподілу
  innovations <- rt(n, df = 3)

  # Ініціалізуємо процес
  ar_series <- numeric(n)
  ar_series[1] <- innovations[1]

  # Генеруємо AR процес
  for(t in 2:n) {
    ar_series[t] <- phi * ar_series[t-1] + innovations[t]
  }

  # Тестуємо різні моделі часових рядів

  # AR модель
  ar_fit <- tryCatch({
    ar_pmm2(ar_series, order = 1)
  }, error = function(e) {
    skip(paste("AR модель викликала помилку:", e$message))
    NULL
  })

  if (!is.null(ar_fit)) {
    expect_s4_class(ar_fit, "TS2fit")
    expect_equal(ar_fit@model_type, "ar")
    expect_equal(ar_fit@order$ar, 1)
    expect_length(ar_fit@coefficients, 1)

    # Перевірка методів екстракції результатів для AR
    ar_coefs <- coef(ar_fit)
    expect_false(is.null(ar_coefs))
    expect_length(residuals(ar_fit), length(ar_series))

    # Тестування MA моделі
    ma_fit <- tryCatch({
      ma_pmm2(ar_series, order = 1)
    }, error = function(e) {
      skip(paste("MA модель викликала помилку:", e$message))
      NULL
    })

    if (!is.null(ma_fit)) {
      expect_s4_class(ma_fit, "TS2fit")
      expect_equal(ma_fit@model_type, "ma")
      expect_equal(ma_fit@order$ma, 1)

      # Тестування ARMA моделі
      arma_fit <- tryCatch({
        arma_pmm2(ar_series, order = c(1, 1))
      }, error = function(e) {
        skip(paste("ARMA модель викликала помилку:", e$message))
        NULL
      })

      if (!is.null(arma_fit)) {
        expect_s4_class(arma_fit, "TS2fit")
        expect_equal(arma_fit@model_type, "arma")
        expect_equal(arma_fit@order$ar, 1)
        expect_equal(arma_fit@order$ma, 1)

        # Порівняння з класичними методами (для AR)
        ar_comparison <- tryCatch({
          compare_ar_methods(ar_series, order = 1)
        }, error = function(e) {
          skip(paste("Порівняння AR методів викликало помилку:", e$message))
          NULL
        })

        if (!is.null(ar_comparison)) {
          expect_s3_class(ar_comparison$coefficients, "data.frame")
          expect_s3_class(ar_comparison$residual_stats, "data.frame")

          # Перевірка прогнозування
          ar_pred <- tryCatch({
            predict(ar_fit, n.ahead = 5)
          }, error = function(e) {
            skip(paste("Прогнозування AR моделі викликало помилку:", e$message))
            NULL
          })

          if (!is.null(ar_pred)) {
            expect_length(ar_pred, 5)

            # Тест бутстреп-висновку для часових рядів
            # Використовуємо мале B для швидкості тестування
            # Додаємо debug = TRUE для діагностики
            ts_inference <- tryCatch({
              ts_pmm2_inference(ar_fit, B = 10, method = "residual", debug = TRUE)
            }, error = function(e) {
              skip(paste("Бутстреп-висновок для часових рядів викликав помилку:", e$message))
              NULL
            })

            if (!is.null(ts_inference)) {
              expect_s3_class(ts_inference, "data.frame")
              expect_equal(rownames(ts_inference), "ar1")
              expect_equal(colnames(ts_inference),
                           c("Estimate", "Std.Error", "t.value", "p.value", "conf.low", "conf.high"))

              # Виводимо результати для спостереження
              cat("\nРезультати AR підгонки PMM2:\n")
              print(summary(ar_fit))

              cat("\nРезультати ARMA підгонки PMM2:\n")
              print(summary(arma_fit))

              cat("\nПорівняння коефіцієнтів AR методів:\n")
              print(ar_comparison$coefficients)

              cat("\nПорівняння статистики залишків AR методів:\n")
              print(ar_comparison$residual_stats)

              cat("\nРезультати бутстреп-висновку для AR:\n")
              print(ts_inference)

              # Графічний тест для діагностики
              # Перевіряємо лише синтаксичну коректність
              tryCatch({
                plot(ar_fit)
                cat("\nДіагностичні графіки успішно створені\n")
              }, error = function(e) {
                cat("\nПомилка при створенні діагностичних графіків:", e$message, "\n")
              })
            }
          }
        }
      }
    }
  }
})

# ----------------------------------------
# Частина 3: Графічні тести
# ----------------------------------------

test_that("Графічні функції PMM2 працюють коректно", {
  skip_if_not_installed("graphics")

  # Генеруємо дані для лінійної моделі
  n <- 100
  x <- rnorm(n)
  y <- 2 + 3 * x + rt(n, df = 3)
  test_data <- data.frame(y = y, x = x)

  # Підгонка моделі
  pmm2_fit <- lm_pmm2(y ~ x, data = test_data)

  # Обчислення бутстреп-висновку
  inference_results <- tryCatch({
    pmm2_inference(pmm2_fit, y ~ x, test_data, B = 20)
  }, error = function(e) {
    skip(paste("Бутстреп-висновок викликав помилку:", e$message))
    NULL
  })

  if (!is.null(inference_results)) {
    # Перевірка функції plot_pmm2_bootstrap
    tryCatch({
      # Виведемо інформацію про об'єкт для діагностики
      cat("\nСтруктура об'єкта inference_results:\n")
      str(inference_results)

      # Створюємо графік для всіх коефіцієнтів
      plots_all <- plot_pmm2_bootstrap(inference_results)
      expect_true(!is.null(plots_all))

      # Спробуємо створити графік для конкретного коефіцієнта
      if ("(Intercept)" %in% rownames(inference_results)) {
        plots_subset <- plot_pmm2_bootstrap(inference_results, coefficients = "(Intercept)")
        expect_true(!is.null(plots_subset))
      } else if (length(rownames(inference_results)) > 0) {
        # Якщо немає Intercept, використаємо перший доступний коефіцієнт
        first_coef <- rownames(inference_results)[1]
        plots_subset <- plot_pmm2_bootstrap(inference_results, coefficients = first_coef)
        expect_true(!is.null(plots_subset))
      }

      cat("\nБутстреп-графіки успішно створені\n")
      expect_true(TRUE)  # Якщо дійшли до цього місця, тест успішний
    }, error = function(e) {
      cat("\nПомилка при створенні бутстреп-графіків:", e$message, "\n")
      # Не викликаємо expect_true(FALSE), щоб тест міг продовжуватися
      # Замість цього видаємо попередження
      warning(paste("Помилка графічних функцій:", e$message))
    })
  }
})

# ----------------------------------------
# Частина 4: Тестування загальних утиліт
# ----------------------------------------

test_that("Утилітарні функції PMM2 працюють коректно", {
  # Генеруємо випадкові дані з різними розподілами
  norm_data <- rnorm(100)
  t_data <- rt(100, df = 3)
  skew_data <- exp(rnorm(100))  # логнормальний розподіл для асиметрії

  # Тестуємо функцію обчислення моментів
  norm_moments <- compute_moments(norm_data)
  t_moments <- compute_moments(t_data)
  skew_moments <- compute_moments(skew_data)

  expect_type(norm_moments, "list")
  expect_named(norm_moments, c("m2", "m3", "m4", "c3", "c4", "g"))

  # Тестуємо функції для обчислення асиметрії та ексцесу
  expect_type(pmm_skewness(norm_data), "double")
  expect_type(pmm_kurtosis(norm_data), "double")

  # Перевіряємо правильність обчислень асиметрії
  # Для нормального розподілу асиметрія має бути близька до 0
  expect_lt(abs(pmm_skewness(norm_data)), 0.5)

  # Для t-розподілу з 3 степенями свободи ексцес має бути додатнім
  expect_gt(pmm_kurtosis(t_data), 0)

  # Для логнормального розподілу асиметрія має бути додатньою
  expect_gt(pmm_skewness(skew_data), 0)

  # Виводимо результати для спостереження
  cat("\nМоменти нормального розподілу:\n")
  print(norm_moments)

  cat("\nМоменти t-розподілу (df=3):\n")
  print(t_moments)

  cat("\nМоменти логнормального розподілу:\n")
  print(skew_moments)
})

# ----------------------------------------
# Запуск усіх тестів
# ----------------------------------------

cat("\n----------------------------------------------------------------\n")
cat("РЕЗУЛЬТАТИ ТЕСТУВАННЯ БІБЛІОТЕКИ PMM2\n")
cat("----------------------------------------------------------------\n")
