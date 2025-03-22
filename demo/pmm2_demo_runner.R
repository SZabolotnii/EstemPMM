# Скрипт для запуску демонстрації PMM2
# Завантажує та запускає файли для демонстрації ефективності PMM2

# Завантаження необхідних пакетів
library(ggplot2)
library(gridExtra)
library(dplyr)
library(parallel)

# Визначення шляхів до файлів
pkg_path <- file.path("~", "R", "EstemPMM")
r_path <- file.path(pkg_path, "R")
demo_path <- file.path(pkg_path, "demo")

# Завантаження необхідних пакетів
library(ggplot2)
library(gridExtra)
library(dplyr)
library(parallel)

# Перевірка наявності всіх необхідних файлів
r_files <- c(
  file.path(r_path, "pmm_classes.R"),
  file.path(r_path, "pmm_utils.R"),
  file.path(r_path, "pmm_main.R"),
  file.path(r_path, "pmm_inference.R")
)

demo_files <- c(
  file.path(demo_path, "pmm2_simulation.R"),
  file.path(demo_path, "pmm2_real_data.R")
)

all_files <- c(r_files, demo_files)
missing_files <- all_files[!file.exists(all_files)]

if (length(missing_files) > 0) {
  stop("Відсутні файли: ", paste(missing_files, collapse = ", "))
}

# Завантаження необхідних файлів пакету EstemPMM
for (file in r_files) {
  source(file)
}

# Завантаження файлів демонстрації
for (file in demo_files) {
  source(file)
}

# Функція для швидкого тестування (з меншим обсягом обчислень)
quick_test <- function() {
  cat("Запуск швидкого тестування PMM2...\n\n")

  # Генерація тестових даних з гамма-розподілом похибок
  set.seed(42)
  n <- 100
  x <- rnorm(n, mean = 0, sd = 1)

  # Гамма-розподіл зі зсувом для нульової середньої
  alpha <- 2
  beta <- 1/sqrt(alpha)
  errors <- rgamma(n, shape = alpha, scale = beta) - alpha*beta

  # Істинні значення коефіцієнтів
  true_a0 <- 2
  true_a1 <- 1.5

  # Генерація даних
  y <- true_a0 + true_a1 * x + errors
  test_data <- data.frame(x = x, y = y, errors = errors)

  # Обчислення моментів
  moments <- compute_moments(errors)

  cat("Характеристики похибок:\n")
  cat("  Асиметрія (c3):", round(moments$c3, 4), "\n")
  cat("  Ексцес (c4):", round(moments$c4, 4), "\n")
  cat("  Теоретичний коефіцієнт g:", round(moments$g, 4), "\n\n")

  # Підгонка моделей
  ols_fit <- lm(y ~ x, data = test_data)
  pmm2_fit <- lm_pmm2(y ~ x, data = test_data)

  # Виведення результатів
  cat("Результати підгонки моделей:\n")
  cat("OLS коефіцієнти:  a0 =", round(coef(ols_fit)[1], 4),
      "a1 =", round(coef(ols_fit)[2], 4), "\n")
  cat("PMM2 коефіцієнти: a0 =", round(pmm2_fit@coefficients[1], 4),
      "a1 =", round(pmm2_fit@coefficients[2], 4), "\n\n")

  # MSE
  ols_mse <- mean(residuals(ols_fit)^2)
  pmm2_mse <- mean(pmm2_fit@residuals^2)
  mse_ratio <- pmm2_mse / ols_mse

  cat("Порівняння якості підгонки:\n")
  cat("  OLS MSE:", round(ols_mse, 4), "\n")
  cat("  PMM2 MSE:", round(pmm2_mse, 4), "\n")
  cat("  Відношення MSE:", round(mse_ratio, 4), "\n")
  cat("  Покращення ефективності:", round((1 - mse_ratio) * 100, 2), "%\n\n")

  # Візуалізація результатів
  grid_data <- data.frame(x = seq(min(x), max(x), length.out = 100))
  grid_data$ols_pred <- predict(ols_fit, newdata = grid_data)
  grid_data$pmm2_pred <- predict(pmm2_fit, newdata = grid_data)

  p <- ggplot(test_data, aes(x = x, y = y)) +
    geom_point(alpha = 0.5) +
    geom_line(data = grid_data, aes(x = x, y = ols_pred, color = "OLS"), size = 1) +
    geom_line(data = grid_data, aes(x = x, y = pmm2_pred, color = "PMM2"), size = 1) +
    labs(title = "Порівняння OLS та PMM2 на тестових даних",
         x = "x", y = "y") +
    scale_color_manual(values = c("OLS" = "blue", "PMM2" = "red"),
                       name = "Метод") +
    theme_minimal()

  print(p)

  cat("Для виконання повної симуляції Монте-Карло використовуйте:\n")
  cat("  source('pmm2_simulation.R')\n")
  cat("  all_results <- run_all_simulations()\n\n")

  cat("Для аналізу реальних даних використовуйте:\n")
  cat("  source('pmm2_real_data.R')\n")
  cat("  results <- apply_to_mpg_data()\n")
}

# Основна функція для запуску демонстрації
run_demo <- function(quick = TRUE) {
  if (quick) {
    quick_test()
  } else {
    cat("Виберіть демонстрацію для запуску:\n")
    cat("1. Симуляції Монте-Карло (займає багато часу)\n")
    cat("2. Аналіз реальних даних Auto MPG\n")
    cat("3. Швидке тестування\n")
    cat("0. Вихід\n\n")

    choice <- as.integer(readline("Ваш вибір (0-3): "))

    switch(choice + 1,
           cat("Вихід з демонстрації\n"),
           {
             cat("Запуск симуляцій Монте-Карло...\n")
             all_results <- run_all_simulations()
             cat("Симуляції завершено!\n")
           },
           {
             cat("Запуск аналізу реальних даних...\n")
             results <- apply_to_mpg_data()
             cat("Аналіз завершено!\n")
           },
           {
             quick_test()
           },
           cat("Некоректний вибір\n")
    )
  }
}

# Запуск демонстрації
# За замовчуванням - швидке тестування
run_demo(quick = TRUE)

# Для інтерактивного вибору демонстрації, розкоментуйте:
# run_demo(quick = FALSE)
