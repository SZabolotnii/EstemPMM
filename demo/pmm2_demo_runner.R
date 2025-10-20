# Перевірка наявності необхідних пакунків
required_pkgs <- c("EstemPMM", "ggplot2", "gridExtra", "dplyr", "parallel", "reshape2")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Перед запуском демо встановіть пакунки: ",
       paste(missing_pkgs, collapse = ", "), call. = FALSE)
}

# Підключення бібліотек після перевірки
library(EstemPMM)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(parallel)
library(reshape2)

# Допоміжна функція підвантаження інших демо-скриптів
source_demo <- function(fname) {
  local_path <- file.path("demo", fname)
  if (file.exists(local_path)) {
    source(local_path, local = FALSE)
    return(invisible(TRUE))
  }
  pkg_path <- system.file("demo", fname, package = "EstemPMM")
  if (nzchar(pkg_path)) {
    source(pkg_path, local = FALSE)
    return(invisible(TRUE))
  }
  warning("Не вдалося знайти демо-файл ", fname)
  invisible(FALSE)
}

# Підвантажуємо функції з інших демонстрацій (без автозапуску)
source_demo("pmm2_simulation.R")
source_demo("pmm2_real_data.R")
source_demo("pmm2_prediction.R")
source_demo("pmm2_simMC_ts.R")

# Швидкий тест: маленька вибірка + AR приклад
quick_test <- function() {
  cat("\n=== Швидкий тест PMM2 ===\n")
  set.seed(123)
  n <- 120
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  errors <- rgamma(n, shape = 2, scale = 1) - 2
  y <- 1.2 + 1.5 * x1 - 0.8 * x2 + errors
  dat <- data.frame(y = y, x1 = x1, x2 = x2)

  ols_fit <- lm(y ~ x1 + x2, data = dat)
  pmm2_fit <- lm_pmm2(y ~ x1 + x2, data = dat)

  cat("\nCoefficients (OLS vs PMM2):\n")
  print(rbind(OLS = coef(ols_fit), PMM2 = coef(pmm2_fit)))

  new_data <- data.frame(x1 = c(-1, 0, 1), x2 = c(0.5, -0.2, 1.1))
  cat("\nPMM2 predictions for new points:\n")
  print(predict(pmm2_fit, newdata = new_data))

  ts_innov <- rgamma(n, shape = 2, scale = 1) - 2
  ts_series <- stats::filter(ts_innov, filter = 0.6, method = "recursive")
  ar_fit <- ar_pmm2(ts_series, order = 1)
  cat("\nAR(1) PMM2 coefficient:\n")
  print(coef(ar_fit))
  cat("\nШвидкий тест завершено.\n")
}

# Основна функція для запуску демонстрації
run_demo <- function(quick = TRUE) {
  if (quick) {
    quick_test()
  } else {
    cat("Виберіть демонстрацію для запуску:\n")
    cat("1. Симуляції Монте-Карло (займає багато часу)\n")
    cat("2. Аналіз реальних даних Auto MPG\n")
    cat("3. Порівняння передбачувальної точності з крос-валідацією\n")
    cat("4. Швидке тестування\n")
    cat("0. Вихід\n\n")

    choice <- as.integer(readline("Ваш вибір (0-4): "))

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
             cat("Запуск аналізу передбачувальної точності...\n")
             prediction_results <- evaluate_prediction_accuracy()
             cat("Аналіз завершено!\n")
           },
           {
             quick_test()
           },
           cat("Некоректний вибір\n")
    )
  }
}

# Запуск демонстрації (за бажанням користувача)
if (interactive()) {
  message("Запуск демо у швидкому режимі. Викличте run_demo(quick = FALSE) для повного меню.")
  run_demo(quick = TRUE)
}
