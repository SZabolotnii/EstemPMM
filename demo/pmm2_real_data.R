# Застосування PMM2 до реальних даних
# Частина 2: Порівняння PMM2 та OLS на наборі даних Auto MPG

# Перевірки наявності пакунків
required_pkgs <- c("EstemPMM", "ggplot2", "gridExtra", "dplyr")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Для цього демо встановіть пакунки: ",
       paste(missing_pkgs, collapse = ", "), call. = FALSE)
}

library(EstemPMM)
library(ggplot2)
library(gridExtra)
library(dplyr)

#############################################################
# Функція для застосування PMM2 до даних Auto MPG
#############################################################

apply_to_mpg_data <- function() {
  # Завантаження даних безпосередньо з UCI Repository
  cat("Завантаження даних Auto MPG з UCI Repository...\n")

  url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/auto-mpg/auto-mpg.data"
  auto_data <- try(read.table(url, sep = "", strip.white = TRUE, na.strings = "?",
                              col.names = c("mpg", "cylinders", "displacement",
                                            "horsepower", "weight", "acceleration",
                                            "model_year", "origin", "car_name")))

  if (inherits(auto_data, "try-error")) {
    # Якщо не вдалося завантажити дані, створимо тестові дані
    cat("Не вдалося завантажити дані з UCI. Створюємо тестові дані...\n")
    set.seed(123)
    n <- 200
    x <- runif(n, 8, 25)  # прискорення
    # Параметри регресії
    a0 <- 5
    a1 <- 1.2
    # Створимо асиметричні помилки
    errors <- rgamma(n, shape = 2, scale = 1) - 2
    y <- a0 + a1 * x + errors
    auto_data <- data.frame(
      mpg = y,
      acceleration = x
    )
  }

  # Попередня обробка
  auto_data <- auto_data[complete.cases(auto_data[c("mpg", "acceleration")]), ]

  # Виберемо змінні для аналізу
  model_data <- data.frame(
    y = auto_data$mpg,
    x = auto_data$acceleration
  )

  cat("Дані обробляються...\n")
  cat("Розмір вибірки:", nrow(model_data), "спостережень\n")

  # Підгонка OLS
  ols_fit <- lm(y ~ x, data = model_data)

  # Підгонка PMM2
  pmm2_fit <- lm_pmm2(y ~ x, data = model_data)

  # Обчислення моментів залишків OLS
  ols_resid <- residuals(ols_fit)
  moments <- EstemPMM::compute_moments(ols_resid)

  # Створення даних для відображення ліній регресії
  x_range <- seq(min(model_data$x), max(model_data$x), length.out = 100)
  # Безпосереднє обчислення значень для ліній регресії
  ols_pred <- coef(ols_fit)[1] + coef(ols_fit)[2] * x_range
  pmm2_pred <- pmm2_fit@coefficients[1] + pmm2_fit@coefficients[2] * x_range

  grid_data <- data.frame(x = x_range, ols_pred = ols_pred, pmm2_pred = pmm2_pred)

  # Графік з даними та лініями регресії
  p1 <- ggplot(model_data, aes(x = x, y = y)) +
    geom_point(alpha = 0.5) +
    geom_line(data = grid_data, aes(x = x, y = ols_pred, color = "OLS"), size = 1) +
    geom_line(data = grid_data, aes(x = x, y = pmm2_pred, color = "PMM2"), size = 1) +
    labs(title = "Залежність MPG від прискорення",
         x = "Прискорення", y = "MPG") +
    scale_color_manual(values = c("OLS" = "blue", "PMM2" = "red"),
                       name = "Метод") +
    theme_minimal()

  # Гістограма залишків OLS
  p2 <- ggplot(data.frame(resid = ols_resid), aes(x = resid)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "skyblue", alpha = 0.7) +
    geom_density(color = "blue", size = 1) +
    labs(title = "Розподіл залишків OLS",
         x = "Залишки", y = "Густина") +
    theme_minimal()

  # QQ-графік залишків
  p3 <- ggplot(data.frame(resid = ols_resid), aes(sample = resid)) +
    geom_qq() +
    geom_qq_line() +
    labs(title = "QQ-графік залишків OLS",
         x = "Теоретичні квантилі", y = "Вибіркові квантилі") +
    theme_minimal()

  # Статистики
  stats_text <- paste(
    paste("OLS коефіцієнти:"),
    paste("  a0 =", round(coef(ols_fit)[1], 4)),
    paste("  a1 =", round(coef(ols_fit)[2], 4)),
    paste("PMM2 коефіцієнти:"),
    paste("  a0 =", round(pmm2_fit@coefficients[1], 4)),
    paste("  a1 =", round(pmm2_fit@coefficients[2], 4)),
    paste(""),
    paste("Статистика помилок:"),
    paste("  Асиметрія (c3) =", round(moments$c3, 4)),
    paste("  Ексцес (c4) =", round(moments$c4, 4)),
    paste("  Теоретичний коефіцієнт g =", round(moments$g, 4)),
    paste(""),
    paste("Критерії оцінки:"),
    paste("  OLS MSE =", round(mean(ols_resid^2), 4)),
    paste("  PMM2 MSE =", round(mean(pmm2_fit@residuals^2), 4)),
    paste("  Відношення MSE =", round(mean(pmm2_fit@residuals^2) / mean(ols_resid^2), 4)),
    paste("  OLS AIC =", round(AIC(ols_fit), 2)),
    paste("  PMM2 AIC =", round(AIC(pmm2_fit), 2)),
    sep = "\n"
  )

  p4 <- ggplot() +
    annotate("text", x = 0, y = 0.5, label = stats_text, hjust = 0) +
    theme_void() +
    xlim(0, 1) + ylim(0, 1) +
    labs(title = "Статистика для Auto MPG даних")

  # Об'єднання графіків - основні результати
  grid.arrange(p1, p2, p3, p4, ncol = 2)

  # Виконання бутстрепного аналізу для оцінки невизначеності
  bootstrap_pmm2 <- function(data, B = 1000) {
    n <- nrow(data)

    # Підгонка початкової моделі
    ols_fit <- lm(y ~ x, data = data)
    pmm2_fit <- lm_pmm2(y ~ x, data = data)

    # Ініціалізація масивів для зберігання результатів
    ols_coefs <- matrix(0, nrow = B, ncol = 2)
    pmm2_coefs <- matrix(0, nrow = B, ncol = 2)

    # Ініціалізація прогрес-бару
    pb <- txtProgressBar(min = 0, max = B, style = 3)

    for (i in 1:B) {
      # Генерація бутстрепної вибірки
      indices <- sample(1:n, n, replace = TRUE)
      boot_data <- data[indices, ]

      # Підгонка моделей
      boot_ols <- lm(y ~ x, data = boot_data)
      boot_pmm2 <- lm_pmm2(y ~ x, data = boot_data)

      # Збереження коефіцієнтів
      ols_coefs[i, ] <- coef(boot_ols)
      pmm2_coefs[i, ] <- boot_pmm2@coefficients

      setTxtProgressBar(pb, i)
    }
    close(pb)

    colnames(ols_coefs) <- c("a0", "a1")
    colnames(pmm2_coefs) <- c("a0", "a1")

    return(list(
      ols = ols_coefs,
      pmm2 = pmm2_coefs
    ))
  }

  # Виконання бутстрепу
  cat("\nВиконується бутстрепний аналіз...\n")
  boot_results <- bootstrap_pmm2(model_data, B = 500)

  cat("\nБутстрепний аналіз завершено!\n")

  # Візуалізація результатів бутстрепу
  boot_data <- data.frame(
    a0_ols = boot_results$ols[, "a0"],
    a0_pmm2 = boot_results$pmm2[, "a0"],
    a1_ols = boot_results$ols[, "a1"],
    a1_pmm2 = boot_results$pmm2[, "a1"]
  )

  p5 <- ggplot(boot_data) +
    geom_density(aes(x = a0_ols, fill = "OLS"), alpha = 0.5) +
    geom_density(aes(x = a0_pmm2, fill = "PMM2"), alpha = 0.5) +
    labs(title = "Бутстрепний розподіл оцінок a0",
         x = "a0", y = "Густина") +
    scale_fill_manual(values = c("OLS" = "blue", "PMM2" = "red"),
                      name = "Метод") +
    theme_minimal()

  p6 <- ggplot(boot_data) +
    geom_density(aes(x = a1_ols, fill = "OLS"), alpha = 0.5) +
    geom_density(aes(x = a1_pmm2, fill = "PMM2"), alpha = 0.5) +
    labs(title = "Бутстрепний розподіл оцінок a1",
         x = "a1", y = "Густина") +
    scale_fill_manual(values = c("OLS" = "blue", "PMM2" = "red"),
                      name = "Метод") +
    theme_minimal()

  # Статистики бутстрепу
  boot_stats_text <- paste(
    paste("Бутстрепні оцінки (95% CI):"),
    paste("OLS a0:", round(mean(boot_data$a0_ols), 4),
          "[", round(quantile(boot_data$a0_ols, 0.025), 4), ",",
          round(quantile(boot_data$a0_ols, 0.975), 4), "]"),
    paste("PMM2 a0:", round(mean(boot_data$a0_pmm2), 4),
          "[", round(quantile(boot_data$a0_pmm2, 0.025), 4), ",",
          round(quantile(boot_data$a0_pmm2, 0.975), 4), "]"),
    paste("OLS a1:", round(mean(boot_data$a1_ols), 4),
          "[", round(quantile(boot_data$a1_ols, 0.025), 4), ",",
          round(quantile(boot_data$a1_ols, 0.975), 4), "]"),
    paste("PMM2 a1:", round(mean(boot_data$a1_pmm2), 4),
          "[", round(quantile(boot_data$a1_pmm2, 0.025), 4), ",",
          round(quantile(boot_data$a1_pmm2, 0.975), 4), "]"),
    paste(""),
    paste("Стандартні відхилення:"),
    paste("SD OLS a0:", round(sd(boot_data$a0_ols), 4)),
    paste("SD PMM2 a0:", round(sd(boot_data$a0_pmm2), 4)),
    paste("SD OLS a1:", round(sd(boot_data$a1_ols), 4)),
    paste("SD PMM2 a1:", round(sd(boot_data$a1_pmm2), 4)),
    paste(""),
    paste("Відношення SD:"),
    paste("SD PMM2 a0 / SD OLS a0:", round(sd(boot_data$a0_pmm2) / sd(boot_data$a0_ols), 4)),
    paste("SD PMM2 a1 / SD OLS a1:", round(sd(boot_data$a1_pmm2) / sd(boot_data$a1_ols), 4)),
    paste(""),
    paste("Емпіричний коефіцієнт g:"),
    paste("g для a0:", round((sd(boot_data$a0_pmm2) / sd(boot_data$a0_ols))^2, 4)),
    paste("g для a1:", round((sd(boot_data$a1_pmm2) / sd(boot_data$a1_ols))^2, 4)),
    sep = "\n"
  )

  p7 <- ggplot() +
    annotate("text", x = 0, y = 0.5, label = boot_stats_text, hjust = 0) +
    theme_void() +
    xlim(0, 1) + ylim(0, 1) +
    labs(title = "Статистика бутстрепу")

  # Відображення результатів бутстрепу
  grid.arrange(p5, p6, p7, ncol = 2)

  # Повернення результатів
  return(list(
    model_data = model_data,
    ols_fit = ols_fit,
    pmm2_fit = pmm2_fit,
    moments = moments,
    boot_results = boot_results
  ))
}

# Щоб запустити аналіз реальних даних, викличте вручну:
# apply_to_mpg_data()
