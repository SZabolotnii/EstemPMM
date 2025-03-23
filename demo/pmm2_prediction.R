# Частина 3: Порівняння передбачувальної точності PMM2 та OLS
# з використанням розділення даних та крос-валідації

#############################################################
# Функція для оцінки передбачувальної точності з розділенням даних
#############################################################

evaluate_prediction_accuracy <- function() {
  # Завантаження даних Auto MPG з UCI Repository
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

  # Встановлення seed для відтворюваності результатів
  set.seed(42)

  #############################################################
  # 1. Оцінка точності з простим розділенням даних (80/20)
  #############################################################

  cat("\n--- Аналіз із простим розділенням даних (80/20) ---\n")

  # Випадкове перемішування даних
  model_data <- model_data[sample(nrow(model_data)), ]

  # Розділення на навчальну (80%) та тестову (20%) вибірки
  train_size <- floor(0.8 * nrow(model_data))
  train_data <- model_data[1:train_size, ]
  test_data <- model_data[(train_size+1):nrow(model_data), ]

  cat("Розмір навчальної вибірки:", nrow(train_data), "\n")
  cat("Розмір тестової вибірки:", nrow(test_data), "\n")

  # Підгонка моделей на навчальній вибірці
  ols_fit <- lm(y ~ x, data = train_data)
  pmm2_fit <- lm_pmm2(y ~ x, data = train_data)

  # Прогнозування на тестовій вибірці
  ols_pred <- predict(ols_fit, newdata = test_data)
  pmm2_pred <- predict(pmm2_fit, newdata = test_data)

  # Обчислення метрик якості прогнозу
  ols_mse <- mean((test_data$y - ols_pred)^2)
  pmm2_mse <- mean((test_data$y - pmm2_pred)^2)

  ols_mae <- mean(abs(test_data$y - ols_pred))
  pmm2_mae <- mean(abs(test_data$y - pmm2_pred))

  # Коефіцієнт детермінації (R²)
  ols_r2 <- 1 - sum((test_data$y - ols_pred)^2) / sum((test_data$y - mean(test_data$y))^2)
  pmm2_r2 <- 1 - sum((test_data$y - pmm2_pred)^2) / sum((test_data$y - mean(test_data$y))^2)

  # Виведення результатів
  cat("\nРезультати на тестовій вибірці:\n")
  cat("MSE (OLS):", round(ols_mse, 4), "\n")
  cat("MSE (PMM2):", round(pmm2_mse, 4), "\n")
  cat("Відношення MSE (PMM2/OLS):", round(pmm2_mse/ols_mse, 4), "\n\n")

  cat("MAE (OLS):", round(ols_mae, 4), "\n")
  cat("MAE (PMM2):", round(pmm2_mae, 4), "\n")
  cat("Відношення MAE (PMM2/OLS):", round(pmm2_mae/ols_mae, 4), "\n\n")

  cat("R² (OLS):", round(ols_r2, 4), "\n")
  cat("R² (PMM2):", round(pmm2_r2, 4), "\n\n")

  # Налаштування для візуалізації
  test_data$ols_pred <- ols_pred
  test_data$pmm2_pred <- pmm2_pred
  test_data$ols_error <- test_data$y - ols_pred
  test_data$pmm2_error <- test_data$y - pmm2_pred

  # Графік прогнозів на тестовій вибірці
  p1 <- ggplot(test_data, aes(x = x, y = y)) +
    geom_point(alpha = 0.5) +
    geom_point(aes(y = ols_pred, color = "OLS"), shape = 1, size = 2) +
    geom_point(aes(y = pmm2_pred, color = "PMM2"), shape = 2, size = 2) +
    labs(title = "Фактичні та прогнозовані значення (тестова вибірка)",
         x = "Прискорення", y = "MPG") +
    scale_color_manual(values = c("OLS" = "blue", "PMM2" = "red"),
                       name = "Метод") +
    theme_minimal()

  # Порівняння похибок прогнозу
  p2 <- ggplot(test_data) +
    geom_histogram(aes(x = ols_error, fill = "OLS"), alpha = 0.5, bins = 20, position = "identity") +
    geom_histogram(aes(x = pmm2_error, fill = "PMM2"), alpha = 0.5, bins = 20, position = "identity") +
    labs(title = "Розподіл похибок прогнозу (тестова вибірка)",
         x = "Похибка прогнозу", y = "Частота") +
    scale_fill_manual(values = c("OLS" = "blue", "PMM2" = "red"),
                      name = "Метод") +
    theme_minimal()

  # Відображення графіків
  grid.arrange(p1, p2, ncol = 2)

  #############################################################
  # 2. Крос-валідація (k-fold)
  #############################################################

  cat("\n--- Аналіз із k-fold крос-валідацією ---\n")

  # Кількість розбиттів для крос-валідації
  k <- 5

  # Випадкове перемішування даних
  set.seed(42)
  model_data <- model_data[sample(nrow(model_data)), ]

  # Розбиття даних на k підвибірок
  folds <- cut(seq(1, nrow(model_data)), breaks = k, labels = FALSE)

  # Ініціалізація масивів для метрик
  ols_mse_cv <- numeric(k)
  pmm2_mse_cv <- numeric(k)
  ols_mae_cv <- numeric(k)
  pmm2_mae_cv <- numeric(k)
  ols_r2_cv <- numeric(k)
  pmm2_r2_cv <- numeric(k)

  # Виконання k-fold крос-валідації
  for (i in 1:k) {
    # Розділення на навчальну і тестову вибірки
    test_indices <- which(folds == i)
    train_indices <- which(folds != i)

    fold_train <- model_data[train_indices, ]
    fold_test <- model_data[test_indices, ]

    # Підгонка моделей
    ols_fit_cv <- lm(y ~ x, data = fold_train)
    pmm2_fit_cv <- lm_pmm2(y ~ x, data = fold_train)

    # Прогнозування
    ols_pred_cv <- predict(ols_fit_cv, newdata = fold_test)
    pmm2_pred_cv <- predict(pmm2_fit_cv, newdata = fold_test)

    # Обчислення метрик
    ols_mse_cv[i] <- mean((fold_test$y - ols_pred_cv)^2)
    pmm2_mse_cv[i] <- mean((fold_test$y - pmm2_pred_cv)^2)

    ols_mae_cv[i] <- mean(abs(fold_test$y - ols_pred_cv))
    pmm2_mae_cv[i] <- mean(abs(fold_test$y - pmm2_pred_cv))

    # Коефіцієнт детермінації
    ols_r2_cv[i] <- 1 - sum((fold_test$y - ols_pred_cv)^2) / sum((fold_test$y - mean(fold_test$y))^2)
    pmm2_r2_cv[i] <- 1 - sum((fold_test$y - pmm2_pred_cv)^2) / sum((fold_test$y - mean(fold_test$y))^2)
  }

  # Обчислення середніх значень метрик по всім розбиттям
  mean_ols_mse <- mean(ols_mse_cv)
  mean_pmm2_mse <- mean(pmm2_mse_cv)

  mean_ols_mae <- mean(ols_mae_cv)
  mean_pmm2_mae <- mean(pmm2_mae_cv)

  mean_ols_r2 <- mean(ols_r2_cv)
  mean_pmm2_r2 <- mean(pmm2_r2_cv)

  # Виведення результатів крос-валідації
  cat("\nРезультати", k, "-fold крос-валідації:\n")

  cat("Середній MSE (OLS):", round(mean_ols_mse, 4),
      "±", round(sd(ols_mse_cv), 4), "\n")
  cat("Середній MSE (PMM2):", round(mean_pmm2_mse, 4),
      "±", round(sd(pmm2_mse_cv), 4), "\n")
  cat("Відношення середніх MSE (PMM2/OLS):",
      round(mean_pmm2_mse/mean_ols_mse, 4), "\n\n")

  cat("Середній MAE (OLS):", round(mean_ols_mae, 4),
      "±", round(sd(ols_mae_cv), 4), "\n")
  cat("Середній MAE (PMM2):", round(mean_pmm2_mae, 4),
      "±", round(sd(pmm2_mae_cv), 4), "\n")
  cat("Відношення середніх MAE (PMM2/OLS):",
      round(mean_pmm2_mae/mean_ols_mae, 4), "\n\n")

  cat("Середній R² (OLS):", round(mean_ols_r2, 4),
      "±", round(sd(ols_r2_cv), 4), "\n")
  cat("Середній R² (PMM2):", round(mean_pmm2_r2, 4),
      "±", round(sd(pmm2_r2_cv), 4), "\n\n")

  # Створення даних для візуалізації результатів
  cv_metrics <- data.frame(
    Fold = rep(1:k, 3),
    Metric = rep(c("MSE", "MAE", "R²"), each = k),
    OLS = c(ols_mse_cv, ols_mae_cv, ols_r2_cv),
    PMM2 = c(pmm2_mse_cv, pmm2_mae_cv, pmm2_r2_cv),
    Ratio = c(pmm2_mse_cv/ols_mse_cv, pmm2_mae_cv/ols_mae_cv, pmm2_r2_cv/ols_r2_cv)
  )

  # Перетворення для кращої візуалізації
  cv_metrics_long <- reshape2::melt(cv_metrics, id.vars = c("Fold", "Metric"),
                                    variable.name = "Method", value.name = "Value")

  # Фільтрація даних для кожної метрики
  mse_data <- subset(cv_metrics_long, Metric == "MSE" & Method != "Ratio")
  mae_data <- subset(cv_metrics_long, Metric == "MAE" & Method != "Ratio")
  r2_data <- subset(cv_metrics_long, Metric == "R²" & Method != "Ratio")

  # Графіки результатів крос-валідації
  p3 <- ggplot(mse_data, aes(x = Method, y = Value, fill = Method)) +
    geom_boxplot() +
    labs(title = "Розподіл MSE за методами",
         y = "MSE") +
    scale_fill_manual(values = c("OLS" = "blue", "PMM2" = "red")) +
    theme_minimal()

  p4 <- ggplot(mae_data, aes(x = Method, y = Value, fill = Method)) +
    geom_boxplot() +
    labs(title = "Розподіл MAE за методами",
         y = "MAE") +
    scale_fill_manual(values = c("OLS" = "blue", "PMM2" = "red")) +
    theme_minimal()

  p5 <- ggplot(r2_data, aes(x = Method, y = Value, fill = Method)) +
    geom_boxplot() +
    labs(title = "Розподіл R² за методами",
         y = "R²") +
    scale_fill_manual(values = c("OLS" = "blue", "PMM2" = "red")) +
    theme_minimal()

  # Статистики порівняння
  comparison_text <- paste(
    paste("Порівняння продуктивності PMM2 відносно OLS:"),
    paste(""),
    paste("Проста перевірка (80/20 розділення):"),
    paste("  Зменшення MSE:", round((1 - pmm2_mse/ols_mse) * 100, 2), "%"),
    paste("  Зменшення MAE:", round((1 - pmm2_mae/ols_mae) * 100, 2), "%"),
    paste("  Покращення R²:", round((pmm2_r2/ols_r2 - 1) * 100, 2), "%"),
    paste(""),
    paste("Крос-валідація (", k, "-fold):", sep=""),
    paste("  Зменшення MSE:", round((1 - mean_pmm2_mse/mean_ols_mse) * 100, 2), "%"),
    paste("  Зменшення MAE:", round((1 - mean_pmm2_mae/mean_ols_mae) * 100, 2), "%"),
    paste("  Покращення R²:", round((mean_pmm2_r2/mean_ols_r2 - 1) * 100, 2), "%"),
    paste(""),
    paste("Моменти розподілу помилок:"),
    sep = "\n"
  )

  # Обчислення моментів залишків OLS на всьому наборі даних
  full_ols_fit <- lm(y ~ x, data = model_data)
  full_ols_resid <- residuals(full_ols_fit)
  moments <- compute_moments(full_ols_resid)

  comparison_text <- paste(
    comparison_text,
    paste("  Асиметрія (c3) =", round(moments$c3, 4)),
    paste("  Ексцес (c4) =", round(moments$c4, 4)),
    paste("  Теоретичний коефіцієнт g =", round(moments$g, 4)),
    paste(""),
    paste("Висновок:"),
    paste("  PMM2", ifelse(mean_pmm2_mse < mean_ols_mse,
                           "перевершує", "не перевершує"),
          "OLS за передбачувальною точністю"),
    paste("  на даному наборі даних."),
    sep = "\n"
  )

  p6 <- ggplot() +
    annotate("text", x = 0, y = 0.5, label = comparison_text, hjust = 0) +
    theme_void() +
    xlim(0, 1) + ylim(0, 1) +
    labs(title = "Порівняльний аналіз")

  # Відображення результатів
  grid.arrange(p3, p4, p5, p6, ncol = 2)

  # Повернення результатів
  return(list(
    simple_split = list(
      ols_mse = ols_mse,
      pmm2_mse = pmm2_mse,
      ols_mae = ols_mae,
      pmm2_mae = pmm2_mae,
      ols_r2 = ols_r2,
      pmm2_r2 = pmm2_r2
    ),
    cross_validation = list(
      k = k,
      ols_mse = ols_mse_cv,
      pmm2_mse = pmm2_mse_cv,
      ols_mae = ols_mae_cv,
      pmm2_mae = pmm2_mae_cv,
      ols_r2 = ols_r2_cv,
      pmm2_r2 = pmm2_r2_cv,
      mean_ols_mse = mean_ols_mse,
      mean_pmm2_mse = mean_pmm2_mse,
      mean_ols_mae = mean_ols_mae,
      mean_pmm2_mae = mean_pmm2_mae,
      mean_ols_r2 = mean_ols_r2,
      mean_pmm2_r2 = mean_pmm2_r2
    ),
    moments = moments
  ))
}

# Запуск функції для порівняння передбачувальної точності
prediction_results <- evaluate_prediction_accuracy()
