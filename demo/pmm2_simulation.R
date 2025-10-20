# Симуляції Монте-Карло для оцінки ефективності PMM2
# Частина 1: Моделювання та порівняння методів на різних розподілах похибок

# Перевірка наявності та підключення необхідних пакунків
required_pkgs <- c("EstemPMM", "ggplot2", "gridExtra", "dplyr", "parallel")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Для цього демо встановіть пакунки: ",
       paste(missing_pkgs, collapse = ", "), call. = FALSE)
}

library(EstemPMM)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(parallel)

#############################################################
# Допоміжні функції для симуляцій Монте-Карло
#############################################################

# Функція для генерації даних з різними розподілами помилок
generate_data <- function(n, distribution, a0, a1, ...) {
  x <- rnorm(n, mean = 0, sd = 1)

  # Генерація помилок з заданим розподілом
  errors <- switch(distribution,
                   "gaussian" = rnorm(n, mean = 0, sd = 1),
                   "t" = rt(n, df = 3),
                   "gamma" = {
                     # Параметризована гамма для нульової середньої
                     alpha <- 2
                     beta <- 1/sqrt(alpha)
                     rgamma(n, shape = alpha, scale = beta) - alpha*beta
                   },
                   "exponential" = {
                     # Експоненційна зі зсувом для нульової середньої
                     lambda <- 1
                     rexp(n, rate = lambda) - 1/lambda
                   },
                   "chi-squared" = {
                     # Хі-квадрат зі зсувом для нульової середньої
                     df <- 3
                     rchisq(n, df = df) - df
                   },
                   "lognormal" = {
                     # Логнормальна зі зсувом для нульової середньої
                     sigma <- 0.5
                     exp(rnorm(n, mean = -sigma^2/2, sd = sigma)) - 1
                   })

  # Обчислення значень y
  y <- a0 + a1 * x + errors

  return(data.frame(x = x, y = y, errors = errors))
}

# Функція для порівняння методів PMM2 та OLS
compare_methods <- function(data, true_a0, true_a1) {
  # Підгонка OLS
  ols_fit <- lm(y ~ x, data = data)

  # Підгонка PMM2
  pmm2_fit <- lm_pmm2(y ~ x, data = data, verbose = FALSE)

  # Обчислення залишків
  ols_resid <- residuals(ols_fit)
  pmm2_resid <- pmm2_fit@residuals

  # Обчислення MSE
  ols_mse <- mean(ols_resid^2)
  pmm2_mse <- mean(pmm2_resid^2)

  # Обчислення AIC
  ols_aic <- AIC(ols_fit)
  pmm2_aic <- AIC(pmm2_fit)

  # Обчислення зміщення оцінок
  ols_bias_a0 <- coef(ols_fit)[1] - true_a0
  ols_bias_a1 <- coef(ols_fit)[2] - true_a1

  pmm2_bias_a0 <- pmm2_fit@coefficients[1] - true_a0
  pmm2_bias_a1 <- pmm2_fit@coefficients[2] - true_a1

  # Моменти розподілу помилок
  moments <- EstemPMM::compute_moments(data$errors)

  return(list(
    ols_coef = coef(ols_fit),
    pmm2_coef = pmm2_fit@coefficients,
    ols_mse = ols_mse,
    pmm2_mse = pmm2_mse,
    ols_aic = ols_aic,
    pmm2_aic = pmm2_aic,
    ols_bias = c(a0 = ols_bias_a0, a1 = ols_bias_a1),
    pmm2_bias = c(a0 = pmm2_bias_a0, a1 = pmm2_bias_a1),
    mse_ratio = pmm2_mse / ols_mse,
    moments = moments
  ))
}

# Функція для проведення симуляцій Монте-Карло
monte_carlo <- function(n_sim, n_samples, distribution, true_a0, true_a1, parallel = FALSE) {

  run_sim <- function(i) {
    data <- generate_data(n_samples, distribution, true_a0, true_a1)
    results <- compare_methods(data, true_a0, true_a1)

    # Зберегти тільки важливі результати для економії пам'яті
    return(list(
      ols_coef = results$ols_coef,
      pmm2_coef = results$pmm2_coef,
      ols_mse = results$ols_mse,
      pmm2_mse = results$pmm2_mse,
      mse_ratio = results$mse_ratio
    ))
  }

  # Використання паралельних обчислень, якщо вказано
  if (parallel && requireNamespace("parallel", quietly = TRUE)) {
    n_cores <- parallel::detectCores() - 1
    results <- parallel::mclapply(1:n_sim, function(i) run_sim(i), mc.cores = n_cores)
  } else {
    # Ініціалізація прогрес-бару
    pb <- txtProgressBar(min = 0, max = n_sim, style = 3)

    results <- vector("list", n_sim)
    for (i in 1:n_sim) {
      results[[i]] <- run_sim(i)
      setTxtProgressBar(pb, i)
    }
    close(pb)
  }

  # Обчислення статистик по результатам
  ols_a0 <- sapply(results, function(x) x$ols_coef[1])
  ols_a1 <- sapply(results, function(x) x$ols_coef[2])

  pmm2_a0 <- sapply(results, function(x) x$pmm2_coef[1])
  pmm2_a1 <- sapply(results, function(x) x$pmm2_coef[2])

  ols_mse <- sapply(results, function(x) x$ols_mse)
  pmm2_mse <- sapply(results, function(x) x$pmm2_mse)

  mse_ratio <- sapply(results, function(x) x$mse_ratio)

  # Генерація одноразового набору даних для обчислення моментів
  data <- generate_data(10000, distribution, true_a0, true_a1)
  moments <- EstemPMM::compute_moments(data$errors)

  return(list(
    ols_a0_mean = mean(ols_a0),
    ols_a0_sd = sd(ols_a0),
    ols_a1_mean = mean(ols_a1),
    ols_a1_sd = sd(ols_a1),

    pmm2_a0_mean = mean(pmm2_a0),
    pmm2_a0_sd = sd(pmm2_a0),
    pmm2_a1_mean = mean(pmm2_a1),
    pmm2_a1_sd = sd(pmm2_a1),

    ols_mse_mean = mean(ols_mse),
    pmm2_mse_mean = mean(pmm2_mse),

    mse_ratio_mean = mean(mse_ratio),

    theoretical_g = moments$g,
    c3 = moments$c3,
    c4 = moments$c4,

    ols_a0 = ols_a0,
    ols_a1 = ols_a1,
    pmm2_a0 = pmm2_a0,
    pmm2_a1 = pmm2_a1
  ))
}

# Функція для візуалізації результатів Монте-Карло
plot_monte_carlo_results <- function(mc_results, distribution_name, true_a0, true_a1) {
  # Розподіл оцінок a0
  p1 <- ggplot() +
    geom_density(aes(x = mc_results$ols_a0, fill = "OLS"), alpha = 0.5) +
    geom_density(aes(x = mc_results$pmm2_a0, fill = "PMM2"), alpha = 0.5) +
    geom_vline(xintercept = true_a0, linetype = "dashed") +
    labs(title = paste("Розподіл оцінок a0 -", distribution_name),
         x = "a0", y = "Густина") +
    scale_fill_manual(values = c("OLS" = "blue", "PMM2" = "red"),
                      name = "Метод") +
    theme_minimal()

  # Розподіл оцінок a1
  p2 <- ggplot() +
    geom_density(aes(x = mc_results$ols_a1, fill = "OLS"), alpha = 0.5) +
    geom_density(aes(x = mc_results$pmm2_a1, fill = "PMM2"), alpha = 0.5) +
    geom_vline(xintercept = true_a1, linetype = "dashed") +
    labs(title = paste("Розподіл оцінок a1 -", distribution_name),
         x = "a1", y = "Густина") +
    scale_fill_manual(values = c("OLS" = "blue", "PMM2" = "red"),
                      name = "Метод") +
    theme_minimal()

  # QQ-графіки
  p3 <- ggplot() +
    geom_qq(aes(sample = mc_results$ols_a0 - true_a0, color = "OLS")) +
    geom_qq(aes(sample = mc_results$pmm2_a0 - true_a0, color = "PMM2")) +
    geom_qq_line(aes(sample = mc_results$ols_a0 - true_a0)) +
    labs(title = paste("QQ-графік для a0 -", distribution_name),
         x = "Теоретичні квантилі", y = "Вибіркові квантилі") +
    scale_color_manual(values = c("OLS" = "blue", "PMM2" = "red"),
                       name = "Метод") +
    theme_minimal()

  p4 <- ggplot() +
    geom_qq(aes(sample = mc_results$ols_a1 - true_a1, color = "OLS")) +
    geom_qq(aes(sample = mc_results$pmm2_a1 - true_a1, color = "PMM2")) +
    geom_qq_line(aes(sample = mc_results$ols_a1 - true_a1)) +
    labs(title = paste("QQ-графік для a1 -", distribution_name),
         x = "Теоретичні квантилі", y = "Вибіркові квантилі") +
    scale_color_manual(values = c("OLS" = "blue", "PMM2" = "red"),
                       name = "Метод") +
    theme_minimal()

  # Відображення статистик
  stats_text <- paste(
    paste("Асиметрія (c3):", round(mc_results$c3, 4)),
    paste("Ексцес (c4):", round(mc_results$c4, 4)),
    paste("Теоретичний коефіцієнт g:", round(mc_results$theoretical_g, 4)),
    paste("Фактичний коефіцієнт (MSE):", round(mc_results$mse_ratio_mean, 4)),
    paste("Середнє значення a0 (OLS):", round(mc_results$ols_a0_mean, 4),
          "±", round(mc_results$ols_a0_sd, 4)),
    paste("Середнє значення a0 (PMM2):", round(mc_results$pmm2_a0_mean, 4),
          "±", round(mc_results$pmm2_a0_sd, 4)),
    paste("Середнє значення a1 (OLS):", round(mc_results$ols_a1_mean, 4),
          "±", round(mc_results$ols_a1_sd, 4)),
    paste("Середнє значення a1 (PMM2):", round(mc_results$pmm2_a1_mean, 4),
          "±", round(mc_results$pmm2_a1_sd, 4)),
    sep = "\n"
  )

  p5 <- ggplot() +
    annotate("text", x = 0, y = 0.5, label = stats_text, hjust = 0) +
    theme_void() +
    xlim(0, 1) + ylim(0, 1) +
    labs(title = paste("Статистика -", distribution_name))

  # Об'єднання графіків
  grid.arrange(p1, p2, p3, p4, p5, ncol = 2)
}

#############################################################
# Налаштування та запуск симуляцій
#############################################################

set.seed(42)

# Параметри симуляцій
n_sim <- 1000       # Кількість симуляцій Монте-Карло
n_samples <- 100    # Розмір вибірки в кожній симуляції
true_a0 <- 2        # Справжнє значення a0
true_a1 <- 1.5      # Справжнє значення a1

# Для швидкого демонстраційного запуску можна використати:
# n_sim <- 100      # Менша кількість симуляцій для швидшого запуску
# n_samples <- 50   # Менший розмір вибірки

# Список розподілів для тестування
distributions <- c("gaussian", "t", "gamma", "exponential", "chi-squared", "lognormal")
distribution_names <- c("Нормальний", "Стьюдента (df=3)", "Гамма (a=2)",
                        "Експоненційний", "Хі-квадрат (df=3)", "Логнормальний")

# Функція для виконання всіх симуляцій
run_all_simulations <- function() {
  results <- list()

  for(i in 1:length(distributions)) {
    cat("\nВиконується симуляція для розподілу:", distribution_names[i], "\n")

    mc_results <- monte_carlo(n_sim, n_samples,
                              distributions[i],
                              true_a0, true_a1,
                              parallel = TRUE)

    results[[distributions[i]]] <- mc_results

    # Візуалізація результатів
    plot_monte_carlo_results(mc_results, distribution_names[i], true_a0, true_a1)

    # Вивід результатів у консоль
    cat("\nРезультати для розподілу:", distribution_names[i], "\n")
    cat("Асиметрія (c3):", round(mc_results$c3, 4), "\n")
    cat("Ексцес (c4):", round(mc_results$c4, 4), "\n")
    cat("Теоретичний коефіцієнт g:", round(mc_results$theoretical_g, 4), "\n")
    cat("Фактичний коефіцієнт (MSE):", round(mc_results$mse_ratio_mean, 4), "\n")
    cat("Покращення ефективності:",
        round((1 - mc_results$mse_ratio_mean) * 100, 2), "%\n")

    cat("a0 (OLS):", round(mc_results$ols_a0_mean, 4),
        "±", round(mc_results$ols_a0_sd, 4), "\n")
    cat("a0 (PMM2):", round(mc_results$pmm2_a0_mean, 4),
        "±", round(mc_results$pmm2_a0_sd, 4), "\n")

    cat("a1 (OLS):", round(mc_results$ols_a1_mean, 4),
        "±", round(mc_results$ols_a1_sd, 4), "\n")
    cat("a1 (PMM2):", round(mc_results$pmm2_a1_mean, 4),
        "±", round(mc_results$pmm2_a1_sd, 4), "\n")
  }

  # Створення підсумкового порівняння всіх розподілів
  summary_table <- data.frame(
    Distribution = distribution_names,
    Skewness = sapply(results, function(x) round(x$c3, 4)),
    Kurtosis = sapply(results, function(x) round(x$c4, 4)),
    Theoretical_g = sapply(results, function(x) round(x$theoretical_g, 4)),
    Actual_g = sapply(results, function(x) round(x$mse_ratio_mean, 4)),
    Improvement = sapply(results, function(x)
      round((1 - x$mse_ratio_mean) * 100, 2))
  )

  cat("\n\nПідсумок всіх симуляцій:\n")
  print(summary_table)

  # Візуалізація порівняння ефективності
  p <- ggplot(summary_table, aes(x = reorder(Distribution, -Improvement))) +
    geom_bar(aes(y = Improvement, fill = "Фактичне"), stat = "identity",
             alpha = 0.7, position = position_dodge()) +
    geom_point(aes(y = (1 - Theoretical_g) * 100, color = "Теоретичне"),
               size = 3) +
    labs(title = "Порівняння ефективності PMM2 відносно OLS",
         x = "Розподіл помилок",
         y = "Покращення ефективності (%)") +
    scale_fill_manual(values = c("Фактичне" = "steelblue"),
                      name = "Покращення") +
    scale_color_manual(values = c("Теоретичне" = "red"),
                       name = "Покращення") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  print(p)

  return(list(results = results, summary = summary_table))
}

# Щоб запустити всі симуляції, викличте вручну:
# run_all_simulations()
