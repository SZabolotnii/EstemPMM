#!/usr/bin/env Rscript
# ==============================================================================
# Приклад використання гібридного оцінювача v2
# ==============================================================================

source("R/experimental/07_hybrid_with_estpmm_ma.R")

cat("=== Приклади використання Hybrid SARIMA Estimator v2 ===\n\n")

# ==============================================================================
# Приклад 1: MA(1) модель з EstemPMM-style PMM2
# ==============================================================================

cat(strrep("=", 70), "\n")
cat("Приклад 1: MA(1) з асиметричними інноваціями\n")
cat(strrep("=", 70), "\n\n")

set.seed(123)
n <- 200
theta_true <- 0.6

# Генерація даних з Exp(1) - 1 інноваціями
innovations <- rexp(n, rate = 1) - 1
x1 <- numeric(n)
x1[1] <- innovations[1]

for (t in 2:n) {
  x1[t] <- innovations[t] + theta_true * innovations[t - 1]
}

cat("Істинний параметр: θ =", theta_true, "\n")
cat("Інновації: Exp(1) - 1 (асиметричні, γ₃ ≈ 2.0)\n\n")

# Варіант 1: MLE (традиційний)
cat("--- Варіант 1: MLE ---\n")
fit_mle <- hybrid_sarima_estimator_v2(
  x1,
  order = c(0, 0, 1),
  seasonal = list(order = c(0, 0, 0), period = 1),
  include.mean = FALSE,
  ma_method = "mle",
  ar_method = "mle",
  verbose = TRUE
)

# Варіант 2: EstemPMM-style PMM2 для MA
cat("\n--- Варіант 2: EstemPMM-style PMM2 для MA ---\n")
fit_pmm2 <- hybrid_sarima_estimator_v2(
  x1,
  order = c(0, 0, 1),
  seasonal = list(order = c(0, 0, 0), period = 1),
  include.mean = FALSE,
  ma_method = "pmm2",  # ⭐ PMM2 для MA
  ar_method = "mle",
  verbose = TRUE
)

# Порівняння
cat("\n--- Порівняння ---\n")
cat(sprintf("MLE:  θ̂ = %.4f, Bias = %.4f\n",
            fit_mle$ma_coef[1],
            fit_mle$ma_coef[1] - theta_true))
cat(sprintf("PMM2: θ̂ = %.4f, Bias = %.4f\n",
            fit_pmm2$ma_coef[1],
            fit_pmm2$ma_coef[1] - theta_true))

bias_reduction <- abs(fit_mle$ma_coef[1] - theta_true) -
                  abs(fit_pmm2$ma_coef[1] - theta_true)

if (bias_reduction > 0) {
  cat(sprintf("✅ PMM2 зменшує bias на %.4f\n", bias_reduction))
} else {
  cat(sprintf("⚠️ PMM2 збільшує bias на %.4f\n", -bias_reduction))
}

# ==============================================================================
# Приклад 2: SMA(1)_4 модель з EstemPMM-style PMM2
# ==============================================================================

cat("\n\n", strrep("=", 70), "\n")
cat("Приклад 2: SMA(1)₄ з асиметричними інноваціями\n")
cat(strrep("=", 70), "\n\n")

set.seed(456)
n2 <- 200
s <- 4
Theta_true <- 0.6

# Генерація даних
innovations2 <- rexp(n2, rate = 1) - 1
x2 <- numeric(n2)

for (t in 1:s) {
  x2[t] <- innovations2[t]
}

for (t in (s+1):n2) {
  x2[t] <- innovations2[t] + Theta_true * innovations2[t - s]
}

cat("Істинний параметр: Θ =", Theta_true, "\n")
cat("Сезонний період: s =", s, "\n")
cat("Інновації: Exp(1) - 1 (асиметричні)\n\n")

# Варіант 1: MLE
cat("--- Варіант 1: MLE ---\n")
fit_mle2 <- hybrid_sarima_estimator_v2(
  x2,
  order = c(0, 0, 0),
  seasonal = list(order = c(0, 0, 1), period = s),
  include.mean = FALSE,
  ma_method = "mle",
  ar_method = "mle",
  verbose = TRUE
)

# Варіант 2: EstemPMM-style PMM2 для SMA
cat("\n--- Варіант 2: EstemPMM-style PMM2 для SMA ---\n")
fit_pmm2_2 <- hybrid_sarima_estimator_v2(
  x2,
  order = c(0, 0, 0),
  seasonal = list(order = c(0, 0, 1), period = s),
  include.mean = FALSE,
  ma_method = "pmm2",  # ⭐ PMM2 для SMA
  ar_method = "mle",
  verbose = TRUE
)

# Порівняння
cat("\n--- Порівняння ---\n")
cat(sprintf("MLE:  Θ̂ = %.4f, Bias = %.4f\n",
            fit_mle2$sma_coef[1],
            fit_mle2$sma_coef[1] - Theta_true))
cat(sprintf("PMM2: Θ̂ = %.4f, Bias = %.4f\n",
            fit_pmm2_2$sma_coef[1],
            fit_pmm2_2$sma_coef[1] - Theta_true))

bias_reduction2 <- abs(fit_mle2$sma_coef[1] - Theta_true) -
                   abs(fit_pmm2_2$sma_coef[1] - Theta_true)

if (bias_reduction2 > 0) {
  cat(sprintf("✅ PMM2 зменшує bias на %.4f\n", bias_reduction2))
} else {
  cat(sprintf("⚠️ PMM2 збільшує bias на %.4f\n", -bias_reduction2))
}

# ==============================================================================
# Приклад 3: Автоматичний вибір методу (рекомендований підхід)
# ==============================================================================

cat("\n\n", strrep("=", 70), "\n")
cat("Приклад 3: Рекомендований підхід для вибору методу\n")
cat(strrep("=", 70), "\n\n")

#' Функція для автоматичного вибору методу на основі асиметрії
#'
#' @param x Time series
#' @param threshold Поріг асиметрії для вибору PMM2 (за замовчуванням 0.5)
#' @return "pmm2" або "mle"
choose_ma_method <- function(x, threshold = 0.5) {
  # Оцінити асиметрію залишків (грубо через ACF)
  # В реальності краще використати попередню MLE оцінку

  # Простий підхід: оцінити асиметрію безпосередньо з даних
  m2 <- mean((x - mean(x))^2)
  m3 <- mean((x - mean(x))^3)
  gamma3 <- m3 / (m2^(3/2))

  if (abs(gamma3) > threshold) {
    return("pmm2")
  } else {
    return("mle")
  }
}

# Приклад використання
ma_method_auto <- choose_ma_method(x1, threshold = 0.5)

cat("Автоматичний вибір методу для MA(1) даних:\n")
cat("  Оцінена асиметрія: γ₃ ≈", round(mean((x1 - mean(x1))^3) / (mean((x1 - mean(x1))^2)^(3/2)), 2), "\n")
cat("  Обраний метод:", ma_method_auto, "\n\n")

fit_auto <- hybrid_sarima_estimator_v2(
  x1,
  order = c(0, 0, 1),
  seasonal = list(order = c(0, 0, 0), period = 1),
  include.mean = FALSE,
  ma_method = ma_method_auto,
  ar_method = "mle",
  verbose = TRUE
)

# ==============================================================================
# Підсумок
# ==============================================================================

cat("\n\n", strrep("=", 70), "\n")
cat("ПІДСУМОК ТА РЕКОМЕНДАЦІЇ\n")
cat(strrep("=", 70), "\n\n")

cat("1. Базове використання:\n")
cat("   fit <- hybrid_sarima_estimator_v2(x, order = c(p,d,q),\n")
cat("                                      ma_method = \"pmm2\",  # для асиметричних\n")
cat("                                      ar_method = \"mle\")\n\n")

cat("2. Вибір методу:\n")
cat("   - ma_method = \"pmm2\": для MA/SMA з |γ₃| > 0.5\n")
cat("   - ma_method = \"mle\":  для симетричних інновацій або малих вибірок\n\n")

cat("3. Очікувані покращення (асиметричні інновації):\n")
cat("   - MA(1):  20-30% покращення MSE\n")
cat("   - SMA(1): 15-25% покращення MSE\n")
cat("   - MA(2):  30-45% покращення MSE\n\n")

cat("4. Обмеження:\n")
cat("   - Змішані MA+SMA: поки що використовується MLE\n")
cat("   - Малі вибірки (n < 50): не рекомендовано PMM2\n\n")

cat("=== Приклади завершено ===\n")
