# Простий сценарій для перевірки основних можливостей EstemPMM

if (!requireNamespace("EstemPMM", quietly = TRUE)) {
  stop("Встановіть пакет EstemPMM перед запуском демо", call. = FALSE)
}

library(EstemPMM)

set.seed(2025)
n <- 150
x1 <- rnorm(n)
x2 <- rnorm(n)
errors <- rgamma(n, shape = 2, scale = 1) - 2
y <- 1.5 + 2 * x1 - 1 * x2 + errors

dat <- data.frame(y = y, x1 = x1, x2 = x2)

cat("\n=== Лінійна регресія: PMM2 vs OLS ===\n")
ols_fit <- lm(y ~ x1 + x2, data = dat)
pmm2_fit <- lm_pmm2(y ~ x1 + x2, data = dat)

print(summary(ols_fit))
print(summary(pmm2_fit))

cat("\nКоефіцієнти PMM2:\n")
print(coef(pmm2_fit))

new_data <- data.frame(x1 = c(-1, 0, 1), x2 = c(0.5, -0.3, 1.2))
cat("\nПрогноз PMM2 для нових даних:\n")
print(predict(pmm2_fit, newdata = new_data))

cat("\n=== Приклад часових рядів: AR(1) з асиметричними похибками ===\n")
ts_innov <- rgamma(n, shape = 2, scale = 1) - 2
ts_series <- filter(ts_innov, filter = 0.6, method = "recursive")
ar_fit <- ar_pmm2(ts_series, order = 1)
print(summary(ar_fit))

cat("\nДемо завершено.\n")
