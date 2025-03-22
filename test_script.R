# Скрипт для тестування пакету PMM2

# 1. Завантажуємо необхідні функції з покращених файлів
source("R/pmm_classes.R")
source("R/pmm_utils.R")
source("R/pmm_main.R")
source("R/pmm_inference.R")

# 2. Запускаємо простий тест з test-pmm2.R
cat("Запускаємо базовий тест PMM2...\n\n")

set.seed(123)
n <- 50
x <- rnorm(n)
# Генеруємо y з негаусовою похибкою:
y <- 1 + 2*x + rt(n, df=4)
dat <- data.frame(x, y)

cat("Підгонка моделі PMM2...\n")
fit <- lm_pmm2(y ~ x, data=dat, max_iter=50, verbose=TRUE)

cat("\n\nРезультати підгонки PMM2:\n")
print(fit)

cat("\n\nСтислий огляд результатів:\n")
summary(fit)

cat("\n\nПроводимо інференцію через бутстреп:\n")
inf_tab <- pmm2_inference(fit, formula=y~x, data=dat, B=50, seed=456)
print(inf_tab)

cat("\n\nПорівнюємо з методом найменших квадратів:\n")
comparison <- compare_with_ols(y ~ x, dat)
print(comparison$coefficients)
print(comparison$residual_stats)

cat("\n\nВізуалізуємо діагностику моделі PMM2:\n")
# Збережемо поточні налаштування графіки
old_par <- par(no.readonly = TRUE)
# Встановимо параметри для графіків
par(mfrow=c(2,2), mar=c(4,4,2,1))
# Побудуємо діагностичні графіки
plot(fit)
# Відновимо налаштування графіки
par(old_par)

cat("\n\nВізуалізуємо розподіл бутстреп-оцінок:\n")
# Встановимо параметри для графіків
par(mfrow=c(1,2), mar=c(4,4,2,1))
# Побудуємо графіки розподілу оцінок
plot_pmm2_bootstrap(inf_tab)
# Відновимо налаштування графіки
par(old_par)

cat("\n\nПерепідгонка зі зміненими параметрами (більша регуляризація):\n")
fit2 <- lm_pmm2(y ~ x, data=dat, max_iter=50, regularize=TRUE, reg_lambda=1e-6, verbose=TRUE)
summary(fit2)

cat("\n\nТестування завершено!\n")
