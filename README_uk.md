# EstemPMM: Метод Поліноміальної Максимізації (PMM2) для регресії та часових рядів

[![R >= 4.0.0](https://img.shields.io/badge/R-%3E%3D%204.0.0-blue)](https://cran.r-project.org/)
[![License: GPL-3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://opensource.org/licenses/GPL-3.0)

## Огляд

`EstemPMM` реалізує Метод Поліноміальної Максимізації (PMM) для оцінювання параметрів лінійних та часових моделей при НЕнормальних похибках (інноваціях), особливо за наявності асиметрії та зміненої крутості (kurtosis). Для степеня полінома `S = 1` PMM збігається з OLS. Для `S = 2` (PMM2) можливе зниження дисперсії оцінок порівняно зі звичайним методом найменших квадратів (OLS) за рахунок використання інформації про асиметрію.

## Теорія

Нехай `c3` — коефіцієнт асиметрії, `c4` — надлишкова крутість (excess kurtosis). Теоретичний коефіцієнт:

```
g = 1 - c3^2 / (2 + c4)
```

Інтерпретація в пакеті:
- `g` ≈ Var(PMM2) / Var(OLS).
- Відсоток покращення ≈ (1 - g) * 100%.

Перевірте узгодженість між емпіричними симуляціями та цією формулою у майбутніх розширеннях.

## Ключові можливості

- Лінійні моделі: `lm_pmm2()`.`
- Часові ряди: `ar_pmm2()`, `ma_pmm2()`, `arma_pmm2()`, `arima_pmm2()`, сезонні `sar_pmm2()`, `sma_pmm2()`, універсальна обгортка `ts_pmm2()`.`
- Адаптивне двоетапне оцінювання (OLS → моменти → PMM2).` 
- Бутстреп для невизначеності: `pmm2_inference()`, `ts_pmm2_inference()`, графіки через `plot_pmm2_bootstrap()`. `
- Оцінка виграшу в дисперсії: `pmm2_variance_factor()`, `pmm2_variance_matrices()`. `
- Порівняння методів: `compare_with_ols()`, `compare_ts_methods()` та спеціалізовані `compare_*`. `
- Діагностика (метод `plot()`): QQ-плоти, залишки тощо. `
- Утиліти моментів: `compute_moments()`, `pmm_skewness()`, `pmm_kurtosis()`. `

## Встановлення

```r
devtools::install_github("SZabolotnii/EstemPMM")
```

## Швидкий приклад (лінійна модель)

```r
library(EstemPMM)
set.seed(1)
n <- 200
x <- rnorm(n)
errors <- rgamma(n, shape = 2, scale = 1) - 2
y <- 2 + 1.5 * x + errors
dat <- data.frame(y, x)

fit <- lm_pmm2(y ~ x, data = dat)
summary(fit)
compare_with_ols(y ~ x, dat)

vf <- pmm2_variance_factor(fit@m2, fit@m3, fit@m4)
vf$g
```

## Приклад для часових рядів

```r
set.seed(42)
n <- 300
innov <- rgamma(n, shape = 2, scale = 1) - 2
y <- numeric(n)
for (t in 3:n) {
  y[t] <- 0.5 * y[t-1] - 0.3 * y[t-2] + innov[t]
}

fit_ar <- ar_pmm2(y, order = 2)
summary(fit_ar)
compare_ar_methods(y, order = 2)
```

## Повний перелік експортованих функцій

Оцінювання:
`lm_pmm2`, `ar_pmm2`, `ma_pmm2`, `arma_pmm2`, `arima_pmm2`, `sar_pmm2`, `sma_pmm2`, `ts_pmm2`

Порівняння:
`compare_with_ols`, `compare_ar_methods`, `compare_ma_methods`, `compare_arma_methods`, `compare_arima_methods`, `compare_sar_methods`, `compare_ts_methods`

Моменти:
`compute_moments`, `pmm_skewness`, `pmm_kurtosis`

Дисперсія:
`pmm2_variance_factor`, `pmm2_variance_matrices`

Симуляції:
`pmm2_monte_carlo_compare`

Бутстреп:
`pmm2_inference`, `ts_pmm2_inference`, `plot_pmm2_bootstrap`

Допоміжні:
`create_sar_matrix`

S4 класи:
`PMM2fit`, `TS2fit`, `ARPMM2`, `MAPMM2`, `ARMAPMM2`, `ARIMAPMM2`, `SARPMM2`, `SMAPMM2`, `BasePMM2`

S4 методи:
`coef`, `fitted`, `predict`, `residuals`, `summary`, `plot`, `AIC`

## S4 Об’єкти

Містять:
- Оцінені коефіцієнти.
- Оцінки моментів (m2, m3, m4).
- Залишки, матриця моделі (для регресії).
- Параметри порядку для часових моделей.

## Адаптивна процедура

1. OLS оцінювання.
2. Оцінка моментів залишків.
3. Побудова поліномів максимізації.
4. PMM2 оцінки з покращеною ефективністю.

## Асиметрія та крутість

```r
mom <- compute_moments(residuals(fit))
mom$skewness
mom$kurtosis

pmm_skewness(residuals(fit))
pmm_kurtosis(residuals(fit))
```

## Бутстреп

```r
inf_res <- pmm2_inference(fit, R = 500)
plot_pmm2_bootstrap(inf_res)
summary(inf_res)

ts_inf <- ts_pmm2_inference(fit_ar, R = 300)
```

## Монте Карло

Скрипти:
- `run_sma_monte_carlo.R`
- `test_sma_quick.R`
- `test_debug_sar.R`

Рекомендовано створити вінієтку “Monte Carlo Efficiency”.

## Ефективність (приклад)

| Розподіл            | Skewness | Kurtosis | Теор. покращення | Емпірично |
|---------------------|----------|----------|------------------|-----------|
| Gamma (shape=0.5)   | 2.83     | 12       | 57%              | ~50%      |
| Exponential         | 2.00     | 6        | 50%              | ~45%      |
| Gamma (shape=2)     | 1.41     | 3        | 40%              | ~35%      |
| Lognormal           | 1.00     | 1.5      | 29%              | ~25%      |

Уточніть, чи подана Kurtosis — надлишкова (excess).

## Застосування

- Фінанси та економіка (асиметричні ризики).
- Біостатистика.
- Технічні вимірювання з негаусівським шумом.
- Будь-які моделі з суттєвою асиметрією.

## Порівняння з альтернативами (для майбутнього)

- Робастні M-оцінки.
- Quantile regression.
- Умови максимальної вигоди (висока асиметрія + помірна надлишкова крутість).

## План розвитку

- [ ] CI та покриття тестів.
- [ ] Вінієтки (вступ, часові ряди, бутстреп).
- [ ] Публікація на CRAN.
- [ ] Оптимізація (векторизація / Rcpp).
- [ ] Експерименти для S > 2.
- [ ] Підтримка факторних змінних.
- [ ] Синтетичні набори даних.
- [ ] Бенчмарки (`bench`).

## Публікації

1. Kunchenko, Y. P., & Lega, Y. G. (1992). *Estimation of Random Variable Parameters by the Polynomial Maximization Method*. Kyiv: Naukova Dumka.
2. Zabolotnii S., Warsza Z.L., Tkachenko O. (2018). Polynomial Estimation of Linear Regression Parameters for the Asymmetric PDF of Errors. …
3. Zabolotnii S., Tkachenko O., Warsza Z.L. (2022). Application … Autoregressive Models …
4. Zabolotnii S., Tkachenko O., Warsza Z.L. (2023). Polynomial Maximization Method … Moving Average Models …

(Додайте DOI або лінки.)

## Автор

- Сергій Заболотній – Cherkasy State Business College

## Ліцензія

GPL-3.

---

Англомовна версія: див. файл [README.md](README.md).