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
# Стабільна версія на CRAN
install.packages("EstemPMM")

# Або останній девелоперський знімок з GitHub
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

### Оцінювання моделей

**Лінійні моделі:**
- `lm_pmm2()` - лінійна регресія з PMM2

**Моделі часових рядів:**
- `ar_pmm2()` - авторегресійні AR(p) моделі
- `ma_pmm2()` - моделі ковзного середнього MA(q)
- `arma_pmm2()` - змішані ARMA(p,q) моделі
- `arima_pmm2()` - ARIMA(p,d,q) моделі з диференціюванням
- `sar_pmm2()` - сезонні авторегресійні SAR(p,P)_s моделі
- `sma_pmm2()` - сезонні моделі ковзного середнього SMA(Q)_s
- `ts_pmm2()` - універсальна обгортка для всіх моделей часових рядів

### Функції порівняння

- `compare_with_ols()` - порівняння PMM2 з OLS для регресії
- `compare_ar_methods()` - порівняння методів оцінювання AR
- `compare_ma_methods()` - порівняння методів оцінювання MA
- `compare_arma_methods()` - порівняння методів оцінювання ARMA
- `compare_arima_methods()` - порівняння методів оцінювання ARIMA
- `compare_sar_methods()` - порівняння методів оцінювання SAR
- `compare_ts_methods()` - універсальне порівняння для часових рядів

### Моменти та дисперсія

- `compute_moments()` - обчислення вибіркових моментів
- `pmm_skewness()` - коефіцієнт асиметрії
- `pmm_kurtosis()` - коефіцієнт ексцесу (крутість)
- `pmm2_variance_factor()` - теоретичний коефіцієнт зниження дисперсії
- `pmm2_variance_matrices()` - порівняння матриць дисперсій

### Статистичний висновок (інференція)

- `pmm2_inference()` - бутстреп-інференція для лінійних моделей
- `ts_pmm2_inference()` - бутстреп-інференція для часових рядів
- `plot_pmm2_bootstrap()` - візуалізація результатів бутстрепу

### Симуляції

- `pmm2_monte_carlo_compare()` - порівняння методів через Монте-Карло

### Допоміжні функції

- `create_sar_matrix()` - створення матриці для SAR моделей

### S4 класи

`PMM2fit`, `BasePMM2`, `ARPMM2`, `MAPMM2`, `ARMAPMM2`, `ARIMAPMM2`, `SARPMM2`, `SMAPMM2`

### S4 методи

`coef()`, `fitted()`, `predict()`, `residuals()`, `summary()`, `plot()`, `AIC()`

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

## Ефективність

### Лінійні моделі

| Розподіл            | Skewness | Kurtosis | Теор. покращення | Емпірично |
|---------------------|----------|----------|------------------|-----------|
| Gamma (shape=0.5)   | 2.83     | 12       | 57%              | ~50%      |
| Exponential         | 2.00     | 6        | 50%              | ~45%      |
| Gamma (shape=2)     | 1.41     | 3        | 40%              | ~35%      |
| Lognormal           | 1.00     | 1.5      | 29%              | ~25%      |

*Примітка: Kurtosis — надлишкова крутість (excess kurtosis).*

### Моделі часових рядів

**Авторегресійні (AR) моделі:** Подібні покращення спостерігаються, особливо для AR(1) та AR(2) з асиметричними інноваціями.

**Моделі ковзного середнього (MA):** PMM2 демонструє суттєві покращення при високій асиметрії інновацій.

**Сезонні моделі (валідовано у листопаді 2025):**
- **SAR моделі:** 20-30% зниження дисперсії з гамма-інноваціями
- **SMA моделі:** До 34.1% зниження дисперсії (емпірично підтверджено на 500 повтореннях), що перевищує теоретичні прогнози для деяких комбінацій параметрів
- Обидві моделі показують стабільну збіжність та обчислювальну ефективність

## Структура пакета

### Основна реалізація (R/)
- `pmm2_classes.R` - визначення S4 класів (PMM2fit, ARPMM2, MAPMM2, ARMAPMM2, ARIMAPMM2, SARPMM2, SMAPMM2)
- `pmm2_main.R` - лінійна регресія з PMM2
- `pmm2_ts_main.R` - моделі часових рядів (AR, MA, ARMA, ARIMA, SAR, SMA)
- `pmm2_ts_design.R` - матриці дизайну та структури лагів
- `pmm2_ts_methods.R` - функції порівняння методів
- `pmm2_common.R` - спільні чисельні процедури
- `pmm2_utils.R` - утиліти для моментів та аналізу дисперсії
- `pmm2_inference.R` - бутстреп-інференція
- `pmm2_monte_carlo.R` - фреймворк для Монте-Карло симуляцій
- `data.R` - документація наборів даних

### Документація
- `vignettes/` - три детальні підручники
- `man/` - більше 30 довідкових файлів
- `docs/` - розширена документація, включаючи звіти по SAR/SMA
- `NEWS.md` - історія версій

### Тестування та валідація
- `tests/testthat/` - 36 модульних тестів
- `test_results/` - звіти валідації методом Монте-Карло
- `demo/` - 8 демонстраційних скриптів

## Документація та віньєтки

| Віньєтка | Опис |
| --- | --- |
| `vignette(\"pmm2-introduction\")` | базові регресійні кейси з PMM2 |
| `vignette(\"pmm2-time-series\")` | AR/MA/ARMA/ARIMA + сезонні SAR/SMA/SARMA приклади |
| `vignette(\"bootstrap-inference\")` | бутстреп-інференція для асиметричних похибок |

Перегенерувати документацію на чистому середовищі:

```r
devtools::document()
devtools::build_vignettes()
```

## Відтворення Монте-Карло експериментів

- **Базовий SMA-бенчмарк (n = 120, γ-інновації):** `Rscript run_sma_monte_carlo.R`.  
  Результати збережено у `test_results/SMA_Monte_Carlo_Report_20251113_500reps.md` (≈34% зменшення дисперсії).
- **Розширене сезонне порівняння (SAR/SMA/SARMA/SARIMA, n = 100/200/500):** `Rscript monte_carlo_seasonal_comparison.R`.  
  Підсумки фіксуються у `test_results/SAR_MONTE_CARLO_REPORT_2025-11-14.md` та файлі `monte_carlo_seasonal_results.rds` (тривалість ≈8 хв на Apple Silicon).

## Перевірки перед публікацією на CRAN

```r
devtools::check()
devtools::document()
devtools::build_vignettes()
system(\"R CMD build .\")
system(\"R CMD check --as-cran EstemPMM_0.1.3.tar.gz\")
```

GitHub Actions (workflow `R-CMD-check.yaml`) автоматично запускає ті ж перевірки на Ubuntu, macOS та Windows для кожного пуша в `main` або `claude/*`.

## Застосування

### Регресійний аналіз
- Економічне та фінансове моделювання з асиметричними похибками
- Аналіз біологічних систем зі скошеними похибками вимірювань
- Технічні вимірювання з негаусівським шумом
- Будь-які задачі регресії з суттєвою асиметрією похибок

### Аналіз часових рядів
- Фінансові часові ряди з асиметричними розподілами доходності
- Економічні індикатори з сезонними патернами (ВВП, інфляція, безробіття)
- Прогнозування споживання енергії з сезонними компонентами
- Аналіз кліматичних даних зі скошеними розподілами температури чи опадів
- Ціни на сировину (нафта, газ) з асиметричними шоками
- Будь-які часові ряди з інноваціями, що суттєво відхиляються від нормальності

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
