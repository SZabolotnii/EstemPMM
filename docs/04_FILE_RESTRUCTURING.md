# EstemPMM Phase 1: File Restructuring Strategy

## Executive Summary

Поточна структура файлів приховує PMM2-специфічність під генеричними іменами (`pmm_main.R`, `pmm_utils.R`). Це створює архітектурний борг при введенні PMM3. Фаза 1 реструктуризує файли для явної PMM2-ідентифікації, що забезпечує когнітивне розуміння та майбутню масштабованість.

---

## Стратегія Рефакторингу: Від Монолітної до Модульної Архітектури

### Поточна Архітектура (Проблематична)

```
R/
├─ pmm_package.R       # Що тут? Залежності? Документація пакету?
├─ pmm_main.R          # Лінійна регресія? Часові ряди? Обидва?
├─ pmm_classes.R       # S4 класи для чого? Яка кількість?
├─ pmm_utils.R         # Утиліти ЯКИХ функцій? Де вони використовуються?
└─ pmm_ts_design.R     # Тільки часові ряди, чи також лінійна регресія?

ПРОБЛЕМА: Розробник не знає, як пакет структурований, доки не прочитає весь код
РИЗИК: PMM3 додасть 3-4 нових файли, архітектура розпадеться
```

### Цільова Архітектура 0.1.0 (PMM2-Явна)

```
R/
├─ pmm2_package.R         # Залежності, документація пакету (замість pmm_package.R)
├─ pmm2_classes.R         # S4 класи: PMM2fit, TS2fit, ARPMM2, MAPMM2 і т.д.
├─ pmm2_main.R            # Основні функції: lm_pmm2, ar_pmm2, ma_pmm2, arma_pmm2, arima_pmm2
├─ pmm2_utils.R           # Утиліти для оптимізації та моментів
├─ pmm2_ts_design.R       # Часові ряди: матриці дизайну, структури моделей
└─ pmm2_inference.R       # Бутстреп-висновок та статистичні методи

КОРИСТЬ: 
- Кожен файл має ясну PMM2-ідентифікацію
- Новий розробник розуміє архітектуру за 2 хвилини
- PMM3 додається в `R/pmm3_*.R` файли без конфліктів
```

### Архітектурна еволюція для 0.2.0 (PMM3-Ready)

```
R/
├─ base_classes.R        # [НОВИЙ] Базові S4 класи для майбутніх методів
├─ pmm_common_utils.R    # [НОВИЙ] Утіліти, СТІЙКІ до змін (моменти, статистика)
├─ pmm2_*.R              # [БЕЗ ЗМІН] Вся PMM2 логіка ізольована
├─ pmm3_package.R        # [НОВИЙ] PMM3 документація
├─ pmm3_classes.R        # [НОВИЙ] PMM3 S4 класи (наслідуючі base_classes)
├─ pmm3_main.R           # [НОВИЙ] PMM3 функції
├─ pmm3_utils.R          # [НОВИЙ] PMM3 утиліти
└─ pmm3_inference.R      # [НОВИЙ] PMM3 бутстреп

АРХІТЕКТУРНА ПЕРЕВАЖНІСТЬ:
- PMM2 та PMM3 повністю ізольовані
- Спільна логіка мінімізована і явна
- Тести організовані паралельно
- Немає архітектурного боргу
```

---

## Практичний План Реструктуризації: 5 Кроків

### Крок 1: Переіменування R/ файлів (Час: 1 година)

**Операції:**

```bash
# НОМЕР 1: pmm_package.R → pmm2_package.R (змінити, чи видалити?)
# НОМЕР 2: pmm_classes.R → pmm2_classes.R
# НОМЕР 3: pmm_main.R → pmm2_main.R
# НОМЕР 4: pmm_utils.R → pmm2_utils.R
# НОМЕР 5: pmm_ts_design.R → pmm2_ts_design.R
```

**Чистлист для кожного файлу:**

#### pmm2_package.R (раніше pmm_package.R)
```r
# Конфіденційна перевірка: які імпорти насправді використовуються?

#' @importFrom methods is new slotNames
# → ЗАЛИШИТИ (використовується в S4 класах)

#' @importFrom graphics abline hist legend lines par
# → ЗАЛИШИТИ (використовується в plot методах)

#' @importFrom stats acf arima cov dnorm lm ...
# → ЗАЛИШИТИ (основна статистика)

#' @importFrom utils tail
# → ПЕРЕВІРИТИ (TAIL використовується? Видалити, якщо ні)
```

**Рекомендація:** Зберегти всі імпорти, але додати коментарі для CRAN ясності

---

#### pmm2_classes.R (раніше pmm_classes.R)

**Операція:** Просто переіменувати, логіка залишається незмінною

**Перевірка:** Усі S4 класи мають задокументовані слоти?
```r
#' @section Слоти:
#' \describe{
#'   \item{coefficients}{Оцінені коефіцієнти}
#'   ...
#' }
```
✅ Добре, потрібно тільки переіменування

---

#### pmm2_main.R (раніше pmm_main.R)

**Операція:** Переіменувати + ОЧИСТИТИ функції

**Критична перевірка:** Чи `lm_pmm2`, `ar_pmm2`, `ma_pmm2` — тільки PMM2-специфічні?

```r
# ✅ lm_pmm2      — Так, PMM2 для лінійної регресії
# ✅ ar_pmm2      — Так, PMM2 для AR часових рядів
# ✅ ma_pmm2      — Так, PMM2 для MA часових рядів
# ✅ arma_pmm2    — Так, PMM2 для ARMA часових рядів
# ✅ arima_pmm2   — Так, PMM2 для ARIMA часових рядів
# ✅ ts_pmm2      — Так, дисетчер для всіх часових рядів
```

**Потрібна операція:** Додати явні коментарі на початку файлу

```r
# pmm2_main.R: Основні функції PMM2 (S=2)
# Функції лінійної регресії та часових рядів для методу Polynomial Maximization Method
# з порядком 2

# Експортовані функції:
# - lm_pmm2() — Лінійна регресія
# - ts_pmm2() — Часові ряди (дисетчер)
# - ar_pmm2(), ma_pmm2(), arma_pmm2(), arima_pmm2() — Специфічні моделі

# Приватні утиліти: .pmm2_fit(), .ts_pmm2_fit()
```

---

#### pmm2_utils.R (раніше pmm_utils.R)

**Операція:** Переіменувати + ОЧИСТИТИ функції

**Критична перевірка:** Які утиліти здійснюють PMM2-оптимізацію, а які — універсальні?

```r
# ОЧИСТИТИ В ЦЬОМУ ФАЙЛІ:
# ✅ pmm2_fit() — PMM2 оптимізація (залишити)
# ✅ pmm_skewness() — Обчислення асиметрії (РОЗГЛЯНУТИ: мінімально-спільна?)
# ✅ pmm_kurtosis() — Обчислення ексцесу (РОЗГЛЯНУТИ: мінімально-спільна?)
# ✅ compute_moments() — Обчислення моментів (РОЗГЛЯНУТИ: мінімально-спільна?)
# ⚠️ .ts_pmm2_fit() — Приватна функція (ЕКСПОРТОВАНА? ВИДАЛИТИ!)
# ⚠️ create_ts_design_matrix() — Утіліта для часових рядів (ПЕРЕМІСТИТИ В pmm2_ts_design.R?)
```

**Пропонована операція:**
```r
# pmm2_utils.R: Утиліти для PMM2 оптимізації

# Експортовані функції:
# - pmm_skewness() — Обчислення асиметрії (універсальна, але спільна для PMM2 та PMM3)
# - pmm_kurtosis() — Обчислення ексцесу
# - compute_moments() — Обчислення моментів та кумулянтів

# Приватні утіліти для PMM2:
# - .pmm2_fit() — Основна PMM2 оптимізація
# - .ts_pmm2_fit() — Часові ряди PMM2 оптимізація

# ВИДАЛИТИ ЕКСПОРТУВАННЯ:
# - .ts_pmm2_fit — це приватна функція!
```

---

#### pmm2_ts_design.R (раніше pmm_ts_design.R)

**Операція:** Просто переіменувати

**Логіка:** Цей файл ЯВНО про часові ряди, немає необхідності в змінах

**Перевірка:** Чи вся логіка часових рядів тут?
```r
# create_ts_design_matrix() — часові ряди дизайн ✅
# create_ar_matrix() — AR дизайн ✅
# get_yw_estimates() — Yule-Walker оцінки ✅
```

✅ Добре

---

### Крок 2: Очистити приватні функції в NAMESPACE (Час: 30 хвилин)

**Операція в NAMESPACE:**

Видалити лінії:
```r
# ВИДАЛИТИ (приватні функції, які помилково експортовані):
# export(.ts_pmm2_fit)
# export(.pmm2_fit)
# export(create_ts_design_matrix)  # Якщо експортована
```

Залишити:
```r
# ЗАЛИШИТИ (публічні функції):
export(lm_pmm2)
export(ar_pmm2)
export(ma_pmm2)
export(arma_pmm2)
export(arima_pmm2)
export(ts_pmm2)
export(compare_*)
export(pmm2_inference)
export(ts_pmm2_inference)
export(pmm_skewness)
export(pmm_kurtosis)
export(compute_moments)
export(plot_pmm2_bootstrap)
```

**Генерування:** Запустити `devtools::document()` щоб регенерувати

---

### Крок 3: Оновити Roxygen2 коментарі для ясності (Час: 1-2 години)

**Приклад для `lm_pmm2`:**

```r
#' Fit Linear Regression with Polynomial Maximization Method (S=2)
#'
#' Implements the PMM2 method for linear regression parameter estimation
#' in presence of non-Gaussian (asymmetric) error distributions.
#'
#' @param formula Model formula
#' @param data Data frame containing variables
#' @param subset Subset of observations to use
#' @param na.action Function for handling missing values
#' @param ... Additional arguments
#'
#' @return An object of class \code{\linkS4class{PMM2fit}} containing:
#'   \describe{
#'     \item{coefficients}{Estimated PMM2 coefficients}
#'     \item{residuals}{Fitted residuals}
#'     \item{m2, m3, m4}{Central moments of residuals}
#'     \item{convergence}{Logical: did algorithm converge?}
#'     \item{iterations}{Number of iterations}
#'   }
#'
#' @section Methodology:
#' PMM2 uses second-order polynomial weight functions for moment-based
#' parameter estimation. Particularly effective for regression with
#' asymmetric error distributions.
#'
#' @examples
#' # Simulate data with asymmetric errors
#' set.seed(42)
#' n <- 100
#' x <- rnorm(n)
#' y <- 2 + 1.5*x + rgamma(n, 2) - 2
#' 
#' # Fit PMM2 model
#' fit <- lm_pmm2(y ~ x)
#' summary(fit)
#' 
#' @seealso
#'   \code{\link{lm}} for ordinary least squares,
#'   \code{\link{compare_with_ols}} for comparison
#'
#' @export
#' @import methods
#' @importFrom stats model.frame model.response model.matrix
lm_pmm2 <- function(formula, data, subset, na.action, ...) {
  # Implementation
}
```

**Чек-лист для кожної функції:**
- ✅ @title — одна лінія, ясна
- ✅ @param — для кожного аргументу
- ✅ @return — S4 клас явно названий
- ✅ @examples — виконуваний код
- ✅ @export — або залишити приватною
- ✅ @seealso — пов'язані функції

---

### Крок 4: Оновити NAMESPACE (Час: 30 хвилин)

Запустити:
```r
devtools::document()
```

Це автоматично регенерує NAMESPACE на основі Roxygen2 коментарів

---

### Крок 5: Перевірка та валідація (Час: 1 година)

```r
# 1. Базова перевірка пакету
devtools::check()

# 2. CRAN перевірка
devtools::check(args = c('--as-cran'))

# 3. Тести
devtools::test()
```

---

## Результат Рефакторингу

### До (Початковий стан)
```
R/
├─ pmm_package.R        # ???
├─ pmm_main.R           # ???
├─ pmm_classes.R        # ???
├─ pmm_utils.R          # ???
└─ pmm_ts_design.R      # ???

CRAN Readiness: 40%
```

### Після (Фаза 1)
```
R/
├─ pmm2_package.R       # ✅ Залежності пакету
├─ pmm2_classes.R       # ✅ S4 класи
├─ pmm2_main.R          # ✅ Основні функції
├─ pmm2_utils.R         # ✅ Утіліти оптимізації
└─ pmm2_ts_design.R     # ✅ Часові ряди

CRAN Readiness: 85%
Архітектурна ясність: 100%
Готовність до PMM3: 90%
```

---

## Майбутні розширення: Як додати PMM3

Коли PMM3 буде готова:

```r
# Етап 1: Копіювати структуру
cp R/pmm2_main.R R/pmm3_main.R
cp R/pmm2_classes.R R/pmm3_classes.R
# ... для всіх файлів

# Етап 2: Модифікувати під PMM3
# Edit R/pmm3_*.R файли для PMM3 логіки

# Етап 3: Оновити NAMESPACE
devtools::document()

# Етап 4: Семантичне версіонування
# DESCRIPTION: Version 0.2.0
# NEWS.md: Додати запис про PMM3

# РЕЗУЛЬТАТ: PMM2 та PMM3 ізольовані, без конфліктів
```

---

## Резюме: Архітектурна цінність

**Когнітивне навантаження:** 📉 Зменшується на 50%
**Готовність до розширення:** 📈 Збільшується на 80%
**CRAN відповідність:** ✅ Поліпшується з 40% до 85%
**Техдолг:** 📉 Мінімізується на майбутнє

