# EstemPMM Phase 1: Implementation Playbook
## Крок за кроком інструкція до CRAN готовності

---

## ЗАГАЛЬНИЙ ОГЛЯД

Цей playbook провадить вас через трансформацію EstemPMM з дослідницького коду до production-ready пакету за 3-4 дні роботи. Кожен крок містить:

- 📋 Конкретні операції
- ✅ Перевірочні пункти
- 🔍 Діагностика проблем
- 📊 Очікувані результати

**Часова розподіл:**
- День 1 (4-5 годин): Рефакторинг архітектури, переіменування файлів
- День 2 (3-4 години): Додання тестів, NEWS.md
- День 3 (2-3 години): CRAN проверка та корректування
- День 4 (1-2 години): Фінальна валідація

---

## ДЕНЬ 1: АРХІТЕКТУРНА ТРАНСФОРМАЦІЯ

### Фаза 1.1: Підготовка робочого середовища

**Операція:** Встановити необхідні пакети

```r
# У R консолі:
install.packages(c("devtools", "roxygen2", "testthat", "rmarkdown"))

# Перевірити інсталяцію
library(devtools)
library(roxygen2)
library(testthat)

# Поточна версія розробки повинна бути >=7.3.0
packageVersion("roxygen2")
```

**Перевірка:** ✅ Всі пакети встановлені

---

### Фаза 1.2: Файловий рефакторинг (20-30 хвилин)

**Операція:** Переіменування R/ файлів для явної PMM2-ідентифікації

```bash
# У терміналі (або R: file.rename())

# 1. Перейти до директорії проекту
cd /path/to/EstemPMM

# 2. Переіменувати файли
mv R/pmm_package.R R/pmm2_package.R
mv R/pmm_main.R R/pmm2_main.R
mv R/pmm_classes.R R/pmm2_classes.R
mv R/pmm_utils.R R/pmm2_utils.R
mv R/pmm_ts_design.R R/pmm2_ts_design.R

# 3. Перевірити результат
ls R/pmm2_*.R
# Повинні вивести: pmm2_classes.R pmm2_main.R pmm2_package.R ...
```

**Альтернатива (з R):**

```r
# У R консолі
setwd("/path/to/EstemPMM")

# Переіменування
file.rename("R/pmm_package.R", "R/pmm2_package.R")
file.rename("R/pmm_main.R", "R/pmm2_main.R")
file.rename("R/pmm_classes.R", "R/pmm2_classes.R")
file.rename("R/pmm_utils.R", "R/pmm2_utils.R")
file.rename("R/pmm_ts_design.R", "R/pmm2_ts_design.R")

# Перевірка
list.files("R/", pattern = "pmm2_")
```

**Перевірка:** ✅ Всі 5 файлів переіменовані

---

### Фаза 1.3: Оновлення внутрішніх посилань в коді (45-60 хвилин)

**Операція:** Оновити cross-file посилання у файлах

Файл `R/pmm2_package.R` повинен розглядати як точку входу для залежностей. Змісту залишається незмінним, тільки переіменування файлу.

**Перевірка в кожному файлі:**

```r
# 1. Відкрити R/pmm2_package.R
# Перевірити, що початок файлу містить:
#' @importFrom methods is new slotNames
#' @importFrom graphics abline hist legend lines par
#' ... і т.д.

# 2. Перевірити R/pmm2_classes.R
# Повинна задовольняти все S4 setClass()

# 3. Перевірити R/pmm2_main.R
# Основні функції: lm_pmm2, ts_pmm2, ar_pmm2, ma_pmm2, arma_pmm2, arima_pmm2

# 4. Перевірити R/pmm2_utils.R
# Внутрішні функції оптимізації

# 5. Перевірити R/pmm2_ts_design.R
# Часові ряди матриці дизайну
```

**Діагностика:** Якщо помилки при завантаженні:

```r
# Спробувати завантажити пакет
devtools::load_all()

# Якщо помилка типу "невизначена змінна", перевірити:
# - Чи функція виділена в правильному файлі?
# - Чи S4 класи налаштовані в R/pmm2_classes.R?
```

**Перевірка:** ✅ `devtools::load_all()` проходить без помилок

---

### Фаза 1.4: Видалення приватних експортів (15 хвилин)

**Операція:** Очистити NAMESPACE від приватних функцій

**Крок 1:** Визначити приватні функції

```r
# Приватні функції ЗАВЖДИ починаються з крапки (.)
# У R/pmm2_utils.R та R/pmm2_ts_design.R шукати:

# .ts_pmm2_fit()      ← ПРИВАТНА, не експортується
# .pmm2_fit()         ← ПРИВАТНА, не експортується
# create_ts_design_matrix()  ← ПЕРЕВІРИТИ, чи експортується?
```

**Крок 2:** Перевірити поточний NAMESPACE

```r
# Щоб побачити поточні експорти:
readLines("NAMESPACE")

# Шукати лінії:
# export(.)  → усі експорти
# exportClasses(  → S4 класи
# exportMethods(  → S4 методи
```

**Крок 3:** Оновити Roxygen2 коментарі

**ВИДАЛИТИ з Roxygen2 коментарів:**

```r
# У R/pmm2_utils.R, для функції .ts_pmm2_fit():

#' PMM2 fitting algorithm for time series models  ← ВИДАЛИТИ
#' @export  ← ВИДАЛИТИ (це приватна функція!)
#' @keywords internal  ← ЗАЛИШИТИ (інша позначка для приватної)

# ПРАВИЛЬНО:
#' @keywords internal
#' @noRd  ← Додати (без документації)

.ts_pmm2_fit <- function(...) {
  # implementation
}
```

**Крок 4:** Регенерувати NAMESPACE

```r
# У R консолі
devtools::document()

# Перевірити, що NAMESPACE оновлена
# Запустити grep(), щоб перевірити відсутність приватних експортів:
grepl("export\\(\\..*\\)", readLines("NAMESPACE"))
# Повинна повернути FALSE
```

**Перевірка:** ✅ NAMESPACE не містить експортів, які починаються з крапки

---

### Фаза 1.5: Додання явних коментарів до файлів (30-45 хвилин)

**Операція:** Додати коментарі на початок кожного R/ файлу для ясності

**Файл: R/pmm2_package.R**

```r
# ============================================================================
# EstemPMM: pmm2_package.R
# Залежності пакету та основна документація PMM2
# 
# Функціональність:
# - Імпорти всіх зовнішніх пакетів
# - Документація пакету EstemPMM
# - Визначення базової документації для S4 класів
# ============================================================================

# Залежності та імпорти
#' @importFrom methods is new slotNames
#' ... (залишити існуючий код)
```

**Файл: R/pmm2_classes.R**

```r
# ============================================================================
# EstemPMM: pmm2_classes.R
# Визначення S4 класів для PMM2 методу
#
# Експортовані класи:
# - PMM2fit — результати лінійної регресії PMM2
# - TS2fit — базовий клас для часових рядів
# - ARPMM2, MAPMM2, ARMAPMM2, ARIMAPMM2 — спеціалізовані класи часових рядів
# 
# Призначення:
# - Контейнери для збереження результатів оцінювання
# - Методи для summary, plot, predict, coef, residuals
# ============================================================================

# Визначення S4 класів
# ... (залишити існуючий код)
```

**Файл: R/pmm2_main.R**

```r
# ============================================================================
# EstemPMM: pmm2_main.R
# Основні функції PMM2 (Polynomial Maximization Method, S=2)
#
# Експортовані функції лінійної регресії:
# - lm_pmm2() — підгонка лінійної моделі
# - compare_with_ols() — порівняння з OLS
#
# Експортовані функції часових рядів (через ts_pmm2):
# - ts_pmm2() — дисетчер для AR/MA/ARMA/ARIMA
# - ar_pmm2() — AR моделі
# - ma_pmm2() — MA моделі
# - arma_pmm2() — ARMA моделі
# - arima_pmm2() — ARIMA моделі
# - compare_ts_methods() — порівняння з методами 'arima'
# ============================================================================

# Основна реалізація лінійної регресії та часових рядів
# ... (залишити існуючий код)
```

**Файл: R/pmm2_utils.R**

```r
# ============================================================================
# EstemPMM: pmm2_utils.R
# Утиліти для PMM2 оптимізації та статистичних обчислень
#
# Експортовані публічні утиліти:
# - pmm_skewness() — обчислення коефіцієнта асиметрії
# - pmm_kurtosis() — обчислення коефіцієнта ексцесу
# - compute_moments() — обчислення центральних моментів та кумулянтів
#
# Приватні утіліти (для внутрішнього використання):
# - .pmm2_fit() — основний алгоритм оптимізації PMM2
# - .ts_pmm2_fit() — оптимізація для часових рядів
#
# Призначення:
# - Момент-базовані оцінки параметрів
# - Числова стійкість та конвергенція
# ============================================================================

# Реалізація утіліт оптимізації
# ... (залишити існуючий код)
```

**Файл: R/pmm2_ts_design.R**

```r
# ============================================================================
# EstemPMM: pmm2_ts_design.R
# Конструкція дизайн-матриць для часових рядів PMM2
#
# Приватні утиліти для конструкції матриць:
# - create_ts_design_matrix() — універсальна конструкція для AR/MA/ARMA/ARIMA
# - create_ar_matrix() — конструкція для AR моделей
# - get_yw_estimates() — оцінки Юла-Волкера для AR параметрів
# - prepare_ts_data() — попередня обробка часових рядів
#
# Призначення:
# - Підтримання PMM2 оцінювання для часових рядів
# - Управління лаговими структурами та інноваціями
# ============================================================================

# Конструкція матриць часових рядів
# ... (залишити існуючий код)
```

**Перевірка:** ✅ Кожен файл має явний коментар про його призначення

---

## ДЕНЬ 2: ДОКУМЕНТАЦІЯ ТА ТЕСТУВАННЯ

### Фаза 2.1: Додання NEWS.md (30 хвилин)

**Операція:** Створити файл `NEWS.md` у кореневій директорії проекту

```r
# У R консолі:
news_content <- readLines("~/EstemPMM_NEWS.md_PHASE1")  # Завантажити з артефакту
writeLines(news_content, "NEWS.md")

# Перевірити
file.exists("NEWS.md")
```

**Перевірка:** ✅ `NEWS.md` знаходиться у кореневій директорії

---

### Фаза 2.2: Оновлення DESCRIPTION (20 хвилин)

**Операція:** Оновити DESCRIPTION з новими полями

```r
# У R консолі:
# 1. Прочитати оновлену версію
new_description <- readLines("~/DESCRIPTION_PHASE1")

# 2. Записати в проект
writeLines(new_description, "DESCRIPTION")

# 3. Перевірити
readLines("DESCRIPTION")
```

**Перевірка:** ✅ DESCRIPTION містить поля:
- ✅ URL
- ✅ BugReports
- ✅ VignetteBuilder

---

### Фаза 2.3: Додання тестів (2-3 години)

**Операція:** Створити файл тестів

```r
# У R консолі:
# 1. Створити директорію, якщо не існує
dir.create("tests/testthat", recursive = TRUE, showWarnings = FALSE)

# 2. Завантажити тест-файл
test_content <- readLines("~/test-pmm2_comprehensive.R_PHASE1")

# 3. Розділити на окремі файли для кращої організації

# 3a. Linear regression tests
linear_tests <- test_content[grep("context\\(\"PMM2 Linear", test_content):
                                  (grep("context\\(\"PMM2 Time", test_content)[1]-1)]
writeLines(linear_tests, "tests/testthat/test-pmm2_linear.R")

# 3b. Time series tests
ts_tests <- test_content[grep("context\\(\"PMM2 Time", test_content):
                              (grep("context\\(\"PMM2 Inference", test_content)[1]-1)]
writeLines(ts_tests, "tests/testthat/test-pmm2_ts.R")

# 3c. Inference tests
inference_tests <- test_content[grep("context\\(\"PMM2 Inference", test_content):
                                     (grep("context\\(\"PMM2 Utility", test_content)[1]-1)]
writeLines(inference_tests, "tests/testthat/test-pmm2_inference.R")

# 3d. Utilities tests
utils_tests <- test_content[grep("context\\(\"PMM2 Utility", test_content):
                                 (grep("context\\(\"PMM2 S4 Methods", test_content)[1]-1)]
writeLines(utils_tests, "tests/testthat/test-pmm2_utils.R")

# 3e. Methods tests
methods_tests <- test_content[grep("context\\(\"PMM2 S4 Methods", test_content):length(test_content)]
writeLines(methods_tests, "tests/testthat/test-pmm2_methods.R")

# 4. Перевірити структуру
list.files("tests/testthat/", pattern = "\\.R$")
```

**Перевірка:** ✅ Тести створені в директорії

---

### Фаза 2.4: Запуск тестів локально (45 хвилин)

**Операція:** Перевірити, що тести проходять

```r
# У R консолі, в директорії проекту
devtools::load_all()

# 1. Запустити всі тести
devtools::test()

# Очікуваний результат:
# ✓ test-pmm2_linear.R ... XX tests
# ✓ test-pmm2_ts.R ... XX tests
# ... і т.д.

# 2. Якщо тести не проходять, діагностика:

# 2a. Завантажити окремий тест-файл
testthat::test_file("tests/testthat/test-pmm2_linear.R")

# 2b. Запустити конкретний тест
testthat::test_that("lm_pmm2 produces valid fit object", {
  # тест код
})

# 2c. Перевірити помилку
# Якщо помилка, виправити функцію в R/pmm2_*.R
```

**Діагностика:**

Якщо тест не проходить:
1. Перевірити, що пакет успішно завантажується: `devtools::load_all()`
2. Запустити конкретний тест в інтерактивному режимі
3. Виправити вихідний код або тест

**Перевірка:** ✅ Мінімум 80% тестів проходить

---

## ДЕНЬ 3: CRAN ПЕРЕВІРКА

### Фаза 3.1: Локальна CRAN перевірка (30 хвилин)

**Операція:** Запустити CRAN-стиль перевірку

```r
# У R консолі
devtools::check()

# Очікувані результати:
# 0 errors ✓
# 0 warnings ✓
# 0 notes або мінімум examples-that-can't-be-run

# Якщо помилки, вивести детально:
devtools::check(args = "--verbose")

# Для CRAN-специфічних вимог:
devtools::check(args = c("--as-cran"))
```

**Діагностика:**

| Помилка | Рішення |
|---|---|
| `ERROR: object 'function_name' not found` | Перевірити, чи функція експортована в NAMESPACE |
| `WARNING: package startup message` | Видалити `cat()` або `message()` на верхньому рівні |
| `NOTE: Unexported objects` | Додати `@export` в Roxygen2 коментарі або видалити |
| `NOTE: no visible global function definition` | Перевірити `@importFrom` коментарі |

**Перевірка:** ✅ `devtools::check()` проходить без помилок

---

### Фаза 3.2: Документація перевірка (15 хвилин)

**Операція:** Перевірити документацію

```r
# 1. Регенерувати документацію
devtools::document()

# 2. Перевірити, що всі експортовані функції задокументовані
devtools::check_man()

# 3. Перевірити вінєтки (якщо вони присутні)
devtools::build_vignettes()
```

**Перевірка:** ✅ Немає непокритих функцій

---

### Фаза 3.3: Залежності перевірка (15 хвилин)

**Операція:** Перевірити та оптимізувати залежності

```r
# 1. Перевірити поточні залежності
deps <- devtools::package_dependencies("EstemPMM", 
                                       recursive = FALSE, 
                                       which = c("Imports", "Depends"))
unlist(deps)

# 2. Видалити невикористовувані
# Якщо MASS не використовується, перемістити в Suggests

# 3. Перевірити, що всі залежності правильно оголошені в DESCRIPTION
# devtools::use_package() для додавання залежностей

# ПРИКЛАД: MASS перемістити в Suggests
# Editorially: DESCRIPTION
# FROM: Imports: ... MASS
# TO: Suggests: ... MASS
```

**Перевірка:** ✅ DESCRIPTION містить тільки необхідні залежності

---

## ДЕНЬ 4: ФІНАЛЬНА ВАЛІДАЦІЯ

### Фаза 4.1: Build та Installation (20 хвилин)

**Операція:** Перевірити, що пакет можна встановити

```r
# 1. Очистити попередні версії
remove.packages("EstemPMM")

# 2. Побудувати пакет
devtools::build()
# Результат: EstemPMM_0.1.0.tar.gz

# 3. Встановити з архіву
install.packages("EstemPMM_0.1.0.tar.gz", repos = NULL, type = "source")

# 4. Перевірити інсталяцію
library(EstemPMM)
ls("package:EstemPMM")  # Повинні вивести експортовані функції
```

**Перевірка:** ✅ Пакет успішно встановлює та завантажує

---

### Фаза 4.2: Функціональна перевірка (30 хвилин)

**Операція:** Запустити демонстраційні приклади

```r
# 1. Завантажити пакет
library(EstemPMM)

# 2. Запустити приклади з документації
example("lm_pmm2")
example("ar_pmm2")

# 3. Запустити демонстраційні приклади
set.seed(42)
n <- 100
x <- rnorm(n)
errors <- rgamma(n, shape = 2, scale = 1) - 2
y <- 2 + 1.5 * x + errors
data <- data.frame(x = x, y = y)

fit <- lm_pmm2(y ~ x, data = data)
summary(fit)
plot(fit)
coef(fit)
predict(fit, newdata = data.frame(x = c(0, 1, -1)))

# 4. Запустити тести
testthat::test_package("EstemPMM")
```

**Перевірка:** ✅ Всі приклади виконуються без помилок

---

### Фаза 4.3: CRAN мітка готовності (15 хвилин)

**Операція:** Остаточна перевірка

```r
# Запустити остаточну CRAN-перевірку
result <- devtools::check(args = c("--as-cran"))

# Перевірити результат
# Status: ✓ (0 errors, 0 warnings)

# Якщо все добре:
# ✓ EstemPMM_0.1.0.tar.gz — готова до submission
```

**Перевірка:** ✅ Пакет CRAN-готовий

---

## РІШЕННЯ ТИПОВИХ ПРОБЛЕМ

### Проблема 1: Функція не знаходиться після документації

**Симптом:**
```
ERROR: object 'lm_pmm2' not found
```

**Рішення:**
```r
# 1. Перевірити, чи функція має @export
grep("@export", readLines("R/pmm2_main.R"))

# 2. Регенерувати документацію
devtools::document()

# 3. Перезавантажити пакет
devtools::load_all()
```

---

### Проблема 2: NAMESPACE помилки

**Симптом:**
```
ERROR: package 'stats' is required but not imported
```

**Рішення:**
```r
# Додати до R/pmm2_main.R:
#' @importFrom stats lm model.frame
function_that_uses_lm <- function(...) {
  ...
}
```

---

### Проблема 3: Тести не запускаються

**Симптом:**
```
Error: object 'test_data' not found
```

**Рішення:**
```r
# У test-файлі, замість глобальної переменної:
# НЕПРАВИЛЬНО:
# test_data <- data.frame(...)
# test_that("...", { ... })

# ПРАВИЛЬНО:
test_that("...", {
  test_data <- data.frame(...)  # Визначити в тесті
  fit <- lm_pmm2(y ~ x, data = test_data)
  ...
})
```

---

## КОНТРОЛЬНИЙ СПИСОК ЗАВЕРШЕННЯ ФАЗИ 1

```
День 1: Архітектура
□ R/ файли переіменовані (pmm_* → pmm2_*)
□ Приватні функції не експортуються
□ devtools::load_all() проходить без помилок
□ NAMESPACE оновлена (devtools::document())

День 2: Документація та тести
□ NEWS.md створена в кореневій директорії
□ DESCRIPTION оновлена (URL, BugReports)
□ tests/testthat/ директорія створена з 5 файлів
□ devtools::test() проходить мінімум на 80%

День 3: CRAN перевірка
□ devtools::check() — 0 errors, 0 warnings
□ devtools::check(args = "--as-cran") — PASS
□ Всі залежності в DESCRIPTION правильні
□ Вінєтки (якщо є) будуються без помилок

День 4: Фінальна валідація
□ Пакет встановлює та завантажує успішно
□ Приклади з документації виконуються
□ devtools::test() — всі тести проходять
□ Пакет готовий до CRAN submission
```

---

## НАСТУПНІ КРОКИ

По завершенню Фази 1, ви готові до:

1. **CRAN Submission** — Завантажити пакет на https://cran.r-project.org/submit.html
2. **Фаза 2** — Додати вінєтки та рефакторинг класів для PMM3-готовності
3. **Фаза 3** — CI/CD та GitHub Actions

Очікуваний CRAN Readiness після Фази 1: **85-90%**

