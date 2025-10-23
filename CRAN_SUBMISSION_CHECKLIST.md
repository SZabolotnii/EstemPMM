# Контрольний список для CRAN Submission
# CRAN Submission Checklist

## ✅ Завершено / Completed

### 1. Структура пакету / Package Structure
- ✅ DESCRIPTION файл коректно заповнений
- ✅ LICENSE файл створено (GPL-3)
- ✅ NAMESPACE згенеровано через roxygen2
- ✅ README.md наявний з описом пакету
- ✅ NEWS.md оновлено з історією версій
- ✅ cran-comments.md підготовлено

### 2. Документація / Documentation
- ✅ Всі функції задокументовані англійською мовою
- ✅ man/ файли згенеровані через roxygen2
- ✅ Вінетки (vignettes/) створено:
  - pmm2_introduction.Rmd
  - pmm2_time_series.Rmd
  - bootstrap_inference.Rmd

### 3. Тести / Tests
- ✅ Структура testthat налаштована
- ✅ Тести для основних функцій:
  - test-pmm2_linear.R
  - test-pmm2_ts.R
  - test-pmm2_inference.R
  - test-pmm2_methods.R
  - test-pmm2_utils.R
  - test-monte-carlo.R

### 4. Демонстрації / Demos
- ✅ Demo файли в demo/ директорії
- ✅ demo/00Index присутній

### 5. Залежності / Dependencies
- ✅ Всі пакети з Imports встановлені зі стандартних джерел
- ✅ Suggests містить опціональні пакети

### 6. .Rbuildignore
- ✅ Виключено непотрібні файли:
  - .git, .github
  - docs/
  - cran-comments.md
  - CRAN_CHECK_INSTRUCTIONS.md

---

## 🔄 Кроки для фінальної перевірки / Final Verification Steps

### Крок 1: Встановлення необхідних інструментів
```r
install.packages(c("devtools", "roxygen2", "testthat", "knitr", "rmarkdown"))
```

### Крок 2: Перегенерувати документацію (якщо потрібно)
```r
library(roxygen2)
roxygenize()
```

### Крок 3: Збудувати пакет
```bash
R CMD build .
```

Очікуваний результат: файл `EstemPMM_0.1.1.tar.gz`

### Крок 4: Запустити CRAN перевірку
```bash
R CMD check --as-cran EstemPMM_0.1.1.tar.gz
```

**Мета:** 0 errors, 0 warnings, мінімум notes

**Допустимі NOTES:**
- "New submission" - це нормально для першої подачі
- "Possibly mis-spelled words in DESCRIPTION" - перевірте, чи термінологія правильна

### Крок 5: Запустити додаткові перевірки
```r
library(devtools)

# Локальна перевірка
check(cran = TRUE)

# Перевірка на різних платформах (опціонально)
check_win_devel()  # Windows
check_rhub()       # rhub службa
```

### Крок 6: Запустити тести
```r
library(testthat)
library(EstemPMM)

test_check("EstemPMM")
```

Всі тести мають пройти успішно.

---

## 📋 Чек-лист перед submission

Перевірте наступні пункти:

- [ ] `R CMD check --as-cran` проходить без ERROR та WARNING
- [ ] Всі тести проходять успішно
- [ ] Вінетки компілюються без помилок
- [ ] Версія в DESCRIPTION актуальна (0.1.1)
- [ ] NEWS.md містить інформацію про поточну версію
- [ ] cran-comments.md заповнений з результатами перевірок
- [ ] Немає конфліктів з іншими CRAN пакетами (перевірте назви функцій)
- [ ] Ліцензія GPL-3 коректно вказана

---

## 📤 Submission процес

### 1. Підготовка tarball
```bash
R CMD build .
```

### 2. Остаточна перевірка
```bash
R CMD check --as-cran EstemPMM_0.1.1.tar.gz
```

### 3. Submission на CRAN

**Метод 1: Через веб-форму (рекомендовано для першого submission)**

1. Перейти на https://cran.r-project.org/submit.html
2. Заповнити форму:
   - Ім'я пакету: EstemPMM
   - Версія: 0.1.1
   - Ім'я та email maintainer'a
3. Завантажити файл `EstemPMM_0.1.1.tar.gz`
4. Вставити вміст файлу `cran-comments.md`
5. Підтвердити submission

**Метод 2: Через devtools**
```r
library(devtools)
release()
```

### 4. Після submission

1. Перевірте email для підтвердження submission
2. Підтвердіть submission через email link протягом 24 годин
3. Очікуйте відповіді від CRAN (зазвичай 2-7 днів)

---

## 🔍 Можливі проблеми та рішення

### Проблема 1: WARNING про розмір vignette
**Рішення:** Переконайтесь, що вінетки не генерують великі графіки. Використовуйте параметри chunk:
```r
knitr::opts_chunk$set(
  fig.width = 7,
  fig.height = 5,
  dpi = 72
)
```

### Проблема 2: NOTE про невідомі слова
**Рішення:** Якщо це технічні терміни (PMM, ARIMA, etc.), додайте їх у \code{} або \acronym{} в документації.

### Проблема 3: WARNING про приклади
**Рішення:** Всі приклади повинні виконуватись < 5 секунд. Довгі приклади обгортайте в `\donttest{}` або `\dontrun{}`.

### Проблема 4: Відсутні тести
**Рішення:** CRAN очікує, що важливі функції мають тести. Переконайтесь, що покриття тестами > 80%.

---

## 📚 Корисні ресурси

- **CRAN Repository Policy:** https://cran.r-project.org/web/packages/policies.html
- **Writing R Extensions:** https://cran.r-project.org/doc/manuals/r-release/R-exts.html
- **R Packages book:** https://r-pkgs.org/
- **CRAN Submission Checklist:** https://cran.r-project.org/web/packages/submission_checklist.html

---

## 📞 Контакти та підтримка

Якщо виникають питання під час submission:
- Email CRAN: cran@r-project.org
- GitHub Issues: https://github.com/SZabolotnii/EstemPMM/issues

---

## 🎯 Очікуваний результат

Після успішного submission та review, пакет EstemPMM буде:
1. Опублікований на CRAN
2. Доступний через `install.packages("EstemPMM")`
3. Індексований на офіційній сторінці CRAN
4. Автоматично тестуватиметься на різних платформах

**Час очікування:** 2-7 днів для першого review

---

*Останнє оновлення: 2025-10-23*
*Версія пакету: 0.1.1*
