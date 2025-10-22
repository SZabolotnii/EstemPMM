# 🚀 PMM3 Implementation Roadmap for EstemPMM
## Метод Максимізації Поліномів третього порядку (S=3)

**Дата створення:** 22 жовтня 2025
**Статус:** Draft для обговорення
**Призначення:** Розширення функціоналу для симетричних негаусових розподілів

---

## 📋 Executive Summary

### Ключові відмінності PMM2 vs PMM3

| Характеристика | PMM2 (S=2) | PMM3 (S=3) |
|----------------|------------|------------|
| **Призначення** | **Асиметричні** розподіли | **Симетричні** розподіли |
| **Порядок моментів** | До 4-го (μ₂, μ₄) | До 6-го (μ₂, μ₄, μ₆) |
| **Кумулянтні коефіцієнти** | γ₃ (асиметрія), γ₄ | γ₄ (ексцес), γ₆ |
| **Найкраща ефективність** | Skewed distributions | Flat-topped (platykurtic) |
| **Обчислювальна складність** | Середня | Вища (нелінійні рівняння) |
| **Розв'язання системи** | Аналітичне (при S=2) | **Чисельне** (Newton-Raphson) |

### Чому PMM3 важливий?

✅ **Для EPD з β>2** (flat-topped): PMM3 **значно точніший** за OLS
✅ **Uniform, triangular, trapezoidal** розподіли: зменшення дисперсії до **5x**
✅ **Симетричні non-Gaussian** дані: універсальний метод
❌ **НЕ використовувати** для асиметричних розподілів (γ₃≠0)

---

## 🎯 Математична Основа PMM3

### 1. Модель лінійної регресії

```
yᵥ = f(θ, Xᵥ) + ξᵥ,  v = 1,…,N
```

де:
- f(θ, Xᵥ) = a₀ + a₁x₁ᵥ + … + aₚxₚᵥ — детермінована компонента
- ξᵥ — випадкова компонента з **симетричним** розподілом
- **Важливо**: E{ξᵥ} = 0, γ₃ = 0 (симетрія!)

### 2. Стохастичний поліном третього порядку

Система рівнянь для знаходження оцінок (для p-ї компоненти θ):

```
Σᵥ₌₁ᴺ [k₁ᵥ⁽ᵖ⁾(yᵥ - αᵥ) + k₂ᵥ⁽ᵖ⁾(yᵥ² - α₂ᵥ) + k₃ᵥ⁽ᵖ⁾(yᵥ³ - α₃ᵥ)] = 0
```

де αᵢᵥ = E{yᵥⁱ} — теоретичні початкові моменти.

### 3. Оптимальні коефіцієнти (формула 7 з статті)

```
k₁ᵥ⁽ᵖ⁾ = (3[f(θ,Xᵥ)]²(μ₄ - 3μ₂²) + 3μ₄μ₂ - μ₆) / Δ₃ · ∂f/∂aₚ

k₂ᵥ⁽ᵖ⁾ = (-3f(θ,Xᵥ)(μ₄ - 3μ₂)) / Δ₃ · ∂f/∂aₚ

k₃ᵥ⁽ᵖ⁾ = (μ₄ - 3μ₂) / Δ₃ · ∂f/∂aₚ
```

де **Δ₃ = μ₂²(μ₄² - μ₂μ₆)**

### 4. Система для знаходження оцінок (формула 8)

Після підстановки коефіцієнтів:

```
Σᵥ₌₁ᴺ xᵥ⁽ᵖ⁻¹⁾ [A(a₀ + a₁xᵥ)³ + B(a₀ + a₁xᵥ)² + C(a₀ + a₁xᵥ) + D] = 0
```

де:
- A = 1
- B = -3ŷ₁
- C = 3ŷ₂ - (μ₆-3μ₄μ₂)/(μ₄-3μ₂²)
- D = ŷ₁(μ₆-3μ₄μ₂)/(μ₄-3μ₂²) - ŷ₃

ŷᵢ = (1/N)Σyᵥⁱ — вибіркові початкові моменти

### 5. Коефіцієнт зменшення дисперсії (формула 12)

**Теоретична ефективність PMM3 відносно OLS:**

```
g₍θₚ₎₃ = σ²₍PMM3₎ / σ²₍OLS₎ = 1 - γ₄² / (6 + 9γ₄ + γ₆)
```

де:
- γ₄ = μ₄/μ₂² - 3 — коефіцієнт ексцесу
- γ₆ = μ₆/μ₂³ - 15μ₄/μ₂² + 30 — кумулянтний коефіцієнт 6-го порядку

**Інтерпретація:**
- g < 1: PMM3 **точніший** за OLS
- g = 1: однакова ефективність (гаусові дані)
- Для flat-topped (γ₄ < 0): g може бути 0.2-0.5 (покращення у 2-5 разів!)

---

## 🏗️ Архітектурна Стратегія Інтеграції

### Phase 2 Architecture (згідно з документацією)

```
R/
├─ base_classes.R          [НОВИЙ] Базові S4 класи
│  ├─ setClass("BasePMM")
│  └─ setClass("BaseTS")
│
├─ pmm_common_utils.R      [НОВИЙ] Спільні утиліти
│  ├─ compute_moments()
│  ├─ pmm_skewness()
│  └─ pmm_kurtosis()
│
├─ pmm2_package.R          [БЕЗ ЗМІН] PMM2 документація
├─ pmm2_classes.R          [МОДИФІКАЦІЯ] PMM2 S4 (наслідування від base)
├─ pmm2_main.R             [БЕЗ ЗМІН] PMM2 функції
├─ pmm2_utils.R            [РЕФАКТОРИНГ] PMM2-специфічні утиліти
├─ pmm2_ts_design.R        [БЕЗ ЗМІН] PMM2 дизайн матриці
│
├─ pmm3_package.R          [НОВИЙ] PMM3 документація
├─ pmm3_classes.R          [НОВИЙ] PMM3 S4 класи
│  ├─ setClass("PMM3fit", contains = "BasePMM")
│  └─ setClass("TS3fit", contains = "BaseTS")
├─ pmm3_main.R             [НОВИЙ] PMM3 функції
│  ├─ lm_pmm3()
│  └─ ts_pmm3()
└─ pmm3_utils.R            [НОВИЙ] PMM3 утиліти
   ├─ .pmm3_fit()
   └─ .ts_pmm3_fit()
```

### Ключові принципи

1. **Isolation Principle**: PMM2 і PMM3 повністю ізольовані
2. **Minimal Common Ground**: Тільки базові математичні функції у `pmm_common_utils.R`
3. **Inheritance**: PMM3fit extends BasePMM
4. **Test Isolation**: Окремі тести для PMM3

---

## 📅 Implementation Roadmap

### ФАЗА 1: Базова інфраструктура (Тиждень 1-2)

#### Крок 1.1: Створення базових класів (2-3 години)

**Файл:** `R/base_classes.R`

```r
# ============================================================================
# EstemPMM: base_classes.R
# Базові S4 класи для PMM методів
# ============================================================================

#' Base PMM Fit Class
#'
#' @slot coefficients Numeric vector of estimated parameters
#' @slot residuals Numeric vector of residuals
#' @slot convergence Logical convergence indicator
#' @slot iterations Integer number of iterations
#' @slot call Original function call
#' @slot method Character PMM method used ("PMM2" or "PMM3")
#'
#' @exportClass BasePMM
setClass("BasePMM",
  slots = c(
    coefficients = "numeric",
    residuals = "numeric",
    convergence = "logical",
    iterations = "integer",
    call = "call",
    method = "character"
  )
)

#' Base Time Series Fit Class
#'
#' @slot model_type Character: "AR", "MA", "ARMA", "ARIMA"
#' @slot order Integer vector (p, d, q)
#'
#' @exportClass BaseTS
setClass("BaseTS",
  contains = "BasePMM",
  slots = c(
    model_type = "character",
    order = "integer"
  )
)
```

#### Крок 1.2: PMM3 S4 класи (2 години)

**Файл:** `R/pmm3_classes.R`

```r
# ============================================================================
# EstemPMM: pmm3_classes.R
# S4 класи для PMM3 (S=3, симетричні розподіли)
# ============================================================================

#' PMM3 Fit Class
#'
#' Результати PMM3 estimation для лінійної регресії
#' з симетричними негаусовими помилками
#'
#' @slot m2 Second central moment (variance)
#' @slot m4 Fourth central moment
#' @slot m6 Sixth central moment
#' @slot gamma4 Kurtosis coefficient
#' @slot gamma6 6th order cumulant coefficient
#' @slot g_coefficient Variance reduction coefficient (vs OLS)
#'
#' @exportClass PMM3fit
setClass("PMM3fit",
  contains = "BasePMM",
  slots = c(
    m2 = "numeric",
    m4 = "numeric",
    m6 = "numeric",
    gamma4 = "numeric",
    gamma6 = "numeric",
    g_coefficient = "numeric"
  ),
  prototype = list(
    method = "PMM3"
  )
)

#' PMM3 Time Series Base Class
#'
#' @exportClass TS3fit
setClass("TS3fit",
  contains = c("PMM3fit", "BaseTS")
)

#' PMM3 AR Model Class
#' @exportClass ARPMM3
setClass("ARPMM3",
  contains = "TS3fit",
  prototype = list(
    model_type = "AR"
  )
)

#' PMM3 MA Model Class
#' @exportClass MAPMM3
setClass("MAPMM3",
  contains = "TS3fit",
  prototype = list(
    model_type = "MA"
  )
)

#' PMM3 ARMA Model Class
#' @exportClass ARMAPMM3
setClass("ARMAPMM3",
  contains = "TS3fit",
  prototype = list(
    model_type = "ARMA"
  )
)

#' PMM3 ARIMA Model Class
#' @exportClass ARIMAPMM3
setClass("ARIMAPMM3",
  contains = "TS3fit",
  prototype = list(
    model_type = "ARIMA"
  )
)
```

#### Крок 1.3: Спільні утиліти (3-4 години)

**Файл:** `R/pmm_common_utils.R`

```r
# ============================================================================
# EstemPMM: pmm_common_utils.R
# Спільні утиліти для PMM2 та PMM3
# ============================================================================

#' Compute Central Moments
#'
#' Обчислює центральні моменти до заданого порядку
#'
#' @param x Numeric vector
#' @param orders Integer vector порядків (default: c(2, 3, 4, 6))
#' @return Named list моментів
#'
#' @examples
#' x <- rnorm(100)
#' moments <- compute_moments(x)
#'
#' @export
compute_moments <- function(x, orders = c(2, 3, 4, 6)) {
  n <- length(x)
  x_mean <- mean(x)
  x_centered <- x - x_mean

  moments <- lapply(orders, function(r) {
    sum(x_centered^r) / n
  })
  names(moments) <- paste0("m", orders)

  return(moments)
}

#' PMM Skewness Coefficient
#'
#' Обчислює коефіцієнт асиметрії (γ₃)
#'
#' @param x Numeric vector or moment list
#' @return Numeric skewness coefficient
#'
#' @export
pmm_skewness <- function(x) {
  if (is.list(x)) {
    # x is moments list
    m2 <- x$m2
    m3 <- x$m3
  } else {
    moments <- compute_moments(x, orders = c(2, 3))
    m2 <- moments$m2
    m3 <- moments$m3
  }

  gamma3 <- m3 / (m2^(3/2))
  return(gamma3)
}

#' PMM Kurtosis Coefficient
#'
#' Обчислює коефіцієнт ексцесу (γ₄)
#'
#' @param x Numeric vector or moment list
#' @return Numeric kurtosis coefficient
#'
#' @export
pmm_kurtosis <- function(x) {
  if (is.list(x)) {
    m2 <- x$m2
    m4 <- x$m4
  } else {
    moments <- compute_moments(x, orders = c(2, 4))
    m2 <- moments$m2
    m4 <- moments$m4
  }

  gamma4 <- (m4 / m2^2) - 3
  return(gamma4)
}

#' PMM 6th Order Cumulant Coefficient
#'
#' Обчислює кумулянтний коефіцієнт 6-го порядку (γ₆)
#'
#' @param x Numeric vector or moment list
#' @return Numeric 6th order cumulant coefficient
#'
#' @export
pmm_gamma6 <- function(x) {
  if (is.list(x)) {
    m2 <- x$m2
    m4 <- x$m4
    m6 <- x$m6
  } else {
    moments <- compute_moments(x, orders = c(2, 4, 6))
    m2 <- moments$m2
    m4 <- moments$m4
    m6 <- moments$m6
  }

  gamma6 <- (m6 / m2^3) - 15 * (m4 / m2^2) + 30
  return(gamma6)
}

#' Variance Reduction Coefficient (PMM3 vs OLS)
#'
#' Обчислює теоретичний коефіцієнт зменшення дисперсії
#'
#' @param gamma4 Kurtosis coefficient
#' @param gamma6 6th order cumulant coefficient
#' @return Numeric g coefficient (< 1 means PMM3 is better)
#'
#' @export
pmm3_variance_reduction <- function(gamma4, gamma6) {
  g <- 1 - (gamma4^2) / (6 + 9*gamma4 + gamma6)
  return(g)
}

#' Test for Symmetry
#'
#' Тестує гіпотезу про симетричність розподілу
#'
#' @param x Numeric vector (residuals)
#' @param alpha Significance level (default: 0.05)
#' @return List with test results
#'
#' @export
test_symmetry <- function(x, alpha = 0.05) {
  gamma3 <- pmm_skewness(x)
  n <- length(x)

  # Asymptotic test: sqrt(n) * gamma3 ~ N(0, 6)
  se_gamma3 <- sqrt(6 / n)
  z_stat <- gamma3 / se_gamma3
  p_value <- 2 * pnorm(-abs(z_stat))

  is_symmetric <- p_value > alpha

  return(list(
    gamma3 = gamma3,
    z_statistic = z_stat,
    p_value = p_value,
    is_symmetric = is_symmetric,
    alpha = alpha
  ))
}
```

---

### ФАЗА 2: Основні PMM3 функції (Тиждень 2-3)

#### Крок 2.1: PMM3 лінійна регресія (6-8 годин)

**Файл:** `R/pmm3_main.R`

```r
# ============================================================================
# EstemPMM: pmm3_main.R
# Основні функції PMM3 (Polynomial Maximization Method, S=3)
# Для СИМЕТРИЧНИХ негаусових розподілів
#
# Експортовані функції:
# - lm_pmm3() — Лінійна регресія
# - compare_pmm2_pmm3_ols() — Порівняння методів
#
# Приватні утіліти:
# - .pmm3_fit() — Основний алгоритм оптимізації
# - .pmm3_newton_raphson() — Newton-Raphson solver
# ============================================================================

#' Linear Regression with PMM3
#'
#' Оцінка параметрів лінійної регресії методом максимізації поліномів
#' третього порядку (S=3) для симетричних негаусових помилок
#'
#' @param formula Formula об'єкт (y ~ x1 + x2 + ...)
#' @param data Data frame (optional)
#' @param max_iter Maximum iterations for Newton-Raphson (default: 100)
#' @param tol Convergence tolerance (default: 1e-6)
#' @param use_ols_start Use OLS as starting values (default: TRUE)
#'
#' @return PMM3fit S4 object
#'
#' @examples
#' # Generate data with flat-topped errors (uniform-like)
#' set.seed(123)
#' x <- rnorm(100)
#' # Uniform errors with zero mean
#' errors <- runif(100, -sqrt(3), sqrt(3))
#' y <- 2 + 3*x + errors
#'
#' fit_pmm3 <- lm_pmm3(y ~ x)
#' summary(fit_pmm3)
#'
#' @export
lm_pmm3 <- function(formula, data = NULL,
                    max_iter = 100, tol = 1e-6,
                    use_ols_start = TRUE) {

  # Parse formula
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  # Extract response and design matrix
  y <- model.response(mf, "numeric")
  X <- model.matrix(attr(mf, "terms"), data = mf)
  n <- length(y)
  p <- ncol(X)

  # Step 1: Get OLS starting values
  if (use_ols_start) {
    ols_fit <- lm(formula, data = data)
    theta_init <- coef(ols_fit)
    residuals_ols <- residuals(ols_fit)
  } else {
    theta_init <- rep(0, p)
    residuals_ols <- y - X %*% theta_init
  }

  # Step 2: Compute moments from OLS residuals
  moments <- compute_moments(residuals_ols, orders = c(2, 4, 6))
  m2 <- moments$m2
  m4 <- moments$m4
  m6 <- moments$m6

  # Step 3: Check for symmetry
  symmetry_test <- test_symmetry(residuals_ols)
  if (!symmetry_test$is_symmetric) {
    warning(paste0(
      "Asymmetry detected (gamma3 = ", round(symmetry_test$gamma3, 3),
      "). PMM3 is designed for symmetric distributions. ",
      "Consider using PMM2 (lm_pmm2) instead."
    ))
  }

  # Step 4: Compute cumulant coefficients
  gamma4 <- pmm_kurtosis(moments)
  gamma6 <- pmm_gamma6(moments)

  # Step 5: Check if Gaussian (gamma4 ≈ 0)
  if (abs(gamma4) < 0.1) {
    message("Distribution is close to Gaussian. OLS estimates are optimal.")
    # Return OLS wrapped as PMM3fit
    return(.pmm3_wrap_ols(ols_fit, moments, gamma4, gamma6, cl))
  }

  # Step 6: PMM3 optimization
  result <- .pmm3_fit(y, X, theta_init, m2, m4, m6, max_iter, tol)

  # Step 7: Create PMM3fit object
  residuals <- y - X %*% result$coefficients
  g_coef <- pmm3_variance_reduction(gamma4, gamma6)

  pmm3_obj <- new("PMM3fit",
    coefficients = as.numeric(result$coefficients),
    residuals = as.numeric(residuals),
    convergence = result$convergence,
    iterations = as.integer(result$iterations),
    call = cl,
    method = "PMM3",
    m2 = m2,
    m4 = m4,
    m6 = m6,
    gamma4 = gamma4,
    gamma6 = gamma6,
    g_coefficient = g_coef
  )

  return(pmm3_obj)
}

#' PMM3 Fitting Algorithm (Private)
#'
#' Основний алгоритм для знаходження PMM3 оцінок
#' Використовує Newton-Raphson для розв'язання нелінійної системи
#'
#' @keywords internal
.pmm3_fit <- function(y, X, theta_init, m2, m4, m6, max_iter, tol) {

  n <- length(y)
  p <- ncol(X)
  theta <- theta_init
  converged <- FALSE

  # Precompute constants
  Delta3 <- m2^2 * (m4^2 - m2*m6)
  A <- 1

  for (iter in 1:max_iter) {
    # Current fitted values
    f_theta <- X %*% theta

    # Sample initial moments
    y_bar_1 <- mean(y)
    y_bar_2 <- mean(y^2)
    y_bar_3 <- mean(y^3)

    # Constants B, C, D (formula 8)
    B <- -3 * y_bar_1
    C <- 3*y_bar_2 - (m6 - 3*m4*m2) / (m4 - 3*m2^2)
    D <- y_bar_1 * (m6 - 3*m4*m2) / (m4 - 3*m2^2) - y_bar_3

    # Build system of equations F(theta) = 0
    F <- numeric(p)
    J <- matrix(0, p, p)

    for (q in 1:p) {
      # Equation for q-th parameter
      for (v in 1:n) {
        x_v_q <- X[v, q]
        f_v <- f_theta[v]

        poly_val <- A*f_v^3 + B*f_v^2 + C*f_v + D
        F[q] <- F[q] + x_v_q * poly_val

        # Jacobian
        dpoly_df <- 3*A*f_v^2 + 2*B*f_v + C
        for (r in 1:p) {
          x_v_r <- X[v, r]
          J[q, r] <- J[q, r] + x_v_q * dpoly_df * x_v_r
        }
      }
    }

    # Newton-Raphson update: theta_new = theta - J^(-1) * F
    tryCatch({
      delta_theta <- solve(J, F)
      theta_new <- theta - delta_theta

      # Check convergence
      if (max(abs(delta_theta)) < tol) {
        converged <- TRUE
        theta <- theta_new
        break
      }

      theta <- theta_new
    }, error = function(e) {
      warning(paste("Newton-Raphson failed at iteration", iter, ":", e$message))
      converged <- FALSE
      break
    })
  }

  if (!converged && iter == max_iter) {
    warning(paste("PMM3 did not converge after", max_iter, "iterations"))
  }

  return(list(
    coefficients = theta,
    convergence = converged,
    iterations = iter
  ))
}

#' Wrap OLS as PMM3fit (Private)
#'
#' @keywords internal
.pmm3_wrap_ols <- function(ols_fit, moments, gamma4, gamma6, call) {
  new("PMM3fit",
    coefficients = coef(ols_fit),
    residuals = residuals(ols_fit),
    convergence = TRUE,
    iterations = 0L,
    call = call,
    method = "PMM3 (OLS fallback)",
    m2 = moments$m2,
    m4 = moments$m4,
    m6 = moments$m6,
    gamma4 = gamma4,
    gamma6 = gamma6,
    g_coefficient = 1.0
  )
}
```

#### Крок 2.2: S4 методи для PMM3fit (3-4 години)

**Файл:** `R/pmm3_methods.R`

```r
# ============================================================================
# EstemPMM: pmm3_methods.R
# S4 методи для PMM3fit об'єктів
# ============================================================================

#' Summary Method for PMM3fit
#'
#' @export
setMethod("summary", "PMM3fit", function(object, ...) {
  cat("PMM3 Linear Regression\n")
  cat("=====================\n\n")

  cat("Call:\n")
  print(object@call)
  cat("\n")

  cat("Coefficients:\n")
  print(object@coefficients)
  cat("\n")

  cat("Convergence:", object@convergence, "\n")
  cat("Iterations:", object@iterations, "\n\n")

  cat("Moments and Cumulants:\n")
  cat(sprintf("  m2 (variance): %.4f\n", object@m2))
  cat(sprintf("  m4: %.4f\n", object@m4))
  cat(sprintf("  m6: %.4f\n", object@m6))
  cat(sprintf("  gamma4 (kurtosis): %.4f\n", object@gamma4))
  cat(sprintf("  gamma6: %.4f\n", object@gamma6))
  cat("\n")

  cat(sprintf("Variance Reduction (vs OLS): g = %.4f\n", object@g_coefficient))
  if (object@g_coefficient < 1) {
    improvement <- (1 - object@g_coefficient) * 100
    cat(sprintf("  --> PMM3 is %.1f%% more efficient than OLS\n", improvement))
  } else {
    cat("  --> OLS is more efficient (Gaussian-like distribution)\n")
  }
})

#' Coef Method for PMM3fit
#'
#' @export
setMethod("coef", "PMM3fit", function(object, ...) {
  return(object@coefficients)
})

#' Residuals Method for PMM3fit
#'
#' @export
setMethod("residuals", "PMM3fit", function(object, ...) {
  return(object@residuals)
})

#' Print Method for PMM3fit
#'
#' @export
setMethod("print", "PMM3fit", function(x, ...) {
  cat("PMM3 Fit\n")
  cat("Coefficients:\n")
  print(x@coefficients)
  invisible(x)
})
```

---

### ФАЗА 3: PMM3 для часових рядів (Тиждень 3-4)

#### Крок 3.1: PMM3 Time Series функції (8-10 годин)

**Файл:** `R/pmm3_ts.R`

```r
# ============================================================================
# EstemPMM: pmm3_ts.R
# PMM3 для часових рядів (AR, MA, ARMA, ARIMA)
# ТІЛЬКИ для симетричних інновацій!
# ============================================================================

#' AR Model with PMM3
#'
#' @param x Time series vector
#' @param p AR order
#' @param ... Additional arguments
#'
#' @return ARPMM3 S4 object
#'
#' @export
ar_pmm3 <- function(x, p, ...) {
  # Implementation similar to ar_pmm2 but using S=3 polynomials
  # TODO: Implement
}

#' MA Model with PMM3
#'
#' @export
ma_pmm3 <- function(x, q, ...) {
  # TODO: Implement
}

#' ARMA Model with PMM3
#'
#' @export
arma_pmm3 <- function(x, p, q, ...) {
  # TODO: Implement
}

#' ARIMA Model with PMM3
#'
#' @export
arima_pmm3 <- function(x, p, d, q, ...) {
  # TODO: Implement
}

#' Generic Time Series Dispatcher for PMM3
#'
#' @export
ts_pmm3 <- function(x, model = c("ar", "ma", "arma", "arima"), ...) {
  model <- match.arg(model)

  switch(model,
    ar = ar_pmm3(x, ...),
    ma = ma_pmm3(x, ...),
    arma = arma_pmm3(x, ...),
    arima = arima_pmm3(x, ...)
  )
}
```

---

### ФАЗА 4: Тестування (Тиждень 4-5)

#### Крок 4.1: Unit тести (6-8 годин)

**Файл:** `tests/testthat/test-pmm3-linear.R`

```r
# ============================================================================
# Test PMM3 Linear Regression
# ============================================================================

library(testthat)
library(EstemPMM)

test_that("lm_pmm3 works with uniform errors", {
  set.seed(123)
  n <- 100
  x <- rnorm(n)
  # Uniform errors: E[ε] = 0, Var[ε] = 1
  errors <- runif(n, -sqrt(3), sqrt(3))
  y <- 2 + 3*x + errors

  fit <- lm_pmm3(y ~ x)

  expect_s4_class(fit, "PMM3fit")
  expect_length(coef(fit), 2)
  expect_true(fit@convergence)

  # Check that gamma4 is negative (flat-topped)
  expect_true(fit@gamma4 < 0)

  # Check that PMM3 is more efficient than OLS
  expect_true(fit@g_coefficient < 1)
})

test_that("lm_pmm3 handles trapezoidal errors", {
  set.seed(456)
  n <- 100
  x <- rnorm(n)

  # Trapezoidal distribution (symmetric)
  u <- runif(n)
  b <- 0.5  # shape parameter
  errors <- ifelse(u < 0.5,
                   -sqrt(2*b*u) + 1,
                   sqrt(2*b*(1-u)) - 1)

  y <- 1 + 2*x + errors

  fit <- lm_pmm3(y ~ x)

  expect_s4_class(fit, "PMM3fit")
  expect_true(abs(coef(fit)[2] - 2) < 0.5)
})

test_that("lm_pmm3 gives warning for asymmetric data", {
  set.seed(789)
  n <- 100
  x <- rnorm(n)
  # Asymmetric errors (exponential)
  errors <- rexp(n, rate = 1) - 1
  y <- 2 + 3*x + errors

  expect_warning(
    lm_pmm3(y ~ x),
    "Asymmetry detected"
  )
})

test_that("lm_pmm3 falls back to OLS for Gaussian data", {
  set.seed(101)
  n <- 100
  x <- rnorm(n)
  y <- 2 + 3*x + rnorm(n)

  fit <- lm_pmm3(y ~ x)

  expect_s4_class(fit, "PMM3fit")
  # For Gaussian: g ≈ 1 (no improvement)
  expect_true(abs(fit@g_coefficient - 1) < 0.1)
})

test_that("lm_pmm3 multivariate regression", {
  set.seed(202)
  n <- 80
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  errors <- runif(n, -sqrt(3), sqrt(3))
  y <- 1 + 2*x1 + 3*x2 + errors

  fit <- lm_pmm3(y ~ x1 + x2)

  expect_s4_class(fit, "PMM3fit")
  expect_length(coef(fit), 3)
  expect_true(fit@convergence)
})
```

---

### ФАЗА 5: Документація та Vignettes (Тиждень 5-6)

#### Крок 5.1: Створення vignette (4-6 годин)

**Файл:** `vignettes/04-pmm3-symmetric-distributions.Rmd`

```rmd
---
title: "PMM3: Estimation for Symmetric Non-Gaussian Distributions"
author: "EstemPMM Development Team"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{PMM3: Estimation for Symmetric Non-Gaussian Distributions}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## Introduction

PMM3 (Polynomial Maximization Method with S=3) is designed for **symmetric non-Gaussian** distributions. This vignette demonstrates:

1. When to use PMM3 vs PMM2
2. How PMM3 works
3. Examples with different symmetric distributions
4. Comparison with OLS and PMM2

## When to Use PMM3?

✅ **Use PMM3 when:**
- Errors are **symmetric** (γ₃ ≈ 0)
- Distribution is **non-Gaussian** (γ₄ ≠ 0)
- Especially for **flat-topped** distributions (γ₄ < 0):
  - Uniform
  - Trapezoidal
  - Triangular
  - Arcsine

❌ **Do NOT use PMM3 when:**
- Errors are **asymmetric** (use PMM2 instead)
- Distribution is Gaussian (use OLS)

## Example 1: Uniform Errors

...

## Example 2: Comparison PMM2 vs PMM3

...

## Monte Carlo Study

...
```

---

## 📊 Тестування та Валідація

### Критерії успіху

1. **Функціональні тести**:
   - ✅ lm_pmm3() працює для uniform, trapezoidal, triangular errors
   - ✅ Коректне визначення симетричності
   - ✅ Warning при асиметричних даних
   - ✅ Convergence для різних розмірів вибірки

2. **Ефективність**:
   - ✅ g < 1 для flat-topped розподілів
   - ✅ g ≈ 1 для Gaussian
   - ✅ PMM3 точніший за OLS для platykurtic (γ₄ < 0)

3. **Інтеграція**:
   - ✅ Сумісність з pmm2_*
   - ✅ Ізоляція коду
   - ✅ Спільні утиліти працюють

---

## 🎯 Пріоритети та Timeline

### HIGH PRIORITY (Критично для базового функціоналу)

1. ✅ **base_classes.R** (Тиждень 1)
2. ✅ **pmm_common_utils.R** (Тиждень 1)
3. ✅ **pmm3_classes.R** (Тиждень 1)
4. ✅ **pmm3_main.R: lm_pmm3()** (Тиждень 2-3)
5. ✅ **Unit tests** (Тиждень 4)

### MEDIUM PRIORITY (Розширений функціонал)

6. ✅ **pmm3_ts.R: ar_pmm3(), ma_pmm3()** (Тиждень 3-4)
7. ✅ **pmm3_methods.R: summary, plot** (Тиждень 3)
8. ✅ **Vignette** (Тиждень 5)

### LOW PRIORITY (Оптимізація та додаткові функції)

9. ⚠️ **pmm3_diagnostics.R**
10. ⚠️ **pmm3_bootstrap.R**
11. ⚠️ **Performance optimization**

---

## 🚨 Потенційні Ризики та Mitigation

### Ризик 1: Newton-Raphson не сходиться

**Mitigation:**
- Використовувати OLS як starting point
- Dodati fallback до OLS якщо convergence fails
- Інформувати користувача через warnings

### Ризик 2: Неправильне визначення симетричності

**Mitigation:**
- Статистичний тест на γ₃ = 0
- Clear warning messages
- Документувати, що PMM3 працює тільки для симетричних

### Ризик 3: Конфлікт між PMM2 та PMM3

**Mitigation:**
- Строга ізоляція файлів (pmm2_* vs pmm3_*)
- Базові класи для спільних властивостей
- Окремі тести

---

## 📚 Список літератури

1. **"Estimation of Linear Regression Parameters of Symmetric Non-Gaussian Errors by Polynomial Maximization Method"** (2019) — Основна стаття про PMM3

2. **"Estimating parameters of linear regression with an exponential power distribution..."** (2021) — EPD та PMM3

3. Існуючий код EstemPMM — PMM2 implementation

---

## ✅ Next Steps

### Immediate Actions (Тиждень 1):

1. ✅ Створити `R/base_classes.R`
2. ✅ Створити `R/pmm_common_utils.R`
3. ✅ Створити `R/pmm3_classes.R`
4. ✅ Написати unit tests для base classes
5. ✅ Code review + git commit

### Week 2-3:

1. Імплементувати `lm_pmm3()`
2. Тестування на реальних даних
3. Benchmarking vs OLS/PMM2

---

**Створено:** 22 жовтня 2025
**Статус:** DRAFT
**Версія:** 1.0

**Автори:**
- SZabolotnii (Architecture, PMM3 Algorithm)
- Claude (Documentation, Planning)
