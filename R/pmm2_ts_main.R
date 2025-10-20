# pmm2_ts_main.R - Уніфікований модуль для моделей часових рядів з PMM2

#' Підігнати модель часового ряду за допомогою методу PMM2
#'
#' @param x Числовий вектор даних часового ряду
#' @param order Специфікація порядку моделі:
#'        - Для AR моделей: одне ціле число (порядок AR)
#'        - Для MA моделей: одне ціле число (порядок MA)
#'        - Для ARMA моделей: вектор c(p, q) (порядки AR та MA)
#'        - Для ARIMA моделей: вектор c(p, d, q) (порядки AR, диференціювання та MA)
#' @param model_type Рядок, що визначає тип моделі: "ar", "ma", "arma", або "arima"
#' @param method Рядок: метод оцінювання, один з "pmm2" (за замовчуванням), "css", "ml", "yw", "ols"
#' @param max_iter Ціле число: максимальна кількість ітерацій для алгоритму
#' @param tol Числове: допуск для збіжності
#' @param include.mean Логічне: чи включати член середнього (перехоплення)
#' @param initial Список або вектор початкових оцінок параметрів (опціонально)
#' @param na.action Функція для обробки відсутніх значень, за замовчуванням - na.fail
#' @param regularize Логічне, додавати малі значення до діагоналі для числової стабільності
#' @param reg_lambda Параметр регуляризації (якщо regularize=TRUE)
#' @param verbose Логічне: чи виводити інформацію про прогрес
#'
#' @details
#' Алгоритм PMM2 працює наступним чином:
#'
#' 1. Підганяє початкову модель за допомогою стандартного методу (OLS, Юла-Волкера, CSS або ML)
#' 2. Обчислює центральні моменти (m2, m3, m4) з початкових залишків/інновацій
#' 3. Використовує ці моменти зі спеціалізованим розв'язувачем (pmm2_algorithm) для знаходження
#'    робастних оцінок параметрів
#'
#' @return Об'єкт S4 \code{TS2fit} відповідного підкласу
#' @export
ts_pmm2 <- function(x, order,
                    model_type = c("ar", "ma", "arma", "arima"),
                    method      = "pmm2",
                    max_iter    = 50,
                    tol         = 1e-6,
                    include.mean= TRUE,
                    initial     = NULL,
                    na.action   = na.fail,
                    regularize  = TRUE,
                    reg_lambda  = 1e-8,
                    verbose     = FALSE) {

  model_type <- match.arg(model_type)
  cl <- match.call()

  if (!is.null(na.action)) {
    x <- na.action(x)
  }

  # 1) Перевірка вхідних даних
  model_params <- validate_ts_parameters(x, order, model_type, include.mean)

  if (model_params$model_type == "ma") {
    q <- model_params$ma_order
    css_fit <- ma_css_fit(model_params$original_x, q, include.mean, verbose)

    if (method == "css") {
      moments <- compute_moments(css_fit$residuals)
      return(new("MAPMM2",
                 coefficients    = as.numeric(css_fit$coefficients),
                 residuals       = as.numeric(css_fit$residuals),
                 m2              = as.numeric(moments$m2),
                 m3              = as.numeric(moments$m3),
                 m4              = as.numeric(moments$m4),
                 convergence     = css_fit$convergence,
                 iterations      = as.numeric(css_fit$iterations),
                 call            = cl,
                 model_type      = "ma",
                 intercept       = if (include.mean) as.numeric(css_fit$intercept) else 0,
                 original_series = as.numeric(model_params$original_x),
                 order           = list(ar = 0L, ma = q, d = 0L)))
    }

    pmm2_fit <- ma_pmm2_fit(model_params$original_x, q, css_fit,
                             max_iter = max_iter, tol = tol,
                             verbose = verbose)
    moments <- compute_moments(pmm2_fit$innovations)

    return(new("MAPMM2",
               coefficients    = as.numeric(pmm2_fit$coefficients),
               residuals       = as.numeric(pmm2_fit$innovations),
               m2              = as.numeric(moments$m2),
               m3              = as.numeric(moments$m3),
               m4              = as.numeric(moments$m4),
               convergence     = pmm2_fit$convergence,
               iterations      = as.numeric(pmm2_fit$iterations),
               call            = cl,
               model_type      = "ma",
               intercept       = if (include.mean) as.numeric(pmm2_fit$intercept) else 0,
               original_series = as.numeric(model_params$original_x),
               order           = list(ar = 0L, ma = q, d = 0L)))
  }

  # 2) Отримання початкових оцінок
  init <- get_initial_estimates(model_params, initial, method, verbose)
  b_init      <- init$b_init
  x_mean      <- init$x_mean
  innovations <- init$innovations
  x_centered  <- init$x_centered
  orig_x      <- init$orig_x
  m2          <- init$m2
  m3          <- init$m3
  m4          <- init$m4

  # 3) Створення дизайн-матриці
  dm <- create_ts_design_matrix(
    x = x_centered,
    model_info = list(
      ar_order = model_params$ar_order,
      ma_order = model_params$ma_order,
      d = model_params$d,
      model_type = model_params$model_type,
      include.mean = model_params$include.mean
    ),
    innovations = innovations
  )

  # 4) Якщо метод == "pmm2", використовувати алгоритм PMM2
  if (method == "pmm2") {
    if (verbose) cat("Початок оптимізації PMM2...\n")

    model_info <- list(
      ar_order = model_params$ar_order,
      ma_order = model_params$ma_order,
      d = model_params$d,
      model_type = model_params$model_type,
      include.mean = model_params$include.mean,
      innovations = innovations,
      x = x_centered,
      x_mean = x_mean,
      original_x = orig_x,
      verbose = verbose
    )

    result <- pmm2_algorithm(
      b_init = b_init,
      X = dm$X,
      y = dm$y,
      m2 = m2,
      m3 = m3,
      m4 = m4,
      max_iter = max_iter,
      tol = tol,
      regularize = regularize,
      reg_lambda = reg_lambda,
      verbose = verbose
    )

    final_coef <- result$b
    converged <- result$convergence
    iterations <- result$iterations

    # Обчислити кінцеві залишки
    final_res <- compute_ts_residuals(final_coef, model_info)
  } else {
    # Для інших методів просто використовувати початкові оцінки
    final_coef <- b_init
    converged <- TRUE
    iterations <- 0
    final_res <- innovations
  }

  # 5) Створити відповідний об'єкт класу
  if (model_type == "ar") {
    result_class <- "ARPMM2"
  } else if (model_type == "ma") {
    result_class <- "MAPMM2"
  } else if (model_type == "arma") {
    result_class <- "ARMAPMM2"
  } else if (model_type == "arima") {
    result_class <- "ARIMAPMM2"
  } else {
    result_class <- "TS2fit"  # Базовий клас за замовчуванням
  }

  # Створити та повернути об'єкт відповідного класу
  new(result_class,
      coefficients    = as.numeric(final_coef),
      residuals       = as.numeric(final_res),
      m2              = as.numeric(m2),
      m3              = as.numeric(m3),
      m4              = as.numeric(m4),
      convergence     = converged,
      iterations      = as.numeric(iterations),
      call            = cl,
      model_type      = model_type,
      intercept       = as.numeric(x_mean),
      original_series = as.numeric(orig_x),
      order           = list(ar = model_params$ar_order,
                             ma = model_params$ma_order,
                             d = model_params$d))
}

#' Підігнати AR модель за допомогою PMM2 (обгортка)
#'
#' @inheritParams ts_pmm2
#' @export
ar_pmm2 <- function(x, order = 1, method = "pmm2", max_iter = 50, tol = 1e-6,
                    include.mean = TRUE, initial = NULL, na.action = na.fail,
                    regularize = TRUE, reg_lambda = 1e-8, verbose = FALSE) {
  ts_pmm2(x, order = order, model_type = "ar", method = method,
          max_iter = max_iter, tol = tol,
          include.mean = include.mean, initial = initial,
          na.action = na.action, regularize = regularize,
          reg_lambda = reg_lambda, verbose = verbose)
}

#' Підігнати MA модель за допомогою PMM2 (обгортка)
#'
#' @inheritParams ts_pmm2
#' @export
ma_pmm2 <- function(x, order = 1, method = "pmm2", max_iter = 50, tol = 1e-6,
                    include.mean = TRUE, initial = NULL, na.action = na.fail,
                    regularize = TRUE, reg_lambda = 1e-8, verbose = FALSE) {
  ts_pmm2(x, order = order, model_type = "ma", method = method,
          max_iter = max_iter, tol = tol,
          include.mean = include.mean, initial = initial,
          na.action = na.action, regularize = regularize,
          reg_lambda = reg_lambda, verbose = verbose)
}

#' Підігнати ARMA модель за допомогою PMM2 (обгортка)
#'
#' @inheritParams ts_pmm2
#' @export
arma_pmm2 <- function(x, order = c(1, 1), method = "pmm2", max_iter = 50, tol = 1e-6,
                      include.mean = TRUE, initial = NULL, na.action = na.fail,
                      regularize = TRUE, reg_lambda = 1e-8, verbose = FALSE) {
  ts_pmm2(x, order = order, model_type = "arma", method = method,
          max_iter = max_iter, tol = tol,
          include.mean = include.mean, initial = initial,
          na.action = na.action, regularize = regularize,
          reg_lambda = reg_lambda, verbose = verbose)
}

#' Підігнати ARIMA модель за допомогою PMM2 (обгортка)
#'
#' @inheritParams ts_pmm2
#' @export
arima_pmm2 <- function(x, order = c(1, 1, 1), method = "pmm2", max_iter = 50, tol = 1e-6,
                       include.mean = TRUE, initial = NULL, na.action = na.fail,
                       regularize = TRUE, reg_lambda = 1e-8, verbose = FALSE) {
  ts_pmm2(x, order = order, model_type = "arima", method = method,
          max_iter = max_iter, tol = tol,
          include.mean = include.mean, initial = initial,
          na.action = na.action, regularize = regularize,
          reg_lambda = reg_lambda, verbose = verbose)
}

# --- MA utilities ---------------------------------------------------------

ma_css_fit <- function(x, q, include_mean = TRUE, verbose = FALSE) {
  fit <- tryCatch(
    stats::arima(x, order = c(0, 0, q), method = "CSS-ML",
                 include.mean = include_mean),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    if (verbose) {
      cat("Не вдалось оцінити MA модель через stats::arima: повертаю нульові коефіцієнти\n")
    }
    coef_css <- rep(0, q)
    intercept <- if (include_mean) mean(x) else 0
    residuals <- x - intercept
    residuals[is.na(residuals)] <- 0
    return(list(
      coefficients = coef_css,
      intercept = intercept,
      residuals = residuals,
      convergence = FALSE,
      iterations = 0L
    ))
  }

  names_coef <- names(fit$coef)
  coef_css <- numeric(q)
  for (j in seq_len(q)) {
    coef_name <- paste0("ma", j)
    if (coef_name %in% names_coef) {
      coef_css[j] <- fit$coef[coef_name]
    } else if (length(fit$coef) >= j) {
      coef_css[j] <- fit$coef[j]
    } else {
      coef_css[j] <- 0
    }
  }

  intercept <- 0
  if (include_mean) {
    intercept_name <- setdiff(names_coef, paste0("ma", seq_len(q)))
    if (length(intercept_name) > 0) {
      intercept <- as.numeric(fit$coef[intercept_name[1]])
    }
  }

  residuals <- as.numeric(fit$residuals)
  residuals[is.na(residuals)] <- 0

  list(
    coefficients = coef_css,
    intercept = intercept,
    residuals = residuals,
    convergence = TRUE,
    iterations = 1L
  )
}

ma_pmm2_fit <- function(x, q, css_fit, max_iter = 50, tol = 1e-6, verbose = FALSE) {
  design <- ma_build_design(css_fit$intercept, css_fit$residuals, x, q)
  moments <- compute_moments(css_fit$residuals)

  b_init <- c(0, css_fit$coefficients)
  solve_res <- ma_solve_pmm2(b_init, design$X, design$y,
                             moments$m2, moments$m3, moments$m4,
                             max_iter = max_iter, tol = tol,
                             verbose = verbose)

  theta <- solve_res$coefficients
  innovations <- ma_compute_innovations(x - css_fit$intercept, theta, q)

  list(
    coefficients = theta,
    intercept = css_fit$intercept,
    innovations = innovations,
    convergence = solve_res$convergence,
    iterations = solve_res$iterations
  )
}

ma_build_design <- function(intercept, residuals, x, q) {
  idx <- seq.int(q + 1L, length(x))
  X <- matrix(1, nrow = length(idx), ncol = q + 1L)
  for (j in seq_len(q)) {
    X[, j + 1L] <- residuals[idx - j]
  }
  y <- x[idx] - intercept
  list(X = X, y = y)
}

ma_solve_pmm2 <- function(b_init, X, Y, m2, m3, m4,
                          max_iter = 50, tol = 1e-6,
                          verbose = FALSE) {
  b <- as.numeric(b_init)
  iterations <- 0L
  converged <- FALSE
  for (iter in seq_len(max_iter)) {
    iterations <- iter
    S <- as.vector(X %*% b)
    Z1 <- m3 * S^2 + (m4 - m2^2 - 2 * m3 * Y) * S +
      (m3 * Y^2 - (m4 - m2^2) * Y - m2 * m3)
    Z <- as.numeric(t(X) %*% Z1)
    JZ11 <- 2 * m3 * S + (m4 - m2^2 - 2 * m3 * Y)
    J <- t(X) %*% (X * JZ11)
    step <- tryCatch(solve(J, Z), error = function(e) NULL)
    if (is.null(step)) {
      if (verbose) cat("Система сингулярна на ітерації", iter, "\n")
      break
    }
    b_new <- b - step
    if (sqrt(sum((b_new - b)^2)) < tol) {
      b <- b_new
      converged <- TRUE
      break
    }
    b <- b_new
  }
  list(coefficients = b[-1], convergence = converged, iterations = iterations)
}

ma_compute_innovations <- function(x, theta, q) {
  n <- length(x)
  innovations <- numeric(n)
  history <- rep(0, q)
  for (t in seq_len(n)) {
    ma_component <- if (q > 0) sum(theta * history) else 0
    innovations[t] <- x[t] - ma_component
    if (q > 0) {
      history <- c(innovations[t], history)[seq_len(q)]
    }
  }
  innovations
}
