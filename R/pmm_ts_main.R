# pmm_ts_main.R - Уніфікований модуль для моделей часових рядів з PMM2

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
      original_x = orig_x
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
