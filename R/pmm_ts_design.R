# pmm_ts_design.R - Функції для роботи з дизайн-матрицями часових рядів

#' Перевірка та підготовка параметрів часового ряду
#'
#' @param x Дані часового ряду
#' @param order Специфікація порядку моделі
#' @param model_type Тип моделі (ar, ma, arma, або arima)
#' @param include.mean Чи включати середнє/перехоплення
#'
#' @return Список перевірених параметрів та інформації про модель
#' @keywords internal
validate_ts_parameters <- function(x, order, model_type, include.mean) {
  # Перевірка вхідних даних
  if (missing(x)) {
    stop("Відсутній аргумент 'x'")
  }

  # Перетворити вхідні дані на числовий вектор
  x <- as.numeric(x)

  if (!is.numeric(x)) {
    stop("'x' має бути числовим вектором")
  }

  # Перевірити на NA та нескінченні значення
  if (any(is.na(x)) || any(is.infinite(x))) {
    warning("Виявлено NA або нескінченні значення у вхідному ряді. Вони будуть видалені.")
    x <- x[!is.na(x) & !is.infinite(x)]
    if (length(x) < 10) {
      stop("Замало валідних спостережень після видалення NA/нескінченних значень")
    }
  }

  if (missing(order)) {
    stop("Відсутній аргумент 'order'")
  }

  # Розбір параметра order в залежності від model_type
  if (model_type == "ar") {
    if (!is.numeric(order) || length(order) != 1)
      stop("Для AR моделей 'order' має бути одним цілим числом")
    ar_order <- as.integer(order)
    ma_order <- 0
    d <- 0
    if (ar_order <= 0)
      stop("Порядок AR має бути додатним")

    # Перевірити, чи достатньо даних
    if (length(x) <= ar_order + 1) {
      stop("Замало спостережень для моделі AR порядку ", ar_order)
    }
  } else if (model_type == "ma") {
    if (!is.numeric(order) || length(order) != 1)
      stop("Для MA моделей 'order' має бути одним цілим числом")
    ar_order <- 0
    ma_order <- as.integer(order)
    d <- 0
    if (ma_order <= 0)
      stop("Порядок MA має бути додатним")

    # Перевірити, чи достатньо даних
    if (length(x) <= ma_order + 1) {
      stop("Замало спостережень для моделі MA порядку ", ma_order)
    }
  } else if (model_type == "arma") {
    if (!is.numeric(order) || length(order) != 2)
      stop("Для ARMA моделей 'order' має бути вектором довжини 2 (порядок AR, порядок MA)")
    ar_order <- as.integer(order[1])
    ma_order <- as.integer(order[2])
    d <- 0
    if (ar_order < 0 || ma_order < 0)
      stop("Порядки AR та MA мають бути невід'ємними")
    if (ar_order == 0 && ma_order == 0)
      stop("Принаймні один з порядків AR або MA має бути додатним")

    # Перевірити, чи достатньо даних
    if (length(x) <= max(ar_order, ma_order) + 1) {
      stop("Замало спостережень для моделі ARMA порядків (", ar_order, ",", ma_order, ")")
    }
  } else if (model_type == "arima") {
    if (!is.numeric(order) || length(order) != 3)
      stop("Для ARIMA моделей 'order' має бути вектором довжини 3 (порядок AR, диференціювання, порядок MA)")
    ar_order <- as.integer(order[1])
    d <- as.integer(order[2])
    ma_order <- as.integer(order[3])
    if (ar_order < 0 || ma_order < 0 || d < 0)
      stop("Порядки AR, диференціювання та MA мають бути невід'ємними")
    if (ar_order == 0 && ma_order == 0 && d == 0)
      stop("Принаймні один з порядків AR, диференціювання або MA має бути додатним")

    # Перевірити, чи достатньо даних після диференціювання
    if (length(x) <= d + max(ar_order, ma_order) + 1) {
      stop("Замало спостережень для моделі ARIMA після диференціювання")
    }
  } else {
    stop("Невідомий тип моделі: ", model_type)
  }

  # Зберегти оригінальний ряд
  orig_x <- as.numeric(x)

  list(
    original_x = orig_x,
    ar_order   = ar_order,
    ma_order   = ma_order,
    d          = d,
    model_type = model_type,
    include.mean = include.mean
  )
}

#' Створення матриці дизайну для AR моделі
#'
#' @param x центрований часовий ряд
#' @param p порядок AR
#' @return Матриця дизайну з лаговими значеннями
#' @keywords internal
create_ar_matrix <- function(x, p) {
  n <- length(x)
  if (n <= p) {
    stop("Недостатньо точок даних для AR порядку p = ", p)
  }

  nr <- n - p
  M <- matrix(0, nr, p)
  for (i in seq_len(p)) {
    M[, i] <- x[(p - i + 1):(n - i)]
  }
  M
}

#' Отримати оцінки Юла-Волкера для AR(p)
#'
#' @param x числовий вектор
#' @param p ціле значення порядку AR
#' @return числовий вектор довжини p (коефіцієнти AR)
#' @keywords internal
get_yw_estimates <- function(x, p) {
  # Це спрощений підхід, який може не обробляти граничні випадки
  r <- numeric(p+1)
  n <- length(x)
  xm <- mean(x)
  for (k in 0:p) {
    r[k+1] <- sum((x[1:(n-k)] - xm)*(x[(k+1):n] - xm))
  }
  R <- matrix(0, p, p)
  for (i in seq_len(p)) {
    for (j in seq_len(p)) {
      R[i,j] <- r[abs(i-j)+1]
    }
  }
  rhs <- r[2:(p+1)]
  phi <- solve(R, rhs)
  phi
}

#' Створення дизайн-матриці для часових рядів
#'
#' @param x Дані часового ряду
#' @param model_info Список з параметрами моделі
#' @param innovations Опціональні інновації/залишки для MA компонентів
#'
#' @return Список з дизайн-матрицею, змінною відгуку та іншими компонентами
#' @keywords internal
create_ts_design_matrix <- function(x, model_info, innovations = NULL) {
  # Витягнення параметрів моделі
  ar_order <- model_info$ar_order
  ma_order <- model_info$ma_order
  d <- model_info$d
  model_type <- model_info$model_type
  include_mean <- model_info$include.mean

  # Перевірка, чи потрібно диференціювати ряд
  if (model_type == "arima" && d > 0) {
    x_diff <- diff(x, differences = d)
  } else {
    x_diff <- x
  }

  # Обробка середнього
  if (include_mean) {
    x_mean <- mean(x_diff, na.rm = TRUE)
    x_centered <- x_diff - x_mean
  } else {
    x_mean <- 0
    x_centered <- x_diff
  }

  # Обчислення максимального лагу та ефективної довжини даних
  max_lag <- max(ar_order, ma_order)
  n_data <- length(x_centered)

  if (n_data <= max_lag) {
    stop("Недостатньо даних для побудови моделі після диференціювання")
  }

  n_rows <- n_data - max_lag
  n_cols <- ar_order + ma_order

  # Створення дизайн-матриці з відповідними розмірами
  X <- matrix(0, nrow = n_rows, ncol = n_cols)
  y <- x_centered[(max_lag + 1):n_data]

  # Додавання AR компонентів
  if (ar_order > 0) {
    col_index <- 1
    for (i in 1:ar_order) {
      X[, col_index] <- x_centered[(max_lag - i + 1):(n_data - i)]
      col_index <- col_index + 1
    }
  }

  # Додавання MA компонентів, якщо потрібно
  if (ma_order > 0) {
    # Якщо інновації не надані, використовуємо нулі
    if (is.null(innovations)) {
      innovations <- rep(0, n_data)
    }

    # Забезпечення правильної довжини інновацій
    if (length(innovations) < n_data) {
      innovations <- c(rep(0, n_data - length(innovations)), innovations)
    } else if (length(innovations) > n_data) {
      innovations <- tail(innovations, n_data)
    }

    # Заповнення MA стовпців
    col_index <- ar_order + 1
    for (j in 1:ma_order) {
      X[, col_index] <- innovations[(max_lag - j + 1):(n_data - j)]
      col_index <- col_index + 1
    }
  }

  # Повернення результатів у вигляді списку
  list(
    X = X,
    y = y,
    x_centered = x_centered,
    x_mean = x_mean,
    n_rows = n_rows,
    n_cols = n_cols,
    effective_length = n_rows,
    innovations = innovations,
    original_x = x
  )
}

#' Отримати початкові оцінки параметрів для моделей часових рядів
#'
#' @param model_params Перевірені параметри моделі з validate_ts_parameters
#' @param initial Опціонально надані користувачем початкові оцінки
#' @param method Метод оцінювання
#' @param verbose Виводити детальну інформацію
#'
#' @return Список, що містить:
#'   \item{b_init}{вектор початкових коефіцієнтів AR/MA}
#'   \item{x_mean}{оцінене середнє (якщо include.mean=TRUE)}
#'   \item{innovations}{початкові залишки/інновації}
#'   \item{x_centered}{центрований (або диференційований + центрований) ряд}
#'   \item{m2}{другий центральний момент початкових залишків}
#'   \item{m3}{третій центральний момент початкових залишків}
#'   \item{m4}{четвертий центральний момент початкових залишків}
#' @keywords internal
get_initial_estimates <- function(model_params,
                                  initial = NULL,
                                  method = "pmm2",
                                  verbose = FALSE) {
  x         <- model_params$original_x
  ar_order  <- model_params$ar_order
  ma_order  <- model_params$ma_order
  d         <- model_params$d
  mtype     <- model_params$model_type
  inc_mean  <- model_params$include.mean

  # Можливо диференціювати для ARIMA
  if (mtype == "arima" && d > 0) {
    x_diff <- diff(x, differences = d)
  } else {
    x_diff <- x
  }

  # Центрувати, якщо потрібно
  if (inc_mean) {
    x_mean <- mean(x_diff, na.rm = TRUE)
    x_centered <- x_diff - x_mean
  } else {
    x_mean <- 0
    x_centered <- x_diff
  }

  if (mtype == "ar") {
    # AR(p): швидкий підхід для початкових значень
    if (is.null(initial)) {
      if (method == "yw") {
        b_init <- get_yw_estimates(x_centered, ar_order)
      } else {
        X <- create_ar_matrix(x_centered, ar_order)
        y <- x_centered[(ar_order + 1):length(x_centered)]
        fit_ols <- lm.fit(x = X, y = y)
        b_init <- fit_ols$coefficients
      }
    } else {
      if (length(initial) != ar_order) {
        stop("Довжина 'initial' має відповідати порядку AR")
      }
      b_init <- initial
    }
    # Інновації з початкової підгонки
    X <- create_ar_matrix(x_centered, ar_order)
    y <- x_centered[(ar_order + 1):length(x_centered)]
    innovations <- as.numeric(y - X %*% b_init)

  } else if (mtype %in% c("ma", "arma", "arima")) {
    # Використовувати stats::arima для початкового припущення або надані користувачем
    arima_order <- c(ar_order, ifelse(mtype=="arima", d, 0), ma_order)
    if (is.null(initial)) {
      init_fit <- NULL
      try_methods <- c("CSS","CSS-ML","ML")
      if(method=="pmm2") try_methods <- c("CSS","CSS-ML","ML")

      for(mm in try_methods) {
        tmp <- tryCatch({
          stats::arima(x, order=arima_order, method=mm,
                       include.mean=inc_mean && (mtype!="arima"))
        }, error=function(e) NULL)
        if(!is.null(tmp)) {
          init_fit <- tmp
          break
        }
      }
      if(is.null(init_fit)) {
        if(verbose) cat("Усі стандартні методи не спрацювали; використовуються спрощені значення.\n")
        init_fit <- list(
          coef = numeric(ar_order + ma_order),
          residuals = if(mtype=="arima") x_diff else x_centered
        )
        if(ar_order>0) names(init_fit$coef)[1:ar_order] <- paste0("ar",1:ar_order)
        if(ma_order>0) names(init_fit$coef)[(ar_order+1):(ar_order+ma_order)] <- paste0("ma",1:ma_order)
      }

      ar_init <- rep(0, ar_order)
      ma_init <- rep(0, ma_order)
      if(ar_order>0) {
        idx <- paste0("ar",1:ar_order)
        ar_init <- if(all(idx %in% names(init_fit$coef))) as.numeric(init_fit$coef[idx]) else rep(0.1, ar_order)
      }
      if(ma_order>0) {
        idx <- paste0("ma",1:ma_order)
        ma_init <- if(all(idx %in% names(init_fit$coef))) as.numeric(init_fit$coef[idx]) else rep(0.1, ma_order)
      }
      if(inc_mean && !is.null(init_fit$coef) && ("intercept" %in% names(init_fit$coef))) {
        x_mean <- init_fit$coef["intercept"]
      }
      innovations <- if(!is.null(init_fit$residuals)) as.numeric(init_fit$residuals) else {
        rnorm(length(x_centered),0,sd(x_centered,na.rm=TRUE))
      }
      b_init <- c(ar_init, ma_init)

    } else {
      # Надано initial
      if(is.list(initial)) {
        if(ar_order>0 && is.null(initial$ar)) {
          stop("Відсутній 'ar' у списку initial, але ar_order>0")
        }
        if(ma_order>0 && is.null(initial$ma)) {
          stop("Відсутній 'ma' у списку initial, але ma_order>0")
        }
        ar_init <- if(ar_order>0) initial$ar else numeric(0)
        ma_init <- if(ma_order>0) initial$ma else numeric(0)
      } else {
        if(length(initial) != (ar_order+ma_order)) {
          stop("Довжина 'initial' має відповідати сумі порядків AR та MA")
        }
        ar_init <- if(ar_order>0) initial[1:ar_order] else numeric(0)
        ma_init <- if(ma_order>0) initial[(ar_order+1):(ar_order+ma_order)] else numeric(0)
      }
      b_init <- c(ar_init, ma_init)

      init_fit <- tryCatch({
        stats::arima(x, order=arima_order,
                     fixed=b_init,
                     include.mean=inc_mean && (mtype!="arima"))
      }, error=function(e) {
        if(verbose) cat("Помилка з наданими користувачем початковими значеннями:",e$message,"\n")
        list(residuals = if(mtype=="arima") x_diff else x_centered)
      })
      innovations <- as.numeric(init_fit$residuals)
    }
  }

  if(anyNA(b_init)) {
    warning("NA в початкових параметрах замінені на 0.")
    b_init[is.na(b_init)] <- 0
  }

  # Обчислення моментів
  moments <- compute_moments(innovations)

  # Повернення результатів
  list(
    b_init      = b_init,
    x_mean      = x_mean,
    innovations = innovations,
    x_centered  = x_centered,
    orig_x      = x,
    m2          = moments$m2,
    m3          = moments$m3,
    m4          = moments$m4
  )
}

#' Оновити інновації MA моделі
#'
#' @param x центрований часовий ряд
#' @param ma_coef вектор коефіцієнтів MA
#' @return вектор інновацій
#' @keywords internal
update_ma_innovations <- function(x, ma_coef) {
  n <- length(x)
  q <- length(ma_coef)

  # Ініціалізувати інновації як нулі
  innovations <- numeric(n)

  # Ітеративно обчислити інновації
  for(t in 1:n) {
    # Обчислити очікуване значення на основі попередніх інновацій
    expected <- 0
    for(j in 1:q) {
      if(t - j > 0) {
        expected <- expected + ma_coef[j] * innovations[t - j]
      }
    }

    # Обчислити поточну інновацію
    innovations[t] <- x[t] - expected
  }

  # Перевірити на неконечні значення
  if(any(is.infinite(innovations)) || any(is.na(innovations))) {
    warning("Виявлено неконечні інновації в update_ma_innovations. Використовуємо регуляризовані значення.")
    # Замінити проблемні значення на середнє або 0
    bad_idx <- is.infinite(innovations) | is.na(innovations)
    if(sum(!bad_idx) > 0) {
      # Якщо є дійсні значення, використовуємо їх середнє
      innovations[bad_idx] <- mean(innovations[!bad_idx])
    } else {
      # Інакше використовуємо 0
      innovations[bad_idx] <- 0
    }
  }

  return(innovations)
}

#' Обчислити кінцеві залишки для моделей часових рядів
#'
#' @param coefs Оцінені коефіцієнти
#' @param model_info Інформація про модель
#' @return Вектор залишків
#' @keywords internal
compute_ts_residuals <- function(coefs, model_info) {
  # Витягнути параметри моделі
  x           <- model_info$x
  ar_order    <- model_info$ar_order
  ma_order    <- model_info$ma_order
  d           <- model_info$d
  model_type  <- model_info$model_type
  include.mean <- model_info$include.mean

  # Розділити коефіцієнти на AR та MA частини
  if (ar_order > 0) {
    ar_coefs <- coefs[1:ar_order]
  } else {
    ar_coefs <- numeric(0)
  }

  if (ma_order > 0) {
    ma_coefs <- coefs[(ar_order+1):(ar_order+ma_order)]
  } else {
    ma_coefs <- numeric(0)
  }

  # Для більш надійного обчислення залишків використовуємо arima з фіксованими параметрами
  arima_order <- c(ar_order, ifelse(model_type == "arima", d, 0), ma_order)

  # Підготувати фіксовані параметри для arima
  fixed_params <- c(ar_coefs, ma_coefs)
  if (include.mean) {
    fixed_params <- c(fixed_params, model_info$x_mean)
    names(fixed_params) <- c(
      if(ar_order > 0) paste0("ar", 1:ar_order) else NULL,
      if(ma_order > 0) paste0("ma", 1:ma_order) else NULL,
      "intercept"
    )
  } else {
    names(fixed_params) <- c(
      if(ar_order > 0) paste0("ar", 1:ar_order) else NULL,
      if(ma_order > 0) paste0("ma", 1:ma_order) else NULL
    )
  }

  # Обчислити залишки з фіксованими параметрами
  final_fit <- tryCatch({
    stats::arima(model_info$original_x,
                 order = arima_order,
                 fixed = fixed_params,
                 include.mean = include.mean)
  }, error = function(e) {
    if (verbose) cat("Помилка при обчисленні кінцевих залишків:", e$message, "\n")
    list(residuals = rep(NA, length(model_info$original_x)))
  })

  # Витягнути залишки та забезпечити правильну довжину
  final_res <- as.numeric(final_fit$residuals)

  if (length(final_res) < length(model_info$original_x)) {
    final_res <- c(rep(NA, length(model_info$original_x) - length(final_res)), final_res)
  }

  return(final_res)
}
