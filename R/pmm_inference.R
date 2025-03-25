# pmm_inference.R - Статистичний висновок для моделей PMM2

#' Бутстреп-висновок для підгонки PMM2
#'
#' @param object об'єкт класу PMM2fit
#' @param formula та сама формула, що використовувалася спочатку
#' @param data фрейм даних, що використовувався спочатку
#' @param B кількість бутстреп-реплікацій
#' @param seed (опціонально) для відтворюваності
#' @param parallel логічне, чи використовувати паралельні обчислення
#' @param cores кількість ядер для використання при паралельних обчисленнях, за замовчуванням - автовизначення
#'
#' @return data.frame з стовпцями: Estimate, Std.Error, t.value, p.value
#' @export
pmm2_inference <- function(object, formula, data, B=200, seed=NULL,
                           parallel=FALSE, cores=NULL) {
  # Встановити зерно для відтворюваності, якщо надано
  if(!is.null(seed)) set.seed(seed)

  # Витягнути коефіцієнти та залишки
  coefs <- object@coefficients
  res   <- object@residuals

  # Перевірка вхідних даних
  if(B < 10) {
    warning("Кількість бутстреп-вибірок (B) дуже мала. Розгляньте використання B >= 100 для більш надійного висновку.")
  }

  if(!inherits(object, "PMM2fit")) {
    stop("Об'єкт має бути класу 'PMM2fit'")
  }

  if(missing(formula) || missing(data)) {
    stop("Обидва 'formula' та 'data' мають бути надані")
  }

  # Побудувати матриці X, y
  mf <- model.frame(formula, data)
  X <- model.matrix(formula, mf)
  y <- model.response(mf)
  n <- nrow(X)

  # Раннє повернення у випадку помилок
  if(is.null(y) || is.null(X)) {
    stop("Не вдалося витягнути відгук або матрицю дизайну з даних")
  }

  # Перевірити, чи слід використовувати паралельні обчислення
  use_parallel <- parallel && requireNamespace("parallel", quietly = TRUE)

  if(use_parallel) {
    if(is.null(cores)) {
      cores <- max(1, parallel::detectCores() - 1)
    }

    boot_results <- parallel::mclapply(seq_len(B), function(b) {
      # 1) Бутстреп залишків
      res_b <- sample(res, size=n, replace=TRUE)

      # 2) Створити новий y
      y_b <- X %*% coefs + res_b

      # 3) Створити нові дані
      data_b <- data
      # Припустити, що ліва сторона є першим терміном у формулі
      lhs <- as.character(formula[[2]])
      data_b[[lhs]] <- as.numeric(y_b)

      # 4) Повторно оцінити модель
      fit_b <- tryCatch({
        lm_pmm2(formula, data_b, max_iter=20, tol=1e-6)
      }, error = function(e) {
        warning("Бутстреп-реплікація ", b, " не вдалася: ", e$message)
        return(NULL)
      })

      if(!is.null(fit_b)) {
        return(fit_b@coefficients)
      } else {
        return(rep(NA, length(coefs)))
      }
    }, mc.cores = cores)

    # Перетворити список на матрицю
    boot_est <- do.call(rbind, boot_results)

  } else {
    # Послідовні обчислення
    # Матриця для зберігання результатів
    boot_est <- matrix(0, nrow=B, ncol=length(coefs))
    colnames(boot_est) <- names(coefs)

    # Відстеження прогресу
    pb <- NULL
    if(interactive() && B > 10) {
      if(requireNamespace("utils", quietly = TRUE)) {
        pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
      }
    }

    for(b in seq_len(B)) {
      # 1) Бутстреп залишків
      res_b <- sample(res, size=n, replace=TRUE)

      # 2) Створити новий y
      y_b <- X %*% coefs + res_b

      # 3) Створити нові дані
      data_b <- data
      # Припустити, що ліва сторона є першим терміном у формулі
      lhs <- as.character(formula[[2]])
      data_b[[lhs]] <- as.numeric(y_b)

      # 4) Повторно оцінити модель
      fit_b <- tryCatch({
        lm_pmm2(formula, data_b, max_iter=20, tol=1e-6)
      }, error = function(e) {
        warning("Бутстреп-реплікація ", b, " не вдалася: ", e$message)
        return(NULL)
      })

      if(!is.null(fit_b)) {
        boot_est[b, ] <- fit_b@coefficients
      } else {
        boot_est[b, ] <- NA
      }

      # Оновити індикатор прогресу
      if(!is.null(pb)) utils::setTxtProgressBar(pb, b)
    }

    # Закрити індикатор прогресу
    if(!is.null(pb)) close(pb)
  }

  # Видалити рядки зі значеннями NA
  na_rows <- apply(boot_est, 1, function(row) any(is.na(row)))
  if(any(na_rows)) {
    warning("Видалено ", sum(na_rows), " бутстреп-реплікацій через помилки оцінювання")
    boot_est <- boot_est[!na_rows, , drop = FALSE]
  }

  # Перевірити, чи маємо достатньо успішних бутстрепів
  if(nrow(boot_est) < 10) {
    stop("Замало успішних бутстреп-реплікацій для обчислення надійного висновку")
  }

  # Обчислити коваріаційну матрицю та стандартні помилки
  cov_mat <- cov(boot_est)
  est <- coefs
  se  <- sqrt(diag(cov_mat))

  # Обчислити t-значення та p-значення
  t_val <- est / se
  # Для великих вибірок використовувати нормальне наближення
  p_val <- 2 * (1 - pnorm(abs(t_val)))

  # Створити вихідний фрейм даних
  out <- data.frame(
    Estimate  = est,
    Std.Error = se,
    t.value   = t_val,
    p.value   = p_val
  )
  rownames(out) <- names(est)

  # Обчислити довірчі інтервали
  ci <- t(apply(boot_est, 2, quantile, probs = c(0.025, 0.975)))
  colnames(ci) <- c("2.5%", "97.5%")

  # Додати довірчі інтервали до виходу
  out$conf.low <- ci[, "2.5%"]
  out$conf.high <- ci[, "97.5%"]

  return(out)
}

#' Побудувати графіки бутстреп-розподілів для підгонки PMM2
#'
#' @param object Результат з pmm2_inference
#' @param coefficients Які коефіцієнти побудувати, за замовчуванням усі
#'
#' @return Невидимо повертає інформацію про гістограму
#' @export
plot_pmm2_bootstrap <- function(object, coefficients = NULL) {
  if(!inherits(object, "data.frame") ||
     !all(c("Estimate", "Std.Error", "conf.low", "conf.high") %in% names(object))) {
    stop("Об'єкт має бути результатом pmm2_inference()")
  }

  # Якщо коефіцієнти не вказані, використовувати всі
  if(is.null(coefficients)) {
    coefficients <- rownames(object)
  }

  # Фільтрувати до запитаних коефіцієнтів
  object_subset <- object[intersect(coefficients, rownames(object)), , drop = FALSE]

  # Перевірка на порожній набір даних
  if(nrow(object_subset) == 0) {
    warning("Жоден із запитаних коефіцієнтів не знайдено в результатах.")
    return(invisible(NULL))
  }

  # Налаштувати компонування графіка
  n_coefs <- nrow(object_subset)
  n_cols <- min(2, n_coefs)
  n_rows <- ceiling(n_coefs / n_cols)

  # Зберегти старі налаштування par і відновити при виході
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  par(mfrow = c(n_rows, n_cols))

  # Створити графік щільності для кожного коефіцієнта
  result <- lapply(seq_len(n_coefs), function(i) {
    coef_name <- rownames(object_subset)[i]
    est <- object_subset[i, "Estimate"]
    ci_low <- object_subset[i, "conf.low"]
    ci_high <- object_subset[i, "conf.high"]
    se <- object_subset[i, "Std.Error"]

    # Перевірка на скінченні значення
    if(!is.finite(est) || !is.finite(ci_low) || !is.finite(ci_high) || !is.finite(se)) {
      warning("Нескінченні або NA значення для коефіцієнта ", coef_name,
              ". Пропускаємо цей графік.")
      return(NULL)
    }

    # Створити заголовок графіка
    main_title <- paste0(coef_name, "\nОцінка: ", round(est, 4))

    # Оцінити діапазон для осі x
    # Використовуємо більш надійний підхід для визначення діапазону
    x_range <- range(c(est, ci_low, ci_high), na.rm = TRUE)
    # Розширити діапазон на 20% в обох напрямках
    x_range_width <- diff(x_range)
    x_range <- x_range + c(-0.2, 0.2) * x_range_width

    # Створити точки для осі x
    x_seq <- seq(x_range[1], x_range[2], length.out = 100)

    # Створити значення щільності для нормального розподілу
    y_seq <- dnorm(x_seq, mean = est, sd = se)

    # Побудувати графік
    plot(x_seq, y_seq, type = "l",
         main = main_title,
         xlab = "Значення",
         ylab = "Щільність")

    # Додати вертикальні лінії для оцінки та CI
    abline(v = est, col = "red", lwd = 2)
    abline(v = ci_low, col = "blue", lty = 2)
    abline(v = ci_high, col = "blue", lty = 2)

    # Додати легенду
    legend("topright",
           legend = c("Оцінка", "95% CI"),
           col = c("red", "blue"),
           lty = c(1, 2),
           lwd = c(2, 1),
           cex = 0.8)

    invisible(list(x = x_seq, y = y_seq, estimate = est,
                   ci_low = ci_low, ci_high = ci_high))
  })

  # Видалити NULL результати
  result <- result[!sapply(result, is.null)]

  # Якщо всі результати NULL, повернути NULL
  if(length(result) == 0) {
    warning("Не вдалося створити жодного графіка.")
    return(invisible(NULL))
  }

  # Додати імена до результатів
  names(result) <- rownames(object_subset)[sapply(seq_len(n_coefs), function(i) {
    !is.null(result[[i]])
  })]

  invisible(result)
}


#' Бутстреп-висновок для моделей часових рядів PMM2
#'
#' @param object об'єкт класу TS2fit
#' @param x (опціонально) оригінальний часовий ряд; якщо NULL, використовує object@original_series
#' @param B кількість бутстреп-реплікацій
#' @param seed (опціонально) для відтворюваності
#' @param block_length довжина блоку для блокового бутстрепу; якщо NULL, використовує евристичне значення
#' @param method тип бутстрепу: "residual" або "block"
#' @param parallel логічне, чи використовувати паралельні обчислення
#' @param cores кількість ядер для паралельних обчислень
#' @param debug логічне, чи виводити додаткову діагностичну інформацію
#'
#' @return data.frame з стовпцями: Estimate, Std.Error, t.value, p.value
#' @export
ts_pmm2_inference <- function(object, x = NULL, B = 200, seed = NULL,
                              block_length = NULL, method = c("residual", "block"),
                              parallel = FALSE, cores = NULL, debug = FALSE) {
  # Перевірити клас об'єкта
  if (!inherits(object, "TS2fit")) {
    stop("Об'єкт має бути класу 'TS2fit'")
  }

  # Вибрати метод бутстрепу
  method <- match.arg(method)

  # Встановити зерно для відтворюваності, якщо надано
  if (!is.null(seed)) set.seed(seed)

  # Витягнути параметри моделі
  model_type <- object@model_type
  ar_order <- object@order$ar
  ma_order <- object@order$ma
  d <- object@order$d
  intercept <- object@intercept
  include_mean <- intercept != 0

  if(debug) {
    cat("Параметри моделі:\n")
    cat("model_type:", model_type, "\n")
    cat("ar_order:", ar_order, "\n")
    cat("ma_order:", ma_order, "\n")
    cat("d:", d, "\n")
    cat("intercept:", intercept, "\n")
    cat("include_mean:", include_mean, "\n")
  }

  # Якщо x не наданий, використовувати оригінальний ряд з об'єкта
  if (is.null(x)) {
    x <- object@original_series
  }

  # Витягнути коефіцієнти та залишки
  coefs <- object@coefficients
  res <- object@residuals

  if(debug) {
    cat("Розмір оригінального ряду:", length(x), "\n")
    cat("Розмір вектора залишків:", length(res), "\n")
    cat("Кількість коефіцієнтів:", length(coefs), "\n")
  }

  # Виконати блоковий бутстреп, якщо вказано
  if (method == "block") {
    # Визначити довжину блоку, якщо не надано
    if (is.null(block_length)) {
      # Використовувати евристику: квадратний корінь з довжини ряду
      block_length <- ceiling(sqrt(length(x)))
    }

    # Перевірити, чи довжина блоку має сенс
    if (block_length < 2) {
      warning("Довжина блоку замала. Встановлюю на 2.")
      block_length <- 2
    }
    if (block_length > length(x) / 4) {
      warning("Довжина блоку завелика. Встановлюю на 1/4 довжини ряду.")
      block_length <- floor(length(x) / 4)
    }

    # Функція для генерації бутстреп-ряду з блокового бутстрепу
    generate_block_bootstrap <- function(x, block_length) {
      n <- length(x)
      blocks_needed <- ceiling(n / block_length)

      # Можливі початкові позиції блоків
      start_positions <- 1:(n - block_length + 1)

      # Вибрати випадкові початкові позиції
      selected_starts <- sample(start_positions, blocks_needed, replace = TRUE)

      # Створити бутстреп-ряд
      boot_series <- numeric(0)
      for (start in selected_starts) {
        boot_series <- c(boot_series, x[start:(start + block_length - 1)])
      }

      # Обрізати до оригінальної довжини
      boot_series[1:n]
    }

    # Логіка бутстрепу відрізняється для блокового методу
    boot_function <- function(b) {
      # Генерувати новий ряд за допомогою блокового бутстрепу
      x_b <- generate_block_bootstrap(x, block_length)

      # Визначаємо правильний формат order в залежності від типу моделі
      if(model_type == "ar") {
        boot_order <- ar_order  # Для AR моделей - одне число
      } else if(model_type == "ma") {
        boot_order <- ma_order  # Для MA моделей - одне число
      } else if(model_type == "arma") {
        boot_order <- c(ar_order, ma_order)  # Для ARMA - вектор довжини 2
      } else if(model_type == "arima") {
        boot_order <- c(ar_order, d, ma_order)  # Для ARIMA - вектор довжини 3
      } else {
        stop("Невідомий тип моделі: ", model_type)
      }

      # Підігнати модель на бутстреп-ряді
      fit_b <- tryCatch({
        ts_pmm2(x_b, order = boot_order,
                model_type = model_type,
                include.mean = include_mean)
      }, error = function(e) {
        warning("Бутстреп-реплікація ", b, " не вдалася: ", e$message)
        return(NULL)
      })

      if (!is.null(fit_b)) {
        return(fit_b@coefficients)
      } else {
        return(rep(NA, length(coefs)))
      }
    }
  } else {
    # Для бутстрепу залишків
    # Функція для генерації нового ряду на основі моделі та бутстреп-залишків
    generate_with_residuals <- function(model, residuals) {
      n <- length(model@original_series)

      # Для ARIMA моделей, потрібно спочатку диференціювати
      if (model_type == "arima" && d > 0) {
        # Тут потрібна складніша логіка для відновлення оригінального ряду
        # після генерації диференційованого ряду
        # Це спрощений підхід:
        diff_x <- numeric(n - d)

        # Бутстреп залишків
        res_b <- sample(residuals[!is.na(residuals)], size = n - d, replace = TRUE)

        # Генерувати диференційований ряд
        if (ar_order > 0) {
          ar_coefs <- model@coefficients[1:ar_order]
        } else {
          ar_coefs <- numeric(0)
        }

        if (ma_order > 0) {
          ma_coefs <- model@coefficients[(ar_order+1):(ar_order+ma_order)]
        } else {
          ma_coefs <- numeric(0)
        }

        # Спрощена симуляція ARIMA процесу
        diff_x <- arima.sim(model = list(
          ar = if(ar_order > 0) ar_coefs else NULL,
          ma = if(ma_order > 0) ma_coefs else NULL
        ), n = n - d, innov = res_b, n.start = max(ar_order, ma_order))

        # Інтегрувати назад
        x_b <- diffinv(diff_x, differences = d)

        # Додати перехоплення, якщо потрібно
        if (include_mean) {
          x_b <- x_b + intercept
        }
      } else {
        # Для AR, MA, ARMA моделей
        res_b <- sample(residuals[!is.na(residuals)], size = n, replace = TRUE)

        # Витягнути коефіцієнти AR і MA
        if (ar_order > 0) {
          ar_coefs <- model@coefficients[1:ar_order]
        } else {
          ar_coefs <- numeric(0)
        }

        if (ma_order > 0) {
          ma_coefs <- model@coefficients[(ar_order+1):(ar_order+ma_order)]
        } else {
          ma_coefs <- numeric(0)
        }

        # Симулювати процес ARMA
        x_b <- arima.sim(model = list(
          ar = if(ar_order > 0) ar_coefs else NULL,
          ma = if(ma_order > 0) ma_coefs else NULL
        ), n = n, innov = res_b, n.start = max(ar_order, ma_order))

        # Додати перехоплення, якщо потрібно
        if (include_mean) {
          x_b <- x_b + intercept
        }
      }

      return(x_b)
    }

    boot_function <- function(b) {
      # Генерувати новий ряд
      x_b <- generate_with_residuals(object, res)

      # Визначаємо правильний формат order в залежності від типу моделі
      if(model_type == "ar") {
        boot_order <- ar_order  # Для AR моделей - одне число
      } else if(model_type == "ma") {
        boot_order <- ma_order  # Для MA моделей - одне число
      } else if(model_type == "arma") {
        boot_order <- c(ar_order, ma_order)  # Для ARMA - вектор довжини 2
      } else if(model_type == "arima") {
        boot_order <- c(ar_order, d, ma_order)  # Для ARIMA - вектор довжини 3
      } else {
        stop("Невідомий тип моделі: ", model_type)
      }

      if(debug && b == 1) {
        cat("Бутстреп реплікація 1:\n")
        cat("Тип моделі:", model_type, "\n")
        cat("Порядок для бутстрепу:", paste(boot_order, collapse=", "), "\n")
      }

      # Підігнати модель на бутстреп-ряді
      fit_b <- tryCatch({
        ts_pmm2(x_b, order = boot_order,
                model_type = model_type,
                include.mean = include_mean)
      }, error = function(e) {
        warning("Бутстреп-реплікація ", b, " не вдалася: ", e$message)
        return(NULL)
      })

      if (!is.null(fit_b)) {
        return(fit_b@coefficients)
      } else {
        return(rep(NA, length(coefs)))
      }
    }
  }

  # Виконати бутстреп: паралельно або послідовно
  use_parallel <- parallel && requireNamespace("parallel", quietly = TRUE)

  if (use_parallel) {
    if (is.null(cores)) {
      cores <- max(1, parallel::detectCores() - 1)
    }

    boot_results <- parallel::mclapply(seq_len(B), function(b) {
      boot_function(b)
    }, mc.cores = cores)

    # Перетворити список на матрицю
    boot_est <- do.call(rbind, boot_results)
  } else {
    # Послідовні обчислення
    boot_est <- matrix(0, nrow = B, ncol = length(coefs))
    colnames(boot_est) <- names(coefs)

    # Відстеження прогресу
    pb <- NULL
    if (interactive() && B > 10) {
      if (requireNamespace("utils", quietly = TRUE)) {
        pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
      }
    }

    for (b in seq_len(B)) {
      boot_est[b, ] <- boot_function(b)

      # Оновити індикатор прогресу
      if (!is.null(pb)) utils::setTxtProgressBar(pb, b)
    }

    # Закрити індикатор прогресу
    if (!is.null(pb)) close(pb)
  }

  # Видалити рядки зі значеннями NA
  na_rows <- apply(boot_est, 1, function(row) any(is.na(row)))
  if (any(na_rows)) {
    warning("Видалено ", sum(na_rows), " бутстреп-реплікацій через помилки оцінювання")
    boot_est <- boot_est[!na_rows, , drop = FALSE]
  }

  # Перевірити, чи маємо достатньо успішних бутстрепів
  if (nrow(boot_est) < 10) {
    stop("Замало успішних бутстреп-реплікацій для обчислення надійного висновку")
  }

  # Перевірити, чи є NaN або Inf значення
  if (any(is.nan(boot_est)) || any(is.infinite(boot_est))) {
    warning("Виявлено NaN або нескінченні значення в бутстреп-реплікаціях. Замінюємо їх на NA.")
    boot_est[is.nan(boot_est) | is.infinite(boot_est)] <- NA
  }

  # Обчислити коваріаційну матрицю та стандартні помилки
  cov_mat <- cov(boot_est, use = "pairwise.complete.obs")
  est <- coefs
  se <- sqrt(diag(cov_mat))

  # Обчислити t-значення та p-значення
  t_val <- est / se
  p_val <- 2 * (1 - pnorm(abs(t_val)))

  # Створити вихідний фрейм даних
  out <- data.frame(
    Estimate = est,
    Std.Error = se,
    t.value = t_val,
    p.value = p_val
  )

  # Додати імена AR і MA параметрів
  param_names <- c()
  if (ar_order > 0) {
    param_names <- c(param_names, paste0("ar", 1:ar_order))
  }
  if (ma_order > 0) {
    param_names <- c(param_names, paste0("ma", 1:ma_order))
  }
  rownames(out) <- param_names

  # Обчислити довірчі інтервали
  ci <- t(apply(boot_est, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE))
  colnames(ci) <- c("2.5%", "97.5%")

  # Додати довірчі інтервали до виходу
  out$conf.low <- ci[, "2.5%"]
  out$conf.high <- ci[, "97.5%"]

  return(out)
}
