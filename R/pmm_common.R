# pmm_common.R - Спільні утиліти для всіх PMM2 моделей


#' Універсальний алгоритм PMM2 для всіх типів моделей
#'
#' @param b_init Початкові оцінки параметрів
#' @param X Матриця дизайну
#' @param y Вектор відгуку
#' @param m2,m3,m4 Центральні моменти
#' @param max_iter Максимальна кількість ітерацій
#' @param tol Допуск для збіжності
#' @param regularize Чи додавати регуляризацію
#' @param reg_lambda Параметр регуляризації
#' @param verbose Чи виводити інформацію про прогрес
#'
#' @return Список з результатами оцінювання
#' @keywords internal
pmm2_algorithm <- function(b_init, X, y, m2, m3, m4,
                           max_iter = 50, tol = 1e-6,
                           regularize = TRUE, reg_lambda = 1e-8,
                           verbose = FALSE) {
  # Поточні оцінки параметрів
  b_cur <- b_init
  converged <- FALSE
  iterations <- 0

  # Обчислити коефіцієнти PMM2 полінома
  A <- m3
  B <- m4 - m2^2 - 2*m3*y
  C <- m3*y^2 - (m4 - m2^2)*y - m2*m3

  # Відстежувати історію збіжності, якщо verbose
  if (verbose) {
    conv_history <- numeric(max_iter)
  }

  # Основний цикл ітерацій
  for (iter in seq_len(max_iter)) {
    iterations <- iter

    # Обчислити прогнозовані значення
    y_pred <- as.vector(X %*% b_cur)

    # Обчислити Z1 = A*y_pred^2 + B*y_pred + C
    Z1 <- A*(y_pred^2) + B*y_pred + C

    # Сформувати вектор Z для кожного параметра
    p <- length(b_cur)
    Z <- numeric(p)
    for (r in 1:p) {
      Z[r] <- sum(Z1 * X[, r])
    }

    # Обчислити похідну JZ11 = 2*A*y_pred + B
    JZ11 <- 2*A*y_pred + B

    # Сформувати матрицю Якобіана
    JZs <- matrix(0, p, p)
    for (i in 1:p) {
      for (j in 1:p) {
        JZs[i, j] <- sum(JZ11 * X[, i] * X[, j])
      }
    }

    # Додати регуляризацію, якщо потрібно
    if (regularize) {
      diag(JZs) <- diag(JZs) + reg_lambda
    }

    # Розв'язати систему JZs * delta = Z
    delta <- tryCatch({
      solve(JZs, Z)
    }, error = function(e) {
      if (verbose) {
        cat("Помилка при розв'язанні лінійної системи:", conditionMessage(e), "\n")
        cat("Додаємо сильнішу регуляризацію\n")
      }
      diag(JZs) <- diag(JZs) + 1e-4
      solve(JZs, Z)
    })

    # Оновити параметри
    b_new <- b_cur - delta
    diff_par <- sqrt(sum((b_new - b_cur)^2))

    # Зберегти історію збіжності, якщо verbose
    if (verbose) {
      conv_history[iter] <- diff_par
      if (iter %% 5 == 0 || iter == 1) {
        cat("Ітерація", iter, ": Зміна параметрів =",
            formatC(diff_par, digits = 8), "\n")
      }
    }

    b_cur <- b_new

    # Перевірити збіжність
    if (diff_par < tol) {
      converged <- TRUE
      if (verbose) cat("Збіжність досягнута після", iter, "ітерацій\n")
      break
    }
  }

  # Попередження, якщо досягнуто максимальну кількість ітерацій без збіжності
  if (!converged && verbose) {
    cat("Попередження: Алгоритм не збігся після", max_iter, "ітерацій\n")
  }

  # Обчислити кінцеві залишки
  final_res <- as.numeric(y - X %*% b_cur)

  # Побудувати історію збіжності, якщо verbose
  if (verbose && iterations > 1) {
    if (requireNamespace("graphics", quietly = TRUE)) {
      graphics::plot(1:iterations, conv_history[1:iterations], type = "b",
                     xlab = "Ітерація", ylab = "Зміна параметрів",
                     main = "Історія збіжності")
      graphics::abline(h = tol, col = "red", lty = 2)
    }
  }

  # Повернути результати
  list(
    b = as.numeric(b_cur),
    convergence = converged,
    iterations = iterations,
    residuals = final_res
  )
}
