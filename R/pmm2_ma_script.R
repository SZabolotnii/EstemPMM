# pmm2_ma_script.R - Альтернативна реалізація PMM2 для MA моделей за шаблоном користувача

#' Оцінка MA(q) моделі з використанням алгоритму, що повторює користувацький скрипт
#'
#' @param x Часовий ряд
#' @param order порядок MA моделі (ціле число)
#' @param method "css" або "pmm2"
#' @param include.mean логічний прапорець, чи оцінювати перехоплення
#' @param max_iter максимальна кількість ітерацій для PMM2
#' @param tol допуск збіжності
#' @param verbose чи друкувати повідомлення
#'
#' @keywords internal
ma_pmm2_script <- function(x, order = 1, method = c("pmm2", "css"),
                           include.mean = TRUE, max_iter = 50, tol = 1e-6,
                           verbose = FALSE) {
  method <- match.arg(method)
  x <- as.numeric(x)
  q <- as.integer(order)
  if (q < 1) stop("order має бути >= 1")
  n <- length(x)
  if (n <= q + 1) stop("Замало спостережень для MA(q)")

  fit_css <- stats::arima(x, order = c(0, 0, q), method = "CSS-ML",
                          include.mean = include.mean)
  coef_css <- as.numeric(coef(fit_css)[seq_len(q)])
  intercept <- if (include.mean && ("intercept" %in% names(coef(fit_css))))
    as.numeric(coef(fit_css)["intercept"]) else 0
  residuals_css <- as.numeric(fit_css$residuals)
  residuals_css[is.na(residuals_css)] <- 0

  if (method == "css") {
    return(list(coefficients = coef_css,
                intercept = intercept,
                innovations = residuals_css,
                convergence = TRUE,
                iterations = 0L))
  }

  design <- build_design_for_ma_script(x, residuals_css, intercept, q)
  moments <- compute_moments(residuals_css)

  b_init <- c(0, coef_css)
  sol <- solve_ma_pmm2_script(b_init, design$X, design$y,
                              moments$m2, moments$m3, moments$m4,
                              max_iter = max_iter, tol = tol,
                              verbose = verbose)

  theta <- sol$coefficients
  final_innov <- compute_ma_innovations(x - intercept, theta, q, backcast = TRUE)

  list(coefficients = theta,
       intercept = intercept,
       innovations = final_innov,
       convergence = sol$convergence,
       iterations = sol$iterations)
}

build_design_for_ma_script <- function(x, residuals, intercept, q) {
  idx <- seq.int(q + 1L, length(x))
  X <- matrix(1, nrow = length(idx), ncol = q + 1L)
  for (j in seq_len(q)) {
    X[, j + 1L] <- residuals[idx - j]
  }
  y <- x[idx] - intercept
  list(X = X, y = y)
}

solve_ma_pmm2_script <- function(b_init, X, Y, m2, m3, m4,
                                 max_iter = 50, tol = 1e-6,
                                 verbose = FALSE) {
  b <- as.numeric(b_init)
  converged <- FALSE
  iterations <- 0L
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

compute_ma_innovations <- function(x, theta, q, backcast = TRUE) {
  n <- length(x)
  innovations <- numeric(n)
  if (q == 0L) return(x)
  history <- rep(0, q)
  for (t in seq_len(n)) {
    ma_component <- sum(theta * history)
    innovations[t] <- x[t] - ma_component
    history <- c(innovations[t], history)[seq_len(q)]
  }
  innovations
}
