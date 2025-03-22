# pmm_utils.R

#' .pmm2_fit: unified implementation for PMM2
#'
#' @param b_init initial parameter estimates (typically from OLS)
#' @param X design matrix
#' @param y response vector
#' @param m2,m3,m4 central moments
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
#' @param regularize add small value to diagonal for numerical stability
#'
#' @keywords internal
.pmm2_fit <- function(b_init, X, y,
                      m2, m3, m4,
                      max_iter=50, tol=1e-6,
                      regularize=TRUE, reg_lambda = 1e-8) {
  n <- nrow(X)
  p <- ncol(X)

  # Обчислюємо A, B, C - для всього набору даних одразу
  A <- m3
  B <- (m4 - m2^2) - 2*m3*y
  C <- m3*(y^2) - ((m4 - m2^2)*y) - m2*m3

  # Поточні оцінки параметрів
  b_cur <- b_init
  converged <- FALSE

  for (iter in seq_len(max_iter)) {
    # Обчислюємо Yx = X %*% b_cur (предиктор)
    Yx <- X[,1]*b_cur[1]
    if (p > 1) {
      for (r in 2:p) {
        Yx <- Yx + b_cur[r]*X[,r]
      }
    }

    # Обчислюємо Z1 = A*Yx^2 + B*Yx + C
    Z1 <- A*(Yx^2) + B*Yx + C

    # Формуємо Z - вектор рівнянь
    Z <- numeric(p)
    Z[1] <- sum(Z1)
    for (r in 2:p) {
      Z[r] <- sum(Z1 * X[,r])
    }

    # Формуємо JZs - матрицю Якобіана
    JZs <- matrix(NA, p, p)
    JZ11 <- 2*A*Yx + B

    JZs[1,1] <- sum(JZ11)
    for (ii in 2:p) {
      tmp <- JZ11*X[,ii]
      JZs[1,ii] <- sum(tmp)
      JZs[ii,1] <- JZs[1,ii]
    }
    for (ii in 2:p) {
      for (jj in 2:p) {
        tmp <- JZ11 * X[,ii] * X[,jj]
        JZs[ii,jj] <- sum(tmp)
      }
    }

    # Невелика регуляризація для уникнення сингулярності якщо потрібно
    if (regularize) {
      diag(JZs) <- diag(JZs) + reg_lambda
    }

    # Розв'язуємо систему JZs * delta = Z
    step <- try(solve(JZs, Z), silent=TRUE)
    if (inherits(step, "try-error")) {
      warning("PMM2: solve(JZs,Z) failed => breaking iteration")
      break
    }

    # Оновлюємо параметри
    b_new <- b_cur - step
    diff_par <- sqrt(sum((b_new - b_cur)^2))
    b_cur <- b_new

    # Перевіряємо збіжність
    if (diff_par < tol) {
      converged <- TRUE
      break
    }
  }

  list(b=b_cur, convergence=converged)
}
