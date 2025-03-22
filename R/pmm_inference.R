#' Bootstrap inference for a PMM2 fit
#'
#' @param object an object of class PMM2fit
#' @param formula,data the same formula and data used originally
#' @param B number of bootstrap replications
#' @param seed (optional) for reproducibility
#'
#' @return A data.frame with columns: Estimate, Std.Error, t.value, p.value
#' @export
pmm2_inference <- function(object, formula, data, B=200, seed=NULL) {
  if(!is.null(seed)) set.seed(seed)

  coefs <- object@coefficients
  res   <- object@residuals

  # будуємо X,y
  mf <- model.frame(formula, data)
  X <- model.matrix(formula, mf)
  y <- model.response(mf)
  n <- nrow(X)

  # матриця для зберігання результатів
  boot_est <- matrix(0, nrow=B, ncol=length(coefs))
  colnames(boot_est) <- names(coefs)

  for(b in seq_len(B)) {
    # 1) бутстреп залишків
    res_b <- sample(res, size=n, replace=TRUE)
    # 2) new y
    y_b <- X %*% coefs + res_b

    # 3) дані
    data_b <- data
    # припускаємо, що ліва змінна - перша в formula
    lhs <- as.character(formula[[2]])
    data_b[[lhs]] <- as.numeric(y_b)

    # 4) повторне оцінювання
    fit_b <- lm_pmm2(formula, data_b, max_iter=20, tol=1e-6)

    boot_est[b, ] <- fit_b@coefficients
  }

  cov_mat <- cov(boot_est)
  est <- coefs
  se  <- sqrt(diag(cov_mat))

  t_val <- est / se
  # Для великої вибірки -> нормальне наближення
  p_val <- 2*(1 - pnorm(abs(t_val)))

  out <- data.frame(
    Estimate  = est,
    Std.Error = se,
    t.value   = t_val,
    p.value   = p_val
  )
  rownames(out) <- names(est)
  out
}
