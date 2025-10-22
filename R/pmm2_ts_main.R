# pmm2_ts_main.R - Unifikovanyi modul dlia modelei chasovykh riadiv z PMM2

#' Pidihnaty model chasovoho riadu za dopomohoiu metodu PMM2
#'
#' @param x Chyslovyi vektor danykh chasovoho riadu
#' @param order Spetsyfikatsiia poriadku modeli:
#'        - Dlia AR modelei: odne tsile chyslo (poriadok AR)
#'        - Dlia MA modelei: odne tsile chyslo (poriadok MA)
#'        - Dlia ARMA modelei: vektor c(p, q) (poriadky AR ta MA)
#'        - Dlia ARIMA modelei: vektor c(p, d, q) (poriadky AR, dyferentsiiuvannia ta MA)
#' @param model_type Riadok, shcho vyznachaie typ modeli: "ar", "ma", "arma", abo "arima"
#' @param method Riadok: metod otsiniuvannia, odyn z "pmm2" (za zamovchuvanniam), "css", "ml", "yw", "ols"
#' @param max_iter Tsile chyslo: maksymalna kilkist iteratsii dlia alhorytmu
#' @param tol Chyslove: dopusk dlia zbizhnosti
#' @param include.mean Lohichne: chy vkliuchaty chlen serednoho (perekhoplennia)
#' @param initial Spysok abo vektor pochatkovykh otsinok parametriv (optsionalno)
#' @param na.action Funktsiia dlia obrobky vidsutnikh znachen, za zamovchuvanniam - na.fail
#' @param regularize Lohichne, dodavaty mali znachennia do diahonali dlia chyslovoi stabilnosti
#' @param reg_lambda Parametr rehuliaryzatsii (iakshcho regularize=TRUE)
#' @param verbose Lohichne: chy vyvodyty informatsiiu pro prohres
#'
#' @details
#' Alhorytm PMM2 pratsiuie nastupnym chynom:
#'
#' 1. Pidhaniaie pochatkovu model za dopomohoiu standartnoho metodu (OLS, Yula-Volkera, CSS abo ML)
#' 2. Obchysliuie tsentralni momenty (m2, m3, m4) z pochatkovykh zalyshkiv/innovatsii
#' 3. Vykorystovuie tsi momenty zi spetsializovanym rozv'iazuvachem (pmm2_algorithm) dlia znakhodzhennia
#'    robastnykh otsinok parametriv
#'
#' @return Ob'iekt S4 \code{TS2fit} vidpovidnoho pidklasu
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

  # 1) Perevirka vkhidnykh danykh
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

  if (model_params$model_type == "arima") {
    p <- model_params$ar_order
    d <- model_params$d
    q <- model_params$ma_order

    css_fit <- tryCatch(
      stats::arima(model_params$original_x,
                   order = c(p, d, q),
                   method = "CSS-ML",
                   include.mean = include.mean),
      error = function(e) NULL
    )

    if (is.null(css_fit)) {
      stop("Ne vdalosia otsinyty ARIMA model klasychnym metodom")
    }

    coef_names <- names(css_fit$coef)
    ar_css <- if (p > 0) as.numeric(css_fit$coef[paste0("ar", seq_len(p))]) else numeric(0)
    ma_css <- if (q > 0) as.numeric(css_fit$coef[paste0("ma", seq_len(q))]) else numeric(0)

    intercept_css <- 0
    if (include.mean && d == 0) {
      intercept_name <- setdiff(coef_names,
                                c(paste0("ar", seq_len(p)), paste0("ma", seq_len(q))))
      if (length(intercept_name) > 0) {
        intercept_css <- as.numeric(css_fit$coef[intercept_name[1]])
      }
    }

    residuals_css <- as.numeric(css_fit$residuals)
    residuals_css[is.na(residuals_css)] <- 0

    x_diff <- if (d > 0) diff(model_params$original_x, differences = d)
              else model_params$original_x
    res_diff <- tail(residuals_css, length(x_diff))

    include_intercept_diff <- include.mean && d == 0

    if (method == "css") {
      res_clean <- res_diff[is.finite(res_diff)]
      moments <- compute_moments(res_clean)
      return(new("ARIMAPMM2",
                 coefficients    = as.numeric(c(ar_css, ma_css)),
                 residuals       = as.numeric(residuals_css),
                 m2              = as.numeric(moments$m2),
                 m3              = as.numeric(moments$m3),
                 m4              = as.numeric(moments$m4),
                 convergence     = TRUE,
                 iterations      = 1L,
                 call            = cl,
                 model_type      = "arima",
                 intercept       = if (include_intercept_diff) intercept_css else 0,
                 original_series = as.numeric(model_params$original_x),
                 order           = list(ar = p, ma = q, d = d)))
    }

    design <- arma_build_design(x_diff, res_diff,
                                 p = p, q = q,
                                 intercept = intercept_css,
                                 include_intercept = include_intercept_diff)

    moments <- compute_moments(res_diff[is.finite(res_diff)])
    b_init <- c(if (include_intercept_diff) 0 else NULL, ar_css, ma_css)

    algo_res <- pmm2_algorithm(
      b_init = b_init,
      X = design$X,
      y = design$y,
      m2 = moments$m2,
      m3 = moments$m3,
      m4 = moments$m4,
      max_iter = max_iter,
      tol = tol,
      regularize = regularize,
      reg_lambda = reg_lambda,
      verbose = verbose
    )

    if (include_intercept_diff) {
      intercept_hat <- algo_res$b[1]
      ar_hat <- if (p > 0) algo_res$b[1 + seq_len(p)] else numeric(0)
      ma_hat <- if (q > 0) algo_res$b[1 + p + seq_len(q)] else numeric(0)
    } else {
      intercept_hat <- 0
      ar_hat <- if (p > 0) algo_res$b[seq_len(p)] else numeric(0)
      ma_hat <- if (q > 0) algo_res$b[p + seq_len(q)] else numeric(0)
    }

    x_mean <- if (include_intercept_diff) intercept_hat else 0
    x_centered <- x_diff - x_mean

    model_info <- list(
      ar_order = p,
      ma_order = q,
      d = d,
      model_type = "arima",
      include.mean = include.mean,
      innovations = res_diff,
      x = x_centered,
      x_mean = x_mean,
      original_x = model_params$original_x,
      verbose = verbose
    )

    final_coef <- c(ar_hat, ma_hat)
    final_res <- compute_ts_residuals(final_coef, model_info)
    res_clean <- final_res[is.finite(final_res)]
    moments_final <- compute_moments(res_clean)

    return(new("ARIMAPMM2",
               coefficients    = as.numeric(final_coef),
               residuals       = as.numeric(final_res),
               m2              = as.numeric(moments_final$m2),
               m3              = as.numeric(moments_final$m3),
               m4              = as.numeric(moments_final$m4),
               convergence     = algo_res$convergence,
               iterations      = as.numeric(algo_res$iterations),
               call            = cl,
               model_type      = "arima",
               intercept       = if (include_intercept_diff) as.numeric(intercept_hat) else 0,
               original_series = as.numeric(model_params$original_x),
               order           = list(ar = p, ma = q, d = d)))
  }

  # 2) Otrymannia pochatkovykh otsinok
  init <- get_initial_estimates(model_params, initial, method, verbose)
  b_init      <- init$b_init
  x_mean      <- init$x_mean
  innovations <- init$innovations
  x_centered  <- init$x_centered
  orig_x      <- init$orig_x
  m2          <- init$m2
  m3          <- init$m3
  m4          <- init$m4

  # 3) Stvorennia dyzain-matrytsi
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

  # 4) Yakshcho metod == "pmm2", vykorystovuvaty alhorytm PMM2
  if (method == "pmm2") {
    if (verbose) cat("Pochatok optymizatsii PMM2...\n")

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

    # Obchyslyty kintsevi zalyshky
    final_res <- compute_ts_residuals(final_coef, model_info)
  } else {
    # Dlia inshykh metodiv prosto vykorystovuvaty pochatkovi otsinky
    final_coef <- b_init
    converged <- TRUE
    iterations <- 0
    final_res <- innovations
  }

  # 5) Stvoryty vidpovidnyi ob'iekt klasu
  if (model_type == "ar") {
    result_class <- "ARPMM2"
  } else if (model_type == "ma") {
    result_class <- "MAPMM2"
  } else if (model_type == "arma") {
    result_class <- "ARMAPMM2"
  } else if (model_type == "arima") {
    result_class <- "ARIMAPMM2"
  } else {
    result_class <- "TS2fit"  # Bazovyi klas za zamovchuvanniam
  }

  # Stvoryty ta povernuty ob'iekt vidpovidnoho klasu
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

#' Pidihnaty AR model za dopomohoiu PMM2 (obhortka)
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

#' Pidihnaty MA model za dopomohoiu PMM2 (obhortka)
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

#' Pidihnaty ARMA model za dopomohoiu PMM2 (obhortka)
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

#' Pidihnaty ARIMA model za dopomohoiu PMM2 (obhortka)
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
      cat("Ne vdalos otsinyty MA model cherez stats::arima: povertaiu nulovi koefitsiienty\n")
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
      if (verbose) cat("Systema synhuliarna na iteratsii", iter, "\n")
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

arma_build_design <- function(x, residuals, p, q, intercept = 0, include_intercept = FALSE) {
  n <- length(x)
  max_lag <- max(p, q)
  if (n <= max_lag) {
    stop("Nedostatno danykh dlia pobudovy ARMA dyzain-matrytsi")
  }

  idx <- seq.int(max_lag + 1L, n)
  columns <- list()
  if (include_intercept) {
    columns <- c(columns, list(rep(1, length(idx))))
  }
  if (p > 0) {
    for (j in seq_len(p)) {
      columns <- c(columns, list(x[idx - j]))
    }
  }
  if (q > 0) {
    for (j in seq_len(q)) {
      columns <- c(columns, list(residuals[idx - j]))
    }
  }

  X <- if (length(columns) > 0) do.call(cbind, columns) else matrix(0, length(idx), 0)
  y <- x[idx] - intercept
  list(X = X, y = y)
}
