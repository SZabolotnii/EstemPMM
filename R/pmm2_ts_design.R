# pmm2_ts_design.R - Funktsii dlia roboty z dyzain-matrytsiamy chasovykh riadiv

#' Perevirka ta pidhotovka parametriv chasovoho riadu
#'
#' @param x Dani chasovoho riadu
#' @param order Spetsyfikatsiia poriadku modeli
#' @param model_type Typ modeli (ar, ma, arma, abo arima)
#' @param include.mean Chy vkliuchaty serednie/perekhoplennia
#'
#' @return Spysok perevirenykh parametriv ta informatsii pro model
#' @keywords internal
validate_ts_parameters <- function(x, order, model_type, include.mean) {
  # Perevirka vkhidnykh danykh
  if (missing(x)) {
    stop("Vidsutnii arhument 'x'")
  }

  # Peretvoryty vkhidni dani na chyslovyi vektor
  x <- as.numeric(x)

  if (!is.numeric(x)) {
    stop("'x' maie buty chyslovym vektorom")
  }

  # Pereviryty na NA ta neskinchenni znachennia
  if (any(is.na(x)) || any(is.infinite(x))) {
    warning("Vyiavleno NA abo neskinchenni znachennia u vkhidnomu riadi. Vony budut vydaleni.")
    x <- x[!is.na(x) & !is.infinite(x)]
    if (length(x) < 10) {
      stop("Zamalo validnykh sposterezhen pislia vydalennia NA/neskinchennykh znachen")
    }
  }

  if (missing(order)) {
    stop("Vidsutnii arhument 'order'")
  }

  # Rozbir parametra order v zalezhnosti vid model_type
  if (model_type == "ar") {
    if (!is.numeric(order) || length(order) != 1)
      stop("Dlia AR modelei 'order' maie buty odnym tsilym chyslom")
    ar_order <- as.integer(order)
    ma_order <- 0
    d <- 0
    if (ar_order <= 0)
      stop("Poriadok AR maie buty dodatnym")

    # Pereviryty, chy dostatno danykh
    if (length(x) <= ar_order + 1) {
      stop("Zamalo sposterezhen dlia modeli AR poriadku ", ar_order)
    }
  } else if (model_type == "ma") {
    if (!is.numeric(order) || length(order) != 1)
      stop("Dlia MA modelei 'order' maie buty odnym tsilym chyslom")
    ar_order <- 0
    ma_order <- as.integer(order)
    d <- 0
    if (ma_order <= 0)
      stop("Poriadok MA maie buty dodatnym")

    # Pereviryty, chy dostatno danykh
    if (length(x) <= ma_order + 1) {
      stop("Zamalo sposterezhen dlia modeli MA poriadku ", ma_order)
    }
  } else if (model_type == "arma") {
    if (!is.numeric(order) || length(order) != 2)
      stop("Dlia ARMA modelei 'order' maie buty vektorom dovzhyny 2 (poriadok AR, poriadok MA)")
    ar_order <- as.integer(order[1])
    ma_order <- as.integer(order[2])
    d <- 0
    if (ar_order < 0 || ma_order < 0)
      stop("Poriadky AR ta MA maiut buty nevid'iemnymy")
    if (ar_order == 0 && ma_order == 0)
      stop("Prynaimni odyn z poriadkiv AR abo MA maie buty dodatnym")

    # Pereviryty, chy dostatno danykh
    if (length(x) <= max(ar_order, ma_order) + 1) {
      stop("Zamalo sposterezhen dlia modeli ARMA poriadkiv (", ar_order, ",", ma_order, ")")
    }
  } else if (model_type == "arima") {
    if (!is.numeric(order) || length(order) != 3)
      stop("Dlia ARIMA modelei 'order' maie buty vektorom dovzhyny 3 (poriadok AR, dyferentsiiuvannia, poriadok MA)")
    ar_order <- as.integer(order[1])
    d <- as.integer(order[2])
    ma_order <- as.integer(order[3])
    if (ar_order < 0 || ma_order < 0 || d < 0)
      stop("Poriadky AR, dyferentsiiuvannia ta MA maiut buty nevid'iemnymy")
    if (ar_order == 0 && ma_order == 0 && d == 0)
      stop("Prynaimni odyn z poriadkiv AR, dyferentsiiuvannia abo MA maie buty dodatnym")

    # Pereviryty, chy dostatno danykh pislia dyferentsiiuvannia
    if (length(x) <= d + max(ar_order, ma_order) + 1) {
      stop("Zamalo sposterezhen dlia modeli ARIMA pislia dyferentsiiuvannia")
    }
  } else {
    stop("Nevidomyi typ modeli: ", model_type)
  }

  # Zberehty oryhinalnyi riad
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

#' Stvorennia matrytsi dyzainu dlia AR modeli
#'
#' @param x tsentrovanyi chasovyi riad
#' @param p poriadok AR
#' @return Matrytsia dyzainu z lahovymy znachenniamy
#' @keywords internal
create_ar_matrix <- function(x, p) {
  n <- length(x)
  if (n <= p) {
    stop("Nedostatno tochok danykh dlia AR poriadku p = ", p)
  }

  nr <- n - p
  M <- matrix(0, nr, p)
  for (i in seq_len(p)) {
    M[, i] <- x[(p - i + 1):(n - i)]
  }
  M
}

#' Otrymaty otsinky Yula-Volkera dlia AR(p)
#'
#' @param x chyslovyi vektor
#' @param p tsile znachennia poriadku AR
#' @return chyslovyi vektor dovzhyny p (koefitsiienty AR)
#' @keywords internal
get_yw_estimates <- function(x, p) {
  # Tse sproshchenyi pidkhid, iakyi mozhe ne obrobliaty hranychni vypadky
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

#' Stvorennia dyzain-matrytsi dlia chasovykh riadiv
#'
#' @param x Dani chasovoho riadu
#' @param model_info Spysok z parametramy modeli
#' @param innovations Optsionalni innovatsii/zalyshky dlia MA komponentiv
#'
#' @return Spysok z dyzain-matrytseiu, zminnoiu vidhuku ta inshymy komponentamy
#' @keywords internal
create_ts_design_matrix <- function(x, model_info, innovations = NULL) {
  # Vytiahnennia parametriv modeli
  ar_order <- model_info$ar_order
  ma_order <- model_info$ma_order
  d <- model_info$d
  model_type <- model_info$model_type
  include_mean <- model_info$include.mean

  # Perevirka, chy potribno dyferentsiiuvaty riad
  if (model_type == "arima" && d > 0) {
    x_diff <- diff(x, differences = d)
  } else {
    x_diff <- x
  }

  # Obrobka serednoho
  if (include_mean) {
    x_mean <- mean(x_diff, na.rm = TRUE)
    x_centered <- x_diff - x_mean
  } else {
    x_mean <- 0
    x_centered <- x_diff
  }

  # Obchyslennia maksymalnoho lahu ta efektyvnoi dovzhyny danykh
  max_lag <- max(ar_order, ma_order)
  n_data <- length(x_centered)

  if (n_data <= max_lag) {
    stop("Nedostatno danykh dlia pobudovy modeli pislia dyferentsiiuvannia")
  }

  n_rows <- n_data - max_lag
  n_cols <- ar_order + ma_order

  # Stvorennia dyzain-matrytsi z vidpovidnymy rozmiramy
  X <- matrix(0, nrow = n_rows, ncol = n_cols)
  y <- x_centered[(max_lag + 1):n_data]

  # Dodavannia AR komponentiv
  if (ar_order > 0) {
    col_index <- 1
    for (i in 1:ar_order) {
      X[, col_index] <- x_centered[(max_lag - i + 1):(n_data - i)]
      col_index <- col_index + 1
    }
  }

  # Dodavannia MA komponentiv, iakshcho potribno
  if (ma_order > 0) {
    # Yakshcho innovatsii ne nadani, vykorystovuiemo nuli
    if (is.null(innovations)) {
      innovations <- rep(0, n_data)
    }

    # Zabezpechennia pravylnoi dovzhyny innovatsii
    if (length(innovations) < n_data) {
      innovations <- c(rep(0, n_data - length(innovations)), innovations)
    } else if (length(innovations) > n_data) {
      innovations <- tail(innovations, n_data)
    }

    # Zapovnennia MA stovptsiv
    col_index <- ar_order + 1
    for (j in 1:ma_order) {
      X[, col_index] <- innovations[(max_lag - j + 1):(n_data - j)]
      col_index <- col_index + 1
    }
  }

  # Povernennia rezultativ u vyhliadi spysku
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

#' Otrymaty pochatkovi otsinky parametriv dlia modelei chasovykh riadiv
#'
#' @param model_params Perevireni parametry modeli z validate_ts_parameters
#' @param initial Optsionalno nadani korystuvachem pochatkovi otsinky
#' @param method Metod otsiniuvannia
#' @param verbose Vyvodyty detalnu informatsiiu
#'
#' @return Spysok, shcho mistyt:
#'   \item{b_init}{vektor pochatkovykh koefitsiientiv AR/MA}
#'   \item{x_mean}{otsinene serednie (iakshcho include.mean=TRUE)}
#'   \item{innovations}{pochatkovi zalyshky/innovatsii}
#'   \item{x_centered}{tsentrovanyi (abo dyferentsiiovanyi + tsentrovanyi) riad}
#'   \item{m2}{druhyi tsentralnyi moment pochatkovykh zalyshkiv}
#'   \item{m3}{tretii tsentralnyi moment pochatkovykh zalyshkiv}
#'   \item{m4}{chetvertyi tsentralnyi moment pochatkovykh zalyshkiv}
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

  # Mozhlyvo dyferentsiiuvaty dlia ARIMA
  if (mtype == "arima" && d > 0) {
    x_diff <- diff(x, differences = d)
  } else {
    x_diff <- x
  }

  # Tsentruvaty, iakshcho potribno
  if (inc_mean) {
    x_mean <- mean(x_diff, na.rm = TRUE)
    x_centered <- x_diff - x_mean
  } else {
    x_mean <- 0
    x_centered <- x_diff
  }

  if (mtype == "ar") {
    # AR(p): shvydkyi pidkhid dlia pochatkovykh znachen
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
        stop("Dovzhyna 'initial' maie vidpovidaty poriadku AR")
      }
      b_init <- initial
    }
    # Innovatsii z pochatkovoi pidhonky
    X <- create_ar_matrix(x_centered, ar_order)
    y <- x_centered[(ar_order + 1):length(x_centered)]
    innovations <- as.numeric(y - X %*% b_init)

  } else if (mtype %in% c("ma", "arma", "arima")) {
    # Vykorystovuvaty stats::arima dlia pochatkovoho prypushchennia abo nadani korystuvachem
    arima_order <- c(ar_order, ifelse(mtype=="arima", d, 0), ma_order)
    if (is.null(initial)) {
      init_fit <- NULL
      try_methods <- c("CSS-ML","ML","CSS")

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
        if(verbose) cat("Usi standartni metody ne spratsiuvaly; vykorystovuiutsia sproshcheni znachennia.\n")
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
      # Nadano initial
      if(is.list(initial)) {
        if(ar_order>0 && is.null(initial$ar)) {
          stop("Vidsutnii 'ar' u spysku initial, ale ar_order>0")
        }
        if(ma_order>0 && is.null(initial$ma)) {
          stop("Vidsutnii 'ma' u spysku initial, ale ma_order>0")
        }
        ar_init <- if(ar_order>0) initial$ar else numeric(0)
        ma_init <- if(ma_order>0) initial$ma else numeric(0)
      } else {
        if(length(initial) != (ar_order+ma_order)) {
          stop("Dovzhyna 'initial' maie vidpovidaty sumi poriadkiv AR ta MA")
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
        if(verbose) cat("Pomylka z nadanymy korystuvachem pochatkovymy znachenniamy:",e$message,"\n")
        list(residuals = if(mtype=="arima") x_diff else x_centered)
      })
      innovations <- as.numeric(init_fit$residuals)
    }
  }

  if(anyNA(b_init)) {
    warning("NA v pochatkovykh parametrakh zamineni na 0.")
    b_init[is.na(b_init)] <- 0
  }

  # Obchyslennia momentiv
  moments <- compute_moments(innovations)

  # Povernennia rezultativ
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

#' Onovyty innovatsii MA modeli
#'
#' @param x tsentrovanyi chasovyi riad
#' @param ma_coef vektor koefitsiientiv MA
#' @return vektor innovatsii
#' @keywords internal
update_ma_innovations <- function(x, ma_coef) {
  n <- length(x)
  q <- length(ma_coef)

  # Initsializuvaty innovatsii iak nuli
  innovations <- numeric(n)

  # Iteratyvno obchyslyty innovatsii
  for(t in 1:n) {
    # Obchyslyty ochikuvane znachennia na osnovi poperednikh innovatsii
    expected <- 0
    for(j in 1:q) {
      if(t - j > 0) {
        expected <- expected + ma_coef[j] * innovations[t - j]
      }
    }

    # Obchyslyty potochnu innovatsiiu
    innovations[t] <- x[t] - expected
  }

  # Pereviryty na nekonechni znachennia
  if(any(is.infinite(innovations)) || any(is.na(innovations))) {
    warning("Vyiavleno nekonechni innovatsii v update_ma_innovations. Vykorystovuiemo rehuliaryzovani znachennia.")
    # Zaminyty problemni znachennia na serednie abo 0
    bad_idx <- is.infinite(innovations) | is.na(innovations)
    if(sum(!bad_idx) > 0) {
      # Yakshcho ie diisni znachennia, vykorystovuiemo ikh serednie
      innovations[bad_idx] <- mean(innovations[!bad_idx])
    } else {
      # Inakshe vykorystovuiemo 0
      innovations[bad_idx] <- 0
    }
  }

  return(innovations)
}

#' Obchyslyty kintsevi zalyshky dlia modelei chasovykh riadiv
#'
#' @param coefs Otsineni koefitsiienty
#' @param model_info Informatsiia pro model
#' @return Vektor zalyshkiv
#' @keywords internal
compute_ts_residuals <- function(coefs, model_info) {
  # Vytiahnuty parametry modeli
  x           <- model_info$x
  ar_order    <- model_info$ar_order
  ma_order    <- model_info$ma_order
  d           <- model_info$d
  model_type  <- model_info$model_type
  include.mean <- model_info$include.mean
  verbose_flag <- isTRUE(model_info$verbose)

  # Rozdilyty koefitsiienty na AR ta MA chastyny
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

  # Dlia bilsh nadiinoho obchyslennia zalyshkiv vykorystovuiemo arima z fiksovanymy parametramy
  arima_order <- c(ar_order, ifelse(model_type == "arima", d, 0), ma_order)

  # Pidhotuvaty fiksovani parametry dlia arima
  fixed_params <- c(ar_coefs, ma_coefs)
  include_intercept <- include.mean && !(model_info$model_type == "arima" && model_info$d > 0)

  if (include_intercept) {
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

  # Obchyslyty zalyshky z fiksovanymy parametramy
  final_fit <- tryCatch({
    stats::arima(model_info$original_x,
                order = arima_order,
                fixed = fixed_params,
                include.mean = include_intercept)
  }, error = function(e) {
    if (verbose_flag) cat("Pomylka pry obchyslenni kintsevykh zalyshkiv:", e$message, "\n")
    list(residuals = rep(NA, length(model_info$original_x)))
  })

  # Vytiahnuty zalyshky ta zabezpechyty pravylnu dovzhynu
  final_res <- as.numeric(final_fit$residuals)

  if (length(final_res) < length(model_info$original_x)) {
    final_res <- c(rep(NA, length(model_info$original_x) - length(final_res)), final_res)
  }

  return(final_res)
}
