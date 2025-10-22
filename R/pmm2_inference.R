# pmm2_inference.R - Statystychnyi vysnovok dlia modelei PMM2

#' Butstrep-vysnovok dlia pidhonky PMM2
#'
#' @param object ob'iekt klasu PMM2fit
#' @param formula ta sama formula, shcho vykorystovuvalasia spochatku
#' @param data freim danykh, shcho vykorystovuvavsia spochatku
#' @param B kilkist butstrep-replikatsii
#' @param seed (optsionalno) dlia vidtvoriuvanosti
#' @param parallel lohichne, chy vykorystovuvaty paralelni obchyslennia
#' @param cores kilkist iader dlia vykorystannia pry paralelnykh obchyslenniakh, za zamovchuvanniam - avtovyznachennia
#'
#' @return data.frame z stovptsiamy: Estimate, Std.Error, t.value, p.value
#' @export
pmm2_inference <- function(object, formula, data, B=200, seed=NULL,
                           parallel=FALSE, cores=NULL) {
  # Vstanovyty zerno dlia vidtvoriuvanosti, iakshcho nadano
  if(!is.null(seed)) set.seed(seed)

  # Vytiahnuty koefitsiienty ta zalyshky
  coefs <- object@coefficients
  res   <- object@residuals

  # Perevirka vkhidnykh danykh
  if(B < 10) {
    warning("Kilkist butstrep-vybirok (B) duzhe mala. Rozhliante vykorystannia B >= 100 dlia bilsh nadiinoho vysnovku.")
  }

  if(!inherits(object, "PMM2fit")) {
    stop("Ob'iekt maie buty klasu 'PMM2fit'")
  }

  if(missing(formula) || missing(data)) {
    stop("Obydva 'formula' ta 'data' maiut buty nadani")
  }

  # Pobuduvaty matrytsi X, y
  mf <- model.frame(formula, data)
  X <- model.matrix(formula, mf)
  y <- model.response(mf)
  n <- nrow(X)

  # Rannie povernennia u vypadku pomylok
  if(is.null(y) || is.null(X)) {
    stop("Ne vdalosia vytiahnuty vidhuk abo matrytsiu dyzainu z danykh")
  }

  # Pereviryty, chy slid vykorystovuvaty paralelni obchyslennia
  use_parallel <- parallel && requireNamespace("parallel", quietly = TRUE)

  if(use_parallel) {
    if(is.null(cores)) {
      cores <- max(1, parallel::detectCores() - 1)
    }

    boot_results <- parallel::mclapply(seq_len(B), function(b) {
      # 1) Butstrep zalyshkiv
      res_b <- sample(res, size=n, replace=TRUE)

      # 2) Stvoryty novyi y
      y_b <- X %*% coefs + res_b

      # 3) Stvoryty novi dani
      data_b <- data
      # Prypustyty, shcho liva storona ie pershym terminom u formuli
      lhs <- as.character(formula[[2]])
      data_b[[lhs]] <- as.numeric(y_b)

      # 4) Povtorno otsinyty model
      fit_b <- tryCatch({
        lm_pmm2(formula, data_b, max_iter=20, tol=1e-6)
      }, error = function(e) {
        warning("Butstrep-replikatsiia ", b, " ne vdalasia: ", e$message)
        return(NULL)
      })

      if(!is.null(fit_b)) {
        return(fit_b@coefficients)
      } else {
        return(rep(NA, length(coefs)))
      }
    }, mc.cores = cores)

    # Peretvoryty spysok na matrytsiu
    boot_est <- do.call(rbind, boot_results)

  } else {
    # Poslidovni obchyslennia
    # Matrytsia dlia zberihannia rezultativ
    boot_est <- matrix(0, nrow=B, ncol=length(coefs))
    colnames(boot_est) <- names(coefs)

    # Vidstezhennia prohresu
    pb <- NULL
    if(interactive() && B > 10) {
      if(requireNamespace("utils", quietly = TRUE)) {
        pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
      }
    }

    for(b in seq_len(B)) {
      # 1) Butstrep zalyshkiv
      res_b <- sample(res, size=n, replace=TRUE)

      # 2) Stvoryty novyi y
      y_b <- X %*% coefs + res_b

      # 3) Stvoryty novi dani
      data_b <- data
      # Prypustyty, shcho liva storona ie pershym terminom u formuli
      lhs <- as.character(formula[[2]])
      data_b[[lhs]] <- as.numeric(y_b)

      # 4) Povtorno otsinyty model
      fit_b <- tryCatch({
        lm_pmm2(formula, data_b, max_iter=20, tol=1e-6)
      }, error = function(e) {
        warning("Butstrep-replikatsiia ", b, " ne vdalasia: ", e$message)
        return(NULL)
      })

      if(!is.null(fit_b)) {
        boot_est[b, ] <- fit_b@coefficients
      } else {
        boot_est[b, ] <- NA
      }

      # Onovyty indykator prohresu
      if(!is.null(pb)) utils::setTxtProgressBar(pb, b)
    }

    # Zakryty indykator prohresu
    if(!is.null(pb)) close(pb)
  }

  # Vydalyty riadky zi znachenniamy NA
  na_rows <- apply(boot_est, 1, function(row) any(is.na(row)))
  if(any(na_rows)) {
    warning("Vydaleno ", sum(na_rows), " butstrep-replikatsii cherez pomylky otsiniuvannia")
    boot_est <- boot_est[!na_rows, , drop = FALSE]
  }

  # Pereviryty, chy maiemo dostatno uspishnykh butstrepiv
  if(nrow(boot_est) < 10) {
    stop("Zamalo uspishnykh butstrep-replikatsii dlia obchyslennia nadiinoho vysnovku")
  }

  # Obchyslyty kovariatsiinu matrytsiu ta standartni pomylky
  cov_mat <- cov(boot_est)
  est <- coefs
  se  <- sqrt(diag(cov_mat))

  # Obchyslyty t-znachennia ta p-znachennia
  t_val <- est / se
  # Dlia velykykh vybirok vykorystovuvaty normalne nablyzhennia
  p_val <- 2 * (1 - pnorm(abs(t_val)))

  # Stvoryty vykhidnyi freim danykh
  out <- data.frame(
    Estimate  = est,
    Std.Error = se,
    t.value   = t_val,
    p.value   = p_val
  )
  rownames(out) <- names(est)

  # Obchyslyty dovirchi intervaly
  ci <- t(apply(boot_est, 2, quantile, probs = c(0.025, 0.975)))
  colnames(ci) <- c("2.5%", "97.5%")

  # Dodaty dovirchi intervaly do vykhodu
  out$conf.low <- ci[, "2.5%"]
  out$conf.high <- ci[, "97.5%"]

  return(out)
}

#' Pobuduvaty hrafiky butstrep-rozpodiliv dlia pidhonky PMM2
#'
#' @param object Rezultat z pmm2_inference
#' @param coefficients Yaki koefitsiienty pobuduvaty, za zamovchuvanniam usi
#'
#' @return Nevydymo povertaie informatsiiu pro histohramu
#' @export
plot_pmm2_bootstrap <- function(object, coefficients = NULL) {
  if(!inherits(object, "data.frame") ||
     !all(c("Estimate", "Std.Error", "conf.low", "conf.high") %in% names(object))) {
    stop("Ob'iekt maie buty rezultatom pmm2_inference()")
  }

  # Yakshcho koefitsiienty ne vkazani, vykorystovuvaty vsi
  if(is.null(coefficients)) {
    coefficients <- rownames(object)
  }

  # Filtruvaty do zapytanykh koefitsiientiv
  object_subset <- object[intersect(coefficients, rownames(object)), , drop = FALSE]

  # Perevirka na porozhnii nabir danykh
  if(nrow(object_subset) == 0) {
    warning("Zhoden iz zapytanykh koefitsiientiv ne znaideno v rezultatakh.")
    return(invisible(NULL))
  }

  # Nalashtuvaty komponuvannia hrafika
  n_coefs <- nrow(object_subset)
  n_cols <- min(2, n_coefs)
  n_rows <- ceiling(n_coefs / n_cols)

  # Zberehty stari nalashtuvannia par i vidnovyty pry vykhodi
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  par(mfrow = c(n_rows, n_cols))

  # Stvoryty hrafik shchilnosti dlia kozhnoho koefitsiienta
  result <- lapply(seq_len(n_coefs), function(i) {
    coef_name <- rownames(object_subset)[i]
    est <- object_subset[i, "Estimate"]
    ci_low <- object_subset[i, "conf.low"]
    ci_high <- object_subset[i, "conf.high"]
    se <- object_subset[i, "Std.Error"]

    # Perevirka na skinchenni znachennia
    if(!is.finite(est) || !is.finite(ci_low) || !is.finite(ci_high) || !is.finite(se)) {
      warning("Neskinchenni abo NA znachennia dlia koefitsiienta ", coef_name,
              ". Propuskaiemo tsei hrafik.")
      return(NULL)
    }

    # Stvoryty zaholovok hrafika
    main_title <- paste0(coef_name, "\nOtsinka: ", round(est, 4))

    # Otsinyty diapazon dlia osi x
    # Vykorystovuiemo bilsh nadiinyi pidkhid dlia vyznachennia diapazonu
    x_range <- range(c(est, ci_low, ci_high), na.rm = TRUE)
    # Rozshyryty diapazon na 20% v obokh napriamkakh
    x_range_width <- diff(x_range)
    x_range <- x_range + c(-0.2, 0.2) * x_range_width

    # Stvoryty tochky dlia osi x
    x_seq <- seq(x_range[1], x_range[2], length.out = 100)

    # Stvoryty znachennia shchilnosti dlia normalnoho rozpodilu
    y_seq <- dnorm(x_seq, mean = est, sd = se)

    # Pobuduvaty hrafik
    plot(x_seq, y_seq, type = "l",
         main = main_title,
         xlab = "Znachennia",
         ylab = "Shchilnist")

    # Dodaty vertykalni linii dlia otsinky ta CI
    abline(v = est, col = "red", lwd = 2)
    abline(v = ci_low, col = "blue", lty = 2)
    abline(v = ci_high, col = "blue", lty = 2)

    # Dodaty lehendu
    legend("topright",
           legend = c("Otsinka", "95% CI"),
           col = c("red", "blue"),
           lty = c(1, 2),
           lwd = c(2, 1),
           cex = 0.8)

    invisible(list(x = x_seq, y = y_seq, estimate = est,
                   ci_low = ci_low, ci_high = ci_high))
  })

  # Vydalyty NULL rezultaty
  result <- result[!sapply(result, is.null)]

  # Yakshcho vsi rezultaty NULL, povernuty NULL
  if(length(result) == 0) {
    warning("Ne vdalosia stvoryty zhodnoho hrafika.")
    return(invisible(NULL))
  }

  # Dodaty imena do rezultativ
  names(result) <- rownames(object_subset)[sapply(seq_len(n_coefs), function(i) {
    !is.null(result[[i]])
  })]

  invisible(result)
}


#' Butstrep-vysnovok dlia modelei chasovykh riadiv PMM2
#'
#' @param object ob'iekt klasu TS2fit
#' @param x (optsionalno) oryhinalnyi chasovyi riad; iakshcho NULL, vykorystovuie object@original_series
#' @param B kilkist butstrep-replikatsii
#' @param seed (optsionalno) dlia vidtvoriuvanosti
#' @param block_length dovzhyna bloku dlia blokovoho butstrepu; iakshcho NULL, vykorystovuie evrystychne znachennia
#' @param method typ butstrepu: "residual" abo "block"
#' @param parallel lohichne, chy vykorystovuvaty paralelni obchyslennia
#' @param cores kilkist iader dlia paralelnykh obchyslen
#' @param debug lohichne, chy vyvodyty dodatkovu diahnostychnu informatsiiu
#'
#' @return data.frame z stovptsiamy: Estimate, Std.Error, t.value, p.value
#' @export
ts_pmm2_inference <- function(object, x = NULL, B = 200, seed = NULL,
                              block_length = NULL, method = c("residual", "block"),
                              parallel = FALSE, cores = NULL, debug = FALSE) {
  # Pereviryty klas ob'iekta
  if (!inherits(object, "TS2fit")) {
    stop("Ob'iekt maie buty klasu 'TS2fit'")
  }

  # Vybraty metod butstrepu
  method <- match.arg(method)

  # Vstanovyty zerno dlia vidtvoriuvanosti, iakshcho nadano
  if (!is.null(seed)) set.seed(seed)

  # Vytiahnuty parametry modeli
  model_type <- object@model_type
  ar_order <- object@order$ar
  ma_order <- object@order$ma
  d <- object@order$d
  intercept <- object@intercept
  include_mean <- intercept != 0

  if(debug) {
    cat("Parametry modeli:\n")
    cat("model_type:", model_type, "\n")
    cat("ar_order:", ar_order, "\n")
    cat("ma_order:", ma_order, "\n")
    cat("d:", d, "\n")
    cat("intercept:", intercept, "\n")
    cat("include_mean:", include_mean, "\n")
  }

  # Yakshcho x ne nadanyi, vykorystovuvaty oryhinalnyi riad z ob'iekta
  if (is.null(x)) {
    x <- object@original_series
  }

  # Vytiahnuty koefitsiienty ta zalyshky
  coefs <- object@coefficients
  res <- object@residuals

  if(debug) {
    cat("Rozmir oryhinalnoho riadu:", length(x), "\n")
    cat("Rozmir vektora zalyshkiv:", length(res), "\n")
    cat("Kilkist koefitsiientiv:", length(coefs), "\n")
  }

  # Vykonaty blokovyi butstrep, iakshcho vkazano
  if (method == "block") {
    # Vyznachyty dovzhynu bloku, iakshcho ne nadano
    if (is.null(block_length)) {
      # Vykorystovuvaty evrystyku: kvadratnyi korin z dovzhyny riadu
      block_length <- ceiling(sqrt(length(x)))
    }

    # Pereviryty, chy dovzhyna bloku maie sens
    if (block_length < 2) {
      warning("Dovzhyna bloku zamala. Vstanovliuiu na 2.")
      block_length <- 2
    }
    if (block_length > length(x) / 4) {
      warning("Dovzhyna bloku zavelyka. Vstanovliuiu na 1/4 dovzhyny riadu.")
      block_length <- floor(length(x) / 4)
    }

    # Funktsiia dlia heneratsii butstrep-riadu z blokovoho butstrepu
    generate_block_bootstrap <- function(x, block_length) {
      n <- length(x)
      blocks_needed <- ceiling(n / block_length)

      # Mozhlyvi pochatkovi pozytsii blokiv
      start_positions <- 1:(n - block_length + 1)

      # Vybraty vypadkovi pochatkovi pozytsii
      selected_starts <- sample(start_positions, blocks_needed, replace = TRUE)

      # Stvoryty butstrep-riad
      boot_series <- numeric(0)
      for (start in selected_starts) {
        boot_series <- c(boot_series, x[start:(start + block_length - 1)])
      }

      # Obrizaty do oryhinalnoi dovzhyny
      boot_series[1:n]
    }

    # Lohika butstrepu vidrizniaietsia dlia blokovoho metodu
    boot_function <- function(b) {
      # Heneruvaty novyi riad za dopomohoiu blokovoho butstrepu
      x_b <- generate_block_bootstrap(x, block_length)

      # Vyznachaiemo pravylnyi format order v zalezhnosti vid typu modeli
      if(model_type == "ar") {
        boot_order <- ar_order  # Dlia AR modelei - odne chyslo
      } else if(model_type == "ma") {
        boot_order <- ma_order  # Dlia MA modelei - odne chyslo
      } else if(model_type == "arma") {
        boot_order <- c(ar_order, ma_order)  # Dlia ARMA - vektor dovzhyny 2
      } else if(model_type == "arima") {
        boot_order <- c(ar_order, d, ma_order)  # Dlia ARIMA - vektor dovzhyny 3
      } else {
        stop("Nevidomyi typ modeli: ", model_type)
      }

      # Pidihnaty model na butstrep-riadi
      fit_b <- tryCatch({
        ts_pmm2(x_b, order = boot_order,
                model_type = model_type,
                include.mean = include_mean)
      }, error = function(e) {
        warning("Butstrep-replikatsiia ", b, " ne vdalasia: ", e$message)
        return(NULL)
      })

      if (!is.null(fit_b)) {
        return(fit_b@coefficients)
      } else {
        return(rep(NA, length(coefs)))
      }
    }
  } else {
    # Dlia butstrepu zalyshkiv
    # Funktsiia dlia heneratsii novoho riadu na osnovi modeli ta butstrep-zalyshkiv
    generate_with_residuals <- function(model, residuals) {
      n <- length(model@original_series)

      # Dlia ARIMA modelei, potribno spochatku dyferentsiiuvaty
      if (model_type == "arima" && d > 0) {
        # Tut potribna skladnisha lohika dlia vidnovlennia oryhinalnoho riadu
        # pislia heneratsii dyferentsiiovanoho riadu
        # Tse sproshchenyi pidkhid:
        diff_x <- numeric(n - d)

        # Butstrep zalyshkiv
        res_b <- sample(residuals[!is.na(residuals)], size = n - d, replace = TRUE)

        # Heneruvaty dyferentsiiovanyi riad
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

        # Sproshchena symuliatsiia ARIMA protsesu
        diff_x <- arima.sim(model = list(
          ar = if(ar_order > 0) ar_coefs else NULL,
          ma = if(ma_order > 0) ma_coefs else NULL
        ), n = n - d, innov = res_b, n.start = max(ar_order, ma_order))

        # Intehruvaty nazad
        x_b <- diffinv(diff_x, differences = d)

        # Dodaty perekhoplennia, iakshcho potribno
        if (include_mean) {
          x_b <- x_b + intercept
        }
      } else {
        # Dlia AR, MA, ARMA modelei
        res_b <- sample(residuals[!is.na(residuals)], size = n, replace = TRUE)

        # Vytiahnuty koefitsiienty AR i MA
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

        # Symuliuvaty protses ARMA
        x_b <- arima.sim(model = list(
          ar = if(ar_order > 0) ar_coefs else NULL,
          ma = if(ma_order > 0) ma_coefs else NULL
        ), n = n, innov = res_b, n.start = max(ar_order, ma_order))

        # Dodaty perekhoplennia, iakshcho potribno
        if (include_mean) {
          x_b <- x_b + intercept
        }
      }

      return(x_b)
    }

    boot_function <- function(b) {
      # Heneruvaty novyi riad
      x_b <- generate_with_residuals(object, res)

      # Vyznachaiemo pravylnyi format order v zalezhnosti vid typu modeli
      if(model_type == "ar") {
        boot_order <- ar_order  # Dlia AR modelei - odne chyslo
      } else if(model_type == "ma") {
        boot_order <- ma_order  # Dlia MA modelei - odne chyslo
      } else if(model_type == "arma") {
        boot_order <- c(ar_order, ma_order)  # Dlia ARMA - vektor dovzhyny 2
      } else if(model_type == "arima") {
        boot_order <- c(ar_order, d, ma_order)  # Dlia ARIMA - vektor dovzhyny 3
      } else {
        stop("Nevidomyi typ modeli: ", model_type)
      }

      if(debug && b == 1) {
        cat("Butstrep replikatsiia 1:\n")
        cat("Typ modeli:", model_type, "\n")
        cat("Poriadok dlia butstrepu:", paste(boot_order, collapse=", "), "\n")
      }

      # Pidihnaty model na butstrep-riadi
      fit_b <- tryCatch({
        ts_pmm2(x_b, order = boot_order,
                model_type = model_type,
                include.mean = include_mean)
      }, error = function(e) {
        warning("Butstrep-replikatsiia ", b, " ne vdalasia: ", e$message)
        return(NULL)
      })

      if (!is.null(fit_b)) {
        return(fit_b@coefficients)
      } else {
        return(rep(NA, length(coefs)))
      }
    }
  }

  # Vykonaty butstrep: paralelno abo poslidovno
  use_parallel <- parallel && requireNamespace("parallel", quietly = TRUE)

  if (use_parallel) {
    if (is.null(cores)) {
      cores <- max(1, parallel::detectCores() - 1)
    }

    boot_results <- parallel::mclapply(seq_len(B), function(b) {
      boot_function(b)
    }, mc.cores = cores)

    # Peretvoryty spysok na matrytsiu
    boot_est <- do.call(rbind, boot_results)
  } else {
    # Poslidovni obchyslennia
    boot_est <- matrix(0, nrow = B, ncol = length(coefs))
    colnames(boot_est) <- names(coefs)

    # Vidstezhennia prohresu
    pb <- NULL
    if (interactive() && B > 10) {
      if (requireNamespace("utils", quietly = TRUE)) {
        pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
      }
    }

    for (b in seq_len(B)) {
      boot_est[b, ] <- boot_function(b)

      # Onovyty indykator prohresu
      if (!is.null(pb)) utils::setTxtProgressBar(pb, b)
    }

    # Zakryty indykator prohresu
    if (!is.null(pb)) close(pb)
  }

  # Vydalyty riadky zi znachenniamy NA
  na_rows <- apply(boot_est, 1, function(row) any(is.na(row)))
  if (any(na_rows)) {
    warning("Vydaleno ", sum(na_rows), " butstrep-replikatsii cherez pomylky otsiniuvannia")
    boot_est <- boot_est[!na_rows, , drop = FALSE]
  }

  # Pereviryty, chy maiemo dostatno uspishnykh butstrepiv
  if (nrow(boot_est) < 10) {
    stop("Zamalo uspishnykh butstrep-replikatsii dlia obchyslennia nadiinoho vysnovku")
  }

  # Pereviryty, chy ie NaN abo Inf znachennia
  if (any(is.nan(boot_est)) || any(is.infinite(boot_est))) {
    warning("Vyiavleno NaN abo neskinchenni znachennia v butstrep-replikatsiiakh. Zaminiuiemo ikh na NA.")
    boot_est[is.nan(boot_est) | is.infinite(boot_est)] <- NA
  }

  # Obchyslyty kovariatsiinu matrytsiu ta standartni pomylky
  cov_mat <- cov(boot_est, use = "pairwise.complete.obs")
  est <- coefs
  se <- sqrt(diag(cov_mat))

  # Obchyslyty t-znachennia ta p-znachennia
  t_val <- est / se
  p_val <- 2 * (1 - pnorm(abs(t_val)))

  # Stvoryty vykhidnyi freim danykh
  out <- data.frame(
    Estimate = est,
    Std.Error = se,
    t.value = t_val,
    p.value = p_val
  )

  # Dodaty imena AR i MA parametriv
  param_names <- c()
  if (ar_order > 0) {
    param_names <- c(param_names, paste0("ar", 1:ar_order))
  }
  if (ma_order > 0) {
    param_names <- c(param_names, paste0("ma", 1:ma_order))
  }
  rownames(out) <- param_names

  # Obchyslyty dovirchi intervaly
  ci <- t(apply(boot_est, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE))
  colnames(ci) <- c("2.5%", "97.5%")

  # Dodaty dovirchi intervaly do vykhodu
  out$conf.low <- ci[, "2.5%"]
  out$conf.high <- ci[, "97.5%"]

  return(out)
}
