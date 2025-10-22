# pmm2_ts_methods.R - Metody dlia roboty z ob'iektamy modelei chasovykh riadiv

#' Vytiahnuty koefitsiienty z ob'iekta TS2fit
#'
#' @param object Ob'iekt TS2fit
#' @param ... Dodatkovi arhumenty (ne vykorystovuiutsia)
#'
#' @return Imenovanyi vektor koefitsiientiv
#' @export
setMethod("coef", "TS2fit",
          function(object, ...) {
            # Otrymaty parametry modeli
            model_type <- object@model_type
            ar_order <- object@order$ar
            ma_order <- object@order$ma

            # Vytiahnuty ta imenuvaty AR koefitsiienty
            if(ar_order > 0) {
              ar_coefs <- object@coefficients[1:ar_order]
              names(ar_coefs) <- paste0("ar", 1:ar_order)
            } else {
              ar_coefs <- numeric(0)
            }

            # Vytiahnuty ta imenuvaty MA koefitsiienty
            if(ma_order > 0) {
              ma_coefs <- object@coefficients[(ar_order+1):(ar_order+ma_order)]
              names(ma_coefs) <- paste0("ma", 1:ma_order)
            } else {
              ma_coefs <- numeric(0)
            }

            # Ob'iednaty koefitsiienty
            result <- c(ar_coefs, ma_coefs)

            # Dodaty perekhoplennia, iakshcho prysutnie
            if(object@intercept != 0) {
              result <- c(intercept = object@intercept, result)
            }

            return(result)
          })

#' Vytiahnuty zalyshky z ob'iekta TS2fit
#'
#' @param object Ob'iekt TS2fit
#' @param ... Dodatkovi arhumenty (ne vykorystovuiutsia)
#'
#' @return Vektor zalyshkiv (innovatsii)
#' @export
setMethod("residuals", "TS2fit",
          function(object, ...) {
            object@residuals
          })

#' Otrymaty pidihnani znachennia dlia AR modeli
#'
#' @param object Ob'iekt TS2fit z model_type="ar"
#' @return Vektor pidihnanykh znachen
#' @keywords internal
get_ar_fitted <- function(object) {
  if(object@model_type != "ar") {
    stop("Tsia funktsiia lyshe dlia AR modelei")
  }

  x <- object@original_series
  ar_order <- object@order$ar
  ar_coef <- object@coefficients[1:ar_order]
  intercept <- object@intercept

  if(intercept != 0) {
    x_centered <- x - intercept
  } else {
    x_centered <- x
  }

  # Stvoryty matrytsiu dyzainu ta obchyslyty pidihnani znachennia
  X <- create_ar_matrix(x_centered, ar_order)
  fitted <- as.vector(X %*% ar_coef) + intercept
  return(fitted)
}

#' Vytiahnuty pidihnani znachennia z ob'iekta TS2fit
#'
#' @param object Ob'iekt TS2fit
#' @param ... Dodatkovi arhumenty (ne vykorystovuiutsia)
#'
#' @return Vektor pidihnanykh znachen
#' @export
setMethod("fitted", "TS2fit",
          function(object, ...) {
            # Otrymaty typ modeli
            model_type <- object@model_type

            # Obchyslyty pidihnani znachennia v zalezhnosti vid typu modeli
            if(model_type == "ar") {
              # Dlia AR modelei vykorystovuvaty priame obchyslennia
              fitted_values <- get_ar_fitted(object)
            } else {
              # Dlia inshykh modelei: pidihnani = oryhinalni minus zalyshky
              orig <- object@original_series
              resid <- object@residuals

              # Vyrivniaty dovzhyny (chasto zalyshky korotshi cherez pochatkovi znachennia)
              len_diff <- length(orig) - length(resid)
              if(len_diff > 0 && !all(is.na(resid))) {
                # Znaity pershe ne-NA znachennia v resid
                first_valid <- min(which(!is.na(resid)))

                # Pobuduvaty vektor fitted z NA v pochatku
                fitted_values <- rep(NA, length(orig))
                valid_indices <- first_valid:length(resid)

                # Vstanovyty diisni znachennia
                fitted_values[(len_diff + valid_indices)] <-
                  orig[(len_diff + valid_indices)] - resid[valid_indices]
              } else {
                fitted_values <- orig - resid
              }
            }

            return(fitted_values)
          })

#' Pobuduvaty diahnostychni hrafiky dlia ob'iektiv TS2fit
#'
#' @param x Ob'iekt TS2fit
#' @param y Ne vykorystovuietsia (dlia sumisnosti metodu S4)
#' @param which Tsilochyselnyi vektor, shcho vkazuie, iaki hrafiky vyrobliaty
#' @param ... dodatkovi arhumenty, peredani funktsiiam pobudovy hrafikiv
#'
#' @return Nevydymo povertaie x
#'
#' @export
setMethod("plot", signature(x = "TS2fit", y = "missing"),
          function(x, y, which = c(1:4), ...) {
            op <- par(no.readonly = TRUE)
            on.exit(par(op))

            # Otrymaty parametry modeli
            model_type <- x@model_type
            ar_order <- x@order$ar
            ma_order <- x@order$ma
            d <- x@order$d

            # Maket hrafika za zamovchuvanniam
            par(mfrow = c(2, 2))

            # Dlia modelei ARIMA my mozhemo zakhotity pobuduvaty oryhinalnyi/dyferentsiiovanyi riad takozh
            if(model_type == "arima" && length(which) > 4) {
              par(mfrow = c(3, 2))
            }

            # Obchyslyty pidihnani znachennia ta zalyshky
            residuals <- as.numeric(x@residuals)
            fitted <- fitted(x)

            # Vyznachyty, iaki hrafiky vidobrazhaty
            plot_idx <- 1
            n_plots <- min(length(which), 6) # Maksymum 6 hrafikiv

            # Dlia modelei ARIMA, my mozhemo khotity inshi hrafiky
            if(model_type == "arima") {
              # Hrafik 1: Oryhinalnyi chasovyi riad (tilky ARIMA)
              if(1 %in% which && plot_idx <= n_plots) {
                plot(x@original_series, type = "l",
                     main = "Oryhinalnyi chasovyi riad",
                     xlab = "Chas",
                     ylab = "Znachennia",
                     ...)
                plot_idx <- plot_idx + 1
              }

              # Hrafik 2: Dyferentsiiovanyi chasovyi riad (tilky ARIMA)
              if(2 %in% which && plot_idx <= n_plots && d > 0) {
                diff_series <- diff(x@original_series, differences = d)
                plot(diff_series, type = "l",
                     main = paste0("Dyferentsiiovanyi riad (d=", d, ")"),
                     xlab = "Chas",
                     ylab = "Znachennia",
                     ...)
                plot_idx <- plot_idx + 1
              }
            }

            # Standartni hrafiky dlia vsikh typiv modelei
            # Hrafik: Zalyshky vs Pidihnani
            if(3 %in% which && plot_idx <= n_plots) {
              plot(fitted, residuals,
                   main = "Zalyshky vs Pidihnani",
                   xlab = "Pidihnani znachennia",
                   ylab = "Zalyshky",
                   ...)
              abline(h = 0, lty = 2)
              lines(lowess(fitted, residuals), col = "red")
              plot_idx <- plot_idx + 1
            }

            # Hrafik: Normalnyi Q-Q hrafik
            if(4 %in% which && plot_idx <= n_plots) {
              qqnorm(residuals, main = "Normalnyi Q-Q hrafik", ...)
              qqline(residuals)
              plot_idx <- plot_idx + 1
            }

            # Hrafik: ACF zalyshkiv
            if(5 %in% which && plot_idx <= n_plots) {
              acf(residuals, main = "ACF zalyshkiv", ...)
              plot_idx <- plot_idx + 1
            }

            # Hrafik: Histohrama zalyshkiv
            if(6 %in% which && plot_idx <= n_plots) {
              hist(residuals,
                   main = "Histohrama zalyshkiv",
                   xlab = "Zalyshky",
                   breaks = "FD",
                   ...)
              plot_idx <- plot_idx + 1
            }

            invisible(x)
          })

#' Metod prohnozuvannia dlia ob'iektiv TS2fit
#'
#' @param object Ob'iekt TS2fit
#' @param n.ahead Kilkist krokiv vpered dlia prohnozuvannia
#' @param ... dodatkovi arhumenty (ne vykorystovuiutsia)
#'
#' @return Vektor abo spysok prohnoziv, zalezhno vid typu modeli
#'
#' @export
setMethod("predict", "TS2fit",
          function(object, n.ahead = 1, ...) {
            # Otrymaty parametry modeli
            model_type <- object@model_type
            ar_order <- object@order$ar
            ma_order <- object@order$ma
            d <- object@order$d
            intercept <- object@intercept

            # Vytiahnuty koefitsiienty
            if(ar_order > 0) {
              ar_coef <- object@coefficients[1:ar_order]
            } else {
              ar_coef <- numeric(0)
            }

            if(ma_order > 0) {
              ma_coef <- object@coefficients[(ar_order+1):(ar_order+ma_order)]
            } else {
              ma_coef <- numeric(0)
            }

            # Dlia AR modelei, realizuvaty priame prohnozuvannia
            if(model_type == "ar") {
              x <- object@original_series
              n <- length(x)
              pred <- numeric(n.ahead)

              # Heneruvaty prohnozy
              for(i in 1:n.ahead) {
                # Vykorystovuvaty oryhinalni dani ta poperedni prohnozy za potreby
                lags <- numeric(ar_order)
                for(j in 1:ar_order) {
                  if(i - j <= 0) {
                    # Vykorystovuvaty oryhinalni dani
                    lags[j] <- x[n - j + i]
                  } else {
                    # Vykorystovuvaty poperedni prohnozy
                    lags[j] <- pred[i - j]
                  }
                }

                # Obchyslyty prohnoz
                pred[i] <- sum(ar_coef * lags) + intercept
              }

              return(pred)

            } else if(model_type == "ma") {
              # Dlia MA modelei, prohnozy za mezhamy poriadku - tse prosto serednie
              innovations <- object@residuals

              # Heneruvaty prohnozy MA
              ma_pred <- function(innovations, ma_coef, n.ahead) {
                n <- length(innovations)
                q <- length(ma_coef)
                pred <- numeric(n.ahead)

                for(i in 1:n.ahead) {
                  for(j in 1:min(i, q)) {
                    if((n - i + j) > 0) {
                      pred[i] <- pred[i] + ma_coef[j] * innovations[n - i + j]
                    }
                  }
                }
                return(pred)
              }

              if(n.ahead > ma_order) {
                return(c(ma_pred(innovations, ma_coef, ma_order),
                         rep(intercept, n.ahead - ma_order)))
              } else {
                return(ma_pred(innovations, ma_coef, n.ahead))
              }

            } else {
              # Dlia ARMA ta ARIMA modelei, vykorystovuvaty prohnozy stats::arima
              # iaki pravylno obrobliaiut obydva komponenty

              # Nalashtuvaty model arima z fiksovanymy parametramy
              arima_order <- c(ar_order, ifelse(model_type == "arima", d, 0), ma_order)

              # Vykorystovuvaty funktsiiu prohnozuvannia z paketu stats
              arima_pred <- stats::predict(
                stats::arima(object@original_series,
                             order = arima_order,
                             include.mean = (intercept != 0),
                             fixed = c(ar_coef, ma_coef, if(intercept != 0) intercept else NULL)),
                n.ahead = n.ahead
              )

              return(arima_pred)
            }
          })

#' Porivniaty PMM2 z klasychnymy metodamy otsiniuvannia chasovykh riadiv
#'
#' @param x Chyslovyi vektor danykh chasovoho riadu
#' @param order Spetsyfikatsiia poriadku modeli (dyv. ts_pmm2 dlia formatu)
#' @param model_type Typ modeli: "ar", "ma", "arma", abo "arima"
#' @param include.mean Lohichne, chy vkliuchaty chlen perekhoplennia
#' @param pmm2_args Spysok dodatkovykh arhumentiv dlia peredachi v ts_pmm2()
#'
#' @return Spysok z pidihnanymy modeliamy ta tablytsiamy porivniannia
#' @export
compare_ts_methods <- function(x, order, model_type = c("ar", "ma", "arma", "arima"),
                               include.mean = TRUE, pmm2_args = list()) {
  # Vybraty arhument model_type
  model_type <- match.arg(model_type)

  # Pidhotuvaty porivniannia modelei na osnovi model_type
  if(model_type == "ar") {
    # Dlia AR modelei
    # Pidihnaty AR model za dopomohoiu metodu Yula-Volkera
    yw_fit <- stats::ar(x, order.max = order, aic = FALSE, method = "yw",
                        demean = include.mean)

    # Pidihnaty AR model za dopomohoiu metodu OLS
    ols_fit <- stats::ar(x, order.max = order, aic = FALSE, method = "ols",
                         demean = include.mean)

    # Pidihnaty AR model za dopomohoiu metodu MLE
    mle_fit <- stats::ar(x, order.max = order, aic = FALSE, method = "mle",
                         demean = include.mean)

    # Pidihnaty AR model za dopomohoiu PMM2
    pmm2_args <- c(list(x = x, order = order, model_type = "ar",
                        include.mean = include.mean), pmm2_args)
    pmm2_fit <- do.call(ts_pmm2, pmm2_args)

    # Vytiahnuty koefitsiienty
    coef_yw <- yw_fit$ar
    coef_ols <- ols_fit$ar
    coef_mle <- mle_fit$ar
    coef_pmm2 <- pmm2_fit@coefficients

    # Obchyslyty zalyshky
    res_yw <- yw_fit$resid[!is.na(yw_fit$resid)]
    res_ols <- ols_fit$resid[!is.na(ols_fit$resid)]
    res_mle <- mle_fit$resid[!is.na(mle_fit$resid)]
    res_pmm2 <- pmm2_fit@residuals

    methods <- c("YW", "OLS", "MLE", "PMM2")

    result_list <- list(
      yw = yw_fit,
      ols = ols_fit,
      mle = mle_fit,
      pmm2 = pmm2_fit
    )

  } else if(model_type %in% c("ma", "arma", "arima")) {
    # Dlia MA, ARMA ta ARIMA modelei

    # Pidhotuvaty poriadok arima na osnovi typu modeli
    if(model_type == "ma") {
      arima_order <- c(0, 0, order)
    } else if(model_type == "arma") {
      arima_order <- c(order[1], 0, order[2])
    } else {
      arima_order <- order
    }

    # Pidihnaty model za dopomohoiu metodu CSS
    css_fit <- arima(x, order = arima_order, method = "CSS", include.mean = include.mean)

    # Pidihnaty model za dopomohoiu metodu ML
    ml_fit <- arima(x, order = arima_order, method = "ML", include.mean = include.mean)

    # Pidihnaty model za dopomohoiu PMM2
    pmm2_args <- c(list(x = x, order = order, model_type = model_type,
                        include.mean = include.mean), pmm2_args)
    pmm2_fit <- do.call(ts_pmm2, pmm2_args)

    # Vytiahnuty nazvy koefitsiientiv AR ta MA na osnovi typu modeli
    if(model_type == "ma") {
      ar_names <- character(0)
      ma_names <- paste0("ma", 1:order)
    } else if(model_type == "arma") {
      ar_names <- paste0("ar", 1:order[1])
      ma_names <- paste0("ma", 1:order[2])
    } else {
      ar_names <- if(order[1] > 0) paste0("ar", 1:order[1]) else character(0)
      ma_names <- if(order[3] > 0) paste0("ma", 1:order[3]) else character(0)
    }

    coef_names <- c(ar_names, ma_names)

    # Vytiahnuty koefitsiienty
    coef_css <- as.numeric(css_fit$coef[coef_names])
    coef_ml <- as.numeric(ml_fit$coef[coef_names])
    coef_pmm2 <- pmm2_fit@coefficients

    # Obchyslyty zalyshky
    res_css <- residuals(css_fit)
    res_ml <- residuals(ml_fit)
    res_pmm2 <- pmm2_fit@residuals

    methods <- c("CSS", "ML", "PMM2")

    result_list <- list(
      css = css_fit,
      ml = ml_fit,
      pmm2 = pmm2_fit
    )
  }

  # Obchyslyty statystyku zalyshkiv dlia vsikh metodiv
  residuals_list <- if(model_type == "ar") {
    list(res_yw, res_ols, res_mle, res_pmm2)
  } else {
    list(res_css, res_ml, res_pmm2)
  }

  compute_res_stats <- function(res) {
    m2 <- mean(res^2, na.rm = TRUE)
    m3 <- mean(res^3, na.rm = TRUE)
    m4 <- mean(res^4, na.rm = TRUE)

    c(RSS = sum(res^2, na.rm = TRUE),
      MAE = mean(abs(res), na.rm = TRUE),
      Skewness = m3 / m2^(3/2),
      Kurtosis = m4 / m2^2)
  }

  res_stats <- data.frame(
    Method = methods,
    do.call(rbind, lapply(residuals_list, compute_res_stats))
  )

  # Stvoryty tablytsiu porivniannia koefitsiientiv
  if(model_type == "ar") {
    coef_names <- paste0("ar", 1:order)
    coef_values <- list(coef_yw, coef_ols, coef_mle, coef_pmm2)
  } else if(model_type == "ma") {
    coef_names <- paste0("ma", 1:order)
    coef_values <- list(coef_css, coef_ml, coef_pmm2)
  } else {
    coef_values <- list(coef_css, coef_ml, coef_pmm2)
  }

  coef_table <- data.frame(
    Coefficient = coef_names,
    do.call(cbind, lapply(seq_along(methods), function(i) {
      result <- coef_values[[i]]
      names(result) <- methods[i]
      return(result)
    }))
  )

  # Povernuty rezultaty
  result_list$coefficients <- coef_table
  result_list$residual_stats <- res_stats

  return(result_list)
}

#' Porivniaty metody AR
#'
#' @inheritParams compare_ts_methods
#' @export
compare_ar_methods <- function(x, order = 1, include.mean = TRUE, pmm2_args = list()) {
  compare_ts_methods(x, order = order, model_type = "ar",
                     include.mean = include.mean, pmm2_args = pmm2_args)
}

#' Porivniaty metody MA
#'
#' @inheritParams compare_ts_methods
#' @export
compare_ma_methods <- function(x, order = 1, include.mean = TRUE, pmm2_args = list()) {
  compare_ts_methods(x, order = order, model_type = "ma",
                     include.mean = include.mean, pmm2_args = pmm2_args)
}

#' Porivniaty metody ARMA
#'
#' @inheritParams compare_ts_methods
#' @export
compare_arma_methods <- function(x, order = c(1, 1), include.mean = TRUE, pmm2_args = list()) {
  compare_ts_methods(x, order = order, model_type = "arma",
                     include.mean = include.mean, pmm2_args = pmm2_args)
}

#' Porivniaty metody ARIMA
#'
#' @inheritParams compare_ts_methods
#' @export
compare_arima_methods <- function(x, order = c(1, 1, 1), include.mean = TRUE, pmm2_args = list()) {
  compare_ts_methods(x, order = order, model_type = "arima",
                     include.mean = include.mean, pmm2_args = pmm2_args)
}
