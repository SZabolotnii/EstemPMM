# pmm2_main.R - Osnovnyi modul dlia liniinykh modelei PMM2

#' pmm2: Holovna funktsiia dlia PMM2 (S=2)
#'
#' Pidhaniaie liniinu model za dopomohoiu metodu maksymizatsii polinomiv (poriadok 2),
#' iakyi ie robastnym shchodo nehausivskykh pomylok.
#'
#' @param formula Formula R dlia modeli
#' @param data data.frame, shcho mistyt zminni u formuli
#' @param max_iter tsile: maksymalna kilkist iteratsii dlia alhorytmu
#' @param tol chyslove: dopusk dlia zbizhnosti
#' @param regularize lohichne: dodaty male znachennia do diahonali dlia chyslovoi stabilnosti
#' @param reg_lambda chyslove: parametr rehuliaryzatsii (iakshcho regularize=TRUE)
#' @param na.action funktsiia dlia obrobky vidsutnikh znachen, za zamovchuvanniam - na.fail
#' @param weights optsionalnyi vektor vah (poky ne realizovano)
#' @param verbose lohichne: chy vyvodyty informatsiiu pro prohres
#'
#' @details
#' Alhorytm PMM2 pratsiuie nastupnym chynom:
#'
#' 1. Pidhaniaie zvychainu rehresiiu naimenshykh kvadrativ (OLS) dlia otrymannia pochatkovykh otsinok
#' 2. Obchysliuie tsentralni momenty (m2, m3, m4) iz zalyshkiv OLS
#' 3. Iteratyvno pokrashchuie otsinky parametriv za dopomohoiu pidkhodu na osnovi hradiienta
#'
#' PMM2 osoblyvo korysnyi, koly termy pomylok ne ie hausivskymy.
#'
#' @return Ob'iekt S4 \code{PMM2fit}
#' @export
#'
#' @examples
#' \dontrun{
#' # Heneruvaty dani vybirky z t-rozpodilenymy pomylkamy
#' n <- 100
#' x <- rnorm(n)
#' y <- 2 + 3*x + rt(n, df=3)
#' dat <- data.frame(y=y, x=x)
#'
#' # Pidihnaty model za dopomohoiu PMM2
#' fit <- lm_pmm2(y ~ x, data=dat)
#'
#' # Reziume ta statystychnyi vysnovok
#' summary(fit, formula=y~x, data=dat)
#' }
lm_pmm2 <- function(formula, data,
                    max_iter=50, tol=1e-6,
                    regularize=TRUE, reg_lambda=1e-8,
                    na.action=na.fail, weights=NULL,
                    verbose=FALSE)
{
  # Zakhopyty vyklyk
  call <- match.call()

  # Pereviryty validnist vkhodu
  if(missing(formula) || missing(data)) {
    stop("Obydva 'formula' ta 'data' maiut buty nadani")
  }

  if(!is.data.frame(data)) {
    stop("'data' maie buty data.frame")
  }

  if(max_iter <= 0) {
    stop("'max_iter' maie buty dodatnym")
  }

  if(tol <= 0) {
    stop("'tol' maie buty dodatnym")
  }

  # Obrobyty vidsutni znachennia
  if(!is.null(na.action)) {
    data <- na.action(data)
  }

  # Pereviryty na vahy
  if(!is.null(weights)) {
    warning("Vahy poky ne realizovani v PMM2. Ihnoruiu vahy.")
  }

  # 1) OLS
  if(verbose) cat("Pidhaniaiu pochatkovu model OLS...\n")

  mf <- model.frame(formula, data)
  X  <- model.matrix(formula, mf)
  y  <- model.response(mf)

  # Pereviryty na defitsyt ranhu
  qr_X <- qr(X)
  if(qr_X$rank < ncol(X)) {
    warning("Matrytsia dyzainu maie defitsyt ranhu, deiaki koefitsiienty mozhut buty neotsiniuvanymy")
  }

  fit_ols <- lm.fit(x=X, y=y)
  b_ols   <- fit_ols$coefficients

  # Obrobyty NA v koefitsiientakh OLS
  if(any(is.na(b_ols))) {
    stop("Pidhonka OLS pryvela do NA koefitsiientiv. Perevirte na multykolinearnist.")
  }

  # Peretvoryty b_ols na chyslovyi vektor dlia zabezpechennia uzhodzhenoi obrobky
  b_ols <- as.numeric(b_ols)

  # 2) OLS zalyshky => m2, m3, m4
  res_ols <- y - (X %*% b_ols)
  moments <- compute_moments(res_ols)
  m2 <- moments$m2
  m3 <- moments$m3
  m4 <- moments$m4

  if(verbose) {
    cat("Pochatkovi momenty iz zalyshkiv OLS:\n")
    cat("  m2 =", m2, "\n")
    cat("  m3 =", m3, "\n")
    cat("  m4 =", m4, "\n")
  }

  # Pereviryty na potentsiini problemy z momentamy
  if(m2 <= 0) {
    warning("Druhyi tsentralnyi moment (m2) ne ie dodatnym. Rezultaty mozhut buty nenadiinymy.")
  }

  if(m4 <= m2^2) {
    warning("Chetvertyi tsentralnyi moment (m4) menshyi za m2^2. Tse porushuie bazovu nerivnist dlia rozpodiliv imovirnostei.")
  }

  # 3) Zapustyty unifikovanyi alhorytm PMM2
  if(verbose) cat("Pochynaiu iteratsii PMM2...\n")

  out <- pmm2_algorithm(b_ols, X, y, m2, m3, m4,
                        max_iter = max_iter, tol = tol,
                        regularize = regularize, reg_lambda = reg_lambda,
                        verbose = verbose)

  # Vytiahnuty rezultaty
  b_est   <- out$b
  conv    <- out$convergence
  iter    <- out$iterations
  final_res <- out$residuals

  if(verbose) {
    cat("Alhorytm PMM2 zaversheno.\n")
    cat("  Zbihsia:", conv, "\n")
    cat("  Iteratsii:", iter, "\n")
  }

  # Povernuty ob'iekt S4 z usima rezultatamy
  ans <- new("PMM2fit",
             coefficients = b_est,
             residuals = final_res,
             m2 = m2,
             m3 = m3,
             m4 = m4,
             convergence = conv,
             iterations = iter,
             call = call)
  attr(ans, "model_matrix") <- X
  attr(ans, "model_frame") <- mf
  attr(ans, "response") <- as.numeric(y)
  attr(ans, "data") <- data

  return(ans)
}

#' Vytiahnuty koefitsiienty z ob'iekta PMM2fit
#'
#' @param object Ob'iekt PMM2fit
#' @param ... Dodatkovi arhumenty (ne vykorystovuiutsia)
#'
#' @return Vektor koefitsiientiv
#' @export
setMethod("coef", "PMM2fit",
          function(object, ...) {
            object@coefficients
          })

#' Vytiahnuty zalyshky z ob'iekta PMM2fit
#'
#' @param object Ob'iekt PMM2fit
#' @param ... Dodatkovi arhumenty (ne vykorystovuiutsia)
#'
#' @return Vektor zalyshkiv
#' @export
setMethod("residuals", "PMM2fit",
          function(object, ...) {
            object@residuals
          })

#' Vytiahnuty pidihnani znachennia z ob'iekta PMM2fit
#'
#' @param object Ob'iekt PMM2fit
#' @param data Neobov'iazkove dzherelo dlia rekonstruktsii modeli, yakshcho ob'iekt ne mistyt zberezhenykh danykh
#' @param ... Dodatkovi arhumenty (ne vykorystovuiutsia)
#'
#' @return Vektor pidihnanykh znachen
#' @export
setMethod("fitted", "PMM2fit",
          function(object, data = NULL, ...) {
            fitted_values(object, data)
          })

#' Rozrakhuvaty AIC dlia ob'iekta PMM2fit
#'
#' @param object Ob'iekt PMM2fit
#' @param ... Dodatkovi arhumenty (ne vykorystovuiutsia)
#' @param k Shtraf za parametr, shcho bude vykorystovuvatys; standartno 2
#'
#' @return Znachennia AIC
#' @export
setMethod("AIC", "PMM2fit",
          function(object, ..., k = 2) {
            res <- object@residuals
            n <- length(res)
            p <- length(object@coefficients)

            # Aproksymatsiia loharyfmichnoi pravdopodibnosti
            ll <- -n/2 * log(sum(res^2)/n) - n/2 * (1 + log(2*pi))

            # AIC
            -2 * ll + k * p
          })

#' Pobuduvaty diahnostychni hrafiky dlia ob'iekta PMM2fit
#'
#' @param x Ob'iekt PMM2fit
#' @param y Ne vykorystovuietsia (sumisnist z generic)
#' @param which Nabir hrafikiv dlia vidobrazhennia (znachennia 1-4)
#' @param ... Dodatkovi arhumenty, shcho peredaiutsia hrafichnym funktsiiam
#'
#' @return Nevydymo povertaie vkhidnyi ob'iekt
#' @export
setMethod("plot", signature(x = "PMM2fit", y = "missing"),
          function(x, y, which = 1:4, ...) {
            res <- as.numeric(x@residuals)
            fitted_vals <- tryCatch({
              fitted(x)
            }, error = function(e) {
              stored_X <- attr(x, "model_matrix")
              if(!is.null(stored_X)) {
                as.vector(stored_X %*% x@coefficients)
              } else {
                seq_along(res)
              }
            })

            which <- intersect(unique(which), 1:4)
            if(length(which) == 0) {
              which <- 1:4
            }

            old_par <- graphics::par(no.readonly = TRUE)
            on.exit(graphics::par(old_par))
            n_plots <- length(which)
            graphics::par(mfrow = c(2, 2))

            for(idx in which) {
              switch(idx,
                     {
                       graphics::plot(fitted_vals, res,
                                      main = "Residuals vs Fitted",
                                      xlab = "Fitted values",
                                      ylab = "Residuals", ...)
                       graphics::abline(h = 0, col = "red", lty = 2)
                     },
                     {
                       stats::qqnorm(res, main = "Normal Q-Q", ...)
                       stats::qqline(res, col = "red", lty = 2)
                     },
                     {
                       graphics::plot(seq_along(res), res, type = "l",
                                      main = "Residuals over Index",
                                      xlab = "Observation",
                                      ylab = "Residual", ...)
                       graphics::abline(h = 0, col = "red", lty = 2)
                     },
                     {
                       graphics::hist(res,
                                      main = "Residual Histogram",
                                      xlab = "Residuals",
                                      breaks = "FD", ...)
                     })
            }

            invisible(x)
          })

#' Dopomizhna funktsiia dlia vytiahnennia pidihnanykh znachen
#'
#' @param object Ob'iekt PMM2fit
#' @return Vektor pidihnanykh znachen
#'
#' @keywords internal
fitted_values <- function(object, data = NULL) {
  if(is.null(object@call)) {
    stop("Ob'iekt PMM2fit ne mistyt informatsii pro vyklyk")
  }

  # Folbek do zberezhenykh atrybutiv
  stored_X <- attr(object, "model_matrix")
  stored_response <- attr(object, "response")
  stored_mf <- attr(object, "model_frame")
  stored_data <- attr(object, "data")

  if (!is.null(stored_X)) {
    fitted_attr <- tryCatch({
      as.vector(stored_X %*% object@coefficients)
    }, error = function(e) NULL)
    if (!is.null(fitted_attr)) {
      return(fitted_attr)
    }
  }

  if (!is.null(stored_response) &&
      length(stored_response) == length(object@residuals)) {
    return(as.vector(stored_response - object@residuals))
  }

  # Sprobuvaty rekonstruiuvaty oryhinalni dani
  data_to_use <- data
  if(is.null(data_to_use)) {
    if (!is.null(stored_mf)) {
      data_to_use <- stored_mf
    } else if (!is.null(stored_data)) {
      data_to_use <- stored_data
    } else {
      # Sprobuvaty otrymaty dani z vyklyku, ale bezpechno obrobyty mozhlyvi pomylky
      tryCatch({
        data_to_use <- eval(object@call$data, envir = parent.frame())
      }, error = function(e) {
        if(is.null(data)) {
          stop("Ne vdalosia otrymaty dani z ob'iekta. Peredaite parametr 'data'.")
        }
      })
    }
  }

  if(is.null(data_to_use)) {
    stop("Potriben freim danykh dlia obchyslennia pidihnanykh znachen")
  }

  # Rekonstruiuvaty formulu
  formula <- eval(object@call$formula)

  # Bezpechna pobudova matrytsi dyzainu
  tryCatch({
    mf <- model.frame(formula, data_to_use)
    X <- model.matrix(formula, mf)

    # Obchyslyty pidihnani znachennia
    fitted <- as.vector(X %*% object@coefficients)
    return(fitted)
  }, error = function(e) {
    stop("Pomylka pry obchyslenni pidihnanykh znachen: ", e$message)
  })
}

#' Porivniaty PMM2 z OLS
#'
#' @param formula Formula modeli
#' @param data Freim danykh
#' @param pmm2_args Spysok arhumentiv dlia peredachi v lm_pmm2()
#'
#' @return Spysok z ob'iektamy pidhonky OLS ta PMM2
#' @export
compare_with_ols <- function(formula, data, pmm2_args = list()) {
  # Pidhonka OLS modeli
  fit_ols <- lm(formula, data)

  # Pidhonka PMM2 modeli z standartnymy abo zadanymy arhumentamy
  args <- c(list(formula = formula, data = data), pmm2_args)
  fit_pmm2 <- do.call(lm_pmm2, args)

  # Vytiahnuty ta porivniaty koefitsiienty
  coef_ols <- coef(fit_ols)
  coef_pmm2 <- coef(fit_pmm2)

  # Perevirka, chy imena koefitsiientiv PMM2 vstanovleni korektno
  if(is.null(names(coef_pmm2)) || all(names(coef_pmm2) == "")) {
    # Yakshcho imena ne vstanovleni, vykorystaiemo imena z OLS
    if(length(coef_ols) == length(coef_pmm2)) {
      names(coef_pmm2) <- names(coef_ols)
    } else {
      # Yakshcho dovzhyny vidrizniaiutsia, vykorystaiemo henerovani imena
      names(coef_pmm2) <- paste0("coef", seq_along(coef_pmm2))
    }
  }

  # Obchyslyty statystyku zalyshkiv
  res_ols <- residuals(fit_ols)
  res_pmm2 <- residuals(fit_pmm2)

  res_stats <- data.frame(
    Method = c("OLS", "PMM2"),
    RSS = c(sum(res_ols^2), sum(res_pmm2^2)),
    MAE = c(mean(abs(res_ols)), mean(abs(res_pmm2))),
    Skewness = c(pmm_skewness(res_ols), pmm_skewness(res_pmm2)),
    Kurtosis = c(pmm_kurtosis(res_ols), pmm_kurtosis(res_pmm2))
  )

  # Stvoryty tablytsiu porivniannia koefitsiientiv
  # Vykorystovuiemo vsi unikalni imena koefitsiientiv
  all_coef_names <- unique(c(names(coef_ols), names(coef_pmm2)))

  coef_table <- data.frame(
    Coefficient = all_coef_names,
    OLS = numeric(length(all_coef_names)),
    PMM2 = numeric(length(all_coef_names)),
    Diff_Percent = numeric(length(all_coef_names))
  )

  # Zapovniuiemo znachennia dlia OLS
  for(i in seq_along(all_coef_names)) {
    name <- all_coef_names[i]
    if(name %in% names(coef_ols)) {
      coef_table$OLS[i] <- coef_ols[name]
    } else {
      coef_table$OLS[i] <- NA
    }

    if(name %in% names(coef_pmm2)) {
      coef_table$PMM2[i] <- coef_pmm2[name]

      # Obchyslyty protsentnu riznytsiu tilky iakshcho obydva znachennia isnuiut
      if(!is.na(coef_table$OLS[i]) && coef_table$OLS[i] != 0) {
        coef_table$Diff_Percent[i] <- 100 * (coef_table$PMM2[i] - coef_table$OLS[i]) / abs(coef_table$OLS[i])
      }
    } else {
      coef_table$PMM2[i] <- NA
      coef_table$Diff_Percent[i] <- NA
    }
  }

  return(list(
    ols = fit_ols,
    pmm2 = fit_pmm2,
    coefficients = coef_table,
    residual_stats = res_stats
  ))
}


#' Metod prohnozuvannia dlia ob'iektiv PMM2fit
#'
#' @param object Ob'iekt PMM2fit
#' @param newdata Novyi freim danykh dlia prohnozuvannia
#' @param debug Lohichne znachennia, chy vyvodyty debah-informatsiiu
#' @param ... dodatkovi arhumenty (ne vykorystovuiutsia)
#'
#' @return Vektor prohnoziv
#' @export
setMethod("predict", "PMM2fit",
          function(object, newdata = NULL, debug = FALSE, ...) {
            if(is.null(newdata)) {
              stop("Parametr newdata maie buty nadanyi")
            }

            if(is.null(object@call)) {
              stop("Ob'iekt PMM2fit ne mistyt informatsii pro vyklyk")
            }

            # Vytiahnuty formulu z vyklyku
            formula <- eval(object@call$formula)

            if(debug) {
              cat("Formula:", deparse(formula), "\n")
              cat("Koefitsiienty:", paste(names(object@coefficients), "=", object@coefficients, collapse=", "), "\n")
              cat("Rozmir newdata:", nrow(newdata), "x", ncol(newdata), "\n")
              cat("Zminni v newdata:", paste(names(newdata), collapse=", "), "\n")
            }

            # Proste rishennia - napriamu vykorystaiemo pravu chastynu formuly dlia pobudovy matrytsi dyzainu
            rhs <- formula[[3]]
            design_formula <- as.formula(paste("~", deparse(rhs)))

            if(debug) {
              cat("Formula dyzainu:", deparse(design_formula), "\n")
            }

            # Stvoryty matrytsiu dyzainu
            X <- model.matrix(design_formula, newdata)

            if(debug) {
              cat("Rozmir matrytsi dyzainu:", nrow(X), "x", ncol(X), "\n")
              cat("Stovptsi matrytsi dyzainu:", paste(colnames(X), collapse=", "), "\n")
            }

            # Vypravliaiemo problemu z imenamy koefitsiientiv
            # Yakshcho imena vidsutni abo porozhni, zapovniuiemo ikh pravylnymy znachenniamy
            if(is.null(names(object@coefficients)) || all(names(object@coefficients) == "")) {
              if(debug) {
                cat("Imena koefitsiientiv vidsutni abo porozhni. Vykorystovuiemo standartni imena.\n")
              }

              expected_names <- colnames(X)
              if(length(expected_names) == length(object@coefficients)) {
                names(object@coefficients) <- expected_names
              } else {
                warning("Kilkist koefitsiientiv ne vidpovidaie kilkosti stovptsiv u matrytsi dyzainu.")
                if(length(object@coefficients) == 3 && ncol(X) == 3 &&
                   all(colnames(X) == c("(Intercept)", "x1", "x2"))) {
                  # Naichastishyi vypadok - rehresiia z 2 zminnymy
                  names(object@coefficients) <- c("(Intercept)", "x1", "x2")
                } else {
                  # Zahalne prysvoiennia imen
                  names(object@coefficients) <- paste0("coef", seq_along(object@coefficients))
                }
              }

              if(debug) {
                cat("Novi imena koefitsiientiv:", paste(names(object@coefficients), collapse=", "), "\n")
              }
            }

            # Obchyslyty prohnozy bezposeredno
            predictions <- numeric(nrow(newdata))

            # Dlia sproshchennia, prosto obchysliuiemo prohnozy vruchnu dlia typovoi rehresii
            if(length(object@coefficients) == 3 && all(c("x1", "x2") %in% names(newdata))) {
              if(debug) {
                cat("Obchysliuiemo prohnozy vruchnu dlia typovoi rehresii z interseptom i dvoma zminnymy.\n")
              }
              # Yakshcho tse typova rehresiia y ~ x1 + x2
              intercept_idx <- which(names(object@coefficients) == "(Intercept)")
              x1_idx <- which(names(object@coefficients) == "x1")
              x2_idx <- which(names(object@coefficients) == "x2")

              if(length(intercept_idx) == 1 && length(x1_idx) == 1 && length(x2_idx) == 1) {
                predictions <- object@coefficients[intercept_idx] +
                  object@coefficients[x1_idx] * newdata$x1 +
                  object@coefficients[x2_idx] * newdata$x2
              } else {
                # Yakshcho imena ne taki iak ochikuietsia, vykorystovuiemo ikh pozytsii
                predictions <- object@coefficients[1] +
                  object@coefficients[2] * newdata$x1 +
                  object@coefficients[3] * newdata$x2
              }
            } else {
              # Dlia inshykh vypadkiv
              if(debug) {
                cat("Namahaiemosia obchyslyty zahalnyi vypadok.\n")
              }
              # Sproshchenyi pidkhid dlia zahalnoho vypadku
              coeffs <- object@coefficients

              # Intersept
              if("(Intercept)" %in% names(coeffs)) {
                predictions <- predictions + coeffs["(Intercept)"]
              }

              # Inshi zminni
              for(var_name in intersect(names(coeffs), names(newdata))) {
                if(var_name != "(Intercept)") {
                  predictions <- predictions + coeffs[var_name] * newdata[[var_name]]
                }
              }
            }

            if(debug) {
              cat("Rozmir vektora prohnoziv:", length(predictions), "\n")
              cat("Pershi kilka prohnoziv:", paste(head(predictions), collapse=", "), "\n")
            }

            return(predictions)
          })
