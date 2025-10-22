# pmm2_classes.R - Iierarkhiia klasiv dlia modelei PMM2

#' Bazovyi klas S4 dlia zberihannia rezultativ PMM2 modelei
#'
#' @slot coefficients chyslovyi vektor otsinenykh parametriv
#' @slot residuals chyslovyi vektor kintsevykh zalyshkiv
#' @slot m2 chyslovyi druhyi tsentralnyi moment pochatkovykh zalyshkiv
#' @slot m3 chyslovyi tretii tsentralnyi moment pochatkovykh zalyshkiv
#' @slot m4 chyslovyi chetvertyi tsentralnyi moment pochatkovykh zalyshkiv
#' @slot convergence lohichnyi abo tsilyi kod, shcho vkazuie, chy alhorytm zbihsia
#' @slot iterations chyslova kilkist vykonanykh iteratsii
#' @slot call oryhinalnyi vyklyk funktsii
#'
#' @exportClass BasePMM2
setClass("BasePMM2",
         slots = c(coefficients = "numeric",
                   residuals    = "numeric",
                   m2           = "numeric",
                   m3           = "numeric",
                   m4           = "numeric",
                   convergence  = "logical",
                   iterations   = "numeric",
                   call         = "call"))

#' Klas S4 dlia zberihannia rezultativ PMM2 rehresiinoi modeli
#'
#' @slot coefficients chyslovyi vektor otsinenykh parametriv
#' @slot residuals chyslovyi vektor kintsevykh zalyshkiv
#' @slot m2 chyslovyi druhyi tsentralnyi moment pochatkovykh zalyshkiv
#' @slot m3 chyslovyi tretii tsentralnyi moment pochatkovykh zalyshkiv
#' @slot m4 chyslovyi chetvertyi tsentralnyi moment pochatkovykh zalyshkiv
#' @slot convergence lohichnyi abo tsilyi kod, shcho vkazuie, chy alhorytm zbihsia
#' @slot iterations chyslova kilkist vykonanykh iteratsii
#' @slot call oryhinalnyi vyklyk funktsii
#'
#' @exportClass PMM2fit
setClass("PMM2fit",
         contains = "BasePMM2")

#' Bazovyi klas S4 dlia zberihannia rezultativ PMM2 modelei chasovykh riadiv
#'
#' @slot coefficients chyslovyi vektor otsinenykh parametriv
#' @slot residuals chyslovyi vektor kintsevykh zalyshkiv
#' @slot m2 chyslovyi druhyi tsentralnyi moment pochatkovykh zalyshkiv
#' @slot m3 chyslovyi tretii tsentralnyi moment pochatkovykh zalyshkiv
#' @slot m4 chyslovyi chetvertyi tsentralnyi moment pochatkovykh zalyshkiv
#' @slot convergence lohichnyi abo tsilyi kod, shcho vkazuie, chy alhorytm zbihsia
#' @slot iterations chyslova kilkist vykonanykh iteratsii
#' @slot call oryhinalnyi vyklyk funktsii
#' @slot model_type symvolnyi riadok, shcho vkazuie typ modeli
#' @slot intercept chyslove znachennia perekhoplennia
#' @slot original_series chyslovyi vektor oryhinalnoho chasovoho riadu
#' @slot order spysok parametriv poriadku
#'
#' @exportClass TS2fit
setClass("TS2fit",
         contains = "BasePMM2",
         slots = c(model_type      = "character",
                   intercept       = "numeric",
                   original_series = "numeric",
                   order           = "list"))

#' Klas S4 dlia zberihannia rezultativ PMM2 modeli AR
#'
#' @exportClass ARPMM2
setClass("ARPMM2", contains = "TS2fit")

#' Klas S4 dlia zberihannia rezultativ PMM2 modeli MA
#'
#' @exportClass MAPMM2
setClass("MAPMM2", contains = "TS2fit")

#' Klas S4 dlia zberihannia rezultativ PMM2 modeli ARMA
#'
#' @exportClass ARMAPMM2
setClass("ARMAPMM2", contains = "TS2fit")

#' Klas S4 dlia zberihannia rezultativ PMM2 modeli ARIMA
#'
#' @exportClass ARIMAPMM2
setClass("ARIMAPMM2", contains = "TS2fit")

#' Uzahalnenyi metod summary dlia ob'iektiv PMM2fit
#'
#' @param object ob'iekt klasu "PMM2fit"
#' @param formula (optsionalno) formula, vykorystana dlia modeli
#' @param data (optsionalno) vykorystani dani
#' @param B kilkist butstrep-replikatsii dlia statystychnoho vysnovku
#' @param ... dodatkovi arhumenty (ne vykorystovuiutsia)
#'
#' @return Vyvodyt reziume na konsol; povertaie ob'iekt (nevydymo).
#'
#' @export
setMethod("summary", "PMM2fit",
          function(object, formula=NULL, data=NULL, B=100, ...) {
            cat("Rezultaty otsiniuvannia metodom PMM2\n")
            if(!is.null(object@call)) {
              cat("Vyklyk:\n")
              print(object@call)
              cat("\n")
            }

            cat("Koefitsiienty:\n")
            print(object@coefficients)

            cat("\nTsentralni momenty pochatkovykh zalyshkiv:\n")
            cat("  m2 =", object@m2, "\n")
            cat("  m3 =", object@m3, "\n")
            cat("  m4 =", object@m4, "\n\n")

            vf <- pmm2_variance_factor(object@m2, object@m3, object@m4)
            if(!is.na(vf$g)) {
              cat("Teoretychni kharakterystyky PMM2 (S = 2):\n")
              cat("  c3 =", vf$c3, "\n")
              cat("  c4 =", vf$c4, "\n")
              cat("  g  =", vf$g, " (ochikuvane vidnoshennia Var[PMM2]/Var[OLS])\n\n")
            }

            cat("Informatsiia pro alhorytm:\n")
            cat("  Status zbizhnosti:", object@convergence, "\n")
            if("iterations" %in% slotNames(object)) {
              cat("  Iteratsii:", object@iterations, "\n\n")
            } else {
              cat("\n")
            }

            # Yakshcho korystuvach khoche pobachyty p-znachennia, vyklykaiemo pmm2_inference:
            if(!is.null(formula) && !is.null(data)) {
              cat("Nablyzhenyi statystychnyi vysnovok cherez butstrep (B=", B, "):\n", sep="")
              inf_tab <- pmm2_inference(object, formula, data, B=B)
              print(inf_tab)
            } else {
              cat("Dlia perehliadu p-znachen, peredaite formula= ta data=\n")
            }
            invisible(object)
          }
)

#' Uzahalnenyi metod summary dlia ob'iektiv TS2fit
#'
#' @param object ob'iekt klasu "TS2fit" abo pidklasu
#' @param ... dodatkovi arhumenty (ne vykorystovuiutsia)
#'
#' @return Vyvodyt reziume na konsol; povertaie ob'iekt (nevydymo).
#'
#' @export
setMethod("summary", "TS2fit",
          function(object, ...) {
            model_type <- object@model_type
            ar_order <- object@order$ar
            ma_order <- object@order$ma
            d <- object@order$d

            cat("Rezultaty otsiniuvannia chasovoho riadu metodom PMM2\n")

            # Vyvesty typ modeli
            cat("Typ modeli: ")
            if (model_type == "ar") {
              cat("AR(", ar_order, ")\n", sep = "")
            } else if (model_type == "ma") {
              cat("MA(", ma_order, ")\n", sep = "")
            } else if (model_type == "arma") {
              cat("ARMA(", ar_order, ",", ma_order, ")\n", sep = "")
            } else if (model_type == "arima") {
              cat("ARIMA(", ar_order, ",", d, ",", ma_order, ")\n", sep = "")
            }

            if(!is.null(object@call)) {
              cat("Vyklyk:\n")
              print(object@call)
              cat("\n")
            }

            # Vyvesty koefitsiienty z rozdilenniam na AR ta MA chastyny
            cat("Koefitsiienty:\n")
            if (ar_order > 0) {
              ar_coefs <- object@coefficients[1:ar_order]
              names(ar_coefs) <- paste0("ar", 1:ar_order)
              cat("AR: ")
              print(ar_coefs)
            }

            if (ma_order > 0) {
              ma_coefs <- object@coefficients[(ar_order+1):(ar_order+ma_order)]
              names(ma_coefs) <- paste0("ma", 1:ma_order)
              cat("MA: ")
              print(ma_coefs)
            }

            # Vyvesty perekhoplennia, iakshcho ie
            if (object@intercept != 0) {
              cat("Perekhoplennia: ", object@intercept, "\n")
            }

            cat("\nTsentralni momenty pochatkovykh zalyshkiv:\n")
            cat("  m2 =", object@m2, "\n")
            cat("  m3 =", object@m3, "\n")
            cat("  m4 =", object@m4, "\n\n")

            vf <- pmm2_variance_factor(object@m2, object@m3, object@m4)
            if(!is.na(vf$g)) {
              cat("Teoretychni kharakterystyky PMM2 (S = 2):\n")
              cat("  c3 =", vf$c3, "\n")
              cat("  c4 =", vf$c4, "\n")
              cat("  g  =", vf$g, " (ochikuvane vidnoshennia Var[PMM2]/Var[OLS])\n\n")
            }

            cat("Informatsiia pro alhorytm:\n")
            cat("  Status zbizhnosti:", object@convergence, "\n")
            if("iterations" %in% slotNames(object)) {
              cat("  Iteratsii:", object@iterations, "\n\n")
            } else {
              cat("\n")
            }

            invisible(object)
          }
)
