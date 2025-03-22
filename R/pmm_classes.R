#' An S4 class to store PMM2 fit results
#'
#' @slot coefficients numeric vector of fitted parameters
#' @slot residuals numeric vector of final residuals
#' @slot m2,m3,m4 numeric central moments from the initial residuals
#' @slot convergence logical or integer code
#'
#' @exportClass PMM2fit
setClass("PMM2fit",
         slots = c(coefficients = "numeric",
                   residuals    = "numeric",
                   m2          = "numeric",
                   m3          = "numeric",
                   m4          = "numeric",
                   convergence = "logical"))

#' Summarize PMM2 results
#'
#' @param object an object of class "PMM2fit"
#' @param formula (optional) the formula used for the model
#' @param data (optional) the data used
#' @param B number of bootstrap replicates for inference
#'
#' @return Prints a summary to console; returns the object (invisibly).
#'
#' @export
setMethod("summary", "PMM2fit",
          function(object, formula=NULL, data=NULL, B=100) {
            cat("Polynomial Maximization Method (PMM2) fit\n")
            cat("Coefficients:\n")
            print(object@coefficients)

            cat("\nCentral moments of initial residuals:\n")
            cat("  m2 =", object@m2, "\n")
            cat("  m3 =", object@m3, "\n")
            cat("  m4 =", object@m4, "\n\n")

            cat("Convergence status:", object@convergence, "\n\n")

            # Якщо користувач хоче побачити p-value, викликаємо pmm2_inference:
            if(!is.null(formula) && !is.null(data)) {
              cat("Approx. inference via bootstrap (B=", B, "):\n", sep="")
              inf_tab <- pmm2_inference(object, formula, data, B=B)
              print(inf_tab)
            } else {
              cat("To see p-values, pass formula= and data=\n")
            }
            invisible(object)
          }
)
