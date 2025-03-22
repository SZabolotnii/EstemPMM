#' An S4 class to store PMM2 fit results
#'
#' @slot coefficients numeric vector of fitted parameters
#' @slot residuals numeric vector of final residuals
#' @slot m2 numeric second central moment from the initial residuals
#' @slot m3 numeric third central moment from the initial residuals
#' @slot m4 numeric fourth central moment from the initial residuals
#' @slot convergence logical or integer code indicating if the algorithm converged
#' @slot iterations integer number of iterations performed
#' @slot call the original function call
#'
#' @exportClass PMM2fit
setClass("PMM2fit",
         slots = c(coefficients = "numeric",
                   residuals    = "numeric",
                   m2          = "numeric",
                   m3          = "numeric",
                   m4          = "numeric",
                   convergence = "logical",
                   iterations  = "numeric",
                   call        = "call"))

#' Summarize PMM2 results
#'
#' @param object an object of class "PMM2fit"
#' @param formula (optional) the formula used for the model
#' @param data (optional) the data used
#' @param B number of bootstrap replicates for inference
#' @param ... additional arguments (not used)
#'
#' @return Prints a summary to console; returns the object (invisibly).
#'
#' @export
setMethod("summary", "PMM2fit",
          function(object, formula=NULL, data=NULL, B=100, ...) {
            cat("Polynomial Maximization Method (PMM2) fit\n")
            if(!is.null(object@call)) {
              cat("Call:\n")
              print(object@call)
              cat("\n")
            }

            cat("Coefficients:\n")
            print(object@coefficients)

            cat("\nCentral moments of initial residuals:\n")
            cat("  m2 =", object@m2, "\n")
            cat("  m3 =", object@m3, "\n")
            cat("  m4 =", object@m4, "\n\n")

            cat("Algorithm information:\n")
            cat("  Convergence status:", object@convergence, "\n")
            # Перевіряємо наявність слоту iterations напряму
            if(is(object, "PMM2fit") && "iterations" %in% slotNames(object)) {
              cat("  Iterations:", object@iterations, "\n\n")
            } else {
              cat("\n")
            }

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

#' Predict method for PMM2 fit objects
#'
#' @param object A PMM2fit object
#' @param newdata A data frame in which to look for variables with which to predict
#' @param ... additional arguments (not used)
#'
#' @return A vector of predictions
#'
#' @export
setMethod("predict", "PMM2fit",
          function(object, newdata = NULL, ...) {
            if(is.null(newdata)) {
              stop("newdata must be provided")
            }

            if(is.null(object@call)) {
              stop("PMM2fit object does not contain call information")
            }

            # Extract the formula from the call
            formula <- eval(object@call$formula)

            # Create model matrix from newdata
            mf <- model.frame(formula, newdata, na.action = na.pass)
            X <- model.matrix(formula, mf)

            # Calculate predictions
            as.vector(X %*% object@coefficients)
          }
)

#' Plot diagnostics for PMM2 fit objects
#'
#' @param x A PMM2fit object
#' @param y Unused (for S4 method compatibility)
#' @param which Integer vector specifying which plots to produce
#' @param ... additional arguments passed to plotting functions
#'
#' @return Invisibly returns x
#'
#' @export
setMethod("plot", signature(x = "PMM2fit", y = "missing"),
          function(x, y, which = c(1:4), ...) {
            op <- par(no.readonly = TRUE)
            on.exit(par(op))

            par(mfrow = c(2, 2))

            residuals <- x@residuals
            fitted <- fitted_values(x)

            # Plot 1: Residuals vs Fitted
            if(1 %in% which) {
              plot(fitted, residuals,
                   main = "Residuals vs Fitted",
                   xlab = "Fitted values",
                   ylab = "Residuals",
                   ...)
              abline(h = 0, lty = 2)
              lines(lowess(fitted, residuals), col = "red")
            }

            # Plot 2: Normal Q-Q Plot
            if(2 %in% which) {
              qqnorm(residuals, main = "Normal Q-Q Plot", ...)
              qqline(residuals)
            }

            # Plot 3: Scale-Location Plot
            if(3 %in% which) {
              plot(fitted, sqrt(abs(residuals)),
                   main = "Scale-Location Plot",
                   xlab = "Fitted values",
                   ylab = expression(sqrt("|Residuals|")),
                   ...)
              lines(lowess(fitted, sqrt(abs(residuals))), col = "red")
            }

            # Plot 4: Histogram of residuals
            if(4 %in% which) {
              hist(residuals,
                   main = "Histogram of Residuals",
                   xlab = "Residuals",
                   breaks = "FD",
                   ...)
            }

            invisible(x)
          }
)

#' Helper function to extract fitted values
#'
#' @param object A PMM2fit object
#' @return A vector of fitted values
#'
#' @keywords internal
fitted_values <- function(object) {
  if(is.null(object@call)) {
    stop("PMM2fit object does not contain call information")
  }

  # Reconstruct the original data
  formula <- eval(object@call$formula)
  data <- eval(object@call$data)

  # Build design matrix
  mf <- model.frame(formula, data)
  X <- model.matrix(formula, mf)

  # Calculate fitted values
  as.vector(X %*% object@coefficients)
}
