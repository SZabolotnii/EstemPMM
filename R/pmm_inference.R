#' Bootstrap inference for a PMM2 fit
#'
#' @param object an object of class PMM2fit
#' @param formula the same formula used originally
#' @param data the data frame used originally
#' @param B number of bootstrap replications
#' @param seed (optional) for reproducibility
#' @param parallel logical, whether to use parallel computation
#' @param cores number of cores to use for parallel computation, default is detected
#'
#' @return A data.frame with columns: Estimate, Std.Error, t.value, p.value
#' @export
pmm2_inference <- function(object, formula, data, B=200, seed=NULL,
                           parallel=FALSE, cores=NULL) {
  # Set seed for reproducibility if provided
  if(!is.null(seed)) set.seed(seed)

  # Extract coefficients and residuals
  coefs <- object@coefficients
  res   <- object@residuals

  # Input validation
  if(B < 10) {
    warning("Number of bootstrap samples (B) is very low. Consider using B >= 100 for more reliable inference.")
  }

  if(!inherits(object, "PMM2fit")) {
    stop("Object must be of class 'PMM2fit'")
  }

  if(missing(formula) || missing(data)) {
    stop("Both 'formula' and 'data' must be provided")
  }

  # Build X, y matrices
  mf <- model.frame(formula, data)
  X <- model.matrix(formula, mf)
  y <- model.response(mf)
  n <- nrow(X)

  # Early return in case of errors
  if(is.null(y) || is.null(X)) {
    stop("Failed to extract response or design matrix from data")
  }

  # Check if we should use parallel computation
  use_parallel <- parallel && requireNamespace("parallel", quietly = TRUE)

  if(use_parallel) {
    if(is.null(cores)) {
      cores <- max(1, parallel::detectCores() - 1)
    }

    boot_results <- parallel::mclapply(seq_len(B), function(b) {
      # 1) Bootstrap residuals
      res_b <- sample(res, size=n, replace=TRUE)

      # 2) Create new y
      y_b <- X %*% coefs + res_b

      # 3) Create new data
      data_b <- data
      # Assume the left-hand side is the first term in formula
      lhs <- as.character(formula[[2]])
      data_b[[lhs]] <- as.numeric(y_b)

      # 4) Re-estimate model
      fit_b <- tryCatch({
        lm_pmm2(formula, data_b, max_iter=20, tol=1e-6)
      }, error = function(e) {
        warning("Bootstrap replicate ", b, " failed: ", e$message)
        return(NULL)
      })

      if(!is.null(fit_b)) {
        return(fit_b@coefficients)
      } else {
        return(rep(NA, length(coefs)))
      }
    }, mc.cores = cores)

    # Convert list to matrix
    boot_est <- do.call(rbind, boot_results)

  } else {
    # Sequential computation
    # Matrix to store results
    boot_est <- matrix(0, nrow=B, ncol=length(coefs))
    colnames(boot_est) <- names(coefs)

    # Progress tracking
    pb <- NULL
    if(interactive() && B > 10) {
      if(requireNamespace("utils", quietly = TRUE)) {
        pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
      }
    }

    for(b in seq_len(B)) {
      # 1) Bootstrap residuals
      res_b <- sample(res, size=n, replace=TRUE)

      # 2) Create new y
      y_b <- X %*% coefs + res_b

      # 3) Create new data
      data_b <- data
      # Assume the left-hand side is the first term in formula
      lhs <- as.character(formula[[2]])
      data_b[[lhs]] <- as.numeric(y_b)

      # 4) Re-estimate model
      fit_b <- tryCatch({
        lm_pmm2(formula, data_b, max_iter=20, tol=1e-6)
      }, error = function(e) {
        warning("Bootstrap replicate ", b, " failed: ", e$message)
        return(NULL)
      })

      if(!is.null(fit_b)) {
        boot_est[b, ] <- fit_b@coefficients
      } else {
        boot_est[b, ] <- NA
      }

      # Update progress bar
      if(!is.null(pb)) utils::setTxtProgressBar(pb, b)
    }

    # Close progress bar
    if(!is.null(pb)) close(pb)
  }

  # Remove rows with NA values
  na_rows <- apply(boot_est, 1, function(row) any(is.na(row)))
  if(any(na_rows)) {
    warning("Removed ", sum(na_rows), " bootstrap replicates due to estimation failures")
    boot_est <- boot_est[!na_rows, , drop = FALSE]
  }

  # Check if we have enough successful bootstraps
  if(nrow(boot_est) < 10) {
    stop("Too few successful bootstrap replicates to compute reliable inference")
  }

  # Compute covariance matrix and standard errors
  cov_mat <- cov(boot_est)
  est <- coefs
  se  <- sqrt(diag(cov_mat))

  # Compute t-values and p-values
  t_val <- est / se
  # For large samples, use normal approximation
  p_val <- 2 * (1 - pnorm(abs(t_val)))

  # Create output data frame
  out <- data.frame(
    Estimate  = est,
    Std.Error = se,
    t.value   = t_val,
    p.value   = p_val
  )
  rownames(out) <- names(est)

  # Calculate confidence intervals
  ci <- t(apply(boot_est, 2, quantile, probs = c(0.025, 0.975)))
  colnames(ci) <- c("2.5%", "97.5%")

  # Add confidence intervals to output
  out$conf.low <- ci[, "2.5%"]
  out$conf.high <- ci[, "97.5%"]

  return(out)
}

#' Plot bootstrap distributions for a PMM2 fit
#'
#' @param object Result from pmm2_inference
#' @param coefficients Which coefficients to plot, defaults to all
#'
#' @return Invisibly returns the histogram information
#' @export
plot_pmm2_bootstrap <- function(object, coefficients = NULL) {
  if(!inherits(object, "data.frame") ||
     !all(c("Estimate", "Std.Error", "conf.low", "conf.high") %in% names(object))) {
    stop("Object must be result from pmm2_inference()")
  }

  # If no coefficients specified, use all
  if(is.null(coefficients)) {
    coefficients <- rownames(object)
  }

  # Filter to requested coefficients
  object_subset <- object[coefficients, , drop = FALSE]

  # Set up plot layout
  n_coefs <- nrow(object_subset)
  n_cols <- min(2, n_coefs)
  n_rows <- ceiling(n_coefs / n_cols)

  # Save old par settings and restore on exit
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  par(mfrow = c(n_rows, n_cols))

  # Create a density plot for each coefficient
  result <- lapply(seq_len(n_coefs), function(i) {
    coef_name <- rownames(object_subset)[i]
    est <- object_subset[i, "Estimate"]
    ci_low <- object_subset[i, "conf.low"]
    ci_high <- object_subset[i, "conf.high"]

    # Create plot title
    main_title <- paste0(coef_name, "\nEst: ", round(est, 4))

    # Estimate range for x-axis
    range_val <- c(ci_low - 0.5 * (est - ci_low),
                   ci_high + 0.5 * (ci_high - est))

    # Create distribution visualization
    x_seq <- seq(range_val[1], range_val[2], length.out = 100)
    se <- object_subset[i, "Std.Error"]
    y_seq <- dnorm(x_seq, mean = est, sd = se)

    # Plot
    plot(x_seq, y_seq, type = "l",
         main = main_title,
         xlab = "Value",
         ylab = "Density")

    # Add vertical lines for estimate and CI
    abline(v = est, col = "red", lwd = 2)
    abline(v = ci_low, col = "blue", lty = 2)
    abline(v = ci_high, col = "blue", lty = 2)

    # Add legend
    legend("topright",
           legend = c("Estimate", "95% CI"),
           col = c("red", "blue"),
           lty = c(1, 2),
           lwd = c(2, 1),
           cex = 0.8)

    invisible(list(x = x_seq, y = y_seq, estimate = est,
                   ci_low = ci_low, ci_high = ci_high))
  })

  names(result) <- rownames(object_subset)
  invisible(result)
}
