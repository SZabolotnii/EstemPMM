# pmm_ts_examples.R

#' Examples of applying the PMM2 method for time series analysis
#'
#' This script demonstrates the use of functions from the EstemPMM package
#' for modeling time series with non-Gaussian distributions

# Load package
library(EstemPMM)

#' Demonstration of PMM2 for analyzing AR, MA, ARMA, and ARIMA models
#' with different types of distributions
run_ts_examples <- function() {
  # Set seed for reproducibility
  set.seed(123)

  cat("\n============ PMM2 DEMONSTRATION FOR TIME SERIES ============\n\n")

  # 1. AR(2) model with t-distribution (heavy tails)
  cat("1. AR(2) model with t-distributed innovations\n")
  ar_coef <- c(0.7, -0.3)
  n <- 300

  tryCatch({
    ar_series <- as.numeric(arima.sim(model = list(ar = ar_coef), n = n,
                                      rand.gen = function(n) rt(n, df=3)))

    # Compare AR model estimation methods
    ar_comparison <- compare_ar_methods(ar_series, order = 2)

    # Output results
    cat("\nComparison of estimated AR(2) model coefficients:\n")
    print(ar_comparison$coefficients)

    cat("\nResidual statistics:\n")
    print(ar_comparison$residual_stats)

    # Create plot
    par(mfrow=c(2,2))
    plot(ar_comparison$pmm2)

    # Store for later prediction
    ar_model <- ar_comparison$pmm2

  }, error = function(e) {
    cat("\nError in AR model analysis:", conditionMessage(e), "\n")
    ar_model <<- NULL
  })

  # 2. MA(2) model with normal distribution (for stability)
  cat("\n\n2. MA(2) model with normally distributed innovations\n")

  tryCatch({
    ma_coef <- c(0.6, -0.2)
    # Use normal innovations for stability
    normal_innov <- rnorm(n)
    ma_series <- as.numeric(arima.sim(model = list(ma = ma_coef), n = n,
                                      innov = normal_innov))

    # Compare MA model estimation methods
    ma_comparison <- compare_ma_methods(ma_series, order = 2)

    # Output results
    cat("\nComparison of estimated MA(2) model coefficients:\n")
    print(ma_comparison$coefficients)

    cat("\nResidual statistics:\n")
    print(ma_comparison$residual_stats)

    # Create plot
    par(mfrow=c(2,2))
    plot(ma_comparison$pmm2)

    # Store for later prediction
    ma_model <- ma_comparison$pmm2

  }, error = function(e) {
    cat("\nError in MA model analysis:", conditionMessage(e), "\n")
    ma_model <<- NULL
  })

  # 3. ARMA(1,1) model with mixture of normal distributions
  cat("\n\n3. ARMA(1,1) model with mixture of normal distributions\n")

  tryCatch({
    # Generate innovations as mixture of two normal distributions
    mix_innov <- numeric(n)
    for(i in 1:n) {
      # 70% from N(0,1) and 30% from N(3,2)
      mix_innov[i] <- ifelse(runif(1) < 0.7, rnorm(1), rnorm(1, mean=3, sd=2))
    }
    # Center for zero mean
    mix_innov <- mix_innov - mean(mix_innov)

    # Create ARMA series
    arma_series <- as.numeric(arima.sim(model = list(ar = 0.7, ma = 0.4), n = n,
                                        innov = mix_innov))

    # Compare ARMA model estimation methods
    arma_comparison <- compare_arma_methods(arma_series, order = c(1, 1))

    # Output results
    cat("\nComparison of estimated ARMA(1,1) model coefficients:\n")
    print(arma_comparison$coefficients)

    cat("\nResidual statistics:\n")
    print(arma_comparison$residual_stats)

    # Create plot
    par(mfrow=c(2,2))
    plot(arma_comparison$pmm2)

    # Store for later prediction
    arma_model <- arma_comparison$pmm2

  }, error = function(e) {
    cat("\nError in ARMA model analysis:", conditionMessage(e), "\n")
    arma_model <<- NULL
  })

  # 4. ARIMA(1,1,1) model with asymmetric distribution
  cat("\n\n4. ARIMA(1,1,1) model with asymmetric distribution\n")

  tryCatch({
    # Generate ARMA series with normal errors (for stability)
    arma_base <- as.numeric(arima.sim(model = list(ar = 0.7, ma = 0.4), n = n))

    # Transform to non-stationary series through integration (cumulative sum)
    arima_series <- as.numeric(cumsum(arma_base))

    # Compare ARIMA model estimation methods
    arima_comparison <- compare_arima_methods(arima_series, order = c(1, 1, 1))

    # Output results
    cat("\nComparison of estimated ARIMA(1,1,1) model coefficients:\n")
    print(arima_comparison$coefficients)

    cat("\nResidual statistics:\n")
    print(arima_comparison$residual_stats)

    # Create plot
    par(mfrow=c(3,2))
    plot(arima_comparison$pmm2, which=1:6)

    # Store for later prediction
    arima_model <- arima_comparison$pmm2

  }, error = function(e) {
    cat("\nError in ARIMA model analysis:", conditionMessage(e), "\n")
    arima_model <<- NULL
  })

  # 5. Forecasting
  cat("\n\n5. Forecasting with PMM2 models\n")

  # AR forecast
  if(exists("ar_model") && !is.null(ar_model)) {
    tryCatch({
      ar_forecast <- predict(ar_model, n.ahead = 10)
      cat("AR(2) forecast for 10 steps ahead:\n")
      print(ar_forecast)
    }, error = function(e) {
      cat("Error in AR forecast:", conditionMessage(e), "\n")
    })
  } else {
    cat("AR model not available for forecasting\n")
  }

  # MA forecast
  if(exists("ma_model") && !is.null(ma_model)) {
    tryCatch({
      ma_forecast <- predict(ma_model, n.ahead = 10)
      cat("\nMA(2) forecast for 10 steps ahead:\n")
      print(ma_forecast)
    }, error = function(e) {
      cat("Error in MA forecast:", conditionMessage(e), "\n")
    })
  } else {
    cat("\nMA model not available for forecasting\n")
  }

  # ARMA forecast
  if(exists("arma_model") && !is.null(arma_model)) {
    tryCatch({
      arma_forecast <- predict(arma_model, n.ahead = 10)
      cat("\nARMA(1,1) forecast for 10 steps ahead:\n")
      print(arma_forecast)
    }, error = function(e) {
      cat("Error in ARMA forecast:", conditionMessage(e), "\n")
    })
  } else {
    cat("\nARMA model not available for forecasting\n")
  }

  # ARIMA forecast
  if(exists("arima_model") && !is.null(arima_model)) {
    tryCatch({
      arima_forecast <- predict(arima_model, n.ahead = 10)
      cat("\nARIMA(1,1,1) forecast for 10 steps ahead:\n")
      print(arima_forecast)
    }, error = function(e) {
      cat("Error in ARIMA forecast:", conditionMessage(e), "\n")
    })
  } else {
    cat("\nARIMA model not available for forecasting\n")
  }

  # 6. PMM2 efficiency plot for different skewness and kurtosis values
  cat("\n\n6. PMM2 efficiency as a function of skewness and kurtosis\n")

  tryCatch({
    # Create grid of skewness and kurtosis values
    skewness_values <- seq(0, 2, by=0.2)  # from 0 to 2
    kurtosis_values <- seq(0, 6, by=0.5)   # from 0 to 6

    # Calculate PMM2 efficiency coefficient
    g2_values <- matrix(0, length(skewness_values), length(kurtosis_values))

    for(i in 1:length(skewness_values)) {
      for(j in 1:length(kurtosis_values)) {
        g2 <- 1 - (skewness_values[i]^2) / (2 + kurtosis_values[j])
        # Limit values between 0 and 1
        g2_values[i, j] <- max(0, min(1, g2))
      }
    }

    # Visualize the relationship
    filled.contour(
      skewness_values,
      kurtosis_values,
      g2_values,
      color.palette = colorRampPalette(c("red", "yellow", "green")),
      xlab = "Skewness (γ₃)",
      ylab = "Kurtosis (γ₄)",
      main = "Ratio of PMM2 to classical estimator variances (g²)",
      key.title = title("g²")
    )

    # Add line corresponding to equality γ₄ + 2 = γ₃²
    contour(skewness_values, kurtosis_values,
            outer(skewness_values^2, rep(1, length(kurtosis_values))) -
              outer(rep(1, length(skewness_values)), kurtosis_values) - 2,
            levels = 0, add = TRUE, lwd = 2, col = "black", lty = 2)

    cat("Plot created. Red color indicates high PMM2 efficiency,\n")
    cat("while blue indicates nearly equivalent efficiency to classical methods.\n")
    cat("Dashed line shows the boundary of permissible values γ₄ + 2 ≥ γ₃²\n")
  }, error = function(e) {
    cat("Error in efficiency plot:", conditionMessage(e), "\n")
  })

  cat("\n============ END OF DEMONSTRATION ============\n")
}

# Run demonstration
run_ts_examples()
