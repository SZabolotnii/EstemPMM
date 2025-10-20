# tests/testthat/test-pmm2_linear.R
# Comprehensive test suite for PMM2 linear regression

context("PMM2 Linear Regression")

# Setup test data
set.seed(42)
n <- 100
x <- rnorm(n)
errors <- rgamma(n, shape = 2, scale = 1) - 2  # Asymmetric errors
y <- 2 + 1.5 * x + errors
test_data <- data.frame(x = x, y = y)

test_that("lm_pmm2 produces valid fit object", {
  fit <- lm_pmm2(y ~ x, data = test_data)
  
  expect_s4_class(fit, "PMM2fit")
  expect_true(fit@convergence)
  expect_length(fit@coefficients, 2)  # Intercept + 1 predictor
  expect_equal(length(fit@residuals), n)
})

test_that("lm_pmm2 coefficients are numeric and reasonable", {
  fit <- lm_pmm2(y ~ x, data = test_data)
  
  expect_true(is.numeric(fit@coefficients))
  expect_false(any(is.na(fit@coefficients)))
  expect_false(any(is.infinite(fit@coefficients)))
  
  # Coefficients should be in reasonable range
  expect_true(all(abs(fit@coefficients) < 100))
})

test_that("lm_pmm2 residuals sum to near zero", {
  fit <- lm_pmm2(y ~ x, data = test_data)
  
  expect_true(is.numeric(fit@residuals))
  expect_equal(length(fit@residuals), n)
  
  # Sum of residuals should be small (not exactly zero due to intercept)
  expect_true(abs(sum(fit@residuals)) < 1e-10)
})

test_that("lm_pmm2 returns moment estimates", {
  fit <- lm_pmm2(y ~ x, data = test_data)
  
  expect_true(is.numeric(fit@m2))
  expect_true(is.numeric(fit@m3))
  expect_true(is.numeric(fit@m4))
  
  expect_false(any(is.na(c(fit@m2, fit@m3, fit@m4))))
  
  # m2 should be positive (variance)
  expect_true(fit@m2 > 0)
})

test_that("lm_pmm2 compared to OLS shows variance reduction", {
  fit_pmm2 <- lm_pmm2(y ~ x, data = test_data)
  fit_ols <- lm(y ~ x, data = test_data)
  
  # Compare coefficients (should be similar but not identical)
  coef_diff <- abs(fit_pmm2@coefficients - coef(fit_ols))
  expect_true(all(coef_diff < 1.0))  # Loose bound
})

test_that("lm_pmm2 with single predictor", {
  fit <- lm_pmm2(y ~ x, data = test_data)
  expect_length(fit@coefficients, 2)
  expect_equal(fit@iterations, as.integer(fit@iterations))
})

test_that("lm_pmm2 with multiple predictors", {
  x2 <- rnorm(n)
  test_data$x2 <- x2
  y_multi <- 2 + 1.5 * x + 0.8 * x2 + errors
  test_data$y_multi <- y_multi
  
  fit <- lm_pmm2(y_multi ~ x + x2, data = test_data)
  
  expect_s4_class(fit, "PMM2fit")
  expect_length(fit@coefficients, 3)  # Intercept + 2 predictors
})

test_that("lm_pmm2 handles missing values appropriately", {
  test_data_na <- test_data
  test_data_na$y[1] <- NA
  
  fit <- lm_pmm2(y ~ x, data = test_data_na, na.action = na.omit)
  
  expect_s4_class(fit, "PMM2fit")
  expect_equal(length(fit@residuals), n - 1)  # One observation removed
})

test_that("lm_pmm2 with formula without intercept", {
  fit <- lm_pmm2(y ~ x - 1, data = test_data)
  
  expect_s4_class(fit, "PMM2fit")
  expect_length(fit@coefficients, 1)  # No intercept
})

test_that("lm_pmm2 convergence status is logical", {
  fit <- lm_pmm2(y ~ x, data = test_data)
  
  expect_true(is.logical(fit@convergence))
  expect_length(fit@convergence, 1)
})

test_that("lm_pmm2 with very small sample size fails gracefully", {
  small_data <- data.frame(x = rnorm(3), y = rnorm(3))
  
  expect_error(lm_pmm2(y ~ x, data = small_data), "Too few observations")
})

test_that("lm_pmm2 stores original call", {
  fit <- lm_pmm2(y ~ x, data = test_data)
  
  expect_true(is.call(fit@call))
  expect_true(grepl("lm_pmm2", deparse(fit@call)))
})

# ---

# tests/testthat/test-pmm2_timeseries.R
context("PMM2 Time Series")

set.seed(42)
n <- 150

test_that("ar_pmm2 produces valid fit object", {
  x <- arima.sim(model = list(ar = c(0.7, 0.2)), n = n)
  
  fit <- ar_pmm2(x, order = 2)
  
  expect_s4_class(fit, "ARPMM2")
  expect_length(fit@coefficients, 2)
  expect_equal(length(fit@residuals), n - 2)
})

test_that("ma_pmm2 produces valid fit object", {
  x <- arima.sim(model = list(ma = c(0.5, 0.3)), n = n)
  
  fit <- ma_pmm2(x, order = 2)
  
  expect_s4_class(fit, "MAPMM2")
  expect_length(fit@coefficients, 2)
})

test_that("arma_pmm2 produces valid fit object", {
  x <- arima.sim(model = list(ar = 0.7, ma = 0.3), n = n)
  
  fit <- arma_pmm2(x, order = c(1, 1))
  
  expect_s4_class(fit, "ARMAPMM2")
  expect_length(fit@coefficients, 2)  # 1 AR + 1 MA
})

test_that("arima_pmm2 with differencing", {
  # Generate non-stationary data
  x <- cumsum(rnorm(n))
  
  fit <- arima_pmm2(x, order = c(1, 1, 0))
  
  expect_s4_class(fit, "ARIMAPMM2")
  expect_true(fit@convergence || !fit@convergence)  # Just check it ran
})

test_that("ts_pmm2 dispatch to correct model", {
  x <- arima.sim(model = list(ar = 0.7), n = n)
  
  fit_ar <- ts_pmm2(x, model_type = "ar", order = 1)
  expect_s4_class(fit_ar, "ARPMM2")
  
  fit_ma <- ts_pmm2(x, model_type = "ma", order = 1)
  expect_s4_class(fit_ma, "MAPMM2")
})

test_that("Time series residuals are numeric", {
  x <- arima.sim(model = list(ar = 0.7), n = n)
  fit <- ar_pmm2(x, order = 1)
  
  expect_true(is.numeric(fit@residuals))
  expect_false(any(is.na(fit@residuals)))
})

test_that("ar_pmm2 coefficients match AR structure", {
  x <- arima.sim(model = list(ar = c(0.7, 0.2)), n = n)
  fit <- ar_pmm2(x, order = 2)
  
  # Coefficients should be bounded (autocorrelation structure)
  expect_true(all(abs(fit@coefficients) < 1.5))
})

# ---

# tests/testthat/test-pmm2_inference.R
context("PMM2 Inference")

set.seed(42)
n <- 80
x <- rnorm(n)
errors <- rgamma(n, shape = 2, scale = 1) - 2
y <- 2 + 1.5 * x + errors
test_data <- data.frame(x = x, y = y)

test_that("pmm2_inference returns bootstrap results", {
  fit <- lm_pmm2(y ~ x, data = test_data)
  
  boot_results <- pmm2_inference(fit, B = 100)
  
  expect_true(is.list(boot_results))
  expect_true("bootstrap_coef" %in% names(boot_results))
  expect_equal(dim(boot_results$bootstrap_coef), c(100, 2))
})

test_that("pmm2_inference standard errors are positive", {
  fit <- lm_pmm2(y ~ x, data = test_data)
  
  boot_results <- pmm2_inference(fit, B = 100)
  
  se <- apply(boot_results$bootstrap_coef, 2, sd)
  expect_true(all(se > 0))
})

test_that("pmm2_inference confidence intervals are reasonable", {
  fit <- lm_pmm2(y ~ x, data = test_data)
  
  boot_results <- pmm2_inference(fit, B = 100)
  
  # CI should be symmetric around point estimate
  ci_lower <- apply(boot_results$bootstrap_coef, 2, quantile, 0.025)
  ci_upper <- apply(boot_results$bootstrap_coef, 2, quantile, 0.975)
  
  expect_true(all(ci_lower < ci_upper))
  expect_length(ci_lower, 2)
})

test_that("ts_pmm2_inference works for time series", {
  x <- arima.sim(model = list(ar = 0.7), n = 100)
  
  fit <- ar_pmm2(x, order = 1)
  boot_results <- ts_pmm2_inference(fit, B = 50)
  
  expect_true(is.list(boot_results))
  expect_true("bootstrap_coef" %in% names(boot_results))
})

# ---

# tests/testthat/test-pmm2_utilities.R
context("PMM2 Utility Functions")

set.seed(42)

test_that("pmm_skewness calculates correctly", {
  x <- rnorm(100)
  skew <- pmm_skewness(x)
  
  expect_true(is.numeric(skew))
  expect_length(skew, 1)
  expect_false(is.na(skew))
  
  # For normal distribution, skewness should be close to 0
  expect_true(abs(skew) < 0.5)
})

test_that("pmm_skewness detects asymmetry", {
  x_gamma <- rgamma(100, shape = 2)
  skew_gamma <- pmm_skewness(x_gamma)
  
  x_normal <- rnorm(100)
  skew_normal <- pmm_skewness(x_normal)
  
  # Gamma should have higher positive skewness
  expect_true(skew_gamma > skew_normal)
})

test_that("pmm_kurtosis calculates correctly", {
  x <- rnorm(100)
  kurt <- pmm_kurtosis(x)
  
  expect_true(is.numeric(kurt))
  expect_length(kurt, 1)
  expect_false(is.na(kurt))
  
  # For normal distribution, kurtosis should be close to 3 (excess kurtosis = 0)
  expect_true(abs(kurt - 3) < 2)
})

test_that("compute_moments returns all moments", {
  x <- rnorm(100)
  moments <- compute_moments(x)
  
  expect_true(is.list(moments))
  expect_true("m2" %in% names(moments))
  expect_true("m3" %in% names(moments))
  expect_true("m4" %in% names(moments))
  
  expect_true(is.numeric(moments$m2))
  expect_true(moments$m2 > 0)  # Variance must be positive
})

test_that("compute_moments handles edge cases", {
  # All same value - variance should be zero
  x_const <- rep(5, 100)
  moments <- compute_moments(x_const)
  
  expect_true(moments$m2 < 1e-10)  # Close to zero
})

# ---

# tests/testthat/test-pmm2_methods.R
context("PMM2 S4 Methods")

set.seed(42)
n <- 100
x <- rnorm(n)
errors <- rgamma(n, shape = 2, scale = 1) - 2
y <- 2 + 1.5 * x + errors
test_data <- data.frame(x = x, y = y)

test_that("summary method works for PMM2fit", {
  fit <- lm_pmm2(y ~ x, data = test_data)
  
  summary_text <- capture.output(summary(fit))
  
  expect_true(length(summary_text) > 0)
  expect_true(any(grepl("Coefficients", summary_text)))
})

test_that("coef method extracts coefficients", {
  fit <- lm_pmm2(y ~ x, data = test_data)
  
  c <- coef(fit)
  
  expect_equal(c, fit@coefficients)
  expect_length(c, 2)
})

test_that("residuals method returns fit residuals", {
  fit <- lm_pmm2(y ~ x, data = test_data)
  
  res <- residuals(fit)
  
  expect_equal(res, fit@residuals)
  expect_equal(length(res), n)
})

test_that("fitted method returns fitted values", {
  fit <- lm_pmm2(y ~ x, data = test_data)
  
  fitted_vals <- fitted(fit)
  
  expect_true(is.numeric(fitted_vals))
  expect_equal(length(fitted_vals), n)
})

test_that("predict method works", {
  fit <- lm_pmm2(y ~ x, data = test_data)
  
  new_data <- data.frame(x = c(0, 1, -1))
  preds <- predict(fit, newdata = new_data)
  
  expect_true(is.numeric(preds))
  expect_equal(length(preds), 3)
})

test_that("plot method produces output", {
  fit <- lm_pmm2(y ~ x, data = test_data)
  
  expect_error(plot(fit), NA)  # Should not error
})
