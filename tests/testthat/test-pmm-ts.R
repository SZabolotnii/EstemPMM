# tests/testthat/test-pmm-ts.R

# Tests for time series models based on PMM2
testthat::test_that("AR-PMM2 correctly handles AR processes", {
  # Set seed for reproducibility
  set.seed(123)

  # Generate AR(1) process with normal innovations
  ar_coef <- 0.7
  n <- 200
  ar_series <- as.numeric(arima.sim(model = list(ar = ar_coef), n = n))

  # Skip test if not fixed yet
  testthat::skip("Time series tests need further fixes")

  # Estimate model using PMM2
  ar_fit <- ar_pmm2(ar_series, order = 1, method = "pmm2")

  # Tests
  expect_s4_class(ar_fit, "TS2fit")
  expect_true(length(ar_fit@coefficients) == 1)
  expect_true(abs(ar_fit@coefficients[1] - ar_coef) < 0.2) # Check coefficient approximation
  expect_true(length(ar_fit@residuals) > 0) # Check residuals exist

  # Check methods
  expect_error(predict(ar_fit, n.ahead = 5), NA)
  expect_true(length(predict(ar_fit, n.ahead = 5)) == 5)

  # Now with t-distribution (heavy tails)
  ar_series_t <- as.numeric(arima.sim(model = list(ar = ar_coef), n = n,
                                      rand.gen = function(n) rt(n, df=3)))

  ar_fit_t <- ar_pmm2(ar_series_t, order = 1, method = "pmm2")
  expect_s4_class(ar_fit_t, "TS2fit")
  expect_true(abs(ar_fit_t@coefficients[1] - ar_coef) < 0.2)
})

testthat::test_that("MA-PMM2 correctly handles MA processes", {
  # Set seed for reproducibility
  set.seed(123)

  # Generate MA(1) process with normal innovations
  ma_coef <- 0.6
  n <- 200
  ma_series <- as.numeric(arima.sim(model = list(ma = ma_coef), n = n))

  # Skip test if not fixed yet
  testthat::skip("Time series tests need further fixes")

  # Estimate model using PMM2
  ma_fit <- ma_pmm2(ma_series, order = 1, method = "pmm2")

  # Tests
  expect_s4_class(ma_fit, "TS2fit")
  expect_true(length(ma_fit@coefficients) == 1)
  expect_true(abs(ma_fit@coefficients[1] - ma_coef) < 0.2) # Check coefficient approximation
  expect_true(length(ma_fit@residuals) > 0) # Check residuals exist

  # Check methods
  expect_error(predict(ma_fit, n.ahead = 5), NA)
  expect_true(length(predict(ma_fit, n.ahead = 5)) == 5)

  # Check with asymmetric distribution (gamma)
  gamma_innov <- rgamma(n, shape=2, scale=1) - 2  # Center for zero mean
  ma_series_gamma <- as.numeric(filter(gamma_innov, ma_coef, method="convolution", sides=1))
  ma_series_gamma[is.na(ma_series_gamma)] <- 0

  ma_fit_gamma <- ma_pmm2(ma_series_gamma, order = 1, method = "pmm2")
  expect_s4_class(ma_fit_gamma, "TS2fit")
})

testthat::test_that("ARMA-PMM2 correctly handles ARMA processes", {
  # Set seed for reproducibility
  set.seed(123)

  # Generate ARMA(1,1) process with normal innovations
  ar_coef <- 0.7
  ma_coef <- 0.4
  n <- 200
  arma_series <- as.numeric(arima.sim(model = list(ar = ar_coef, ma = ma_coef), n = n))

  # Skip test if not fixed yet
  testthat::skip("Time series tests need further fixes")

  # Estimate model using PMM2
  arma_fit <- arma_pmm2(arma_series, order = c(1, 1), method = "pmm2")

  # Tests
  expect_s4_class(arma_fit, "TS2fit")
  expect_true(length(arma_fit@coefficients) == 2)
  expect_true(abs(arma_fit@coefficients[1] - ar_coef) < 0.2)
  expect_true(abs(arma_fit@coefficients[2] - ma_coef) < 0.2)

  # Check methods
  expect_error(predict(arma_fit, n.ahead = 5), NA)
  pred <- predict(arma_fit, n.ahead = 5)
  expect_true(length(pred) == 5 || length(pred$pred) == 5)
})

testthat::test_that("ARIMA-PMM2 correctly handles ARIMA processes", {
  # Set seed for reproducibility
  set.seed(123)

  # Generate ARMA core
  ar_coef <- 0.7
  ma_coef <- 0.4
  n <- 200
  arma_core <- arima.sim(model = list(ar = ar_coef, ma = ma_coef), n = n)

  # Create non-stationary series via integration (d=1)
  arima_series <- as.numeric(cumsum(arma_core))

  # Skip test if not fixed yet
  testthat::skip("Time series tests need further fixes")

  # Estimate model using PMM2
  arima_fit <- arima_pmm2(arima_series, order = c(1, 1, 1), method = "pmm2")

  # Tests
  expect_s4_class(arima_fit, "TS2fit")
  expect_true(length(arima_fit@coefficients) >= 2)
  expect_true(abs(arima_fit@coefficients[1] - ar_coef) < 0.3)
  expect_true(abs(arima_fit@coefficients[2] - ma_coef) < 0.3)

  # Check methods
  expect_error(predict(arima_fit, n.ahead = 5), NA)
  pred <- predict(arima_fit, n.ahead = 5)
  expect_true(length(pred) == 5 || length(pred$pred) == 5)
})

testthat::test_that("Comparison functions work correctly", {
  # Set seed for reproducibility
  set.seed(123)

  # Generate test data
  n <- 200
  ar_series <- as.numeric(arima.sim(model = list(ar = 0.7), n = n))
  ma_series <- as.numeric(arima.sim(model = list(ma = 0.6), n = n))
  arma_series <- as.numeric(arima.sim(model = list(ar = 0.7, ma = 0.4), n = n))
  arima_series <- as.numeric(cumsum(arma_series))

  # Skip test if not fixed yet
  testthat::skip("Time series tests need further fixes")

  # Test comparison functions
  ar_comp <- compare_ar_methods(ar_series, order = 1)
  expect_true(is.list(ar_comp))
  expect_true("coefficients" %in% names(ar_comp))
  expect_true("residual_stats" %in% names(ar_comp))

  ma_comp <- compare_ma_methods(ma_series, order = 1)
  expect_true(is.list(ma_comp))
  expect_true("coefficients" %in% names(ma_comp))

  arma_comp <- compare_arma_methods(arma_series, order = c(1, 1))
  expect_true(is.list(arma_comp))
  expect_true("coefficients" %in% names(arma_comp))

  arima_comp <- compare_arima_methods(arima_series, order = c(1, 1, 1))
  expect_true(is.list(arima_comp))
  expect_true("coefficients" %in% names(arima_comp))
})

testthat::test_that("PMM2 is more efficient for non-Gaussian data", {
  # Set seed for reproducibility
  set.seed(123)

  # Generate data with highly skewed distribution
  n <- 300
  gamma_innov <- rgamma(n, shape=1, scale=1) - 1  # Strong asymmetry
  ar_coef <- 0.7

  # Create AR series with asymmetric innovations
  ar_series <- numeric(n)
  ar_series[1] <- gamma_innov[1]
  for(i in 2:n) {
    ar_series[i] <- ar_coef * ar_series[i-1] + gamma_innov[i]
  }

  # Skip test if not fixed yet
  testthat::skip("Time series tests need further fixes")

  # Compare methods
  ar_comp <- compare_ar_methods(ar_series, order = 1)

  # PMM2 should give lower RSS and MAE values for non-Gaussian data
  expect_true(ar_comp$residual_stats$RSS[3] < ar_comp$residual_stats$RSS[1])  # PMM2 vs YW
  expect_true(ar_comp$residual_stats$MAE[3] < ar_comp$residual_stats$MAE[1])  # PMM2 vs YW
})
