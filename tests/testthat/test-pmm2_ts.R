set.seed(42)
n <- 150

test_that("ar_pmm2 produces valid fit object", {
  x <- arima.sim(model = list(ar = c(0.7, 0.2)), n = n)

  fit <- ar_pmm2(x, order = 2)

  expect_s4_class(fit, "ARPMM2")
  expect_length(fit@coefficients, 2)
  expect_gt(length(fit@residuals), 0)
})

test_that("ma_pmm2 produces valid fit object", {
  x <- arima.sim(model = list(ma = c(0.5, -0.3)), n = n)

  fit <- ma_pmm2(x, order = 2)

  expect_s4_class(fit, "MAPMM2")
  expect_length(fit@coefficients, 2)
})

test_that("arma_pmm2 produces valid fit object", {
  x <- arima.sim(model = list(ar = 0.7, ma = 0.3), n = n)

  fit <- arma_pmm2(x, order = c(1, 1))

  expect_s4_class(fit, "ARMAPMM2")
  expect_length(fit@coefficients, 2)
})

test_that("arima_pmm2 handles differencing", {
  x <- cumsum(rnorm(n))

  fit <- arima_pmm2(x, order = c(1, 1, 0))

  expect_s4_class(fit, "ARIMAPMM2")
  expect_length(fit@coefficients, 1)
})

test_that("ts_pmm2 dispatches to correct classes", {
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

test_that("ar_pmm2 coefficients are bounded", {
  x <- arima.sim(model = list(ar = c(0.7, 0.2)), n = n)
  fit <- ar_pmm2(x, order = 2)

  expect_true(all(abs(fit@coefficients) < 2))
})

test_that("arima_pmm2 works with real oil price data", {
  # Load real data from DCOILWTICO.csv
  data_path <- system.file("extdata", "DCOILWTICO.csv", package = "EstemPMM")

  # If not found in package, try relative path
  if (!file.exists(data_path) || data_path == "") {
    data_path <- "../../data/DCOILWTICO.csv"
  }

  # Skip test if data file doesn't exist
  skip_if_not(file.exists(data_path), "Data file DCOILWTICO.csv not found")

  # Read the data
  oil_data <- read.csv(data_path, stringsAsFactors = FALSE)

  # Extract price column and remove NA values
  prices <- as.numeric(oil_data$DCOILWTICO)
  prices_clean <- prices[!is.na(prices)]

  # Check we have enough data
  expect_gt(length(prices_clean), 100)

  # Fit ARIMA(1,1,1) model
  fit <- arima_pmm2(prices_clean, order = c(1, 1, 1))

  # Test 1: Check object class
  expect_s4_class(fit, "ARIMAPMM2")

  # Test 2: Check coefficients length (AR + MA = 2)
  expect_length(fit@coefficients, 2)

  # Test 3: Check convergence
  expect_true(fit@convergence)

  # Test 4: Check residuals are numeric and finite
  expect_true(is.numeric(fit@residuals))
  residuals_clean <- fit@residuals[is.finite(fit@residuals)]
  expect_gt(length(residuals_clean), 0)

  # Test 5: Check moments are computed
  expect_true(is.finite(fit@m2))
  expect_true(is.finite(fit@m3))
  expect_true(is.finite(fit@m4))

  # Test 6: Check order is correct
  expect_equal(fit@order$ar, 1)
  expect_equal(fit@order$ma, 1)
  expect_equal(fit@order$d, 1)
})

test_that("arima_pmm2 handles different orders on oil data", {
  # Load data
  data_path <- system.file("extdata", "DCOILWTICO.csv", package = "EstemPMM")
  if (!file.exists(data_path) || data_path == "") {
    data_path <- "../../data/DCOILWTICO.csv"
  }
  skip_if_not(file.exists(data_path), "Data file DCOILWTICO.csv not found")

  oil_data <- read.csv(data_path, stringsAsFactors = FALSE)
  prices <- as.numeric(oil_data$DCOILWTICO)
  prices_clean <- prices[!is.na(prices)]

  # Use subset for faster testing
  prices_subset <- head(prices_clean, 200)

  # Test ARIMA(0,1,1) - IMA model
  fit_011 <- arima_pmm2(prices_subset, order = c(0, 1, 1))
  expect_s4_class(fit_011, "ARIMAPMM2")
  expect_length(fit_011@coefficients, 1)  # Only MA coefficient

  # Test ARIMA(1,1,0) - ARI model
  fit_110 <- arima_pmm2(prices_subset, order = c(1, 1, 0))
  expect_s4_class(fit_110, "ARIMAPMM2")
  expect_length(fit_110@coefficients, 1)  # Only AR coefficient

  # Test ARIMA(2,1,1)
  fit_211 <- arima_pmm2(prices_subset, order = c(2, 1, 1))
  expect_s4_class(fit_211, "ARIMAPMM2")
  expect_length(fit_211@coefficients, 3)  # AR1 + AR2 + MA1
})
