test_that("estpmm_style_ma works for MA(1) models", {
    set.seed(123)
    n <- 200
    theta <- 0.6

    # Generate MA(1) with exponential errors (asymmetric)
    innovations <- rexp(n, rate = 1) - 1
    x <- arima.sim(n = n, list(ma = theta), innov = innovations)

    # Fit using new PMM2 method
    fit <- EstemPMM:::estpmm_style_ma(x, q = 1, include.mean = FALSE, verbose = FALSE)

    expect_type(fit, "list")
    expect_true(fit$convergence)
    expect_length(fit$ma_coef, 1)
    expect_equal(length(fit$innovations), n)

    # Check that estimates are reasonable
    expect_equal(fit$ma_coef[[1]], theta, tolerance = 0.2)

    # Check method string
    expect_equal(fit$method, "EstemPMM-style PMM2")
})

test_that("estpmm_style_sma works for SMA(1) models", {
    set.seed(456)
    n <- 200
    Theta <- 0.6
    s <- 4

    # Generate SMA(1)_4
    innovations <- rexp(n, rate = 1) - 1
    # Manual generation for SMA
    x <- numeric(n)
    for (t in 1:n) {
        sma_term <- if (t > s) Theta * innovations[t - s] else 0
        x[t] <- innovations[t] + sma_term
    }

    # Fit using new PMM2 method
    fit <- EstemPMM:::estpmm_style_sma(x, Q = 1, s = s, include.mean = FALSE, verbose = FALSE)

    expect_type(fit, "list")
    expect_true(fit$convergence)
    expect_length(fit$sma_coef, 1)

    # Check estimates
    expect_equal(fit$sma_coef[[1]], Theta, tolerance = 0.2)
})

test_that("sarima_pmm2 integrates ma_method='pmm2' correctly", {
    set.seed(789)
    n <- 200

    # MA(1) case
    innovations <- rexp(n, rate = 1) - 1
    x <- arima.sim(n = n, list(ma = 0.5), innov = innovations)

    # Fit with ma_method = "pmm2"
    fit <- sarima_pmm2(x,
        order = c(0, 0, 1, 0), seasonal = list(order = c(0, 0), period = 1),
        ma_method = "pmm2", include.mean = FALSE
    )

    expect_s4_class(fit, "SARIMAPMM2")
    expect_true(fit@convergence)
    expect_equal(length(fit@coefficients), 1)

    # Fit with ma_method = "mle" (default)
    fit_mle <- sarima_pmm2(x,
        order = c(0, 0, 1, 0), seasonal = list(order = c(0, 0), period = 1),
        ma_method = "mle", include.mean = FALSE
    )

    expect_s4_class(fit_mle, "SARIMAPMM2")

    # Ensure results are slightly different (since methods differ)
    # Note: For normal data they might be close, but for exp data they should differ
    expect_false(isTRUE(all.equal(fit@coefficients, fit_mle@coefficients)))
})
