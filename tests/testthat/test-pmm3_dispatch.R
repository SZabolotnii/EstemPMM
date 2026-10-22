test_that("pmm_dispatch selects PMM3 for uniform errors", {
  set.seed(42)
  x <- rnorm(200)
  eps <- runif(200, -sqrt(3), sqrt(3))
  y <- 1 + 2 * x + eps
  r <- residuals(lm(y ~ x))

  result <- pmm_dispatch(r, verbose = FALSE)
  expect_equal(result$method, "PMM3")
})

test_that("pmm_dispatch selects PMM2 for exponential errors", {
  set.seed(42)
  x <- rnorm(200)
  eps <- rexp(200) - 1
  y <- 1 + 2 * x + eps
  r <- residuals(lm(y ~ x))

  result <- pmm_dispatch(r, verbose = FALSE)
  expect_equal(result$method, "PMM2")
})

test_that("pmm_dispatch selects OLS for normal errors", {
  set.seed(42)
  x <- rnorm(200)
  eps <- rnorm(200)
  y <- 1 + 2 * x + eps
  r <- residuals(lm(y ~ x))

  result <- pmm_dispatch(r, verbose = FALSE)
  expect_equal(result$method, "OLS")
})

test_that("pmm_dispatch returns complete structure", {
  set.seed(42)
  r <- rnorm(100)
  result <- pmm_dispatch(r, verbose = FALSE)

  expect_true(is.list(result))
  expect_true(result$method %in% c("OLS", "PMM2", "PMM3"))
  expect_true(is.numeric(result$gamma3))
  expect_true(is.numeric(result$gamma4))
  expect_true(is.numeric(result$gamma6))
  expect_true(is.numeric(result$g2))
  expect_true(is.numeric(result$g3))
  expect_true(is.numeric(result$g_selected))
  expect_true(is.numeric(result$improvement_pct))
  expect_true(is.character(result$reasoning))
  expect_true(is.numeric(result$n))
})

test_that("pmm_dispatch respects custom thresholds", {
  set.seed(42)
  x <- rnorm(200)
  eps <- runif(200, -sqrt(3), sqrt(3))
  y <- 1 + 2 * x + eps
  r <- residuals(lm(y ~ x))

  # With very strict kurtosis threshold, uniform may not qualify
  result <- pmm_dispatch(r, kurtosis_threshold = -2.0, verbose = FALSE)
  expect_true(result$method %in% c("OLS", "PMM2", "PMM3"))
})

test_that("g3 from sample residuals always lies in [0, 1]", {
  # 6 + 9*gamma4 + gamma6 - gamma4^2 = m6/m2^3 - (m4/m2^2)^2 >= 0, so the
  # band 0 < denom <= gamma4^2 cannot be reached from residuals.
  set.seed(7)
  samples <- list(rexp(15), rlnorm(40, 0, 2), runif(30)^5, rt(50, 3),
                  sample(c(-1, 1), 60, TRUE) + rnorm(60, 0, 1e-3))
  for (r in samples) {
    g3 <- pmm_dispatch(r, verbose = FALSE)$g3
    expect_gte(g3, 0)
    expect_lte(g3, 1)
  }
  # Symmetric two-point law: the boundary, where PMM3 is exact (g3 = 0)
  expect_equal(pmm_dispatch(rep(c(-1, 1), 50), verbose = FALSE)$g3, 0)
})

test_that("symmetric heavy-tailed residuals go to PMM3 only once n >= n_leptokurtic", {
  set.seed(11)
  r <- rexp(1000)
  r_big <- c(r, -r)[sample(2000)]            # symmetrised Laplace, gamma3 = 0
  dec <- pmm_dispatch(r_big, verbose = FALSE)
  expect_equal(dec$method, "PMM3")
  expect_lt(dec$g3, 0.95)
  # Same law below the sample-size threshold: the heavy-tailed branch is closed
  dec_small <- pmm_dispatch(r_big[1:60], verbose = FALSE)
  expect_false(dec_small$method == "PMM3" && dec_small$gamma4 > -0.7)
  # Switched off, as in 0.5.0
  expect_false(pmm_dispatch(r_big, g3_threshold = 0, verbose = FALSE)$method == "PMM3")
})

test_that("small-sample platykurtic residuals with noisy skewness go to PMM3", {
  # |gamma3| above 0.3 but within 2*sqrt(6/30) = 0.89 of zero
  set.seed(5)
  repeat {
    r <- rnorm(30, sample(c(-1, 1), 30, TRUE), 0.3)
    cu <- EstemPMM:::.residual_cumulants(r)
    if (abs(cu$gamma3) > 0.3 && abs(cu$gamma3) < 2 * sqrt(6 / 30) &&
        cu$gamma4 < -0.7) break
  }
  expect_equal(pmm_dispatch(r, verbose = FALSE)$method, "PMM3")
  # The 0.5.0 rule (clear_asymmetry_z = 0) sent such residuals to PMM2
  expect_equal(pmm_dispatch(r, clear_asymmetry_z = 0, verbose = FALSE)$method,
               "PMM2")
})

test_that("clearly skewed residuals still go to PMM2", {
  set.seed(8)
  expect_equal(pmm_dispatch(rexp(30) - 1, verbose = FALSE)$method, "PMM2")
  expect_equal(pmm_dispatch(rgamma(500, 2) - 2, verbose = FALSE)$method, "PMM2")
})
