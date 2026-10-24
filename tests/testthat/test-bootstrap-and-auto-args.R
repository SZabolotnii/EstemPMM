# Residual bootstrap of ARMA/ARIMA fits and argument forwarding of
# pmm_lm(method = "auto").

test_that("residual bootstrap runs for ARMA and ARIMA fits", {
  # arima.sim() needs a burn-in of at least p + q; max(p, q) failed here
  set.seed(1)
  x <- arima.sim(list(ar = 0.5, ma = 0.3), n = 200, innov = rexp(200) - 1)
  fit <- arma_pmm2(x, order = c(1, 1), ma_solver = "recursive")
  inf <- ts_pmm2_inference(fit, method = "residual", B = 20, seed = 3)
  expect_equal(nrow(inf), 2L)
  expect_true(all(is.finite(inf$Std.Error)))

  xi <- cumsum(x)
  fit_i <- arima_pmm2(xi, order = c(1, 1, 1), ma_solver = "recursive")
  inf_i <- ts_pmm2_inference(fit_i, method = "residual", B = 20, seed = 3)
  expect_true(all(is.finite(inf_i$Std.Error)))
})

test_that("pmm_lm(auto) drops arguments the selected fitter does not take", {
  set.seed(2)
  d <- data.frame(x = rnorm(200))
  d$y <- 1 + d$x + runif(200, -1, 1)
  expect_message(fit <- pmm_lm(y ~ x, d, reg_lambda = 1e-6),
                 "ignoring argument\\(s\\) it does not take: reg_lambda")
  expect_identical(attr(fit, "dispatch")$method, "PMM3")
  expect_s4_class(fit, "PMM3fit")
  # An explicit method still rejects a foreign argument
  expect_error(pmm_lm(y ~ x, d, method = "pmm3", reg_lambda = 1e-6),
               "unused argument")
})
