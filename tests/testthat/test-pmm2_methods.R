set.seed(42)
n <- 100
x <- rnorm(n)
errors <- rgamma(n, shape = 2, scale = 1) - 2
y <- 2 + 1.5 * x + errors
test_data <- data.frame(x = x, y = y)

test_that("summary method runs without error", {
  fit <- lm_pmm2(y ~ x, data = test_data)

  expect_error(summary(fit), NA)
})

test_that("summary reports the asymptotic standard errors from vcov", {
  fit <- lm_pmm2(y ~ x, data = test_data)
  expect_false(attr(fit, "fallback"))
  out <- capture.output(summary(fit))
  expect_true(any(grepl("Std. Error", out, fixed = TRUE)))
  se <- sqrt(diag(vcov(fit)))
  vf <- pmm2_variance_factor(fit@m2, fit@m3, fit@m4)
  ols_se <- sqrt(diag(fit@m2 * solve(crossprod(cbind(1, x)))))
  # Slope: sqrt(g2) times the OLS SE. Intercept: plus (1 - g2) m2 / n from
  # the moments estimated on the OLS residuals
  expect_equal(unname(se[2]), unname(sqrt(vf$g) * ols_se[2]))
  expect_equal(unname(se[1]^2),
               unname(vf$g * ols_se[1]^2 + (1 - vf$g) * fit@m2 / n))
})

test_that("coef method extracts coefficients", {
  fit <- lm_pmm2(y ~ x, data = test_data)

  c <- coef(fit)

  expect_equal(c, fit@coefficients)
  expect_length(c, 2)
})

test_that("residuals method returns residual vector", {
  fit <- lm_pmm2(y ~ x, data = test_data)

  res <- residuals(fit)

  expect_true(is.numeric(res))
  expect_equal(length(res), n)
})

test_that("fitted method returns fitted values", {
  fit <- lm_pmm2(y ~ x, data = test_data)

  fitted_vals <- fitted(fit)

  expect_true(is.numeric(fitted_vals))
  expect_equal(length(fitted_vals), n)
})

test_that("predict method supports new data", {
  fit <- lm_pmm2(y ~ x, data = test_data)

  new_data <- data.frame(x = c(0, 1, -1))
  preds <- predict(fit, newdata = new_data)

  expect_true(is.numeric(preds))
  expect_equal(length(preds), 3)
})

test_that("plot method works for PMM2fit", {
  fit <- lm_pmm2(y ~ x, data = test_data)

  expect_error(plot(fit), NA)
})
