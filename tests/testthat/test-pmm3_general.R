test_that("lm_pmm3(symmetric = FALSE) reproduces the AISM revision solver", {
  ref <- readRDS(test_path("fixtures", "reference_aism_pmm3.rds"))
  for (nm in names(ref$cases)) {
    cs <- ref$cases[[nm]]
    fit <- lm_pmm3(y ~ x, data = cs$data, symmetric = FALSE)
    # pmm_reg(outer = 1): moments from the OLS residuals, the same equation
    expect_equal(unname(coef(fit)), cs$outer1, tolerance = 1e-10, info = nm)
    # pmm_reg(outer = 3) re-estimates the moments; adaptive = TRUE iterates
    # that to the fixed point, so the two are close but not identical
    fa <- lm_pmm3(y ~ x, data = cs$data, symmetric = FALSE, adaptive = TRUE)
    expect_lt(max(abs(coef(fa) - cs$outer3)), 1e-2)
  }
})

test_that("the general degree-3 normal system reduces to the symmetric g3", {
  # Uniform(-1, 1): odd moments vanish, so h2 = 0 and g equals the
  # closed-form g3 = 1 - gamma4^2 / (6 + 9 gamma4 + gamma6)
  alpha <- c(0, 1/3, 0, 1/5, 0, 1/7)
  pc <- EstemPMM:::.pmm_poly_coef(alpha, 3)
  expect_equal(pc$h[2], 0)
  vf <- pmm3_variance_factor(alpha[2], alpha[4], alpha[6])
  expect_equal(pc$g, vf$g3)
  # Exponential(1): g2 = 1 - 4/(2 + 6) = 0.5 and g3 = 1/3 (exact cumulants)
  ex <- c(0, 1, 2, 9, 44, 265)
  expect_equal(EstemPMM:::.pmm_poly_coef(ex, 2)$g, 0.5)
  expect_equal(EstemPMM:::.pmm_poly_coef(ex, 3)$g, 1 / 3)
})

test_that("general PMM3 fit object: slots, vcov and summary", {
  set.seed(31)
  n <- 300
  x <- runif(n, 0, 5)
  d <- data.frame(x = x, y = 2 + 1.5 * x + rgamma(n, 2, 1) - 2)
  fit <- lm_pmm3(y ~ x, data = d, symmetric = FALSE)
  expect_s4_class(fit, "PMM3fit")
  expect_false(fit@symmetric)
  expect_length(fit@poly_coef, 3)
  expect_true(is.na(fit@kappa))
  expect_true(fit@convergence)
  expect_lt(fit@g_coefficient, 1)

  X <- cbind(1, x)
  g <- fit@g_coefficient
  a <- solve(crossprod(X), colSums(X))
  expect_equal(unname(vcov(fit)),
               unname(g * fit@m2 * solve(crossprod(X)) +
                        (1 - g) * fit@m2 / nrow(X) * tcrossprod(a)))
  expect_output(summary(fit), "general polynomial")

  sym <- suppressWarnings(lm_pmm3(y ~ x, data = d))
  expect_true(sym@symmetric)
  expect_length(sym@poly_coef, 0)

  # pmm_lm forwards the argument
  f2 <- pmm_lm(y ~ x, data = d, method = "pmm3", symmetric = FALSE)
  expect_equal(coef(f2), coef(fit))
})

test_that("symmetric PMM3 warning points to the general estimator", {
  set.seed(32)
  d <- data.frame(x = rnorm(200))
  d$y <- 1 + d$x + rexp(200) - 1
  expect_warning(lm_pmm3(y ~ x, data = d), "symmetric = FALSE")
  expect_no_warning(lm_pmm3(y ~ x, data = d, symmetric = FALSE))
})

test_that("general PMM3 falls back to OLS when the iteration fails", {
  d <- data.frame(x = 1:6, y = c(0, 0, 0, 0, 0, 10))
  fit <- suppressWarnings(lm_pmm3(y ~ x, data = d, symmetric = FALSE))
  b_ols <- unname(coef(lm(y ~ x, data = d)))
  if (!fit@convergence) {
    expect_equal(unname(coef(fit)), b_ols)
    expect_true(attr(fit, "fallback"))
    X <- cbind(1, d$x)
    expect_equal(unname(vcov(fit)), unname(fit@m2 * solve(crossprod(X))))
  } else {
    expect_true(all(is.finite(coef(fit))))
  }
})
