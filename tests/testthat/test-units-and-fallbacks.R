# The linear-model PMM fits must not depend on how the regressors and the
# response are measured, and must say so when they return the OLS values.

test_that("PMM2 converges with an uncentred regressor (calendar years)", {
  set.seed(4)
  yr <- 2000:2019
  y <- 1 + 0.05 * (yr - 2000) + rexp(20) - 1
  d <- data.frame(y = y, yr = yr, t = yr - 2000)
  raw <- lm_pmm2(y ~ yr, data = d)
  cen <- lm_pmm2(y ~ t, data = d)
  expect_true(raw@convergence)
  expect_false(attr(raw, "fallback"))
  expect_equal(unname(coef(raw)[2]), unname(coef(cen)[2]), tolerance = 1e-6)
  expect_equal(unname(coef(raw)[1] + 2000 * coef(raw)[2]),
               unname(coef(cen)[1]), tolerance = 1e-6)
})

test_that("general PMM3 converges with an uncentred regressor", {
  # The plain score norm failed here (28 of 200 such samples)
  set.seed(18)
  yr <- 2000:2029
  y <- 1 + 0.05 * (yr - 2000) + rgamma(30, 2) - 2
  d <- data.frame(y = y, yr = yr, t = yr - 2000)
  raw <- lm_pmm3(y ~ yr, data = d, symmetric = FALSE)
  cen <- lm_pmm3(y ~ t, data = d, symmetric = FALSE)
  expect_true(raw@convergence)
  expect_equal(unname(coef(raw)[2]), unname(coef(cen)[2]), tolerance = 1e-6)
})

test_that("PMM2 is unaffected by a large offset in the response", {
  set.seed(11)
  x <- rnorm(80)
  y <- 2 + x + rexp(80) - 1
  base <- lm_pmm2(y ~ x, data = data.frame(x, y))
  big <- lm_pmm2(y ~ x, data = data.frame(x, y = y + 1e7))
  expect_true(big@convergence)
  expect_false(attr(big, "fallback"))
  expect_equal(unname(coef(big)[2]), unname(coef(base)[2]), tolerance = 1e-6)
})

test_that("PMM3 fits are equivariant to the units of the response", {
  set.seed(12)
  x <- rnorm(200)
  e_skew <- rgamma(200, 2) - 2
  e_flat <- runif(200, -1, 1)
  # The default tol is measured in residual-SD units, so it needs no scaling
  # (an absolute 1e-6 stopped after one step at k = 1e-5)
  for (k in c(1e-5, 1e4)) {
    g1 <- lm_pmm3(y ~ x, data = data.frame(x, y = 1 + x + e_skew), symmetric = FALSE)
    gk <- lm_pmm3(y ~ x, data = data.frame(x, y = k * (1 + x + e_skew)),
                  symmetric = FALSE)
    expect_false(attr(gk, "fallback"))
    expect_equal(unname(coef(gk)), k * unname(coef(g1)), tolerance = 1e-6)
    expect_equal(gk@g_coefficient, g1@g_coefficient, tolerance = 1e-8)

    s1 <- lm_pmm3(y ~ x, data = data.frame(x, y = 1 + x + e_flat))
    sk <- lm_pmm3(y ~ x, data = data.frame(x, y = k * (1 + x + e_flat)))
    expect_false(attr(sk, "fallback"))
    expect_equal(unname(coef(sk)), k * unname(coef(s1)), tolerance = 1e-6)
  }
})

test_that("PMM2 fits are equivariant to the units of the response", {
  set.seed(13)
  x <- rnorm(200)
  y <- 1 + x + rexp(200) - 1
  p1 <- lm_pmm2(y ~ x, data = data.frame(x, y))
  for (k in c(1e-5, 1e4)) {
    pk <- lm_pmm2(y ~ x, data = data.frame(x, y = k * y))
    expect_true(pk@convergence)
    expect_equal(unname(coef(pk)), k * unname(coef(p1)), tolerance = 1e-6)
  }
})

test_that("standardised normal system gives the raw-moment solution", {
  alpha <- c(0, 2, 1.5, 14, 25, 160)
  i <- 1:3
  a <- function(k) if (k == 0) 1 else alpha[k]
  F <- outer(i, i, Vectorize(function(p, q) a(p + q) - a(p) * a(q)))
  b <- i * vapply(i - 1, a, numeric(1))
  pc <- EstemPMM:::.pmm_poly_coef(alpha, 3)
  expect_equal(pc$h, as.numeric(solve(F, b)), tolerance = 1e-12)
  expect_equal(pc$g, 1 / (alpha[2] * sum(b * solve(F, b))), tolerance = 1e-12)
})

test_that("near-Gaussian lm_pmm3 is flagged as an OLS fall-back", {
  set.seed(13)
  d <- data.frame(x = rnorm(60))
  d$y <- 1 + d$x + runif(60, -1, 1)
  orig <- compute_moments_pmm3
  testthat::local_mocked_bindings(compute_moments_pmm3 = function(residuals) {
    m <- orig(residuals)
    m$kappa <- NA_real_
    m
  })
  fit <- lm_pmm3(y ~ x, data = d)
  expect_true(attr(fit, "fallback"))
  X <- attr(fit, "model_matrix")
  expect_equal(unname(vcov(fit)), unname(fit@m2 * solve(crossprod(X))))
})

test_that("a singular Jacobian falls back to the start instead of erroring", {
  X <- cbind(1, 1:10)
  y <- as.numeric(1:10)
  expect_warning(
    res <- EstemPMM:::pmm2_algorithm(c(0, 1), X, y, m2 = 1, m3 = 1, m4 = NaN),
    "singular Jacobian")
  expect_true(res$fallback)
  expect_false(res$convergence)
  expect_equal(res$b, c(0, 1))
})

test_that("intercept covariance includes the plug-in moment term", {
  set.seed(14)
  x <- rnorm(150)
  d <- data.frame(x = x, y = 1 + x + rexp(150) - 1)
  fit <- lm_pmm2(y ~ x, data = d)
  X <- attr(fit, "model_matrix")
  g <- pmm2_variance_factor(fit@m2, fit@m3, fit@m4)$g
  a <- solve(crossprod(X), colSums(X))
  V <- g * fit@m2 * solve(crossprod(X)) +
    (1 - g) * fit@m2 / nrow(X) * tcrossprod(a)
  expect_equal(unname(vcov(fit)), unname(V))
  # a is the intercept's unit vector, so the slope block is unchanged
  expect_equal(unname(vcov(fit)[2, 2]),
               unname(g * fit@m2 * solve(crossprod(X))[2, 2]))
})

test_that("predict(WPMMfit) codes factors with the training levels", {
  set.seed(15)
  d <- data.frame(x = rnorm(90), f = factor(rep(c("a", "b", "c"), 30)))
  d$y <- 1 + d$x + c(a = 0, b = 2, c = -1)[as.character(d$f)] + rexp(90) - 1
  fit <- lm_wpmm2(y ~ x + f, data = d)
  cf <- coef(fit)
  nd <- data.frame(x = c(0.5, -1, 2), f = c("c", "c", "b"))
  expected <- cf[1] + cf["x"] * nd$x + ifelse(nd$f == "b", cf["fb"], cf["fc"])
  expect_equal(predict(fit, newdata = nd), unname(expected))
  nd$f <- factor(nd$f)
  expect_equal(predict(fit, newdata = nd), unname(expected))
  nd$x[2] <- NA
  p <- predict(fit, newdata = nd)
  expect_length(p, 3)
  expect_true(is.na(p[2]))
})
