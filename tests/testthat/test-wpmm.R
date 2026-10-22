test_that("lm_wpmm2/lm_wpmm3 reproduce the gmdhpmm weak-moment estimators", {
  ref <- readRDS(test_path("fixtures", "reference_gmdhpmm.rds"))
  for (nm in names(ref$cases)) {
    cs <- ref$cases[[nm]]
    f2 <- suppressWarnings(lm_wpmm2(ref$formula, cs$data))
    f3 <- suppressWarnings(lm_wpmm3(ref$formula, cs$data))
    # gmdhpmm adds a 1e-8 ridge penalty; the observed gap is below 2e-9
    expect_equal(unname(coef(f2)), unname(cs$wpmm2), tolerance = 1e-7, info = nm)
    expect_equal(unname(coef(f3)), unname(cs$wpmm3), tolerance = 1e-7, info = nm)
  }
})

make_wpmm_data <- function(n, err, seed) {
  set.seed(seed)
  x <- rnorm(n)
  data.frame(x = x, y = 1 + 2 * x + err(n))
}

test_that("WPMMfit object has the expected structure and methods", {
  d <- make_wpmm_data(150, function(n) rt(n, 3), 1)
  fit <- lm_wpmm3(y ~ x, data = d)
  expect_s4_class(fit, "WPMMfit")
  expect_true(is(fit, "PMMfit"))
  expect_true(fit@convergence)
  expect_equal(fit@degree, 3)
  expect_named(coef(fit), c("(Intercept)", "x"))
  expect_equal(fitted(fit) + residuals(fit), d$y)
  expect_equal(predict(fit), fitted(fit))
  expect_equal(predict(fit, newdata = data.frame(x = c(0, 1))),
               unname(coef(fit)[1] + c(0, 1) * coef(fit)[2]))
  expect_equal(nobs(fit), 150L)
  expect_named(fit@weak_moments, c("m2", "m4", "m6"))
  expect_true(all(fit@weights >= 0))

  V <- vcov(fit)
  expect_equal(dim(V), c(2L, 2L))
  expect_equal(dimnames(V), list(names(coef(fit)), names(coef(fit))))
  expect_true(all(diag(V) > 0))
  ci <- confint(fit, level = 0.9)
  expect_equal(colnames(ci), c("5 %", "95 %"))
  expect_true(all(ci[, 1] < coef(fit) & coef(fit) < ci[, 2]))
  expect_output(summary(fit), "sandwich standard errors")
  expect_output(show(fit), "WPMM3 algorithm: converged")

  fit2 <- lm_wpmm2(y ~ x, data = d)
  expect_equal(fit2@degree, 2)
  expect_named(fit2@weak_moments, c("c2", "c3", "c4"))
  expect_output(summary(fit2), "a2 =")
})

test_that("weak estimators resist gross outliers that break OLS", {
  d <- make_wpmm_data(200, function(n) rnorm(n), 2)
  d$y[1:10] <- d$y[1:10] + 50 * sign(d$x[1:10])
  b_ols <- coef(lm(y ~ x, data = d))[2]
  expect_gt(abs(b_ols - 2), 1)
  expect_lt(abs(coef(lm_wpmm2(y ~ x, data = d))[2] - 2), 0.3)
  expect_lt(abs(coef(lm_wpmm3(y ~ x, data = d))[2] - 2), 0.3)
  expect_gte(sum(lm_wpmm3(y ~ x, data = d)@weights == 0), 10)
})

test_that("sandwich bread includes the window slope", {
  # Hand-built A^{-1} B A^{-1} with psi' from a numerical derivative of W(r) r,
  # with the window centre, scale and weak moments held fixed
  d <- make_wpmm_data(300, function(n) rexp(n) - 1, 3)
  for (deg in 2:3) {
    fit <- if (deg == 2) lm_wpmm2(y ~ x, data = d) else lm_wpmm3(y ~ x, data = d)
    r <- residuals(fit)
    X <- cbind(1, d$x)
    ww <- EstemPMM:::.wpmm_weights(r, deg, 2.5)
    mu <- ww$center; sw <- ww$scale
    w <- exp(-0.5 * ((r - mu) / sw)^2)
    m <- sum(w * r) / sum(w)
    psi <- function(e) {
      wv <- exp(-0.5 * ((e - mu) / sw)^2)
      rc <- e - m
      W <- if (deg == 2) wv * (1 + ww$correction * rc) else
        wv * (ww$correction - rc^2) / abs(ww$correction)
      pmax(W, 0) * e
    }
    h <- 1e-6
    dpsi <- (psi(r + h) - psi(r - h)) / (2 * h)
    A <- crossprod(X, X * dpsi)
    B <- crossprod(X * psi(r))
    V <- solve(A) %*% B %*% solve(A)
    expect_equal(unname(vcov(fit)), V, tolerance = 1e-5, info = deg)
  }
})

test_that("input validation and pmm_lm routing", {
  d <- make_wpmm_data(80, function(n) rt(n, 4), 5)
  expect_error(lm_wpmm2(y ~ x, data = as.list(d)), "data.frame")
  expect_error(lm_wpmm3(y ~ x, data = d, sigma_mult = 0), "sigma_mult")
  expect_error(lm_wpmm3(y ~ x, data = d, max_iter = 0), "max_iter")
  expect_warning(lm_wpmm2(y ~ x, data = d, max_iter = 1), "did not converge")

  f <- pmm_lm(y ~ x, data = d, method = "wpmm3", sigma_mult = 3)
  expect_s4_class(f, "WPMMfit")
  expect_equal(f@sigma_mult, 3)
  expect_equal(coef(f), coef(lm_wpmm3(y ~ x, data = d, sigma_mult = 3)))
  expect_equal(deparse(f@call$formula), "y ~ x")
  expect_s4_class(pmm_lm(y ~ x, data = d, method = "wpmm2"), "WPMMfit")
})
