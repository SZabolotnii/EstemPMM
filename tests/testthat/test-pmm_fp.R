# Smoke tests for pmm_fp() scaffold.

set.seed(20260511)

# Hill-type DGP with gamma errors (asymmetric, gamma3 > 0)
n  <- 50
x  <- runif(n, 0.5, 5)
mu <- 1 + 2 * x ^ 0.5
err <- rgamma(n, shape = 2, scale = 1) - 2
y  <- mu + err
dat <- data.frame(x = x, y = y)

test_that("pmm_fp returns a PMM_FP_fit object on track='pos'", {
    suppressWarnings(
        fit <- pmm_fp(y ~ x, data = dat, track = "pos", max_terms = 2)
    )
    expect_s4_class(fit, "PMM_FP_fit")
    expect_s4_class(fit, "PMM2fit")
    expect_true(length(fit@selected_powers) >= 1L)
    expect_equal(fit@track, "pos")
    expect_equal(fit@criterion, "BIC")
})

test_that("pmm_fp accepts track='full'", {
    suppressWarnings(
        fit <- pmm_fp(y ~ x, data = dat, track = "full", max_terms = 2)
    )
    expect_s4_class(fit, "PMM_FP_fit")
    expect_equal(fit@track, "full")
})

test_that("pmm_fp shifts non-positive x and records the shift", {
    dat2 <- dat
    dat2$x <- dat2$x - 1   # introduce non-positive values
    suppressWarnings(
        fit <- pmm_fp(y ~ x, data = dat2, track = "pos", max_terms = 1)
    )
    expect_gt(fit@shift_eps, 0)
})

test_that("pmm_fp validates arguments", {
    expect_error(pmm_fp(y ~ x, data = dat, order = 5),
                 "order")
    expect_error(pmm_fp(y ~ x, data = dat, max_terms = 0),
                 "max_terms")
})
