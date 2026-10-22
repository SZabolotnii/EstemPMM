# vcov() of time-series fits whose coefficients are the classical start rather
# than a PMM root: the PMM variance factor (g2 or g3) must not be applied.

skew_innov <- function(n) rexp(n) - 1
lap_innov <- function(n) rexp(n) * sample(c(-1, 1), n, TRUE)

ar_ols_vcov <- function(f, p) {
  X <- create_ar_matrix(f@original_series - f@intercept, p)
  f@m2 * solve(crossprod(X))
}

test_that("AR PMM2 fall-back reports the classical covariance", {
  set.seed(70)
  x <- as.numeric(arima.sim(list(ar = c(0.5, -0.3, 0.2)), 30,
                            innov = skew_innov(30), n.start = 50,
                            start.innov = skew_innov(50)))
  expect_warning(f <- ar_pmm2(x, order = 3), "returning the initial estimate")
  expect_true(attr(f, "fallback"))
  expect_false(f@convergence)
  expect_equal(unname(vcov(f)), unname(ar_ols_vcov(f, 3)))

  set.seed(71)
  x <- as.numeric(arima.sim(list(ar = 0.5), 300, innov = skew_innov(300)))
  f <- ar_pmm2(x, order = 1)
  expect_false(attr(f, "fallback"))
  g2 <- pmm2_variance_factor(f@m2, f@m3, f@m4)$g
  expect_equal(unname(vcov(f)), unname(g2 * ar_ols_vcov(f, 1)))
})

test_that("AR PMM3 fall-back reports the classical covariance", {
  found <- FALSE
  for (s in 1:100) {
    set.seed(s)
    x <- as.numeric(arima.sim(list(ar = c(0.5, -0.3, 0.2, 0.1, -0.1)), 30,
                              innov = lap_innov(30), n.start = 50,
                              start.innov = lap_innov(50)))
    f <- suppressWarnings(ar_pmm3(x, order = 5))
    if (isTRUE(attr(f, "fallback"))) { found <- TRUE; break }
  }
  expect_true(found)
  expect_equal(unname(vcov(f)), unname(ar_ols_vcov(f, 5)))
})

test_that("recursive MA fall-back and CSS fits report the classical covariance", {
  set.seed(35)
  w <- as.numeric(arima.sim(list(ma = 0.6), 300, innov = skew_innov(300)))
  classical <- function(f) {
    des <- pmm2_recursive_design(f)
    pmm2_variance_matrices(des$X, f@m2, f@m3, f@m4)
  }

  f_css <- ma_pmm2(w, order = 1, method = "css", ma_solver = "recursive",
                   include.mean = FALSE)
  expect_true(attr(f_css, "fallback"))
  expect_equal(unname(vcov(f_css)), unname(classical(f_css)$ols))

  f_pmm <- ma_pmm2(w, order = 1, ma_solver = "recursive", include.mean = FALSE)
  expect_false(attr(f_pmm, "fallback"))
  expect_equal(unname(vcov(f_pmm)), unname(classical(f_pmm)$pmm2))

  testthat::local_mocked_bindings(
    pmm2_recursive_solve = function(w, par_init, ...) {
      list(par = as.numeric(par_init), eps = NULL, convergence = FALSE,
           iterations = 0L, gnorm = NA_real_, fallback = TRUE)
    })
  f_ma <- ma_pmm2(w, order = 1, ma_solver = "recursive", include.mean = FALSE)
  expect_true(attr(f_ma, "fallback"))
  expect_equal(unname(vcov(f_ma)), unname(classical(f_ma)$ols))

  f_arima <- arima_pmm2(cumsum(w), order = c(0, 1, 1), ma_solver = "recursive")
  expect_true(attr(f_arima, "fallback"))
  expect_equal(unname(vcov(f_arima)), unname(classical(f_arima)$ols))
})
