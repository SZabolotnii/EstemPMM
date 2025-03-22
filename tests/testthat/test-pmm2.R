test_that("pmm2 works on a simple linear example", {
  set.seed(123)
  n <- 50
  x <- rnorm(n)
  # Генеруємо y з негаусовою похибкою:
  y <- 1 + 2*x + rt(n, df=4)
  dat <- data.frame(x, y)

  fit <- lm_pmm2(y ~ x, data=dat, max_iter=50)

  expect_s4_class(fit, "PMM2fit")
  expect_true(length(fit@coefficients) == 2)
  expect_true(fit@convergence)

  # Перевіряємо, чи немає NA:
  expect_false(anyNA(fit@coefficients))

  # Спробуємо p-value:
  inf_tab <- pmm2_inference(fit, formula=y~x, data=dat, B=50)
  expect_equal(nrow(inf_tab), 2)
})
