# Regression guard against EstemPMM 0.5.0 (fixture: dev_scripts/make_reference_050.R).
# Every case that 0.5.0 reported as converged must be reproduced by the
# safeguarded PMM2 iteration, which takes the full Newton step whenever it
# reduces the norm of the estimating function. A case that reaches the same
# root through one damped step can differ in the last bit, hence the 1e-12
# tolerance; 23 of the 24 cases checked here are in fact bit-identical.

ref <- readRDS(test_path("fixtures", "reference_0.5.0.rds"))

# Changed on purpose in 0.5.1. With lognormal errors, n = 30 and p = 12, 0.5.0
# reached by undamped Newton steps a root where half of the Jacobian weights
# are negative (a second branch of the estimating equation, error 2.6 against
# 0.95 for OLS). The line search refuses the step that leads there, stalls,
# and the fit returns the OLS start with convergence = FALSE.
changed <- "lm_lnorm_n30_p12"

refit <- function(cs) {
  if (cs$kind == "lm") {
    fitter <- get(cs$fitter, envir = asNamespace("EstemPMM"))
    suppressWarnings(fitter(cs$formula, data = cs$data))
  } else {
    suppressWarnings(do.call(ts_pmm2, cs$args))
  }
}

# The WTI ARIMA(1,1,1) start comes from stats::arima(method = "CSS-ML"), whose
# optimiser stops at a platform-dependent point: the AR and MA roots nearly
# cancel (-0.26 vs 0.29), so the likelihood is flat. The fixture was made on
# macOS; Linux and Windows differ by 2e-7 (linearised) and 5e-11 (recursive),
# and changing only the optimiser's reltol moves the fit by 1e-6.
start_sensitive <- c("wti_arima111_lin", "wti_arima111_rec")

for (nm in setdiff(names(ref$cases), changed)) {
  cs <- ref$cases[[nm]]
  if (!isTRUE(cs$converged)) next
  tol <- if (nm %in% start_sensitive) 1e-5 else 1e-12
  test_that(paste("0.5.0 estimate reproduced:", nm), {
    fit <- refit(cs)
    expect_true(fit@convergence)
    expect_equal(unname(fit@coefficients), unname(cs$coef), tolerance = tol)
  })
}

test_that("a stalled PMM2 iteration returns the OLS start and says so", {
  cs <- ref$cases[[changed]]
  expect_warning(fit <- lm_pmm2(cs$formula, data = cs$data),
                 "line search stalled")
  expect_false(fit@convergence)
  ols <- stats::coef(stats::lm(cs$formula, data = cs$data))
  expect_equal(unname(fit@coefficients), unname(ols), tolerance = 1e-10)
  # The returned estimate is OLS, so its covariance is the OLS one, not g2 times it
  expect_true(attr(fit, "fallback"))
  X <- attr(fit, "model_matrix")
  expect_equal(unname(vcov(fit)), unname(fit@m2 * solve(crossprod(X))))
  expect_output(summary(fit), "OLS fall-back")
})
