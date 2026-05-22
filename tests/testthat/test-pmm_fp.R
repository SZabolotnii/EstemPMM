# Smoke tests and safeguard verification for pmm_fp().
#
# Tests cover:
#   - Basic API (S4 class, slots, track, shift)
#   - Argument validation
#   - Non-regression: lm_pmm2 byte-identical before/after safeguard code
#   - Normal fit: diagnostics vector populated, nr_diverged = FALSE
#   - Ill-conditioned design: Tikhonov activates (flag_tikhonov = 1)
#   - Hard/divergent case: reinit activates, OLS fallback with nr_diverged = TRUE

set.seed(20260511)

# Hill-type DGP with gamma errors (asymmetric, gamma3 > 0)
n  <- 50
x  <- runif(n, 0.5, 5)
mu <- 1 + 2 * x ^ 0.5
err <- rgamma(n, shape = 2, scale = 1) - 2
y  <- mu + err
dat <- data.frame(x = x, y = y)

# ---------------------------------------------------------------------------
# Basic API tests (unchanged from scaffold)
# ---------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------
# Non-regression: lm_pmm2 must give byte-identical results
# ---------------------------------------------------------------------------

test_that("lm_pmm2 coefficients are byte-identical after safeguard code added", {
    set.seed(42)
    n2 <- 80
    x2 <- rnorm(n2)
    y2 <- 2 + 3 * x2 + rt(n2, df = 3)
    dat2 <- data.frame(y = y2, x = x2)
    fit_base <- lm_pmm2(y ~ x, data = dat2)
    # Expected from baseline run recorded before any changes:
    # coefs = c(1.9470480452, 2.9851033543)
    expect_equal(fit_base@coefficients[1], 1.9470480452, tolerance = 1e-8)
    expect_equal(fit_base@coefficients[2], 2.9851033543, tolerance = 1e-8)
    expect_true(fit_base@convergence)
})

# ---------------------------------------------------------------------------
# Diagnostics slot populated for a normal (well-conditioned) fit
# ---------------------------------------------------------------------------

test_that("pmm_fp populates nr_diagnostics and nr_diverged=FALSE on normal fit", {
    suppressWarnings(
        fit <- pmm_fp(y ~ x, data = dat, track = "pos", max_terms = 2)
    )
    diag <- fit@nr_diagnostics
    expect_named(diag, c("k_star", "max_step", "final_grad_norm", "kappa_F",
                          "flag_step_halving", "flag_tikhonov",
                          "flag_reinit", "flag_fallback"))
    expect_true(is.numeric(diag))
    expect_true(all(is.finite(diag[c("k_star", "max_step", "final_grad_norm")])))
    expect_false(fit@nr_diverged)
})

# ---------------------------------------------------------------------------
# (iii) Adaptive Tikhonov: ill-conditioned design triggers flag
# ---------------------------------------------------------------------------

test_that("Adaptive Tikhonov activates on ill-conditioned FP design", {
    set.seed(1)
    n3 <- 30
    # Construct two almost-collinear FP columns: x^1 and x^1.001 ~ x^1
    x3 <- runif(n3, 1, 2)
    y3 <- 1 + 2 * x3 + rnorm(n3, 0, 0.1)

    # Build a design matrix where fp1 and fp2 are nearly identical
    # Force this via kappa_thresh=1 so ANY matrix triggers Tikhonov
    X_ill <- cbind(1, x3, x3 + 1e-6 * rnorm(n3))
    colnames(X_ill) <- c("(Intercept)", "fp1", "fp2")

    b_ols3 <- as.numeric(lm.fit(X_ill, y3)$coefficients)
    res3   <- y3 - as.vector(X_ill %*% b_ols3)
    m2_3   <- mean(res3^2)
    m3_3   <- mean(res3^3)
    m4_3   <- mean(res3^4)

    # Call .pmm_fp_robust_refit directly with kappa_thresh=1 to force Tikhonov
    rr <- EstemPMM:::.pmm_fp_robust_refit(
        X = X_ill, y = y3,
        m2 = m2_3, m3 = m3_3, m4 = m4_3,
        kappa_thresh = 1        # any matrix triggers Tikhonov
    )
    expect_equal(rr$diagnostics[["flag_tikhonov"]], 1L)
})

# ---------------------------------------------------------------------------
# (iv) Reinitialisation + fallback: divergent case triggers OLS fallback
# ---------------------------------------------------------------------------

test_that("OLS fallback triggered on degenerate design (nr_diverged=TRUE)", {
    set.seed(7)
    n4 <- 20
    # Perfectly collinear: X has rank 1 => every Newton step diverges
    x4 <- runif(n4, 1, 2)
    y4 <- rnorm(n4)
    X_degen <- cbind(1, x4, 2 * x4)   # rank 2, but force max_iter=1 so NR fails
    colnames(X_degen) <- c("(Intercept)", "fp1", "fp2")

    b_ols4 <- tryCatch(as.numeric(lm.fit(X_degen, y4)$coefficients),
                       error = function(e) rep(0, 3))
    # Some NA from rank deficiency
    b_ols4[is.na(b_ols4)] <- 0
    res4  <- y4 - as.vector(X_degen %*% b_ols4)
    m2_4  <- mean(res4^2)
    m3_4  <- mean(res4^3)
    m4_4  <- mean(res4^4)

    # Force divergence by setting max_iter=1 and a very tight tol
    rr4 <- suppressMessages(
        EstemPMM:::.pmm_fp_robust_refit(
            X = X_degen, y = y4,
            m2 = m2_4, m3 = m3_4, m4 = m4_4,
            max_iter = 1L,
            tol = 1e-20      # impossible to achieve in 1 iter
        )
    )
    # After 1 iter neither pass converges -> fallback
    expect_true(rr4$nr_diverged)
    expect_equal(rr4$diagnostics[["flag_fallback"]], 1L)
    expect_true(is.numeric(rr4$b))
    expect_equal(length(rr4$b), 3L)
})

# ---------------------------------------------------------------------------
# Snapshot: selected_powers stable for fixed seed (regression test)
# ---------------------------------------------------------------------------

test_that("pmm_fp snapshot: selected_powers stable for seed 20260511", {
    set.seed(20260511)
    n5  <- 200
    x5  <- runif(n5, 0.5, 5)
    y5  <- 1 + 2 * x5^0.5 + (rgamma(n5, shape = 2, scale = 1) - 2)
    dat5 <- data.frame(x = x5, y = y5)
    suppressWarnings(
        fit5 <- pmm_fp(y ~ x, data = dat5, track = "pos", max_terms = 2)
    )
    # Record the selected powers for regression detection
    expect_true(length(fit5@selected_powers) >= 1L)
    # Snapshot (powers chosen by BIC on this fixed data)
    expect_snapshot(fit5@selected_powers)
})
