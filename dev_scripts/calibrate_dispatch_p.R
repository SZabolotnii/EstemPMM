# calibrate_dispatch_p.R -- pmm_dispatch() calibration with several regressors.
#
#   Rscript dev_scripts/calibrate_dispatch_p.R [reps] [cores]
#
# calibrate_dispatch.R used one regressor and so never saw the wrong-root
# failures of the symmetric PMM3 solver at p >= 6 (fixed in 37e57e9). This
# run repeats the calibration on a grid of model sizes with the current
# solvers: y = 1 + 2 x1 + 0 x2 + ... + 0 x_{p-1} + e, x ~ N(0, 1), p columns
# including the intercept; loss = squared error of the x1 slope. Every
# replication stores the three estimates' errors, the solver fall-backs, the
# default decision and the residual cumulants, so analyze_dispatch_p.R can
# re-tune the rule offline.

args  <- commandArgs(TRUE)
reps  <- if (length(args) >= 1) as.integer(args[1]) else 2000L
cores <- if (length(args) >= 2) as.integer(args[2]) else 8L

pkgload::load_all(".", quiet = TRUE)
RNGkind("L'Ecuyer-CMRG")

# Same 20 laws as calibrate_dispatch.R
laws <- list(
  uniform     = function(n) stats::runif(n, -1, 1),
  bimodal     = function(n) stats::rnorm(n, sample(c(-1, 1), n, TRUE), 0.3),
  beta22      = function(n) stats::rbeta(n, 2, 2) - 0.5,
  triangular  = function(n) stats::runif(n) + stats::runif(n) - 1,
  normal      = function(n) stats::rnorm(n),
  logistic    = function(n) stats::rlogis(n),
  laplace     = function(n) stats::rexp(n) * sample(c(-1, 1), n, TRUE),
  t5          = function(n) stats::rt(n, 5),
  t3          = function(n) stats::rt(n, 3),
  beta25      = function(n) stats::rbeta(n, 2, 5) - 2 / 7,
  weibull36   = function(n) stats::rweibull(n, 3.6) - gamma(1 + 1 / 3.6),
  skewnorm    = function(n) { z <- abs(stats::rnorm(n)); u <- stats::rnorm(n)
                              d <- 0.95; d * z + sqrt(1 - d^2) * u - d * sqrt(2 / pi) },
  gamma10     = function(n) stats::rgamma(n, 10) - 10,
  gamma4      = function(n) stats::rgamma(n, 4) - 4,
  gamma2      = function(n) stats::rgamma(n, 2) - 2,
  exp         = function(n) stats::rexp(n) - 1,
  gamma05     = function(n) stats::rgamma(n, 0.5) - 0.5,
  lnorm025    = function(n) stats::rlnorm(n, 0, 0.25) - exp(0.25^2 / 2),
  lnorm05     = function(n) stats::rlnorm(n, 0, 0.5) - exp(0.5^2 / 2),
  lnorm1      = function(n) stats::rlnorm(n, 0, 1) - exp(0.5)
)
cells <- expand.grid(law = names(laws), n = c(30L, 50L, 100L, 200L, 500L, 1000L),
                     p = c(2L, 6L, 11L), stringsAsFactors = FALSE)

run_cell <- function(i) {
  law <- cells$law[i]; n <- cells$n[i]; p <- cells$p[i]
  out <- vector("list", reps)
  for (r in seq_len(reps)) {
    X <- matrix(stats::rnorm(n * (p - 1)), n)
    d <- data.frame(y = 1 + 2 * X[, 1] + laws[[law]](n), X)
    ols <- stats::lm(y ~ ., data = d)
    f2 <- suppressWarnings(lm_pmm2(y ~ ., data = d))
    f3 <- suppressWarnings(lm_pmm3(y ~ ., data = d))
    dec <- pmm_dispatch(stats::residuals(ols), verbose = FALSE)
    out[[r]] <- data.frame(
      law = law, n = n, p = p, rep = r,
      se_ols = (stats::coef(ols)[2] - 2)^2,
      se_pmm2 = (f2@coefficients[2] - 2)^2,
      se_pmm3 = (f3@coefficients[2] - 2)^2,
      fb_pmm2 = isTRUE(attr(f2, "fallback")),
      fb_pmm3 = isTRUE(attr(f3, "fallback")),
      pick = dec$method, gamma3 = dec$gamma3, gamma4 = dec$gamma4,
      gamma6 = dec$gamma6, g2 = dec$g2, g3 = dec$g3,
      row.names = NULL)
  }
  do.call(rbind, out)
}

t0 <- Sys.time()
set.seed(20261023)
res <- parallel::mclapply(seq_len(nrow(cells)), run_cell, mc.cores = cores,
                          mc.set.seed = TRUE)
bad <- vapply(res, inherits, logical(1), "try-error")
if (any(bad)) stop("cells failed: ", paste(which(bad), collapse = ", "))
res <- do.call(rbind, res)

dir.create("mc_results_dispatch", showWarnings = FALSE)
saveRDS(list(reps = res, reps_per_cell = reps,
             estempmm = as.character(utils::packageVersion("EstemPMM")),
             commit = system("git rev-parse --short HEAD", intern = TRUE),
             elapsed = Sys.time() - t0),
        "mc_results_dispatch/calibrate_dispatch_p.rds", version = 2)
cat(sprintf("%d rows, %d cells, %.1f min\n", nrow(res), nrow(cells),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
