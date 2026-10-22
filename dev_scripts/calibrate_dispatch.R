# calibrate_dispatch.R -- Monte Carlo calibration of pmm_dispatch() (OPEN_ISSUES B-5, D-1).
#
#   Rscript dev_scripts/calibrate_dispatch.R [reps] [cores]
#
# Question: how much efficiency does the default dispatcher lose against the
# best of OLS / PMM2 / PMM3 in each (error law, n) cell, and do other
# thresholds lose less? The design is y = 1 + 2 x + e, x ~ N(0, 1); the loss
# is the squared error of the slope. Every replication stores the three
# estimates' errors, the default decision and the residual cumulants, so the
# thresholds can be re-tuned offline (analyze_dispatch.R) without refitting.
#
# Runs on the source tree (pkgload), not on the installed package.

args  <- commandArgs(TRUE)
reps  <- if (length(args) >= 1) as.integer(args[1]) else 2000L
cores <- if (length(args) >= 2) as.integer(args[2]) else 8L

pkgload::load_all(".", quiet = TRUE)
RNGkind("L'Ecuyer-CMRG")

# Error laws with mean 0. Grouped: symmetric platykurtic, symmetric
# mesokurtic/leptokurtic, skewed.
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
ns <- c(30L, 50L, 100L, 200L, 500L, 1000L)
cells <- expand.grid(law = names(laws), n = ns, stringsAsFactors = FALSE)

run_cell <- function(i) {
  law <- cells$law[i]; n <- cells$n[i]
  out <- vector("list", reps)
  for (r in seq_len(reps)) {
    x <- stats::rnorm(n)
    d <- data.frame(x = x, y = 1 + 2 * x + laws[[law]](n))
    ols <- stats::lm(y ~ x, data = d)
    f2 <- suppressWarnings(lm_pmm2(y ~ x, data = d))
    f3 <- suppressWarnings(lm_pmm3(y ~ x, data = d))
    dec <- pmm_dispatch(stats::residuals(ols), verbose = FALSE)
    out[[r]] <- data.frame(
      law = law, n = n, rep = r,
      se_ols = (stats::coef(ols)[2] - 2)^2,
      se_pmm2 = (f2@coefficients[2] - 2)^2,
      se_pmm3 = (f3@coefficients[2] - 2)^2,
      conv_pmm2 = f2@convergence, conv_pmm3 = f3@convergence,
      pick = dec$method, gamma3 = dec$gamma3, gamma4 = dec$gamma4,
      gamma6 = dec$gamma6, g2 = dec$g2, g3 = dec$g3,
      row.names = NULL)
  }
  do.call(rbind, out)
}

t0 <- Sys.time()
set.seed(20261022)
res <- parallel::mclapply(seq_len(nrow(cells)), run_cell, mc.cores = cores,
                          mc.set.seed = TRUE)
bad <- vapply(res, inherits, logical(1), "try-error")
if (any(bad)) stop("cells failed: ", paste(which(bad), collapse = ", "))
res <- do.call(rbind, res)

# Population cumulants of each law (2e6 draws). t3 has no 4th moment and t5
# no 6th, so their g2 / g3 are not defined.
pop <- do.call(rbind, lapply(names(laws), function(l) {
  set.seed(1); e <- laws[[l]](2e6); cu <- .residual_cumulants(e)
  data.frame(law = l, gamma3 = cu$gamma3, gamma4 = cu$gamma4,
             gamma6 = cu$gamma6,
             g2 = .g2_factor(cu$gamma3, cu$gamma4),
             g3 = .g3_factor(cu$gamma4, cu$gamma6))
}))
pop[pop$law == "t3", c("gamma4", "gamma6", "g2", "g3")] <- NA
pop[pop$law == "t5", c("gamma6", "g3")] <- NA

dir.create("mc_results_dispatch", showWarnings = FALSE)
saveRDS(list(reps = res, pop = pop, reps_per_cell = reps,
             estempmm = as.character(utils::packageVersion("EstemPMM")),
             commit = system("git rev-parse --short HEAD", intern = TRUE),
             elapsed = Sys.time() - t0),
        "mc_results_dispatch/calibrate_dispatch.rds", version = 2)
cat(sprintf("%d rows, %d cells, %.1f min\n", nrow(res), nrow(cells),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
