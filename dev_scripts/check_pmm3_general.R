# Monte Carlo check of lm_pmm3(symmetric = FALSE): slope MSE against OLS,
# PMM2 and the symmetric PMM3, and coverage of the 95% interval from vcov().
# Design of the PMM-SLS-BRIDGE revision: y = 2 + 1.5 x + e, x ~ U(0, 5).
#   Rscript dev_scripts/check_pmm3_general.R [R] [cores]

suppressMessages(pkgload::load_all(".", quiet = TRUE))
library(parallel)
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 1000L
cores <- if (length(args) >= 2) as.integer(args[2]) else 8L
RNGkind("L'Ecuyer-CMRG"); set.seed(20261022)

gen_err <- function(n, type) switch(type,
  chisq3  = rchisq(n, 3) - 3,
  gamma1  = rexp(n) - 1,
  gamma2  = rgamma(n, 2, 1) - 2,
  lnorm05 = rlnorm(n, 0, 0.5) - exp(0.125),
  unif    = runif(n, -sqrt(3), sqrt(3)),
  gauss   = rnorm(n),
  t7      = rt(n, 7))

one <- function(i, tp, n) {
  x <- runif(n, 0, 5)
  d <- data.frame(x = x, y = 2 + 1.5 * x + gen_err(n, tp))
  fg <- suppressWarnings(lm_pmm3(y ~ x, d, symmetric = FALSE))
  se <- sqrt(vcov(fg)[2, 2])
  c(ols  = unname(coef(lm(y ~ x, d))[2]),
    pmm2 = unname(coef(suppressWarnings(lm_pmm2(y ~ x, d)))[2]),
    sym  = unname(coef(suppressWarnings(lm_pmm3(y ~ x, d)))[2]),
    gen  = unname(coef(fg)[2]),
    se   = se,
    cover = unname(abs(coef(fg)[2] - 1.5) < 1.96 * se),
    conv = fg@convergence)
}

tab <- NULL
for (tp in c("chisq3", "gamma2", "gamma1", "lnorm05", "unif", "gauss", "t7"))
  for (n in c(50, 200, 1000)) {
    m <- do.call(rbind, mclapply(seq_len(R), one, tp = tp, n = n, mc.cores = cores))
    mse <- colMeans((m[, c("ols", "pmm2", "sym", "gen")] - 1.5)^2)
    tab <- rbind(tab, data.frame(err = tp, n = n,
      RE_pmm2 = mse["ols"] / mse["pmm2"], RE_sym = mse["ols"] / mse["sym"],
      RE_gen = mse["ols"] / mse["gen"],
      se_over_sd = mean(m[, "se"]) / sd(m[, "gen"]),
      cover = mean(m[, "cover"]), conv = mean(m[, "conv"]), row.names = NULL))
  }
print(format(tab, digits = 3), row.names = FALSE)
