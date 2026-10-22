# Monte Carlo check of the symmetric PMM3 solver (lm_pmm3): slope MSE relative
# to OLS, worst error, non-convergence and OLS fall-back rates on a grid of
# error laws, sample sizes and numbers of regressors (500 replications).
#   Rscript dev_scripts/check_pmm3_solver.R [package path]
args <- commandArgs(trailingOnly = TRUE)
pkg <- if (length(args)) args[1] else "."
suppressMessages(pkgload::load_all(pkg, quiet = TRUE))
library(parallel); RNGkind("L'Ecuyer-CMRG"); set.seed(7)
gen <- function(n, tp) switch(tp,
  unif = runif(n, -1, 1), bimod = rnorm(n, sample(c(-1, 1), n, TRUE), 0.3),
  laplace = rexp(n) * sample(c(-1, 1), n, TRUE), t3 = rt(n, 3), t5 = rt(n, 5),
  gamma2 = rgamma(n, 2) - 2, gauss = rnorm(n))
one <- function(i, tp, n, p) {
  X <- matrix(rnorm(n * (p - 1)), n)
  d <- data.frame(y = 1 + 2 * X[, 1] + gen(n, tp), X)
  f <- suppressWarnings(lm_pmm3(y ~ ., d))
  c(b = unname(coef(f)[2]), ols = unname(coef(lm(y ~ ., d))[2]),
    conv = f@convergence, fb = isTRUE(attr(f, "fallback")))
}
tab <- NULL
for (tp in c("unif", "bimod", "laplace", "t3", "t5", "gamma2", "gauss"))
  for (np in list(c(30, 2), c(30, 6), c(30, 12), c(100, 2), c(100, 10))) {
    m <- do.call(rbind, mclapply(1:500, one, tp = tp, n = np[1], p = np[2], mc.cores = 9))
    tab <- rbind(tab, data.frame(tp, n = np[1], p = np[2],
      ratio = mean((m[, "b"] - 2)^2) / mean((m[, "ols"] - 2)^2),
      max_err = max(abs(m[, "b"] - 2)), nonconv = mean(!m[, "conv"]), fb = mean(m[, "fb"])))
  }
print(format(tab, digits = 3), row.names = FALSE)
