# Monte Carlo check of the asymptotic PMM2 covariance returned by vcov():
# mean standard error against the empirical SD of the slope, coverage of the
# 95% interval, and how often lm_pmm2() fell back to OLS.
# y = 1 + 2 x1 + 0 x2 + ... + e, x ~ N(0, 1), p columns including intercept.
#   Rscript dev_scripts/check_pmm2_se.R [R] [cores]

suppressMessages(pkgload::load_all(".", quiet = TRUE))
library(parallel)
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 1000L
cores <- if (length(args) >= 2) as.integer(args[2]) else 8L
RNGkind("L'Ecuyer-CMRG"); set.seed(20261022)

gen_err <- function(n, type) switch(type,
  chisq3  = rchisq(n, 3) - 3,
  gamma2  = rgamma(n, 2, 1) - 2,
  exp     = rexp(n) - 1,
  lnorm05 = rlnorm(n, 0, 0.5) - exp(0.125),
  lnorm1  = rlnorm(n, 0, 1) - exp(0.5),
  gauss   = rnorm(n))

one <- function(i, tp, n, p) {
  X <- matrix(rnorm(n * (p - 1)), n)
  d <- data.frame(y = 1 + 2 * X[, 1] + gen_err(n, tp), X)
  f <- suppressWarnings(lm_pmm2(y ~ ., data = d))
  se <- sqrt(vcov(f)[2, 2])
  b <- unname(coef(f)[2])
  c(b = b, ols = unname(coef(lm(y ~ ., d))[2]), se = se,
    cover = abs(b - 2) < 1.96 * se, fallback = isTRUE(attr(f, "fallback")))
}

tab <- NULL
for (tp in c("chisq3", "gamma2", "exp", "lnorm05", "lnorm1", "gauss"))
  for (np in list(c(30, 2), c(30, 6), c(50, 2), c(200, 2), c(1000, 2))) {
    m <- do.call(rbind, mclapply(seq_len(R), one, tp = tp, n = np[1], p = np[2],
                                 mc.cores = cores))
    tab <- rbind(tab, data.frame(err = tp, n = np[1], p = np[2],
      RE = mean((m[, "ols"] - 2)^2) / mean((m[, "b"] - 2)^2),
      se_over_sd = mean(m[, "se"]) / sd(m[, "b"]),
      cover = mean(m[, "cover"]), fallback = mean(m[, "fallback"])))
  }
print(format(tab, digits = 3), row.names = FALSE)
