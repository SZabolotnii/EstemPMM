suppressMessages(pkgload::load_all("~/Project/Research/EstemPMM", quiet = TRUE))
library(parallel); RNGkind("L'Ecuyer-CMRG"); set.seed(2026)
errs <- list(
  t3      = function(n) rt(n, 3),
  laplace = function(n) rexp(n) * sample(c(-1, 1), n, TRUE),
  gauss   = function(n) rnorm(n),
  exp     = function(n) rexp(n) - 1,
  cont    = function(n) rnorm(n) + ifelse(runif(n) < 0.1, rnorm(n, 0, 10), 0),
  cauchy  = function(n) rcauchy(n)
)
R <- 1000
one <- function(i, nm, n) {
  x <- rnorm(n); d <- data.frame(x = x, y = 1 + 2 * x + errs[[nm]](n))
  res <- c(ols = unname(coef(lm(y ~ x, d))[2]))
  for (deg in 2:3) {
    f <- suppressWarnings(if (deg == 2) lm_wpmm2(y ~ x, d) else lm_wpmm3(y ~ x, d))
    se <- sqrt(vcov(f)[2, 2])
    res[paste0("b", deg)] <- coef(f)[2]; res[paste0("se", deg)] <- se
    res[paste0("cv", deg)] <- abs(coef(f)[2] - 2) < 1.96 * se
    res[paste0("conv", deg)] <- f@convergence
  }
  res
}
tab <- NULL
for (nm in names(errs)) for (n in c(100, 400)) {
  m <- do.call(rbind, mclapply(1:R, one, nm = nm, n = n, mc.cores = 9))
  mse <- function(b) mean((b - 2)^2)
  tab <- rbind(tab, data.frame(err = nm, n = n,
    sd2 = sd(m[, "b2"]), se2 = mean(m[, "se2"]), cov2 = mean(m[, "cv2"]),
    sd3 = sd(m[, "b3"]), se3 = mean(m[, "se3"]), cov3 = mean(m[, "cv3"]),
    mse2_ols = mse(m[, "b2"]) / mse(m[, "ols"]), mse3_ols = mse(m[, "b3"]) / mse(m[, "ols"]),
    conv = mean(m[, "conv2"] & m[, "conv3"])))
}
print(format(tab, digits = 3), row.names = FALSE)
