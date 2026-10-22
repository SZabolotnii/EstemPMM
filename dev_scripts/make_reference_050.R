# make_reference_050.R -- freeze EstemPMM 0.5.0 estimates as a regression fixture.
#
# Run against the installed CRAN 0.5.0 build, never against the source tree:
#   Rscript dev_scripts/make_reference_050.R
# Writes tests/testthat/fixtures/reference_0.5.0.rds. The inputs are stored in
# the fixture, so the test does not depend on the RNG of the R version that
# replays it. test-reference-050.R requires that every case 0.5.0 reported as
# converged is reproduced by the current sources.

library(EstemPMM)
stopifnot(as.character(packageVersion("EstemPMM")) == "0.5.0")

set.seed(20261022)
cases <- list()

add_lm <- function(name, dat, formula, fitter = lm_pmm2) {
  fit <- suppressWarnings(fitter(formula, data = dat))
  cases[[name]] <<- list(kind = "lm", fitter = deparse(substitute(fitter)),
                         formula = formula, data = dat,
                         coef = fit@coefficients, converged = fit@convergence)
}

add_ts <- function(name, x, order, model_type, ma_solver = "linearized") {
  args <- list(x = x, order = order, model_type = model_type)
  if (model_type != "ar") args$ma_solver <- ma_solver
  fit <- suppressWarnings(do.call(ts_pmm2, args))
  cases[[name]] <<- list(kind = "ts", args = args,
                         coef = fit@coefficients, converged = fit@convergence)
}

# Linear regression: bundled data and simulated laws
data(auto_mpg)
am <- stats::na.omit(auto_mpg[, c("mpg", "horsepower", "weight")])
add_lm("auto_mpg_hp", am, mpg ~ horsepower)
add_lm("auto_mpg_hp_wt", am, mpg ~ horsepower + weight)

sim_lm <- function(n, p, err) {
  X <- matrix(stats::rnorm(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  y <- 1 + X %*% rep(c(1, -0.5), length.out = p) + err(n)
  data.frame(y = as.numeric(y), X)
}
laws <- list(
  gauss = function(n) stats::rnorm(n),
  gamma2 = function(n) stats::rgamma(n, 2, 1) - 2,
  exp1 = function(n) stats::rexp(n) - 1,
  lnorm = function(n) stats::rlnorm(n, 0, 0.8) - exp(0.32),
  lnorm_heavy = function(n) stats::rlnorm(n, 0, 1.5) - exp(1.125)
)
for (law in names(laws)) {
  for (np in list(c(200, 2), c(60, 10), c(30, 12))) {
    dat <- sim_lm(np[1], np[2], laws[[law]])
    f <- stats::as.formula(paste("y ~", paste(names(dat)[-1], collapse = " + ")))
    add_lm(sprintf("lm_%s_n%d_p%d", law, np[1], np[2]), dat, f)
  }
}
# PMM3 on a platykurtic law (its solver is not touched in 0.5.1; guard only)
dat <- sim_lm(200, 2, function(n) stats::runif(n, -1, 1))
add_lm("lm3_unif_n200_p2", dat, y ~ x1 + x2, fitter = lm_pmm3)

# Time series: simulated skewed innovations and WTI
innov <- function(n) stats::rgamma(n, 2, 1) - 2
x_ar <- as.numeric(stats::arima.sim(list(ar = 0.6), 300, innov = innov(300)))
x_ma <- as.numeric(stats::arima.sim(list(ma = -0.5), 300, innov = innov(300)))
x_arma <- as.numeric(stats::arima.sim(list(ar = 0.5, ma = 0.3), 300, innov = innov(300)))
x_arima <- cumsum(stats::arima.sim(list(ma = -0.5), 300, innov = innov(300)))
add_ts("ar1", x_ar, 1, "ar")
add_ts("ma1_lin", x_ma, 1, "ma")
add_ts("ma1_rec", x_ma, 1, "ma", "recursive")
add_ts("arma11_lin", x_arma, c(1, 1), "arma")
add_ts("arma11_rec", x_arma, c(1, 1), "arma", "recursive")
add_ts("arima011_lin", x_arima, c(0, 1, 1), "arima")
add_ts("arima011_rec", x_arima, c(0, 1, 1), "arima", "recursive")

data(DCOILWTICO)
wti <- stats::na.omit(DCOILWTICO$DCOILWTICO)
add_ts("wti_arima111_lin", wti, c(1, 1, 1), "arima")
add_ts("wti_arima111_rec", wti, c(1, 1, 1), "arima", "recursive")

out <- "tests/testthat/fixtures/reference_0.5.0.rds"
saveRDS(list(version = "0.5.0", cases = cases), out, version = 2)
conv <- vapply(cases, `[[`, logical(1), "converged")
cat(sprintf("%d cases written to %s; %d converged in 0.5.0\n",
            length(cases), out, sum(conv)))
print(names(cases)[!conv])
