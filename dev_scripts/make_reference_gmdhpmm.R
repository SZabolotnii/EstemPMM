# Reference estimates for lm_wpmm2() / lm_wpmm3() from the gmdhpmm package
# (EstemPMM-GMDH, R/weak_pmm.R), whose weak-moment estimators were ported here.
# Writes tests/testthat/fixtures/reference_gmdhpmm.rds; test-wpmm.R replays it.
# Run from the package root with gmdhpmm installed:
#   Rscript dev_scripts/make_reference_gmdhpmm.R

library(gmdhpmm)

errs <- list(
  t3     = function(n) rt(n, 3) / sqrt(3),
  asym   = function(n) { e <- rnorm(n); h <- runif(n) < 0.12
                         e[h] <- -abs(rt(sum(h), 2)) * 2.5; e },
  gauss  = function(n) rnorm(n),
  exp    = function(n) rexp(n) - 1,
  cauchy = function(n) rcauchy(n)
)

cases <- list()
for (nm in names(errs)) for (n in c(60, 200)) {
  set.seed(1000 + n)
  v1 <- rnorm(n); v2 <- rnorm(n)
  mu <- 0.5 + v1 - 0.8 * v2 + 0.4 * v1 * v2 - 0.3 * v1^2 + 0.2 * v2^2
  d <- kg2_design(v1, v2, mu + errs[[nm]](n))
  cases[[paste0(nm, "_n", n)]] <- list(
    data  = d,
    wpmm2 = fit_kg2_weak_pmm2(d),
    wpmm3 = fit_kg2_weak_pmm3(d))
}

saveRDS(list(gmdhpmm_version = as.character(utils::packageVersion("gmdhpmm")),
             formula = y ~ b1 + b2 + b12 + b11 + b22,
             cases = cases),
        "tests/testthat/fixtures/reference_gmdhpmm.rds")
cat("Saved", length(cases), "cases (gmdhpmm",
    as.character(utils::packageVersion("gmdhpmm")), ")\n")
