# Reference estimates for lm_pmm3(symmetric = FALSE) from the general degree-3
# solver of the PMM-SLS-BRIDGE revision (AISM-D-26-00101),
# experiments/revision_2026-09/common.R::pmm_reg. Writes
# tests/testthat/fixtures/reference_aism_pmm3.rds; test-pmm3_general.R replays it.
#   Rscript dev_scripts/make_reference_aism.R [path/to/common.R]
#
# pmm_reg(outer = 1) solves the same equation as lm_pmm3(symmetric = FALSE):
# moments from the OLS residuals, held fixed. pmm_reg(outer = 3), the variant
# behind the revision tables, re-estimates the moments twice more and is stored
# for comparison with adaptive = TRUE.

args <- commandArgs(trailingOnly = TRUE)
common <- if (length(args)) args[1] else
  "~/Project/Research/PMM-SLS-BRIDGE/experiments/revision_2026-09/common.R"
source(common)

cases <- list()
for (tp in c("chisq3", "gamma2", "gamma1", "lnorm05", "unif", "gauss", "t7"))
  for (n in c(50, 200)) {
    set.seed(20260915L + 7L * n + sum(utf8ToInt(tp)))
    x <- runif(n, 0, 5)
    d <- data.frame(x = x, y = BETA[1] + BETA[2] * x + gen_err(n, tp))
    X <- cbind(1, x)
    cases[[paste0(tp, "_n", n)]] <- list(
      data   = d,
      outer1 = pmm_reg(d$y, X, 3, outer = 1),
      outer3 = pmm_reg(d$y, X, 3))
  }

saveRDS(list(source = "PMM-SLS-BRIDGE experiments/revision_2026-09/common.R::pmm_reg",
             cases = cases),
        "tests/testthat/fixtures/reference_aism_pmm3.rds")
cat("Saved", length(cases), "cases\n")
