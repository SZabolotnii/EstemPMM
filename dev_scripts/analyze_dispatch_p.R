# analyze_dispatch_p.R -- tables for the multi-regressor pmm_dispatch() calibration.
#
#   Rscript dev_scripts/analyze_dispatch_p.R
#
# Reads mc_results_dispatch/calibrate_dispatch_p.rds (calibrate_dispatch_p.R),
# writes mc_results_dispatch/p_grid/*.csv and prints a summary. Regret of a
# rule in a cell = MSE(rule) / min(MSE of OLS, PMM2, PMM3) for the x1 slope.
# Rules are tuned on odd replications and judged on even ones.

x <- readRDS("mc_results_dispatch/calibrate_dispatch_p.rds")
d <- x$reps
key <- paste(d$law, d$n, d$p)
cells <- unique(d[, c("law", "n", "p")])
cells$key <- paste(cells$law, cells$n, cells$p)
law_order <- unique(d$law)
out_dir <- "mc_results_dispatch/p_grid"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

se_of <- function(pick) {
  ifelse(pick == "PMM2", d$se_pmm2, ifelse(pick == "PMM3", d$se_pmm3, d$se_ols))
}
cell_mse <- function(se, rows = rep(TRUE, nrow(d))) {
  tapply(se[rows], key[rows], mean)[cells$key]
}

# Rule family. With nm = "n" and r_min = 0 it is pmm_dispatch(); nm = "n-p"
# uses the residual degrees of freedom in both sample-size thresholds; r_min
# requires n / p >= r_min before the leptokurtic PMM3 arm fires.
rule <- function(z = 2, g3_thr = 0.95, n3 = 100, nm = "n", r_min = 0) {
  a3 <- abs(d$gamma3)
  m <- if (nm == "n") d$n else d$n - d$p
  asym <- a3 > 0.3
  platy <- d$gamma4 < -0.7
  clear <- a3 > z * sqrt(6 / m)
  ifelse(asym & platy & !clear, "PMM3",
  ifelse(a3 > 1.0, "PMM2",
  ifelse(asym & platy, "PMM2",
  ifelse(asym & d$g2 < 0.95, "PMM2",
  ifelse(!asym & platy, "PMM3",
  ifelse(!asym & d$g3 < g3_thr & m >= n3 & d$n / d$p >= r_min, "PMM3", "OLS"))))))
}
rule_050 <- function() rule(z = 0, g3_thr = 0)

def_pick <- rule()
stopifnot(identical(def_pick, d$pick))   # the reimplementation is pmm_dispatch()

# ---- 1. Cell table under the current default --------------------------------
mse_o <- cell_mse(d$se_ols); mse_2 <- cell_mse(d$se_pmm2); mse_3 <- cell_mse(d$se_pmm3)
best <- pmin(mse_o, mse_2, mse_3)
share <- function(m) tapply(d$pick == m, key, mean)[cells$key]
tab <- data.frame(cells[, c("law", "n", "p")],
  pmm2_vs_ols = mse_2 / mse_o, pmm3_vs_ols = mse_3 / mse_o,
  best = c("OLS", "PMM2", "PMM3")[apply(cbind(mse_o, mse_2, mse_3), 1, which.min)],
  pick_pmm2 = share("PMM2"), pick_pmm3 = share("PMM3"),
  default_vs_ols = cell_mse(se_of(def_pick)) / mse_o,
  rule050_vs_ols = cell_mse(se_of(rule_050())) / mse_o,
  regret_default = cell_mse(se_of(def_pick)) / best,
  fb_pmm2 = tapply(d$fb_pmm2, key, mean)[cells$key],
  fb_pmm3 = tapply(d$fb_pmm3, key, mean)[cells$key], row.names = NULL)
tab <- tab[order(tab$p, match(tab$law, law_order), tab$n), ]
utils::write.csv(tab, file.path(out_dir, "cells_default.csv"), row.names = FALSE)

cat(sprintf("EstemPMM %s @ %s, %d reps per cell, %d cells\n\n", x$estempmm,
            x$commit, x$reps_per_cell, nrow(cells)))
cat("== Default rule by p and n: regret and ratio to OLS ==\n")
print(do.call(rbind, lapply(split(tab, list(tab$n, tab$p)), function(t) data.frame(
  p = t$p[1], n = t$n[1], mean_regret = exp(mean(log(t$regret_default))),
  worst_regret = max(t$regret_default), worst_cell = t$law[which.max(t$regret_default)],
  worst_vs_ols = max(t$default_vs_ols), worst_vs_ols_cell = t$law[which.max(t$default_vs_ols)],
  worst_050_vs_ols = max(t$rule050_vs_ols)))), digits = 3, row.names = FALSE)

cat("\n== Symmetric leptokurtic laws: PMM3 vs OLS and the default's choice ==\n")
print(tab[tab$law %in% c("logistic", "laplace", "t5", "t3"),
          c("law", "n", "p", "pmm2_vs_ols", "pmm3_vs_ols", "pick_pmm3",
            "default_vs_ols", "fb_pmm3")], digits = 3, row.names = FALSE)

cat("\n== Cells where the default is more than 5% worse than OLS ==\n")
print(tab[tab$default_vs_ols > 1.05, c("law", "n", "p", "pmm2_vs_ols", "pmm3_vs_ols",
          "pick_pmm2", "pick_pmm3", "default_vs_ols", "rule050_vs_ols")],
      digits = 3, row.names = FALSE)

# ---- 2. Variants, tuned on odd reps and judged on even reps ----------------
odd <- d$rep %% 2 == 1
ev <- !odd
score <- function(pick, rows) {
  m <- cell_mse(se_of(pick), rows)
  b <- pmin(cell_mse(d$se_ols, rows), cell_mse(d$se_pmm2, rows), cell_mse(d$se_pmm3, rows))
  o <- cell_mse(d$se_ols, rows)
  r <- m / b
  c(mean_log_regret = mean(log(r)), worst_regret = max(r), worst_vs_ols = max(m / o),
    p2 = mean(log(r[cells$p == 2])), p6 = mean(log(r[cells$p == 6])),
    p11 = mean(log(r[cells$p == 11])))
}
grid <- expand.grid(z = c(1.5, 2), g3_thr = c(0.9, 0.95), n3 = c(100, 200),
                    nm = c("n", "n-p"), r_min = c(0, 10, 20),
                    stringsAsFactors = FALSE)
sc <- t(vapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]
  c(score(rule(g$z, g$g3_thr, g$n3, g$nm, g$r_min), odd),
    score(rule(g$z, g$g3_thr, g$n3, g$nm, g$r_min), ev))
}, numeric(12)))
colnames(sc) <- c(paste0("tune.", c("mlr", "worst", "worst_ols", "p2", "p6", "p11")),
                  paste0("test.", c("mlr", "worst", "worst_ols", "p2", "p6", "p11")))
grid <- cbind(grid, sc)
grid <- grid[order(grid$tune.mlr), ]
utils::write.csv(grid, file.path(out_dir, "variants.csv"), row.names = FALSE)
ref <- rbind(rule_050 = c(score(rule_050(), odd), score(rule_050(), ev)),
             default = c(score(def_pick, odd), score(def_pick, ev)))
colnames(ref) <- colnames(sc)
utils::write.csv(ref, file.path(out_dir, "reference_rules.csv"))

cat("\n== Reference rules (tune = odd reps, test = even reps; mlr = mean log regret) ==\n")
print(round(ref, 4))
cat("\n== Variants: top 12 by tuning mean log regret ==\n")
print(head(grid, 12), digits = 3, row.names = FALSE)

# Per-cell view of the best variant against the default on the held-out half
b <- grid[1, ]
bp <- rule(b$z, b$g3_thr, b$n3, b$nm, b$r_min)
best_ev <- pmin(cell_mse(d$se_ols, ev), cell_mse(d$se_pmm2, ev), cell_mse(d$se_pmm3, ev))
cmp <- data.frame(cells[, c("law", "n", "p")],
  default = cell_mse(se_of(def_pick), ev) / best_ev,
  variant = cell_mse(se_of(bp), ev) / best_ev,
  variant_vs_ols = cell_mse(se_of(bp), ev) / cell_mse(d$se_ols, ev), row.names = NULL)
cmp <- cmp[order(cmp$p, match(cmp$law, law_order), cmp$n), ]
utils::write.csv(cmp, file.path(out_dir, "heldout_default_vs_best_variant.csv"), row.names = FALSE)
cat("\n== Best variant:", paste(names(b)[1:5], unlist(b[1:5]), sep = "=", collapse = ", "), "==\n")
cat(sprintf("cells better by >3%%: %d, worse by >3%%: %d\n",
            sum(cmp$variant < cmp$default / 1.03), sum(cmp$variant > cmp$default * 1.03)))
print(cmp[abs(cmp$variant / cmp$default - 1) > 0.03, ], digits = 3, row.names = FALSE)
