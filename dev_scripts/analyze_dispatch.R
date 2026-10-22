# analyze_dispatch.R -- tables for the pmm_dispatch() calibration (B-5, D-1).
#
#   Rscript dev_scripts/analyze_dispatch.R
#
# Reads mc_results_dispatch/calibrate_dispatch.rds (calibrate_dispatch.R) and
# writes mc_results_dispatch/*.csv plus a printed summary. "Regret" of a rule
# in a cell is MSE(rule) / min(MSE(OLS), MSE(PMM2), MSE(PMM3)) for the slope:
# what the rule loses against the best single method chosen with hindsight.
# Thresholds are tuned on the odd replications and judged on the even ones.

x <- readRDS("mc_results_dispatch/calibrate_dispatch.rds")
d <- x$reps
pop <- x$pop
cell_key <- paste(d$law, d$n)
cells <- unique(d[, c("law", "n")])
cells$key <- paste(cells$law, cells$n)

# Squared error of the estimator a rule picks, replication by replication
se_of <- function(pick) {
  ifelse(pick == "PMM2", d$se_pmm2, ifelse(pick == "PMM3", d$se_pmm3, d$se_ols))
}

cell_mse <- function(se, rows = rep(TRUE, nrow(d))) {
  m <- tapply(se[rows], cell_key[rows], mean)
  m[cells$key]
}

# Decision rule in the shape of pmm_dispatch(), plus one optional arm:
# symmetric residuals with estimated g3 below g3_thr go to PMM3 whatever the
# sign of the excess kurtosis.
rule <- function(s_thr = 0.3, strong = 1.0, k_thr = -0.7, g2_thr = 0.95,
                 g3_thr = -Inf) {
  a3 <- abs(d$gamma3)
  asym <- a3 > s_thr
  platy <- d$gamma4 < k_thr
  ifelse(a3 > strong, "PMM2",
  ifelse(asym & d$g2 < g2_thr, "PMM2",
  ifelse(!asym & platy, "PMM3",
  ifelse(asym & platy, "PMM2",
  ifelse(!asym & d$g3 < g3_thr, "PMM3", "OLS")))))
}

# ---- 1. Cell table under the default rule ----------------------------------
mse_ols <- cell_mse(d$se_ols); mse_2 <- cell_mse(d$se_pmm2)
mse_3 <- cell_mse(d$se_pmm3)
best <- pmin(mse_ols, mse_2, mse_3)
def_pick <- rule()
stopifnot(identical(def_pick, d$pick))   # the reimplementation is pmm_dispatch()
mse_def <- cell_mse(se_of(def_pick))
share <- function(m) tapply(d$pick == m, cell_key, mean)[cells$key]

tab <- data.frame(
  law = cells$law, n = cells$n,
  g2_theory = pop$g2[match(cells$law, pop$law)],
  g3_theory = pop$g3[match(cells$law, pop$law)],
  pmm2_vs_ols = mse_2 / mse_ols, pmm3_vs_ols = mse_3 / mse_ols,
  best = c("OLS", "PMM2", "PMM3")[apply(cbind(mse_ols, mse_2, mse_3), 1, which.min)],
  pick_ols = share("OLS"), pick_pmm2 = share("PMM2"), pick_pmm3 = share("PMM3"),
  dispatch_vs_ols = mse_def / mse_ols, regret = mse_def / best,
  nonconv_pmm2 = tapply(!d$conv_pmm2, cell_key, mean)[cells$key],
  row.names = NULL)
tab <- tab[order(match(tab$law, pop$law), tab$n), ]
dir.create("mc_results_dispatch", showWarnings = FALSE)
utils::write.csv(tab, "mc_results_dispatch/cells_default.csv", row.names = FALSE)

cat(sprintf("EstemPMM %s @ %s, %d reps per cell\n\n", x$estempmm, x$commit,
            x$reps_per_cell))
cat("== Cells (default rule) ==\n")
print(tab, digits = 3, row.names = FALSE)

cat("\n== Default rule by n: regret and loss against OLS ==\n")
by_n <- do.call(rbind, lapply(split(tab, tab$n), function(t) data.frame(
  n = t$n[1], mean_regret = exp(mean(log(t$regret))), worst_regret = max(t$regret),
  worst_cell = t$law[which.max(t$regret)], cells_worse_than_ols = sum(t$dispatch_vs_ols > 1.02))))
print(by_n, digits = 3, row.names = FALSE)

cat("\n== D-1: symmetric heavy tails, MSE ratio to OLS ==\n")
print(tab[tab$law %in% c("t3", "t5", "laplace", "logistic"),
          c("law", "n", "pmm2_vs_ols", "pmm3_vs_ols", "pick_pmm2", "pick_pmm3",
            "dispatch_vs_ols")], digits = 3, row.names = FALSE)

# ---- 2. Threshold search: tune on odd reps, judge on even reps -------------
odd <- d$rep %% 2 == 1
score <- function(pick, rows) {
  m <- cell_mse(se_of(pick), rows)
  b <- pmin(cell_mse(d$se_ols, rows), cell_mse(d$se_pmm2, rows),
            cell_mse(d$se_pmm3, rows))
  o <- cell_mse(d$se_ols, rows)
  c(mean_log_regret = mean(log(m / b)), worst_regret = max(m / b),
    worst_vs_ols = max(m / o))
}
grid <- expand.grid(s_thr = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.7),
                    strong = c(0.7, 1, 1.5, Inf),
                    k_thr = c(-0.3, -0.5, -0.7, -0.9),
                    g2_thr = c(0.85, 0.9, 0.95, 0.98, 1),
                    g3_thr = c(-Inf, 0.8, 0.9, 0.95))
sc <- t(vapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]
  score(rule(g$s_thr, g$strong, g$k_thr, g$g2_thr, g$g3_thr), odd)
}, numeric(3)))
grid <- cbind(grid, tune = sc)
grid <- grid[order(grid$tune.mean_log_regret), ]
top <- head(grid, 10)
held <- t(vapply(seq_len(nrow(top)), function(i) {
  g <- top[i, ]
  score(rule(g$s_thr, g$strong, g$k_thr, g$g2_thr, g$g3_thr), !odd)
}, numeric(3)))
top <- cbind(top, test = held)
def_test <- score(def_pick, !odd)
utils::write.csv(grid, "mc_results_dispatch/threshold_grid_tune.csv", row.names = FALSE)
utils::write.csv(top, "mc_results_dispatch/threshold_top10_heldout.csv", row.names = FALSE)

cat("\n== Default rule on held-out (even) reps ==\n")
print(round(def_test, 4))
cat("\n== Top 10 rules by tuning-half mean log regret, with held-out scores ==\n")
print(top, digits = 4, row.names = FALSE)

# Per-cell view of the best held-out rule against the default
bi <- which.min(top$test.mean_log_regret)
bg <- top[bi, ]
bp <- rule(bg$s_thr, bg$strong, bg$k_thr, bg$g2_thr, bg$g3_thr)
ev <- !odd
cmp <- data.frame(cells[, c("law", "n")],
  default = cell_mse(se_of(def_pick), ev) / pmin(cell_mse(d$se_ols, ev), cell_mse(d$se_pmm2, ev), cell_mse(d$se_pmm3, ev)),
  best_rule = cell_mse(se_of(bp), ev) / pmin(cell_mse(d$se_ols, ev), cell_mse(d$se_pmm2, ev), cell_mse(d$se_pmm3, ev)),
  row.names = NULL)
cmp <- cmp[order(match(cmp$law, pop$law), cmp$n), ]
utils::write.csv(cmp, "mc_results_dispatch/heldout_regret_default_vs_best.csv", row.names = FALSE)
cat("\n== Held-out regret per cell: default vs best rule ==\n")
print(bg[, c("s_thr", "strong", "k_thr", "g2_thr", "g3_thr")], row.names = FALSE)
print(cmp[abs(cmp$default - cmp$best_rule) > 0.02, ], digits = 3, row.names = FALSE)

# ---- 3. Sample-size-aware rules, with a safety constraint ------------------
# The sampling SE of skewness is about sqrt(6/n) (0.45 at n = 30, 0.08 at
# n = 1000), so a fixed |gamma3| threshold is too loose in small samples and
# too strict in large ones. Family: asymmetric iff |gamma3| > z * sqrt(6/n);
# symmetric residuals go to PMM3 when the estimated g3 < g3_thr and n >= n3,
# or when gamma4 < k_thr (the platykurtic arm, any n). Constraint: on the
# tuning half no cell may be more than 10% worse than OLS.
rule_n <- function(z, g2_thr, g3_thr, n3, k_thr = -0.7) {
  asym <- abs(d$gamma3) > z * sqrt(6 / d$n)
  ifelse(asym & d$g2 < g2_thr, "PMM2",
  ifelse(asym, "OLS",
  ifelse(d$gamma4 < k_thr | (d$g3 < g3_thr & d$n >= n3), "PMM3", "OLS")))
}
grid_n <- expand.grid(z = c(1, 1.5, 2, 2.5, 3), g2_thr = c(0.9, 0.95, 0.98, 1),
                      g3_thr = c(-Inf, 0.85, 0.9, 0.95, 0.98),
                      n3 = c(50, 100, 200), k_thr = c(-0.5, -0.7, -0.9))
sc_n <- t(vapply(seq_len(nrow(grid_n)), function(i) {
  g <- grid_n[i, ]
  score(rule_n(g$z, g$g2_thr, g$g3_thr, g$n3, g$k_thr), odd)
}, numeric(3)))
grid_n <- cbind(grid_n, tune = sc_n)
safe <- grid_n[grid_n$tune.worst_vs_ols <= 1.10, ]
safe <- safe[order(safe$tune.mean_log_regret), ]
top_n <- head(safe, 10)
held_n <- t(vapply(seq_len(nrow(top_n)), function(i) {
  g <- top_n[i, ]
  score(rule_n(g$z, g$g2_thr, g$g3_thr, g$n3, g$k_thr), !odd)
}, numeric(3)))
top_n <- cbind(top_n, test = held_n)
utils::write.csv(grid_n, "mc_results_dispatch/threshold_grid_n_tune.csv", row.names = FALSE)
utils::write.csv(top_n, "mc_results_dispatch/threshold_n_top10_heldout.csv", row.names = FALSE)
cat("\n== n-aware rules, safe on the tuning half (worst vs OLS <= 1.10): top 10 ==\n")
cat(sprintf("%d of %d rules satisfy the constraint\n", nrow(safe), nrow(grid_n)))
print(top_n, digits = 4, row.names = FALSE)

bn <- top_n[1, ]
pn <- rule_n(bn$z, bn$g2_thr, bn$g3_thr, bn$n3, bn$k_thr)
best_ev <- pmin(cell_mse(d$se_ols, ev), cell_mse(d$se_pmm2, ev), cell_mse(d$se_pmm3, ev))
cmp_n <- data.frame(cells[, c("law", "n")],
  default = cell_mse(se_of(def_pick), ev) / best_ev,
  n_rule = cell_mse(se_of(pn), ev) / best_ev,
  n_rule_vs_ols = cell_mse(se_of(pn), ev) / cell_mse(d$se_ols, ev),
  row.names = NULL)
cmp_n <- cmp_n[order(match(cmp_n$law, pop$law), cmp_n$n), ]
utils::write.csv(cmp_n, "mc_results_dispatch/heldout_regret_default_vs_nrule.csv", row.names = FALSE)
cat("\n== Held-out regret by n: default vs best safe n-aware rule ==\n")
print(do.call(rbind, lapply(split(cmp_n, cmp_n$n), function(t) data.frame(
  n = t$n[1], default_mean = exp(mean(log(t$default))), n_rule_mean = exp(mean(log(t$n_rule))),
  default_worst = max(t$default), n_rule_worst = max(t$n_rule),
  n_rule_worst_vs_ols = max(t$n_rule_vs_ols)))), digits = 3, row.names = FALSE)
cat("\n== Cells where the two differ by more than 3% ==\n")
print(cmp_n[abs(cmp_n$default - cmp_n$n_rule) > 0.03, ], digits = 3, row.names = FALSE)

# ---- 4. Hybrid: clear asymmetry first, then platykurtosis ------------------
# Section 3 fixed the symmetric cells but lost PMM2 on skewed laws at n <= 50,
# because z * sqrt(6/n) is ~0.9 at n = 30. The default rule's small-n losses
# come from its "asymmetric AND platykurtic -> PMM2" arm instead. Hybrid:
#   |gamma3| > z * sqrt(6/n) and g2 < g2_thr   -> PMM2  (clear asymmetry)
#   gamma4 < k_thr                              -> PMM3  (platykurtic)
#   |gamma3| > s_fix and g2 < g2_thr            -> PMM2  (moderate asymmetry)
#   g3 < g3_thr and n >= n3                     -> PMM3  (leptokurtic arm)
#   otherwise                                   -> OLS
rule_h <- function(z, s_fix, k_thr, g2_thr, g3_thr, n3) {
  a3 <- abs(d$gamma3)
  ifelse(a3 > z * sqrt(6 / d$n) & d$g2 < g2_thr, "PMM2",
  ifelse(d$gamma4 < k_thr, "PMM3",
  ifelse(a3 > s_fix & d$g2 < g2_thr, "PMM2",
  ifelse(d$g3 < g3_thr & d$n >= n3, "PMM3", "OLS"))))
}
grid_h <- expand.grid(z = c(1.5, 2, 2.5, 3), s_fix = c(0.3, 0.5, Inf),
                      k_thr = c(-0.5, -0.7), g2_thr = c(0.95, 1),
                      g3_thr = c(-Inf, 0.9, 0.95, 0.98), n3 = c(100, 200))
sc_h <- t(vapply(seq_len(nrow(grid_h)), function(i) {
  g <- grid_h[i, ]
  score(rule_h(g$z, g$s_fix, g$k_thr, g$g2_thr, g$g3_thr, g$n3), odd)
}, numeric(3)))
grid_h <- cbind(grid_h, tune = sc_h)
safe_h <- grid_h[grid_h$tune.worst_vs_ols <= 1.10, ]
safe_h <- safe_h[order(safe_h$tune.mean_log_regret), ]
top_h <- head(safe_h, 8)
held_h <- t(vapply(seq_len(nrow(top_h)), function(i) {
  g <- top_h[i, ]
  score(rule_h(g$z, g$s_fix, g$k_thr, g$g2_thr, g$g3_thr, g$n3), !odd)
}, numeric(3)))
top_h <- cbind(top_h, test = held_h)
utils::write.csv(grid_h, "mc_results_dispatch/threshold_grid_hybrid_tune.csv", row.names = FALSE)
utils::write.csv(top_h, "mc_results_dispatch/threshold_hybrid_top_heldout.csv", row.names = FALSE)
cat("\n== Hybrid rules, safe on the tuning half: top 8 ==\n")
cat(sprintf("%d of %d rules satisfy the constraint\n", nrow(safe_h), nrow(grid_h)))
print(top_h, digits = 4, row.names = FALSE)

bh <- top_h[1, ]
ph <- rule_h(bh$z, bh$s_fix, bh$k_thr, bh$g2_thr, bh$g3_thr, bh$n3)
cmp_h <- data.frame(cells[, c("law", "n")],
  default = cell_mse(se_of(def_pick), ev) / best_ev,
  hybrid = cell_mse(se_of(ph), ev) / best_ev,
  hybrid_vs_ols = cell_mse(se_of(ph), ev) / cell_mse(d$se_ols, ev),
  row.names = NULL)
cmp_h <- cmp_h[order(match(cmp_h$law, pop$law), cmp_h$n), ]
utils::write.csv(cmp_h, "mc_results_dispatch/heldout_regret_default_vs_hybrid.csv", row.names = FALSE)
cat("\n== Held-out regret by n: default vs best safe hybrid ==\n")
print(do.call(rbind, lapply(split(cmp_h, cmp_h$n), function(t) data.frame(
  n = t$n[1], default_mean = exp(mean(log(t$default))), hybrid_mean = exp(mean(log(t$hybrid))),
  default_worst = max(t$default), hybrid_worst = max(t$hybrid),
  hybrid_worst_cell = t$law[which.max(t$hybrid)],
  hybrid_worst_vs_ols = max(t$hybrid_vs_ols)))), digits = 3, row.names = FALSE)
cat("\n== Cells where hybrid is worse than default by more than 3% ==\n")
print(cmp_h[cmp_h$hybrid > cmp_h$default * 1.03, ], digits = 3, row.names = FALSE)

# ---- 5. Targeted amendments to the default rule ---------------------------
# Sections 2-4 show no rule that beats the default everywhere: at n = 30 a
# skewed law cannot be told from a symmetric one reliably, so any rule that
# stops sending mild skewness to PMM2 loses up to 25% on skewed laws. Keep
# the default's PMM2 arms and change only what the cell table blames:
#   M1: symmetric residuals with g3 < g3_thr and n >= n3 -> PMM3
#       (leptokurtic laws, where the default never tries PMM3);
#   M2: "asymmetric AND platykurtic -> PMM2" only for |gamma3| > z*sqrt(6/n),
#       otherwise PMM3 (small-n bimodal/uniform residuals look skewed).
rule_m <- function(m1 = FALSE, m2 = FALSE, g3_thr = 0.95, n3 = 200, z = 2) {
  a3 <- abs(d$gamma3)
  asym <- a3 > 0.3
  platy <- d$gamma4 < -0.7
  clear <- a3 > z * sqrt(6 / d$n)
  ifelse(a3 > 1.0 & !(m2 & platy & !clear), "PMM2",
  ifelse(asym & platy, if (m2) ifelse(clear, "PMM2", "PMM3") else "PMM2",
  ifelse(asym & d$g2 < 0.95, "PMM2",
  ifelse(!asym & platy, "PMM3",
  ifelse(m1 & !asym & d$g3 < g3_thr & d$n >= n3, "PMM3", "OLS")))))
}
variants <- list(
  default = rule_m(),
  M1_n200_g95 = rule_m(m1 = TRUE),
  M1_n100_g95 = rule_m(m1 = TRUE, n3 = 100),
  M1_n200_g98 = rule_m(m1 = TRUE, g3_thr = 0.98),
  M2_z2 = rule_m(m2 = TRUE),
  M2_z15 = rule_m(m2 = TRUE, z = 1.5),
  M1M2 = rule_m(m1 = TRUE, m2 = TRUE),
  M1M2_n100 = rule_m(m1 = TRUE, m2 = TRUE, n3 = 100))
stopifnot(identical(variants$default, def_pick))
reg_ev <- function(p) cell_mse(se_of(p), ev) / best_ev
r_def <- reg_ev(def_pick)
summ <- do.call(rbind, lapply(names(variants), function(v) {
  p <- variants[[v]]
  r <- reg_ev(p)
  o <- cell_mse(se_of(p), ev) / cell_mse(d$se_ols, ev)
  tn <- score(p, odd)
  data.frame(rule = v, tune_mean_log_regret = tn[1],
             test_mean_log_regret = mean(log(r)), test_worst_regret = max(r),
             test_worst_vs_ols = max(o),
             cells_worse_than_default = sum(r > r_def * 1.03),
             cells_better_than_default = sum(r < r_def / 1.03))
}))
utils::write.csv(summ, "mc_results_dispatch/amendments_heldout.csv", row.names = FALSE)
cat("\n== Targeted amendments, held-out (even reps) ==\n")
print(summ, digits = 4, row.names = FALSE)
for (v in c("M1M2", "M1M2_n100")) {
  r <- reg_ev(variants[[v]])
  ch <- data.frame(cells[, c("law", "n")], default = r_def, amended = r,
                   amended_vs_ols = cell_mse(se_of(variants[[v]]), ev) / cell_mse(d$se_ols, ev))
  ch <- ch[abs(ch$amended / ch$default - 1) > 0.03, ]
  ch <- ch[order(match(ch$law, pop$law), ch$n), ]
  cat(sprintf("\n-- %s: cells changed by more than 3%% --\n", v))
  print(ch, digits = 3, row.names = FALSE)
}
