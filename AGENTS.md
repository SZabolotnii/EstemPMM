# AGENTS.md

This file provides guidance to AI coding agents working with code in this repository. It is the single source of instructions: Codex reads it directly, and Claude Code reads it through `CLAUDE.md`, which is a symlink to this file. Edit this file, never `CLAUDE.md`.

## What This Package Does

EstemPMM is an R package implementing the Polynomial Maximization Method (PMM) for parameter estimation in linear regression and time series models with non-Gaussian errors. Two estimator variants:

- **PMM2** (S=2): For asymmetric (skewed) errors. Efficiency factor `g2 = 1 - gamma3^2 / (2 + gamma4)`.
- **PMM3** (S=3): For symmetric platykurtic errors (negative excess kurtosis). Efficiency factor `g3 = 1 - gamma4^2 / (6 + 9*gamma4 + gamma6)`.

`pmm_dispatch()` auto-selects between OLS, PMM2, and PMM3 based on residual cumulants.

This package is the portfolio's reference PMM engine: downstream projects (`EstemPMM-GMDH`, `PMM2-ARIMA`, `PMM3_vs_MLE-TN`, `sdr-kunchenko` via rpy2) depend on it, and the `JSS` software paper calls the exported functions directly. Treat exported names, argument names and defaults as a public API.

## Release state

Version 0.5.0 was submitted to CRAN (see `cran-comments.md`, `NEWS.md`). Do not upload to CRAN without explicit authorization. Release tarballs `EstemPMM_*.tar.gz` are kept in the repo root and are build-ignored.

## Build & Test Commands

```bash
# Check package (CRAN-level check is run against the built tarball)
R CMD build . && R CMD check --as-cran EstemPMM_<version>.tar.gz

# Install from source
R CMD INSTALL .

# Run all tests against the source tree
Rscript -e 'devtools::test()'

# Run a single test file
Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-pmm2_linear.R")'

# Regenerate roxygen docs and NAMESPACE (NAMESPACE and man/ are generated: never hand-edit)
Rscript -e 'roxygen2::roxygenise()'
```

CI: `.github/workflows/R-CMD-check.yaml` runs R CMD check on macOS/Windows/Ubuntu (release, devel, oldrel-1).

**`DESCRIPTION` has an explicit `Collate:` field.** A new file in `R/` must be added there, or it is not loaded. Class files declare ordering via roxygen `@include`.

## Architecture

### Two API layers

- **Unified API** (`pmm_unified_api.R`, since 0.4.0): `pmm_lm()`, `pmm_ar()`, `pmm_ma()`, `pmm_arma()`, `pmm_arima()`, `pmm_sarima()`. These are thin wrappers that choose the order through a `method` argument. `pmm_lm(method = "auto")` consults `pmm_dispatch()`.
- **Low-level fitters**: `lm_pmm2/3()`, `ar_/ma_/arma_/arima_pmm2/3()`, `ts_pmm2()`, and the seasonal `sar_/sma_/sarma_/sarima_pmm2()` (`sarimax_wrapper.R`). They are still exported, emit no deprecation warning, and the JSS manuscript uses them. Do not deprecate or remove them without flagging it.

### Source layout (`R/`)

Files split into PMM2 (`pmm2_*`), PMM3 (`pmm3_*`) and shared (`pmm_*`) code. Non-obvious points:
- Core math: `pmm2_common.R` (PMM2 Newton iteration, `g2`) and `pmm3_solver.R` (6th-order cumulants, `g3`).
- `pmm3_dispatch.R` holds the shared `pmm_dispatch()` (S3 class `PMMdispatch`), despite its `pmm3_` prefix.
- TS design matrices: `pmm2_ts_design.R`. The MA innovation recursion: `pmm2_ma_estimator.R` and `pmm2_ts_recursive.R`.
- Inference: `pmm2_inference.R` (bootstrap), plus S4 `vcov`/`confint`/`logLik`/`nobs` methods. The PMM3 vcov derivation is in `inst/notes/pmm3_vcov_derivation.md`.
- `pmm2_unified.R`: generic nonlinear PMM2 (iterative and one-step), which needs `numDeriv`. `optimized_direct_pmm2.R` is **deprecated**: it lost badly in Monte Carlo, so do not revive it.

### MA solver (0.5.0)

`ts_pmm2()`, `ma_pmm2()`, `arma_pmm2()` and `arima_pmm2()` accept `ma_solver = c("linearized", "recursive")`.
- `"linearized"` (the **default**, kept so 0.4.x results reproduce bit-for-bit) builds the MA design matrix once from the CSS residuals. That costs a factor of `1 - theta^2` in asymptotic variance.
- `"recursive"` recomputes the innovations and score regressors by the exact ARMA recursion at every candidate parameter. It has a closed-form sandwich `vcov`/`confint`.
- The seasonal family silently keeps the linearised path. `sarima_pmm2()` has an **unrelated** `ma_method` argument (`"mle"`/`"pmm2"`). Do not confuse it with `ma_solver`.

### S4 class hierarchy (`pmm_base_classes.R`)

```
PMMfit (virtual root: coefficients, residuals, convergence, iterations, call)
├── BasePMM2 (virtual; m2, m3, m4)  ── PMM2fit (linear)
├── BasePMM3 (virtual; m2, m4, m6, gamma4, gamma6, g_coefficient, kappa) ── PMM3fit
└── PMMtsfit (virtual; model_type, intercept, original_series, order)
    TS2fit ⊂ {BasePMM2, PMMtsfit} → ARPMM2, MAPMM2, ARMAPMM2, ARIMAPMM2
                                    (+ SARPMM2, SMAPMM2, SARMAPMM2, SARIMAPMM2)
    TS3fit ⊂ {BasePMM3, PMMtsfit} → ARPMM3, MAPMM3, ARMAPMM3, ARIMAPMM3
```

`show` methods are defined once in `pmm_show_methods.R` and inherited by every subclass. This hierarchy is Figure 1 of the JSS paper, so changes must stay in sync with it.

### Key mathematical pattern

Both PMM2 and PMM3 follow the same workflow:
1. Fit an initial OLS model (or CSS/MLE for time series).
2. Compute the residual cumulants: skewness for PMM2; kurtosis and 6th-order for PMM3.
3. Calculate the efficiency factor (`g2` or `g3`). If it is >= 1, fall back to OLS.
4. Adjust the coefficients using the PMM correction (Newton iteration).
5. Return a fit object with both the initial and the PMM estimates.

### Tests, data, dev scripts

- Tests use `testthat` edition 3 (`tests/testthat/`). The `test-coverage-boost*.R` files are coverage fillers. `test-v040-new-features.R` covers the unified API, `show` and `PMMdispatch`.
- Bundled datasets: `auto_mpg`, `djia2002` and `DCOILWTICO` (`data/`, documented in `R/data.R`). `data/DCOILWTICO.csv` must be a regular file, not a symlink, or a clean checkout will not build.
- `dev_scripts/`, `mc_results_*`, `test_results/`, `articles/`, `docs/` and `notes/` are build-ignored research and Monte Carlo material, not part of the package. When you add a top-level file that must not ship, also add it to `.Rbuildignore`.
