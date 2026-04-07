# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Package Does

EstemPMM is an R package implementing the Polynomial Maximization Method (PMM) for parameter estimation in linear regression and time series models with non-Gaussian errors. Two estimator variants:

- **PMM2** (S=2): For asymmetric (skewed) errors. Efficiency factor `g2 = 1 - gamma3^2 / (2 + gamma4)`.
- **PMM3** (S=3): For symmetric platykurtic errors (negative excess kurtosis). Efficiency factor `g3 = 1 - gamma4^2 / (6 + 9*gamma4 + gamma6)`.

`pmm_dispatch()` auto-selects between OLS, PMM2, and PMM3 based on residual cumulants.

## Build & Test Commands

```bash
# Check package (R CMD check equivalent)
R CMD check .

# Build tarball
R CMD build .

# Install from source
R CMD INSTALL .

# Run all tests
Rscript -e 'testthat::test_check("EstemPMM")'

# Run a single test file
Rscript -e 'testthat::test_file("tests/testthat/test-pmm2_linear.R")'

# Regenerate roxygen docs and NAMESPACE
Rscript -e 'roxygen2::roxygenise()'

# Load package for interactive dev (without installing)
Rscript -e 'devtools::load_all()'
```

## Architecture

### R source organization (`R/`)

The code splits into three layers: **PMM2**, **PMM3**, and **shared infrastructure**.

**PMM2 layer:**
- `pmm2_main.R` — `lm_pmm2()` linear regression entry point
- `pmm2_ts_main.R` — Time series entry points: `ar_pmm2()`, `ma_pmm2()`, `arma_pmm2()`, `arima_pmm2()`, `ts_pmm2()`
- `pmm2_common.R` — Core PMM2 math: moment computation, variance factor `g2`, coefficient adjustment
- `pmm2_classes.R` — S4 classes: `PMM2fit`, `BasePMM2`, `ARPMM2`, `MAPMM2`, `ARMAPMM2`, `ARIMAPMM2`
- `pmm2_ts_methods.R` — S4 method implementations (`summary`, `coef`, `residuals`, `plot`, `predict`, `AIC`) for time series PMM2 classes
- `pmm2_inference.R` — Bootstrap-based inference (`pmm2_inference()`, `ts_pmm2_inference()`)
- `pmm2_unified.R` — Unified interface for PMM2 nonlinear estimators
- `sarimax_wrapper.R` — Seasonal models: `sar_pmm2()`, `sma_pmm2()`, `sarma_pmm2()`, `sarima_pmm2()`

**PMM3 layer:**
- `pmm3_main.R` — `lm_pmm3()` linear regression entry point
- `pmm3_ts_main.R` — Time series entry points: `ar_pmm3()`, `ma_pmm3()`, `arma_pmm3()`, `arima_pmm3()`
- `pmm3_solver.R` — Core PMM3 math: 6th-order cumulant computation, variance factor `g3`
- `pmm3_utils.R` — `pmm3_variance_factor()`, `pmm_gamma6()`, `pmm_kurtosis()`
- `pmm3_classes.R` / `pmm3_ts_classes.R` — S4 classes: `PMM3fit`, `TS3fit` and subclasses
- `pmm3_ts_methods.R` — S4 methods for PMM3 time series classes

**Shared:**
- `pmm3_dispatch.R` — `pmm_dispatch()` auto-selector; computes cumulants, g2/g3 factors, recommends method
- `pmm2_utils.R` — `compute_moments()`, `pmm_skewness()`, `pmm_kurtosis()`, `test_symmetry()`
- `pmm2_monte_carlo.R` — `pmm2_monte_carlo_compare()` for simulation studies
- `data.R` — Bundled dataset documentation (`auto_mpg`)

### S4 class hierarchy

All model fits use S4 classes with standard generics (`summary`, `coef`, `residuals`, `fitted`, `plot`, `predict`, `AIC`):
- Linear: `PMM2fit`, `PMM3fit`
- Time series PMM2: `BasePMM2` -> `ARPMM2`, `MAPMM2`, `ARMAPMM2`, `ARIMAPMM2` (also `SARPMM2`, `SMAPMM2`, `SARMAPMM2`, `SARIMAPMM2`)
- Time series PMM3: `TS3fit` -> `ARPMM3`, `MAPMM3`, `ARMAPMM3`, `ARIMAPMM3`

### Key mathematical pattern

Both PMM2 and PMM3 follow the same workflow:
1. Fit initial OLS/MLE model
2. Compute residual cumulants (skewness for PMM2, kurtosis+6th-order for PMM3)
3. Calculate efficiency factor (g2 or g3); if >= 1, fall back to OLS
4. Adjust coefficients using the PMM correction formula
5. Return fit object with both OLS and PMM estimates

### Tests

Tests use `testthat` (edition 3). Key test files mirror the source structure: `test-pmm2_linear.R`, `test-pmm2_ts.R`, `test-pmm3_linear.R`, `test-pmm3_ts.R`, `test-pmm3_dispatch.R`, `test-seasonal-models.R`.

### Dev scripts

`dev_scripts/` contains Monte Carlo simulation runners and diagnostic tools (not part of the package).
