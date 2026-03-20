# CRAN Submission Comments — EstemPMM 0.3.0

## Summary
- **Submission type:** Major feature release (update from 0.1.1 on CRAN)
- **New features:** PMM3 (S=3) estimation for symmetric platykurtic errors,
  automatic method selection via `pmm_dispatch()`, PMM3 time series models
  (AR/MA/ARMA/ARIMA), two bundled real-world datasets (`auto_mpg`, `djia2002`)

## Test environments
| Platform | R version | Notes |
| --- | --- | --- |
| macOS Tahoe 26.3.1 (arm64) | 4.5.2 | `R CMD check --as-cran` |
| macOS Tahoe 26.3.1 (arm64) | 4.5.2 | Local development and testing |

## R CMD check results
```
Status: 0 ERROR, 0 WARNING, 1 NOTE
```

### Notes explained
1. **Non-standard top-level files/directories:** `CRAN_CHECK_INSTRUCTIONS.md`,
   `CRAN_SUBMISSION_CHECKLIST.md`, `README_uk.md` — supplementary documentation
   files for contributors. Does not affect package functionality.

## Reverse dependencies
None. This package has no reverse dependencies on CRAN.

## What's New in 0.3.0

### PMM3 (S=3) for Symmetric Platykurtic Errors
- `lm_pmm3()` — linear regression with PMM3 Newton-Raphson solver
- `ar_pmm3()`, `ma_pmm3()`, `arma_pmm3()`, `arima_pmm3()` — time series models
- `pmm_dispatch()` — automatic OLS/PMM2/PMM3 method selection
- `compute_moments_pmm3()`, `pmm3_variance_factor()`, `pmm_gamma6()`, `test_symmetry()`
- Standalone S4 classes: `PMM3fit`, `TS3fit`, `ARPMM3`, `MAPMM3`, `ARMAPMM3`, `ARIMAPMM3`

### New Datasets
- `auto_mpg` — UCI Auto MPG (N=398), used in 3 published PMM papers
- `djia2002` — Dow Jones daily data (N=127), used in PMM2 AR paper

### Monte Carlo Validation (M=500, n=100/200/500)
PMM3 achieves 30-70% MSE reduction for platykurtic innovations across all model types:
- Linear regression: g3 = 0.35-0.70 (vs OLS)
- AR models: MSE ratio = 0.30-0.68 (vs OLS)
- ARIMA models: MSE ratio = 0.32-0.70 (vs MLE)
- 100% convergence rate across all configurations

### Other Changes
- New vignette: "PMM3: Estimation for Symmetric Platykurtic Errors"
- `numDeriv` added to Suggests (for ARIMA PMM3 numerical Jacobian)
- Backward compatible: all existing PMM2 API unchanged
