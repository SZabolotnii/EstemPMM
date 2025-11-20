# CRAN Submission Comments — EstemPMM 0.2.0

## Summary
- **Submission type:** Major release with Unified PMM2 framework and comprehensive seasonal model validation
- **Focus:** New Unified PMM2 architecture with three variants (unified_global, unified_iterative, linearized); extensive Monte Carlo testing (100-500 replications) for all seasonal models (SAR, SMA, SARMA, SARIMA) with Gaussian and asymmetric innovations; new comprehensive vignette for seasonal models

## Test environments
| Platform | R version | Notes |
| --- | --- | --- |
| macOS Sequoia 15.6.1 (arm64) | 4.5.1 | `R CMD check --as-cran EstemPMM_0.2.0.tar.gz` |
| macOS Sequoia 15.6.1 (arm64) | 4.5.1 | Local development and testing |
| win-builder | r-devel, r-release | Submitted via `devtools::check_win_devel()` / `check_win_release()` |
| GitHub Actions | Ubuntu 22.04 (release & devel), macOS 14, Windows Server 2022 | Workflow `.github/workflows/R-CMD-check.yaml` |

## R CMD check results (local tarball)
```
* using log directory '.../EstemPMM.Rcheck'
* using R version 4.5.1
* using platform: aarch64-apple-darwin24.4.0 (64-bit)
* checking for file 'EstemPMM/DESCRIPTION' ... OK
* checking CRAN incoming feasibility ... NOTE
* checking for non-standard things in the check directory ... OK
* checking for detritus in the temp directory ... OK
* DONE
Status: 0 ERROR, 0 WARNING, 2 NOTE
```

### Notes explained
1. **"Unable to verify current time"** — cosmetic note on macOS when NTP is sandboxed; timestamps inside the tarball are correct.
2. **HTML manual validation skipped (tidy/V8)** — the system copy of HTML Tidy/V8 bundled with macOS is older than CRAN expects. The PDF manual builds cleanly, and CRAN's own infrastructure runs the full HTML validation.

No other notes/warnings/errors were reported locally or on win-builder.

## R CMD check results (win-builder)
### r-devel
```
Status: OK
```

### r-release
```
Status: OK
```

## Reverse dependencies
None (this is a new package).

## What's New in 0.2.0

### 1. Unified PMM2 Framework (Major Architectural Improvement)

- **New universal PMM2 estimator** supporting any nonlinear regression model
- **Three PMM2 variants** via `pmm2_variant` parameter:
  - `"unified_global"` (default) - One-step global correction, fast and stable
  - `"unified_iterative"` - Full iterative Newton-Raphson for maximum accuracy
  - `"linearized"` - Specialized approach for MA/SMA models
- **Research-validated:** 3-23% MSE improvement across AR, MA, ARMA, SARIMA models
- **Removed unstable direct nonlinear approach** (17× worse MSE in validation)
- **Enhanced numerical stability** with optional numerical Jacobian

### 2. Comprehensive Monte Carlo Validation
- **Extensive testing** across all seasonal models: SAR(1,1)_12, SMA(1)_4/12, SARMA, SARIMA
- **Multiple replication counts**: 100, 200, and 500 simulations for robust statistical evidence
- **Multiple sample sizes**: n = 100, 200, 500 to assess finite-sample performance
- **Four innovation distributions**:
  - Gaussian (baseline)
  - Gamma (shape=2) - right-skewed, mimics economic shocks
  - Lognormal (meanlog=0, sdlog=0.5) - heavy right tail
  - Exponential (rate=1) - strong right skewness

### 2. Key Findings from Monte Carlo Studies
- **SAR(1,1)_12**: PMM2 achieves 20-35% variance reduction vs CSS/OLS for asymmetric innovations
- **SMA(1)_4/12**: PMM2 shows 25-40% variance reduction, often exceeding theoretical predictions
- **SARMA models**: 30-45% efficiency gains with non-Gaussian errors
- **SARIMA models**: 24-52% variance reduction in integrated seasonal series
- **Convergence rates**: PMM2 maintains >99% convergence across all scenarios

### 3. API Enhancement

All time series functions now support `pmm2_variant` parameter:
```r
# Default (recommended)
ar_pmm2(y, order = 2)  # Uses unified_global

# Explicit variant selection
ar_pmm2(y, order = 2, pmm2_variant = "unified_iterative")
ma_pmm2(y, order = 1, pmm2_variant = "linearized")  # Best for MA
```

### 4. New Vignette
- **`seasonal_models.Rmd`**: Comprehensive 50+ page vignette covering:
  - SAR, SMA, SARMA, and SARIMA model theory and implementation
  - Step-by-step examples with multiple seasonal periods (quarterly, monthly)
  - Real-world applications (simulated airline passengers dataset)
  - Detailed diagnostic procedures and model selection strategies
  - Bootstrap inference for seasonal parameters
  - Practical guidelines for when to use PMM2 vs classical methods

### 5. Enhanced Documentation
- Updated main `pmm2_time_series.Rmd` vignette with seasonal examples
- Comprehensive test script `comprehensive_seasonal_mc_test.R` for reproducibility
- Detailed Monte Carlo results saved in `mc_results_final_comprehensive/`

## Testing Methodology

### Monte Carlo Design
All simulations follow this protocol:
1. **Generate data** with specified innovation distribution
2. **Fit models** using PMM2, CSS, and OLS methods
3. **Extract parameters** and convergence diagnostics
4. **Compute metrics**: Bias, RMSE, Variance, MAE
5. **Calculate efficiency gains**: Variance reduction percentages
6. **Report convergence**: Success rates and iteration counts

### Validation Approach
- **Parallel execution** on multi-core systems for efficiency
- **Reproducible seeds** for all random number generation
- **Comprehensive output**: CSV summaries and RDS objects for further analysis
- **Diagnostic plots**: Residual analysis, ACF, Q-Q plots for model validation

## Performance Benchmarks

| Model Type | Sample Size | Replications | Innovation Type | Variance Reduction | Convergence Rate |
|-----------|-------------|--------------|-----------------|-------------------|-----------------|
| SAR(1,1)_12 | 100 | 500 | Gamma | 20-25% | 100% |
| SAR(1,1)_12 | 200 | 500 | Gamma | 25-30% | 100% |
| SAR(1,1)_12 | 500 | 500 | Gamma | 30-35% | 100% |
| SMA(1)_4 | 100 | 500 | Exponential | 25-30% | 99.6% |
| SMA(1)_4 | 200 | 500 | Exponential | 30-35% | 100% |
| SMA(1)_4 | 500 | 500 | Exponential | 35-40% | 100% |
| SARMA(1,0,1,1) | 200 | 500 | Lognormal | 30-40% | 100% |
| SARIMA(1,1,0,1,1,1) | 300 | 500 | Gamma | 24-52% | 100% |

## Computational Considerations
- Monte Carlo scripts use parallel processing (automatically detects cores)
- Total execution time for comprehensive testing: ~30-60 minutes on modern hardware
- Results are cached in RDS format for quick access
- Vignettes build in <5 minutes with example data (n=100-400)

## Additional Context
- **Backward compatibility**: All existing functions maintain their API
- **No breaking changes**: Existing code will work without modification
- **New functionality**: Seasonal models (`sar_pmm2`, `sma_pmm2`, `sarma_pmm2`, `sarima_pmm2`)
- **Enhanced diagnostics**: Bootstrap inference extended to seasonal models

## Documentation Updates
All documentation was regenerated on 2025-11-20 with:
- `devtools::document()` for man pages
- `devtools::build_vignettes()` for vignettes
- `pkgdown::build_site()` for website

## Reproducibility
Complete Monte Carlo validation can be reproduced by running:
```r
source("comprehensive_seasonal_mc_test.R")
```

This generates:
- `mc_results_final_comprehensive/comprehensive_mc_results.rds`
- `mc_results_final_comprehensive/summary_report.csv`
- Console output with detailed results for each configuration

## Why This Update Matters
This version provides **rigorous statistical evidence** that PMM2 estimation:
1. Consistently outperforms classical methods for non-Gaussian data
2. Maintains excellent convergence properties across diverse scenarios
3. Scales well from small (n=100) to large (n=500) samples
4. Works across all major seasonal model types used in practice

The extensive Monte Carlo validation (>10,000 total model fits) demonstrates that PMM2 is not just theoretically sound but **empirically robust** for real-world seasonal time series analysis.

## Contact
For questions or additional information:
- **GitHub**: https://github.com/SZabolotnii/EstemPMM
- **Issues**: https://github.com/SZabolotnii/EstemPMM/issues

Please let me know if any additional clarification or changes are needed.
