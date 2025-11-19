# EstemPMM News

## Version 0.1.4 (Development)

### New Features

- **EstemPMM-style PMM2 Estimator for MA/SMA Models** - Advanced parameter estimation for moving average components
  - New `ma_method` argument in `sarima_pmm2()` with options `"mle"` (default) and `"pmm2"`
  - `estpmm_style_ma()` - PMM2 estimator for pure MA(q) models using CSS residuals as fixed regressors
  - `estpmm_style_sma()` - PMM2 estimator for pure SMA(Q) models
  - **`estpmm_style_ma_sma()` - PMM2 estimator for mixed MA+SMA models** ⭐ NEW
  - Full support for MA(q)+SMA(Q) combinations in `sarima_pmm2()` with `ma_method="pmm2"`
  - Expected 20-45% MSE reduction for MA/SMA parameters under asymmetric innovation distributions
  - Implemented in `R/pmm2_ma_estimator.R` module with complete helper functions
  - Comprehensive unit tests (35 total) in `tests/testthat/test-ma-pmm2.R`
  - Addresses limitations identified in Monte Carlo simulations for MA parameter estimation
  - Full backward compatibility - default behavior unchanged

### Bug Fixes

- **Fixed Function Name Conflicts** - Removed obsolete versions of `ma_solve_pmm2`, `ma_compute_innovations`, `sma_compute_innovations`, and `sma_build_design` from `pmm2_ts_main.R` that were overwriting new implementations in `pmm2_ma_estimator.R`
- **Fixed Seasonal Period Validation** - `sarima_pmm2()` now correctly allows `s=1` when no seasonal components (P=0, D=0, Q=0) are specified
- **Corrected ts Object Handling** - MA/SMA estimators now properly convert `ts` objects to numeric vectors before arithmetic operations

## Version 0.1.3 (2025-11-13)

### Documentation

- Expanded both `README.md` and `README_uk.md` with organized function tables, seasonal SAR/SMA workflow examples, and refreshed Monte Carlo efficiency results so new users can discover the seasonal functionality faster.
- Added Part 8 to `vignettes/pmm2_time_series.Rmd`, walking through `sar_pmm2()`/`sma_pmm2()` usage, convergence tips, and practical guidance for seasonal datasets.
- Captured the seasonal-model release summary directly in `NEWS.md`, keeping the changelog aligned with the refreshed documentation.
- Added CRAN-facing housekeeping: refreshed `cran-comments.md`, `CRAN_CHECK_INSTRUCTIONS.md`, and `CRAN_SUBMISSION_CHECKLIST.md`, plus README sections on rebuilding docs, reproducing Monte Carlo studies, and running `R CMD check --as-cran`.


## Version 0.1.2 (2025-11-13)

### New Features

- **Seasonal Autoregressive Models (`sar_pmm2()`)** - Full implementation of SAR(p,P)_s models for seasonal time series
  - Supports arbitrary seasonal periods (e.g., 12 for monthly, 4 for quarterly data)
  - Multiple estimation methods: PMM2, OLS
  - Demonstrated 20-30% variance reduction with asymmetric innovations
  - Full integration with S4 class system (`SARPMM2` class)

- **Seasonal Moving Average Models (`sma_pmm2()`)** - Complete SMA(Q)_s implementation
  - Flexible seasonal lag specification
  - CSS and PMM2 estimation methods
  - Empirically validated with 500 Monte Carlo replications
  - Achieved 34.1% variance reduction (exceeding theoretical predictions)
  - Robust convergence and computational efficiency

- **Enhanced Comparison Functions**
  - `compare_sar_methods()` - Compare SAR estimation approaches
  - `compare_ts_methods()` - Universal wrapper now supports SAR and SMA models

- **Documentation and Validation**
  - Added comprehensive SAR/SMA documentation in `docs/` directory
  - Monte Carlo validation reports with detailed efficiency metrics
  - Ukrainian language analysis reports
  - Updated both English and Ukrainian READMEs with seasonal models

### Bug Fixes

- **Fixed `predict()` method for PMM2fit class** - The prediction method now correctly handles arbitrary variable names instead of requiring hardcoded "x1", "x2" names. The method now uses general matrix multiplication approach (`X %*% coefficients`) that works with any variable naming convention.
- **Improved coefficient name matching** - Enhanced logic to ensure coefficient names always match design matrix columns, with automatic reordering when necessary.
- **Fixed SAR mean iterations display** - Corrected `sprintf()` call to properly show mean iteration count in comparison output
- **Fixed Seasonal Model Residuals** - `sar_pmm2`, `sarma_pmm2`, and `sarima_pmm2` now correctly pad residuals with zeros (instead of `NA`) to match the original series length, ensuring compatibility with standard diagnostic tools.
- **Fixed S4 Class Definitions** - Reordered class and method definitions in `pmm2_classes.R` to prevent "no definition for class" warnings during package loading.
- **Corrected Multiplicative SAR Specification** - Updated tests to correctly expect 3 coefficients (AR, SAR, Interaction) for multiplicative SAR(1)x(1) models.

### Improvements

- **More robust prediction algorithm** - Simplified prediction code by removing hardcoded special cases and using a unified matrix multiplication approach for all scenarios.
- **Better error messages** - Added clearer error messages when design matrix and coefficients don't match.
- **Enhanced `.gitignore`** - Added `test_results/` directory to version control exclusions

## Version 0.1.1 (2025-10-23)

### Maintenance

- Updated `DESCRIPTION` (latest release date, Suggests list for packages used in the demos).
- Verified the package with `R CMD check --as-cran` (now warning-free after installing `qpdf`).
- Regenerated vignettes (HTML and tangled `.R`) and included them in `inst/doc` for distribution.
- Updated `.Rbuildignore` and `.gitignore`, keeping only files required for CRAN.

## Version 0.1.0 (2025-01-15)

### Initial Release: PMM2 Foundation

**New Features:**
- `lm_pmm2()` - Linear regression estimation using Polynomial Maximization Method (S=2)
- `ar_pmm2()` - Autoregressive (AR) time series modeling with PMM2
- `ma_pmm2()` - Moving Average (MA) time series modeling with PMM2
- `arma_pmm2()` - ARMA time series modeling with PMM2
- `arima_pmm2()` - ARIMA time series modeling with PMM2
- `pmm2_inference()` - Bootstrap inference for linear models
- `ts_pmm2_inference()` - Bootstrap inference for time series models
- Statistical utilities: `pmm_skewness()`, `pmm_kurtosis()`, `compute_moments()`
- Comparison functions: `compare_with_ols()`, `compare_ts_methods()`, `compare_ar_methods()`, `compare_ma_methods()`, `compare_arma_methods()`, `compare_arima_methods()`

**S4 Classes:**
- `PMM2fit` - Results container for linear regression models
- `TS2fit` - Base class for time series results
- `ARPMM2`, `MAPMM2`, `ARMAPMM2`, `ARIMAPMM2` - Specialized time series result classes

**Methods:**
- `summary()` - Model summary statistics
- `coef()` - Extract coefficients
- `fitted()` - Fitted values
- `predict()` - Predictions for new data
- `residuals()` - Model residuals
- `plot()` - Diagnostic plots

**Documentation:**
- Comprehensive Roxygen2 documentation for all exported functions
- README with theoretical background and basic usage examples
- Demonstration script `pmm2_demo_runner.R` showing practical applications

### Package Architecture

**Module Organization:**
- `R/pmm2_main.R` - Primary PMM2 fitting functions
- `R/pmm2_classes.R` - S4 class definitions
- `R/pmm2_utils.R` - Utility functions for moment computation and optimization
- `R/pmm2_ts_design.R` - Time series design matrix construction

**Dependencies:**
- Core: `methods`, `stats`, `graphics`, `utils`
- Optional: `MASS` (for advanced statistical functions, available in Suggests)

**Quality Assurance:**
- Unit tests covering core PMM2 functionality
- Edge case handling for numerical stability
- Convergence diagnostics and warnings

### Known Limitations

- PMM2 only (S=2 order polynomial) - higher orders not yet implemented
- Single-stage estimation (no multi-stage procedures)
- Time series models assume stationarity for AR, MA components
- ARIMA differencing handled via preprocessing, not integrated into core algorithm

### Roadmap for Future Versions

**0.2.0 (PMM3 Ready Architecture):**
- PMM3 implementation (S=3 polynomial methods)
- Refactored base classes supporting method extensibility
- Vignette documentation with practical use cases
- Enhanced bootstrap procedures for small samples
- GitHub Actions CI/CD integration

**0.3.0 (Advanced Methods):**
- Adaptive PMM order selection
- Robust variance estimation
- Model selection criteria (AIC/BIC for PMM)
- Generalized Linear Models (GLM) with PMM

**1.0.0 (Stable API):**
- API stabilization and backward compatibility guarantee
- Extended performance benchmarks
- Specialized applications (econometrics, biostatistics)

### Citation

If you use EstemPMM in your research, please cite the relevant publications:

**For Linear Regression (lm_pmm2):**
Zabolotnii S., Warsza Z.L., Tkachenko O. (2018) Polynomial Estimation of Linear
Regression Parameters for the Asymmetric PDF of Errors. In: Szewczyk R.,
Zieliński C., Kaliczyńska M. (eds) Automation 2018. AUTOMATION 2018. Advances in
Intelligent Systems and Computing, vol 743. Springer, Cham.
https://doi.org/10.1007/978-3-319-77179-3_75

**For Autoregressive Models (ar_pmm2):**
Zabolotnii S., Tkachenko O., Warsza Z.L. (2022) Application of the Polynomial
Maximization Method for Estimation Parameters of Autoregressive Models with
Asymmetric Innovations. In: Szewczyk R., Zieliński C., Kaliczyńska M. (eds)
Automation 2022. AUTOMATION 2022. Advances in Intelligent Systems and Computing,
vol 1427. Springer, Cham. https://doi.org/10.1007/978-3-031-03502-9_37

**For Moving Average Models (ma_pmm2):**
Zabolotnii S., Tkachenko O., Warsza Z.L. (2023) Polynomial Maximization Method
for Estimation Parameters of Asymmetric Non-gaussian Moving Average Models. In:
Szewczyk R., et al. (eds) Automation 2023. AUTOMATION 2023. Lecture Notes in
Networks and Systems, vol 630. Springer, Cham.

### Technical Notes

**Algorithm Stability:**
- Regularization parameter automatically adjusted for ill-conditioned systems
- Step size limiting prevents divergence in optimization
- Convergence history tracking for diagnostics

**Numerical Considerations:**
- Moment estimation uses robust methods to handle outliers
- Design matrices constructed with numerical stability in mind
- NA/Inf values detected and handled appropriately
