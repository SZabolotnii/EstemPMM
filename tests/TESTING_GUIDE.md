# EstemPMM Testing Guide

## Overview

This guide describes the testing infrastructure for the EstemPMM package, with special focus on the seasonal models introduced in version 0.1.2+.

## Test Structure

### Test Files

```
tests/
├── testthat.R                      # Main test runner
└── testthat/
    ├── test-pmm2_linear.R         # Linear regression tests (98 lines)
    ├── test-pmm2_ts.R             # Basic time series tests (64 lines)
    ├── test-pmm2_inference.R      # Bootstrap inference tests (35 lines)
    ├── test-pmm2_methods.R        # S4 methods tests (55 lines)
    ├── test-pmm2_utils.R          # Utility function tests (46 lines)
    ├── test-monte-carlo.R         # Monte Carlo tests (59 lines)
    └── test-seasonal-models.R     # **NEW** Seasonal models tests (500+ lines)
```

### Coverage by Functionality

| Functionality | Test File | Lines | Coverage Target |
|--------------|-----------|-------|-----------------|
| Linear Regression (`lm_pmm2`) | `test-pmm2_linear.R` | 98 | >90% |
| AR/MA/ARMA Models | `test-pmm2_ts.R` | 64 | >85% |
| **SAR Models** | `test-seasonal-models.R` | ~150 | >85% |
| **SMA Models** | `test-seasonal-models.R` | ~150 | >85% |
| **SARMA Models** | `test-seasonal-models.R` | ~100 | >80% |
| **SARIMA Models** | `test-seasonal-models.R` | ~100 | >80% |
| Bootstrap Methods | `test-pmm2_inference.R` | 35 | >75% |
| S4 Methods | `test-pmm2_methods.R` | 55 | >90% |
| Monte Carlo | `test-monte-carlo.R` | 59 | >70% |

## Running Tests

### Quick Test Run

```r
# From R console in package directory
library(testthat)
library(EstemPMM)

# Run all tests
test_check("EstemPMM")
```

### Run Specific Test File

```r
# Run only seasonal model tests
test_file("tests/testthat/test-seasonal-models.R")

# Run only linear regression tests
test_file("tests/testthat/test-pmm2_linear.R")
```

### Run with devtools

```r
library(devtools)

# Run all tests
test()

# Run with coverage
test_coverage()
```

### Command Line

```bash
# From package root directory
R CMD check --no-manual EstemPMM_*.tar.gz

# Or with Rscript
Rscript -e "devtools::test()"
```

## Coverage Analysis

### Running Coverage Analysis

#### Method 1: Using covr package

```r
library(covr)

# Full package coverage
cov <- package_coverage()
print(cov)

# Generate HTML report
report(cov)

# Get percentage
percent_coverage(cov)
```

#### Method 2: Using provided script

```bash
# Make script executable
chmod +x run_coverage_analysis.R

# Run analysis
Rscript run_coverage_analysis.R
```

This will generate:
- Console summary of coverage by file
- `coverage_report.html` - Interactive HTML report
- `coverage_data.rds` - Raw coverage data

### Interpreting Coverage Results

**Coverage Targets:**
- ✅ **Excellent**: >90% - Comprehensive testing
- ✓ **Good**: 80-90% - Adequate coverage, minor gaps
- ⚠️ **Fair**: 70-80% - Needs improvement
- ❌ **Poor**: <70% - Significant gaps

**Focus Areas for Coverage:**
1. **Core functionality** (exported functions): Target >90%
2. **S4 methods** (coef, summary, plot, etc.): Target >90%
3. **Internal helpers**: Target >70%
4. **Error handling**: Target >80%

## New Seasonal Model Tests

### test-seasonal-models.R

This comprehensive test file covers all new seasonal functionality:

#### SAR (Seasonal AutoRegressive) Tests
- ✅ Basic SAR(0,P)_s model fitting
- ✅ Multiplicative SAR(p,P)_s models
- ✅ Coefficient bounds and validity
- ✅ Different seasonal periods (4, 12, etc.)
- ✅ Mean/intercept handling
- ✅ S4 class structure (SARPMM2)

#### SMA (Seasonal Moving Average) Tests
- ✅ Basic SMA(Q)_s model fitting
- ✅ Higher-order models (Q>1)
- ✅ Method comparison (PMM2 vs CSS)
- ✅ Convergence parameters
- ✅ Innovation slot validation
- ✅ S4 class structure (SMAPMM2)

#### SARMA (Combined Seasonal) Tests
- ✅ Full SARMA(p,P,q,Q)_s specification
- ✅ Pure seasonal models
- ✅ Order slot validation
- ✅ Mean parameter handling
- ✅ Residual properties
- ✅ S4 class structure (SARMAPMM2)

#### SARIMA (With Differencing) Tests
- ✅ Differencing orders (d, D)
- ✅ Full model specification
- ✅ Non-stationary data handling
- ✅ Mean parameter auto-determination
- ✅ S4 class structure (SARIMAPMM2)

#### S4 Methods Tests
- ✅ `coef()` extraction
- ✅ `residuals()` computation
- ✅ `fitted()` values
- ✅ `summary()` display
- ✅ Inheritance from TS2fit

#### Edge Cases
- ✅ Short time series
- ✅ Large seasonal periods
- ✅ Constant series
- ✅ Convergence issues

#### Integration Tests
- ✅ Comparison with `stats::arima()`
- ✅ Class inheritance verification
- ✅ Original series storage
- ✅ Moment statistics computation

## Test Data Generation

### Simulating Seasonal Data

```r
# SAR(1)_12 data
simulate_sar <- function(n = 120, phi_s = 0.6, s = 12) {
  x <- numeric(n)
  innovations <- rnorm(n, sd = 1)
  for (t in (s+1):n) {
    x[t] <- phi_s * x[t - s] + innovations[t]
  }
  return(x)
}

# SMA(1)_12 data
simulate_sma <- function(n = 120, theta_s = 0.5, s = 12) {
  innovations <- rnorm(n, sd = 1)
  x <- numeric(n)
  for (t in 1:n) {
    if (t > s) {
      x[t] <- innovations[t] + theta_s * innovations[t - s]
    } else {
      x[t] <- innovations[t]
    }
  }
  return(x)
}
```

## Adding New Tests

### Template for New Test

```r
test_that("descriptive test name", {
  # Setup
  set.seed(123)
  n <- 100
  x <- rnorm(n)

  # Execute
  fit <- your_function(x, ...)

  # Assert
  expect_s4_class(fit, "ExpectedClass")
  expect_length(fit@coefficients, expected_length)
  expect_true(fit@convergence)
  expect_true(is.numeric(fit@residuals))
})
```

### Best Practices

1. **Use descriptive test names** that explain what is being tested
2. **Set seeds** for reproducibility (`set.seed()`)
3. **Test one thing per test** - atomic tests are easier to debug
4. **Use appropriate expectations**:
   - `expect_s4_class()` for class validation
   - `expect_length()` for vector lengths
   - `expect_true()/expect_false()` for logical assertions
   - `expect_equal()` for numeric comparisons (with tolerance)
   - `expect_error()` for error cases
5. **Comment complex test logic**
6. **Keep tests fast** (<1 second per test when possible)

## Continuous Integration

### GitHub Actions

The package includes `.github/workflows/R-CMD-check.yaml` which:
- Runs on push to `main`, `master`, or `claude/**` branches
- Tests on multiple platforms (Ubuntu, Windows, macOS)
- Tests on multiple R versions (release, devel, oldrel)
- Generates coverage reports

### Local Pre-commit Checks

Before committing:

```bash
# Run tests
Rscript -e "devtools::test()"

# Check coverage
Rscript run_coverage_analysis.R

# Run R CMD check
R CMD build .
R CMD check --as-cran EstemPMM_*.tar.gz
```

## Coverage Goals by Version

| Version | Overall Coverage | Seasonal Models | Core Functions |
|---------|------------------|-----------------|----------------|
| 0.1.0 | 75% | N/A | 85% |
| 0.1.2 | 70% | 60% (new) | 85% |
| **0.1.3 (current)** | **>80%** | **>85%** | **>90%** |
| 0.2.0 (planned) | >85% | >90% | >90% |

## Troubleshooting

### Tests Fail Locally but Pass on CI

- Check R version compatibility
- Ensure all Suggests packages are installed
- Check for platform-specific issues

### Coverage Report Generation Fails

```r
# Install required packages
install.packages(c("covr", "DT", "htmltools"))

# Try generating simpler report
covr::report(package_coverage(), browse = FALSE)
```

### Tests Take Too Long

- Reduce Monte Carlo replications in tests
- Use `skip_on_cran()` for slow tests
- Consider using `test_that("...", { skip("Slow test") })`

## Resources

- **testthat documentation**: https://testthat.r-lib.org/
- **covr documentation**: https://covr.r-lib.org/
- **R Packages book (Testing)**: https://r-pkgs.org/testing-basics.html
- **CRAN Testing Policies**: https://cran.r-project.org/web/packages/policies.html

## Contact

For questions about testing:
- Open an issue: https://github.com/SZabolotnii/EstemPMM/issues
- Email: zabolotniua@gmail.com

---

*Last updated: 2025-11-14*
*Package version: 0.1.3*
