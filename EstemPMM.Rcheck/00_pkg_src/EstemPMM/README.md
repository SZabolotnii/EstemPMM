# EstemPMM: Polynomial Maximization Method for Regression Analysis

[![R](https://img.shields.io/badge/R-%3E%3D%204.0.0-blue)](https://cran.r-project.org/)
[![License: GPL-3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://opensource.org/licenses/GPL-3.0)

## Overview

The EstemPMM package implements the Polynomial Maximization Method (PMM) for estimating linear regression parameters in cases where the error distribution differs from normal, particularly when it has an asymmetric character.

PMM allows obtaining parameter estimates with lower variance compared to the classical Ordinary Least Squares (OLS) method, especially when the error distribution has significant asymmetry.

## Theoretical Background

PMM uses polynomials of degree S for parameter estimation. When S=1, PMM estimates coincide with OLS estimates. When S=2 and there is asymmetry in the error distribution, PMM can provide reduced variance of estimates.

The PMM framework was originally developed by Yu. P. Kunchenko, who formulated the polynomial maximization approach for estimating distribution parameters [Kunchenko & Lega, 1992].

The theoretical coefficient of variance reduction for S=2 is calculated by the formula:

```
g = 1 - c3^2 / (2 + c4)
```

where `c3` is the skewness coefficient and `c4` is the kurtosis coefficient.

## PMM2 Variants (version 0.2.0)

EstemPMM offers **three implementation variants of PMM2**, optimized for different scenarios:

### 1. Unified Global (default, recommended)

```r
ar_pmm2(y, order = 2, pmm2_variant = "unified_global")  # default
```

**Characteristics:**
- ✅ One-step correction after classical estimation (MLE/CSS)
- ✅ Fast (~50% faster than iterative)
- ✅ Stable for all model types
- ✅ 3-23% MSE improvement (Monte Carlo R=50)

**When to use:** For most applications requiring the best speed/accuracy tradeoff.

### 2. Unified Iterative (maximum accuracy)

```r
arima_pmm2(y, order = c(1,0,1), pmm2_variant = "unified_iterative")
```

**Characteristics:**
- ✅ Full Newton-Raphson iterative procedure
- ✅ Highest accuracy (up to 23% MSE improvement)
- ✅ Automatic numerical Jacobian (via `numDeriv`)
- ⚠️ Slower than one-step

**When to use:** When accuracy is critical and computation time is not a constraint.

### 3. Linearized (specialist for MA/SMA)

```r
ma_pmm2(y, order = 1, pmm2_variant = "linearized")
```

**Characteristics:**
- ✅ Specialized linear approach for MA/SMA models
- ✅ Optimal efficiency for pure MA models (21-44% MSE improvement)
- ✅ Fast convergence (3-5 iterations)
- ⚠️ Only works for MA/SMA models

**When to use:** For pure MA(q) or SMA(Q) models with asymmetric innovations.

### Variant Comparison (Monte Carlo R=50, n=200)

| Model | Unified Global | Unified Iterative | Linearized | Best |
|-------|----------------|-------------------|------------|------|
| **AR(1)** | -2.2% MSE | **-2.9% MSE** | N/A | Iterative |
| **MA(1)** | **-23.0% MSE** | -19.9% MSE | **-21.6% MSE** | Global |
| **SARIMA** | -15.6% MSE | **-16.4% MSE** | N/A | Iterative |

**Conclusion:** 
- Unified Global — best universal choice (default)
- Unified Iterative — for maximum accuracy
- Linearized — optimal for MA/SMA models


## Installation

```r
# Install the released version from CRAN
install.packages("EstemPMM")

# Or install the development snapshot from GitHub
devtools::install_github("SZabolotnii/EstemPMM")
```

## Basic Usage

```r
library(EstemPMM)

# Create data with asymmetric errors
n <- 100
x <- rnorm(n)
errors <- rgamma(n, shape = 2, scale = 1) - 2  # Shift for zero mean
y <- 2 + 1.5 * x + errors
data <- data.frame(x = x, y = y)

# Fit the model using PMM2
fit <- lm_pmm2(y ~ x, data = data)

# Review results
summary(fit)

# Compare with OLS
ols_fit <- lm(y ~ x, data = data)
compare_with_ols(y ~ x, data)
```

## Time Series Example

```r
library(EstemPMM)

# Simulate AR(2) process with asymmetric innovations
set.seed(42)
n <- 300
innov <- rgamma(n, shape = 2, scale = 1) - 2  # Centered gamma errors
y <- numeric(n)
for (t in 3:n) {
  y[t] <- 0.5 * y[t-1] - 0.3 * y[t-2] + innov[t]
}

# Fit AR(2) model using PMM2
fit_ar <- ar_pmm2(y, order = 2)
summary(fit_ar)

# Compare with other methods
compare_ar_methods(y, order = 2)

# Seasonal models example
# Fit SAR(1,1)_12 model (seasonal lag 12, e.g., monthly data)
fit_sar <- sar_pmm2(y, order = c(1, 1), season = list(period = 12))

# Fit SMA(1)_12 model
fit_sma <- sma_pmm2(y, order = 1, season = list(period = 12))

# Unified SARMA(1,1)×(1,1)_12 model (combines AR, SAR, MA, SMA)
fit_sarma <- sarma_pmm2(y, order = c(1, 1, 1, 1), season = list(period = 12))

# Unified SARIMA(1,1,1)×(1,1,1)_12 model with differencing
fit_sarima <- sarima_pmm2(y, order = c(1, 1, 1, 1),
                          seasonal = list(order = c(1, 1), period = 12))
```

## Main Functions

### Linear Models
- `lm_pmm2()`: Fit linear regression using PMM for S=2
- `compare_with_ols()`: Compare PMM2 with OLS estimates

### Time Series Models
- `ar_pmm2()`: Autoregressive AR(p) models
- `ma_pmm2()`: Moving Average MA(q) models
- `arma_pmm2()`: ARMA(p,q) models
- `arima_pmm2()`: ARIMA(p,d,q) models with differencing
- `sar_pmm2()`: Seasonal Autoregressive SAR(p,P)_s models
- `sma_pmm2()`: Seasonal Moving Average SMA(Q)_s models
- `sarma_pmm2()`: Unified Seasonal ARMA SARMA(p,q)×(P,Q)_s models
- `sarima_pmm2()`: Unified Seasonal ARIMA SARIMA(p,d,q)×(P,D,Q)_s models
- `ts_pmm2()`: Universal wrapper for all time series models

### Comparison Functions
- `compare_ar_methods()`: Compare AR estimation methods
- `compare_ma_methods()`: Compare MA estimation methods
- `compare_arma_methods()`: Compare ARMA estimation methods
- `compare_arima_methods()`: Compare ARIMA estimation methods
- `compare_sar_methods()`: Compare SAR estimation methods
- `compare_ts_methods()`: Universal time series comparison

### Statistical Inference
- `pmm2_inference()`: Bootstrap inference for linear models
- `ts_pmm2_inference()`: Bootstrap inference for time series models

## Documentation & Vignettes

The package ships with three vignettes that walk through typical workflows:

| Vignette | Focus |
| --- | --- |
| `vignette("pmm2-introduction", package = "EstemPMM")` | Linear PMM2 estimation and comparison with OLS |
| `vignette("pmm2-time-series", package = "EstemPMM")` | AR/MA/ARMA/ARIMA plus seasonal SAR/SMA examples |
| `vignette("bootstrap-inference", package = "EstemPMM")` | Resampling-based inference for asymmetric errors |

To rebuild documentation locally:

```r
devtools::document()
devtools::build_vignettes()
```

## Reproducing Monte Carlo Benchmarks

- **Canonical SMA benchmark (n = 120, γ-innovations):** `Rscript run_sma_monte_carlo.R`  
  This script reproduces the 34% variance reduction highlighted in `test_results/SMA_Monte_Carlo_Report_20251113_500reps.md`.
- **Full seasonal comparison (SAR / SMA / SARMA / SARIMA, n ∈ {100, 200, 500}):** `Rscript monte_carlo_seasonal_comparison.R`  
  Results are saved to `monte_carlo_seasonal_results.rds` and summarized in `test_results/SAR_MONTE_CARLO_REPORT_2025-11-14.md`.

The seasonal script takes ~8 minutes on Apple Silicon (M4) because it rebuilds vignettes and runs 500 replications per scenario. The canonical SMA script finishes in < 10 seconds.

## CRAN Readiness & Testing

Before submitting to CRAN, run:

```r
devtools::check()             # quick local sanity check
devtools::document()
devtools::build_vignettes()
system("R CMD build .")
system("R CMD check --as-cran EstemPMM_0.1.3.tar.gz")
```

Continuous integration via `.github/workflows/R-CMD-check.yaml` exercises Ubuntu, macOS, and Windows matrix builds on every push to `main` or `claude/*`.
- `plot_pmm2_bootstrap()`: Visualize bootstrap results

### Utilities
- `compute_moments()`: Calculate sample moments
- `pmm_skewness()`, `pmm_kurtosis()`: Distribution shape measures
- `pmm2_variance_factor()`: Theoretical variance reduction factor
- `pmm2_variance_matrices()`: Compare variance matrices
- `pmm2_monte_carlo_compare()`: Monte Carlo simulations

### S4 Methods
- `summary()`: Display fitting results
- `predict()`: Make predictions based on a PMM model
- `plot()`: Diagnostic plots for PMM models
- `coef()`, `fitted()`, `residuals()`: Extract model components
- `AIC()`: Akaike Information Criterion

## Project Structure

The package consists of several key R files:

### Core Implementation (R/)
- `pmm2_classes.R`: S4 class definitions (PMM2fit, ARPMM2, MAPMM2, ARMAPMM2, ARIMAPMM2, SARPMM2, SMAPMM2, SARMAPMM2, SARIMAPMM2)
- `pmm2_main.R`: Linear regression with PMM2 (`lm_pmm2`)
- `pmm2_ts_main.R`: Time series models (AR, MA, ARMA, ARIMA, SAR, SMA, universal wrapper)
- `pmm2_ts_design.R`: Design matrices and lag structures for time series
- `pmm2_ts_methods.R`: Comparison functions for all time series models
- `pmm2_common.R`: Shared numerical routines and optimization
- `pmm2_utils.R`: Utility functions for moments, variance analysis
- `pmm2_inference.R`: Bootstrap inference for linear and time series models
- `pmm2_monte_carlo.R`: Monte Carlo simulation framework
- `data.R`: Dataset documentation (DCOILWTICO - WTI crude oil prices)

### Documentation
- `vignettes/`: Three comprehensive tutorials (introduction, time series, bootstrap)
- `man/`: 30+ help files with examples
- `docs/`: Extended documentation including SAR/SMA analysis reports
- `NEWS.md`: Version history and changelog

### Testing & Validation
- `tests/testthat/`: 36 unit tests across 7 test files
- `test_results/`: Monte Carlo validation reports with empirical results
- `demo/`: 8 demonstration scripts for all model types

### Development Scripts
- `run_sma_monte_carlo.R`: SMA model validation (500 replications)
- `test_debug_sar.R`: SAR model debugging and verification
- `test_sma_quick.R`: Quick SMA tests

## Variance Diagnostics

You can inspect the theoretical skewness, kurtosis, and expected variance reduction
obtained from the PMM2 algorithm using the helper utilities:

```r
library(EstemPMM)

set.seed(42)
x <- rnorm(200)
errors <- rgamma(200, shape = 2, scale = 1) - 2
y <- 1 + 0.7 * x + errors
dat <- data.frame(y, x)

fit <- lm_pmm2(y ~ x, data = dat)

# Retrieve theoretical cumulants and variance ratio g
vf <- pmm2_variance_factor(fit@m2, fit@m3, fit@m4)
vf$g  # Expected Var(PMM2) / Var(OLS)

# Compare variance matrices for OLS and PMM2
vm <- pmm2_variance_matrices(attr(fit, "model_matrix"),
                             fit@m2, fit@m3, fit@m4)
vm$pmm2
```

## Demo Script

The package includes a detailed demonstration script `pmm2_demo_runner.R` that shows:

1. Comparison of PMM2 and OLS on data with different error distributions
2. Monte Carlo simulations for efficiency evaluation
3. Application to real data (Auto MPG dataset)
4. Bootstrap analysis for uncertainty estimation

To run the demonstration, execute:

```r
source("pmm2_demo_runner.R")
all_results <- run_all_simulations()  # For Monte Carlo simulations
results <- apply_to_mpg_data()  # For real data analysis
```

## Results and Efficiency

PMM2 is particularly effective for distributions with high asymmetry:

### Linear Models
| Distribution    | Skewness | Kurtosis | Theoretical Improvement | Actual Improvement |
|-----------------|----------|----------|------------------------|-------------------|
| Gamma (a=0.5)   | 2.83     | 12       | 57%                    | ~50%              |
| Exponential     | 2.00     | 6        | 50%                    | ~45%              |
| Gamma (a=2)     | 1.41     | 3        | 40%                    | ~35%              |
| Lognormal       | 1.00     | 1.5      | 29%                    | ~25%              |

### Time Series Models

**Autoregressive (AR) Models:** Similar improvements observed, especially for AR(1) and AR(2) with asymmetric innovations.

**Moving Average (MA) Models:** PMM2 shows substantial improvements when innovation distribution has high skewness.

**Seasonal Models (Validated November 2025):**
- **SAR Models:** 20-30% variance reduction with gamma innovations
- **SMA Models:** Up to 34.1% variance reduction (empirically validated with 500 replications), exceeding theoretical predictions for certain parameter combinations
- Both models show stable convergence and computational efficiency

## Adaptive Estimation

The package implements an adaptive procedure for PMM estimation:
1. Find initial OLS estimates and calculate residuals
2. Estimate moments and cumulants of the OLS residuals
3. Calculate refined PMM estimates using these moment estimates

This approach doesn't require prior knowledge of the error distribution properties.

## Applications

The method is particularly useful in:

### Regression Analysis
- Economic and financial modeling with asymmetric error distributions
- Biological systems analysis with skewed measurement errors
- Technical measurements with non-Gaussian noise
- Any regression problem where error distributions exhibit significant skewness

### Time Series Analysis
- Financial time series with asymmetric return distributions
- Economic indicators with seasonal patterns (GDP, inflation, unemployment)
- Energy consumption forecasting with seasonal components
- Climate data analysis with skewed temperature or precipitation distributions
- Commodity prices (e.g., oil, gas) with asymmetric shocks
- Any time series where innovations deviate substantially from normality

## Authors

- Serhii Zabolotnii - Cherkasy State Business College


## Scientific Publications

The Polynomial Maximization Method and its applications are described in the following peer-reviewed publications:

### Foundational Reference
Kunchenko, Y. P., & Lega, Y. G. (1992). *Estimation of Random Variable Parameters by the Polynomial Maximization Method*. Kyiv: Naukova Dumka. 180 pp.

### Linear Regression (PMM2 for lm_pmm2)
Zabolotnii S., Warsza Z.L., Tkachenko O. (2018) Polynomial Estimation of Linear Regression Parameters for the Asymmetric PDF of Errors. In: Szewczyk R., Zieliński C., Kaliczyńska M. (eds) Automation 2018. AUTOMATION 2018. Advances in Intelligent Systems and Computing, vol 743. Springer, Cham. https://doi.org/10.1007/978-3-319-77179-3_75

### Autoregressive Models (PMM2 for ar_pmm2)
Zabolotnii S., Tkachenko O., Warsza Z.L. (2022) Application of the Polynomial Maximization Method for Estimation Parameters of Autoregressive Models with Asymmetric Innovations. In: Szewczyk R., Zieliński C., Kaliczyńska M. (eds) Automation 2022. AUTOMATION 2022. Advances in Intelligent Systems and Computing, vol 1427. Springer, Cham. https://doi.org/10.1007/978-3-031-03502-9_37

### Moving Average Models (PMM2 for ma_pmm2)
Zabolotnii S., Tkachenko O., Warsza Z.L. (2023) Polynomial Maximization Method for Estimation Parameters of Asymmetric Non-gaussian Moving Average Models. In: Szewczyk R., et al. (eds) Automation 2023. AUTOMATION 2023. Lecture Notes in Networks and Systems, vol 630. Springer, Cham. https://doi.org/10.1007/978-3-031-25844-2_21

## License

This project is distributed under the GPL-3 License.
