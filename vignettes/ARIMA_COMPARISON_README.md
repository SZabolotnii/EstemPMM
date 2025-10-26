# ARIMA Model Comparison: PMM2 vs Classical Methods

This directory contains scripts and reports for comparing classical ARIMA estimation methods (CSS and ML) with the Polynomial Maximization Method order 2 (PMM2) on real-world crude oil price data.

## Files

- **`arima_oil_comparison.Rmd`**: R Markdown document with comprehensive analysis including:
  - Data loading and preprocessing
  - Exploratory data analysis (EDA)
  - Model fitting using CSS-ML and PMM2 methods
  - Comparison using AIC/BIC information criteria
  - Residual diagnostics and validation
  - Coefficient comparison across methods
  - Conclusions and recommendations

- **`run_arima_comparison.R`**: Quick script to render the full HTML/PDF report

- **`../demo/arima_oil_quick_demo.R`**: Simplified demo for quick testing without generating plots

## Data

The analysis uses WTI (West Texas Intermediate) crude oil price data from the Federal Reserve Economic Data (FRED) database:
- **File**: `data/DCOILWTICO.csv`
- **Period**: October 2020 - October 2025
- **Frequency**: Daily
- **Source**: https://fred.stlouisfed.org/series/DCOILWTICO

## Running the Analysis

### Option 1: Quick Demo (Console Output Only)

```r
# Load the package
library(EstemPMM)

# Run quick demo
demo("arima_oil_quick_demo", package = "EstemPMM")

# Or run directly
source("demo/arima_oil_quick_demo.R")
```

This will:
- Fit ARIMA(1,1,1), ARIMA(0,1,1), ARIMA(1,1,0), and ARIMA(2,1,1) models
- Compare CSS-ML and PMM2 methods
- Display AIC/BIC values in console
- Identify best models by information criteria

**Expected output**:
```
Model         Method  AIC     BIC     LogLik   RSS
ARIMA(1,1,1)  CSS-ML  xxxx.xx xxxx.xx -xxx.xx  xxxxx.xx
ARIMA(1,1,1)  PMM2    xxxx.xx xxxx.xx -xxx.xx  xxxxx.xx
...

Best by AIC: ARIMA(x,x,x) (METHOD)
Best by BIC: ARIMA(x,x,x) (METHOD)
```

### Option 2: Full Report with Visualizations

```r
# Navigate to vignettes directory
setwd("vignettes")

# Render full HTML report
source("run_arima_comparison.R")

# Or use rmarkdown directly
library(rmarkdown)
render("arima_oil_comparison.Rmd")
```

This generates:
- `arima_oil_comparison.html`: Interactive HTML report with plots
- `arima_oil_comparison.pdf`: PDF version (if LaTeX available)

**Report includes**:
- Time series plots of oil prices
- ACF/PACF analysis
- Model comparison tables
- Information criteria plots
- Residual diagnostic plots
- Ljung-Box tests
- Coefficient comparison across methods

### Option 3: Custom Analysis

```r
library(EstemPMM)

# Load data
data_path <- "data/DCOILWTICO.csv"
oil_data <- read.csv(data_path)
prices <- na.omit(as.numeric(oil_data$DCOILWTICO))

# Fit ARIMA(1,1,1) using both methods
fit_css <- arima(prices, order = c(1, 1, 1), method = "CSS-ML")
fit_pmm2 <- arima_pmm2(prices, order = c(1, 1, 1))

# Compare information criteria
cat("CSS-ML: AIC =", AIC(fit_css), "BIC =", BIC(fit_css), "\n")
cat("PMM2:   AIC =", AIC(fit_pmm2), "BIC =", BIC(fit_pmm2), "\n")

# Compare coefficients
print(coef(fit_css))
print(coef(fit_pmm2))

# Plot diagnostics
par(mfrow = c(2, 2))
plot(fit_pmm2)
```

## Models Tested

The analysis compares the following ARIMA specifications:

1. **ARIMA(1,1,1)**: Standard integrated model
   - 1 AR term, 1 difference, 1 MA term
   - Good balance between fit and parsimony

2. **ARIMA(0,1,1)**: Integrated Moving Average (IMA)
   - Random walk with drift and MA component
   - Common for financial time series

3. **ARIMA(1,1,0)**: Autoregressive Integrated (ARI)
   - Random walk with AR component
   - Captures momentum effects

4. **ARIMA(2,1,1)**: Extended model
   - 2 AR terms, 1 difference, 1 MA term
   - More flexible but risks overfitting

## Evaluation Criteria

### Information Criteria

- **AIC (Akaike Information Criterion)**:
  - AIC = -2*log(L) + 2*k
  - Balances fit quality with model complexity
  - Lower values indicate better models

- **BIC (Bayesian Information Criterion)**:
  - BIC = -2*log(L) + log(n)*k
  - Stronger penalty for complex models
  - Preferred for model selection with large samples

### Residual Diagnostics

- **RSS (Residual Sum of Squares)**: Measures fit quality
- **MAE (Mean Absolute Error)**: Average prediction error
- **Skewness**: Measures asymmetry in residuals
- **Kurtosis**: Measures tail heaviness
- **Ljung-Box Test**: Tests for residual autocorrelation

## Expected Results

Based on typical financial time series characteristics:

1. **Non-stationarity**: Oil prices require differencing (d=1)
2. **Moderate autocorrelation**: AR and MA terms capture dynamics
3. **Non-Gaussian residuals**: Fat tails and skewness common
4. **PMM2 robustness**: Should perform well when residuals deviate from normality

## Interpretation Guidelines

### AIC/BIC Values
- Lower is better
- Difference of 2-3 points is marginally significant
- Difference of 10+ points is strongly significant

### Method Selection
- Use **BIC** when sample size is large (>500) - stronger penalty
- Use **AIC** for smaller samples or when prediction is primary goal
- Consider **PMM2** when residual kurtosis > 5 (fat tails)

### Model Adequacy
- Ljung-Box test p-value > 0.05: residuals are white noise (good)
- Residual ACF within confidence bands: no remaining structure
- Q-Q plot close to diagonal: normality assumption reasonable

## Dependencies

Required R packages:
```r
install.packages(c("rmarkdown", "knitr", "ggplot2", "moments", "tseries"))
```

## References

1. Kunchenko, Y., & Lega, Y. (1992). *Polynomial parameter estimations of close to Gaussian random variables*. Kiev: Kyiv Polytechnic Institute.

2. Akaike, H. (1974). A new look at the statistical model identification. *IEEE Transactions on Automatic Control*, 19(6), 716-723.

3. Schwarz, G. (1978). Estimating the dimension of a model. *The Annals of Statistics*, 6(2), 461-464.

## Contact

For questions or issues, please open an issue on the GitHub repository.

## License

This analysis is part of the EstemPMM package and follows the same license.
