# EstemPMM Demonstrations

This directory contains interactive demonstrations of the **Polynomial Maximization Method (PMM2)** for statistical estimation.

---

## Quick Start

To run a demo, use:

```r
library(EstemPMM)
demo("demo_name", package = "EstemPMM")
```

**Example:**
```r
demo("pmm2_comparison_boxplots", package = "EstemPMM")
```

---

## Available Demonstrations

### 🌟 Recommended First Demo

**`pmm2_comparison_boxplots`** (2-3 minutes)
Monte Carlo comparison of PMM2 vs OLS with comprehensive boxplot visualizations.

- **What it shows:**
  - Parameter estimate distributions across 500 Monte Carlo runs
  - PMM2 efficiency gains with skewed and heavy-tailed errors
  - Variance reduction: 10-30% improvement over OLS
  - Visual comparisons through boxplots and density plots

- **Run with:**
  ```r
  demo("pmm2_comparison_boxplots")
  ```

---

### Quick Demos (< 1 minute)

#### `test_pmm`
Quick sanity check demonstrating basic PMM2 functionality.

- Linear regression example
- AR(1) time series example
- Coefficient comparisons: OLS vs PMM2

**Run with:**
```r
demo("test_pmm")
```

---

### Intermediate Demos (1-5 minutes)

#### `pmm_ts_examples` (SIMPLIFIED & INTERACTIVE)
Comprehensive time series modeling with PMM2.

- **Models covered:**
  - AR(1) with skewed innovations (gamma distribution)
  - MA(1) with heavy-tailed innovations (Student-t)
  - ARMA(1,1) with contaminated innovations (outliers)
  - ARIMA(1,1,1) with integrated series (chi-squared)

- **Features:**
  - Side-by-side comparison: CSS vs PMM2
  - Interactive progression through examples (press Enter)
  - Comprehensive visualizations (ACF, PACF, Q-Q plots, boxplots)
  - True vs estimated coefficients
  - Residual diagnostics

- **NEW: Simplified version**
  - No external dependencies (base R graphics only)
  - Clear educational progression
  - Self-contained examples
  - Duration: 2-3 minutes (interactive)

**Run with:**
```r
demo("pmm_ts_examples")
```

---

#### `pmm2_real_data` (SIMPLIFIED)
Application to real-world Auto MPG dataset.

- **Analysis includes:**
  - Linear regression: MPG ~ acceleration
  - Comparison of OLS vs PMM2
  - Residual diagnostics and moment analysis
  - Visual comparison with base R graphics
  - Interpretation guide for real-world applications

- **NEW: Simplified version**
  - No external dependencies (no ggplot2)
  - Faster execution (1-2 minutes)
  - Focus on core analysis
  - Refers to vignettes for advanced topics (bootstrap)

**Run with:**
```r
demo("pmm2_real_data")
```

---

#### `pmm2_prediction` (SIMPLIFIED)
Prediction accuracy comparison using train/test split.

- **Demonstrates:**
  - 80/20 train-test data split
  - Model training on training set
  - Out-of-sample prediction on test data
  - Performance metrics: MSE, MAE, R²

- **Visualizations:**
  - Training vs test data
  - Prediction errors boxplot
  - Predicted vs actual scatter plot
  - Model comparison summary table

- **NEW: Simplified version**
  - Removed k-fold cross-validation (too complex)
  - Removed multiple experiments
  - Faster execution (1-2 minutes)
  - Focus on core prediction concepts
  - No external dependencies

**Run with:**
```r
demo("pmm2_prediction")
```

---

### Advanced Demos (5-10 minutes)

#### `pmm2_simulation`
Monte Carlo simulation studies for linear regression.

- **Error distributions:**
  - Gaussian (baseline)
  - Student-t (heavy tails)
  - Gamma (right-skewed)
  - Exponential (asymmetric)
  - Chi-squared (highly skewed)
  - Log-normal (extreme skewness)

- **Outputs:**
  - Bias and variance comparisons
  - MSE ratios
  - ggplot2 visualizations
  - Theoretical vs empirical efficiency gains

**Run with:**
```r
demo("pmm2_simulation")
```

---

#### `pmm2_simMC_ts`
Monte Carlo simulations for time series models.

- **Models:**
  - AR(p), MA(q), ARMA(p,q), ARIMA(p,d,q)

- **Methods compared:**
  - MLE, CSS, OLS, Yule-Walker, PMM2

- **Features:**
  - Customizable innovations (gamma, t-distribution)
  - Coefficient estimation accuracy
  - Bias-variance decomposition

**Run with:**
```r
demo("pmm2_simMC_ts")
```

---


---

## Recommended Learning Path

For new users, we recommend the following sequence:

1. **Start here:** `test_pmm` (30 seconds)
   Get familiar with basic syntax - NO DEPENDENCIES ✨

2. **See the benefits:** `pmm2_comparison_boxplots` (2-3 minutes) ⭐
   Understand PMM2's efficiency gains visually - NO DEPENDENCIES ✨

3. **Real data:** `pmm2_real_data` (1-2 minutes)
   See PMM2 applied to Auto MPG dataset - SIMPLIFIED, NO DEPENDENCIES ✨

4. **Prediction:** `pmm2_prediction` (1-2 minutes)
   Learn train/test validation - SIMPLIFIED, NO DEPENDENCIES ✨

5. **Time series:** `pmm_ts_examples` (2-3 minutes) ⭐
   Learn AR, MA, ARMA, ARIMA modeling - INTERACTIVE, NO DEPENDENCIES ✨

6. **Deep dive:** `pmm2_simulation` (5-10 minutes)
   Explore Monte Carlo evidence (requires ggplot2)

---

## Dependencies

### Core Requirement
- `EstemPMM` package (obviously!)

### External Dependencies

**Most demos now use BASE R GRAPHICS ONLY!** ✨

The following simplified demos have **NO external dependencies**:
- `test_pmm` - Base R only
- `pmm2_comparison_boxplots` - Base R only ⭐
- `pmm_ts_examples` - Base R only (SIMPLIFIED, INTERACTIVE) ✨
- `pmm2_real_data` - Base R only (SIMPLIFIED)
- `pmm2_prediction` - Base R only (SIMPLIFIED)

**Advanced demos** may require optional packages:
- `pmm2_simulation` - Uses `ggplot2`, `gridExtra`, `dplyr`, `parallel`
- `pmm2_simMC_ts` - Uses `dplyr`, `ggplot2`

**Installation (if needed for advanced demos):**
```r
install.packages(c("ggplot2", "gridExtra", "dplyr"))
```

**Note:** Demos will check for these packages and provide clear error messages if missing.

---

## Tips for Running Demos

### Controlling Simulation Size

For faster execution of Monte Carlo demos, you can modify parameters:

```r
# In your R session before running demo:
n_sim <- 100  # Instead of default 500 or 1000
```

### Saving Demo Output

To save plots or results:

```r
# Redirect graphics to PDF
pdf("pmm2_comparison.pdf")
demo("pmm2_comparison_boxplots")
dev.off()
```

### Parallel Computing

Some demos support parallel processing:

```r
# Set number of cores (default: detectCores() - 1)
options(mc.cores = 4)
demo("pmm2_simulation")
```

---

## Troubleshooting

### "Package not found" Error

Make sure EstemPMM is loaded:
```r
library(EstemPMM)
demo("demo_name")
```

### Missing Dependencies

Install required packages:
```r
install.packages(c("ggplot2", "gridExtra", "dplyr"))
```

### Slow Performance

- Reduce simulation size (edit demo file)
- Use parallel computing (if available)
- Run simpler demos first

---

## Contributing

Have ideas for new demonstrations? Found a bug? Please open an issue on GitHub:

https://github.com/SZabolotnii/EstemPMM/issues

---

## Additional Resources

- **Vignettes:** Run `browseVignettes("EstemPMM")` for detailed tutorials
- **Documentation:** `?lm_pmm2`, `?ar_pmm2`, `?pmm2_inference`
- **GitHub:** https://github.com/SZabolotnii/EstemPMM
- **Paper:** Zabolotnii et al. (2018), Springer AUTOMATION 2018

---

**Happy exploring with PMM2!** 🎯
