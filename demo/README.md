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

#### `pmm_ts_examples`
Comprehensive time series modeling with PMM2.

- **Models covered:**
  - AR(2) with heavy-tailed errors (Student-t)
  - MA(2) with normal innovations
  - ARMA(1,1) with mixture distributions
  - ARIMA(1,1,1) with asymmetric errors

- **Features:**
  - Comparison of CSS vs PMM2 estimates
  - Forecasting with PMM2 models
  - Efficiency visualization (skewness vs kurtosis)

**Run with:**
```r
demo("pmm_ts_examples")
```

---

#### `pmm2_real_data`
Application to real-world Auto MPG dataset.

- **Analysis includes:**
  - Linear regression: MPG ~ acceleration
  - Comparison of OLS vs PMM2
  - Residual diagnostics
  - Bootstrap confidence intervals
  - Visualization with ggplot2

**Run with:**
```r
demo("pmm2_real_data")
```

---

#### `pmm2_prediction`
Prediction accuracy comparison with train/test splits and cross-validation.

- **Experiments:**
  1. Simple 80/20 train-test split
  2. k-fold cross-validation
  3. Linear and quadratic models

- **Metrics:**
  - MSE, MAE, R²
  - Out-of-sample prediction errors
  - Visual comparisons

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

#### `pmm2_demo_runner`
Interactive menu-driven demonstration system.

- **Options:**
  1. Monte Carlo simulations
  2. Real data analysis
  3. Prediction accuracy evaluation
  4. Quick test

**Run with:**
```r
demo("pmm2_demo_runner")
```

**Note:** This demo sources other files and may not work correctly after package installation. Use individual demos instead.

---

## Recommended Learning Path

For new users, we recommend the following sequence:

1. **Start here:** `test_pmm` (30 seconds)
   Get familiar with basic syntax

2. **See the benefits:** `pmm2_comparison_boxplots` (2-3 minutes) ⭐
   Understand PMM2's efficiency gains visually

3. **Linear models:** `pmm2_real_data` (2-3 minutes)
   See PMM2 applied to real data

4. **Time series:** `pmm_ts_examples` (1-2 minutes)
   Learn AR, MA, ARMA, ARIMA modeling

5. **Deep dive:** `pmm2_simulation` (5-10 minutes)
   Explore Monte Carlo evidence

---

## Dependencies

### Core Requirement
- `EstemPMM` package (obviously!)

### Optional Packages (for enhanced visualizations)

Some demos use additional packages for better graphics:

- `ggplot2` - Advanced plotting
- `gridExtra` - Multi-panel layouts
- `dplyr` - Data manipulation
- `reshape2` - Data reshaping
- `parallel` - Parallel computing for faster simulations

**Installation:**
```r
install.packages(c("ggplot2", "gridExtra", "dplyr", "reshape2"))
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
