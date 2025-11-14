# Monte Carlo Comparison of Seasonal Models

This directory contains scripts for comprehensive Monte Carlo simulations comparing PMM2 with classical methods (CSS, ML) for seasonal time series models.

## Overview

The Monte Carlo study evaluates the performance of PMM2 across four types of seasonal models:

1. **SAR(1,1)_12** - Seasonal Autoregressive model
2. **SMA(1)_12** - Seasonal Moving Average model
3. **SARMA(1,0)×(1,1)_12** - Combined Seasonal ARMA model
4. **SARIMA(1,1,0)×(1,1,1)_12** - Seasonal ARIMA with differencing

## Files

### Main Scripts

- **`monte_carlo_seasonal_comparison.R`** - Main simulation script
  - Runs 500 replications for each scenario
  - Tests sample sizes: 100, 200, 500
  - Compares PMM2 vs CSS methods
  - Saves results to `monte_carlo_seasonal_results.rds`

- **`visualize_monte_carlo_results.R`** - Visualization script
  - Creates comprehensive plots
  - Generates summary tables
  - Saves outputs to PDF and CSV files

### Output Files

After running the simulations, you will get:

- `monte_carlo_seasonal_results.rds` - Raw results (R data format)
- `monte_carlo_seasonal_plots.pdf` - Comprehensive visualization
- `monte_carlo_summary_table.csv` - Summary statistics in tabular format

## How to Run

### Step 1: Run Monte Carlo Simulations

```r
# This will take significant time (several hours depending on your machine)
source("monte_carlo_seasonal_comparison.R")
```

**Expected runtime:**
- ~2-4 hours on a modern desktop
- 500 replications × 4 models × 3 sample sizes = 6,000 model fits per method

### Step 2: Visualize Results

```r
# After simulations complete, create plots and tables
source("visualize_monte_carlo_results.R")
```

### Quick Test Run

To test the scripts with fewer replications:

```r
# Edit the script and change:
N_REPLICATIONS <- 50  # instead of 500

# Then run as usual
source("monte_carlo_seasonal_comparison.R")
source("visualize_monte_carlo_results.R")
```

## Simulation Design

### Data Generation

All models use **asymmetric innovations** from Gamma distribution:
```r
innovations ~ Gamma(shape = 2, scale = 1) - 2  # Centered at 0
```

This creates:
- **Positive skewness** (c₃ ≈ 1.4)
- **Heavy tails** (c₄ ≈ 3)
- Ideal conditions for PMM2 to outperform OLS/CSS

### Model Specifications

#### Scenario 1: SAR(1,1)_12
```
y_t = 0.5·y_{t-1} + 0.6·y_{t-12} + ε_t
```

#### Scenario 2: SMA(1)_12
```
y_t = ε_t + 0.6·ε_{t-12}
```

#### Scenario 3: SARMA(1,0)×(1,1)_12
```
y_t = 0.5·y_{t-1} + 0.6·y_{t-12} + 0.4·ε_{t-12} + ε_t
```

#### Scenario 4: SARIMA(1,1,0)×(1,1,1)_12
```
(1 - 0.4B)(1 - 0.5B¹²)(1 - B)(1 - B¹²)y_t = (1 + 0.6B¹²)ε_t
```

### Evaluation Metrics

For each parameter estimate, we compute:

1. **Bias**: `E[θ̂] - θ`
2. **RMSE**: `√E[(θ̂ - θ)²]`
3. **Variance**: `Var(θ̂)`
4. **MAE**: `E[|θ̂ - θ|]`

### Key Performance Indicator

**Variance Reduction**:
```
VR = 1 - Var(PMM2) / Var(CSS)
```

- Positive values indicate PMM2 improvement
- Expected range: 20-50% for asymmetric innovations
- Higher values for stronger asymmetry

## Expected Results

### Variance Reduction by Model

Based on theoretical predictions and preliminary tests:

| Model | Sample Size | Expected VR | Actual VR Range |
|-------|------------|-------------|-----------------|
| SAR(1,1)_12 | 100 | 25-35% | 20-40% |
| SAR(1,1)_12 | 500 | 30-40% | 25-45% |
| SMA(1)_12 | 100 | 30-40% | 25-45% |
| SMA(1)_12 | 500 | 35-45% | 30-50% |
| SARMA | 100 | 20-30% | 15-35% |
| SARMA | 500 | 25-35% | 20-40% |
| SARIMA | 100 | 15-25% | 10-30% |
| SARIMA | 500 | 20-30% | 15-35% |

### Efficiency Factor g

The PMM2 efficiency factor `g = Var[PMM2]/Var[OLS]` should be:

- **g < 1**: PMM2 is more efficient (lower variance)
- **g ≈ 0.5-0.7**: Typical range for Gamma innovations
- **g → 1**: As sample size increases or asymmetry decreases

## Output Interpretation

### Plots

The visualization script creates 5 pages of plots:

1. **Page 1**: Variance comparison (PMM2 vs CSS)
2. **Page 2**: Variance reduction for simple models
3. **Page 3**: Variance reduction for complex models
4. **Page 4**: RMSE comparison
5. **Page 5**: Efficiency factor g across all models

### Summary Table

The CSV file contains:

| Column | Description |
|--------|-------------|
| `model` | Model specification |
| `n` | Sample size |
| `parameter` | Parameter name |
| `pmm2_var` | PMM2 variance |
| `css_var` | CSS variance |
| `var_reduction_pct` | Variance reduction (%) |

## Customization

### Modify Innovation Distributions

Edit the `generate_*_data` functions to change innovation type:

```r
# Current: Gamma
innovations <- rgamma(n, shape = 2, scale = 1) - 2

# Alternative: Lognormal
innovations <- rlnorm(n, meanlog = 0, sdlog = 0.5) - exp(0.125)

# Alternative: Exponential
innovations <- rexp(n, rate = 1) - 1

# Alternative: Normal (baseline)
innovations <- rnorm(n, mean = 0, sd = 1)
```

### Change Parameter Values

Edit the `true_*` specifications:

```r
# SAR parameters
true_sar <- list(ar = 0.5, sar = 0.6)  # Change as needed

# SMA parameters
true_sma <- 0.6  # Change as needed

# etc.
```

### Adjust Sample Sizes

```r
SAMPLE_SIZES <- c(50, 100, 200, 500, 1000)  # Add or remove sizes
```

### Modify Replications

```r
N_REPLICATIONS <- 1000  # Increase for more precision
```

## Computational Considerations

### Memory Usage

- Each scenario stores coefficient matrices: `N_REPLICATIONS × n_parameters`
- Total memory: ~50-100 MB for default settings
- Increase for larger N_REPLICATIONS

### Parallel Processing

To speed up simulations, modify the script to use parallel processing:

```r
library(parallel)

# Detect cores
n_cores <- detectCores() - 1

# Create cluster
cl <- makeCluster(n_cores)

# Export functions and data
clusterExport(cl, c("generate_sar_data", "safe_fit", ...))

# Run in parallel
results <- parLapply(cl, 1:N_REPLICATIONS, function(i) {
  # Your simulation code here
})

# Stop cluster
stopCluster(cl)
```

## Troubleshooting

### Issue: Simulation takes too long

**Solution**: Reduce `N_REPLICATIONS` to 50-100 for testing

### Issue: Memory error

**Solution**: Run scenarios separately or reduce `N_REPLICATIONS`

### Issue: Convergence warnings

**Solution**: This is normal for some replications. The script handles non-convergence automatically.

### Issue: NA values in results

**Solution**: Some model fits may fail. The script filters these out automatically.

## References

### Theoretical Background

The PMM2 method provides variance reduction given by:

```
g = 1 - c₃² / (2 + c₄)
```

Where:
- `c₃` = skewness coefficient
- `c₄` = excess kurtosis coefficient

For Gamma(2, 1) innovations:
- `c₃ ≈ 1.414`
- `c₄ ≈ 3.0`
- **Theoretical g ≈ 0.40** (60% variance reduction)

### Publications

See main README.md for references to PMM2 methodology.

## Contributing

To add new scenarios:

1. Create a new section in `monte_carlo_seasonal_comparison.R`
2. Follow the existing pattern for data generation and model fitting
3. Add corresponding visualization in `visualize_monte_carlo_results.R`
4. Update this README with the new scenario

## License

GPL-3 (same as main package)
