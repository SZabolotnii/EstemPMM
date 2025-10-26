# Comprehensive ARIMA Study: PMM2 vs Classical Methods

This directory contains a complete research pipeline for comparing PMM2 and classical ARIMA estimation methods on WTI crude oil price data.

## Quick Start

```r
# Run the complete study (3-5 minutes)
source("analysis/run_full_study.R")

# View the analytical report
cat analysis/results/ANALYTICAL_REPORT.md

# Or convert to HTML for better viewing
pandoc analysis/results/ANALYTICAL_REPORT.md \
  -o analysis/results/ANALYTICAL_REPORT.html \
  -s --toc --metadata title="ARIMA Comparison Study"
```

## Pipeline Overview

The study consists of three main steps:

### Step 1: Model Fitting and Data Analysis
**Script:** `comprehensive_study.R`

Performs:
- Data loading and preprocessing (DCOILWTICO.csv)
- Descriptive statistics and stationarity tests
- Fitting 6 ARIMA specifications using both CSS-ML and PMM2 methods
- Collection of comprehensive performance metrics:
  - Information criteria (AIC, BIC, Log-likelihood)
  - Error metrics (RSS, RMSE, MAE, MAPE)
  - Residual diagnostics (skewness, kurtosis)
  - Computation time
  - Convergence status

**Models tested:**
- ARIMA(0,1,1) - Integrated Moving Average
- ARIMA(1,1,0) - Autoregressive Integrated
- ARIMA(1,1,1) - Standard model
- ARIMA(2,1,1) - Extended AR
- ARIMA(1,1,2) - Extended MA
- ARIMA(2,1,2) - Most flexible

**Output:**
- `results/full_results.csv` - Complete results table
- `results/method_comparison.csv` - Method-wise comparison
- `results/descriptive_stats.csv` - Data characteristics
- `results/*.rds` - R objects for further analysis

### Step 2: Visualization Creation
**Script:** `create_visualizations.R`

Generates 10 publication-quality plots:

1. **AIC Comparison** - Bar chart comparing AIC across models
2. **BIC Comparison** - Bar chart comparing BIC across models
3. **RMSE Comparison** - Error comparison by method
4. **Computation Time** - Performance efficiency analysis
5. **Kurtosis Comparison** - Residual tail behavior
6. **Skewness Comparison** - Residual asymmetry
7. **Performance Heatmap** - Normalized metrics across all models
8. **Method Differences** - Direct PMM2 vs CSS-ML comparison
9. **Best Model Diagnostics** - Residual analysis of top model
10. **Summary Statistics** - Aggregate performance by method

**Output:**
- `results/plots/*.png` - All visualizations (300 DPI, publication-ready)

### Step 3: Report Generation
**Script:** `generate_report.R`

Creates a comprehensive analytical report in Markdown format with:

- **Executive Summary** - Key findings and recommendations
- **Data Characteristics** - Descriptive statistics and stationarity analysis
- **Model Specifications** - Detailed description of tested models
- **Comprehensive Results** - Full results table and rankings
- **Method Comparison** - Head-to-head analysis with win rates
- **Residual Analysis** - Distribution characteristics
- **Computational Efficiency** - Speed comparison
- **Conclusions** - Main findings and practical recommendations
- **Limitations** - Study constraints and future work
- **References** - Academic citations

**Output:**
- `results/ANALYTICAL_REPORT.md` - Complete research report

## Directory Structure

```
analysis/
├── README.md                      # This file
├── run_full_study.R               # Master script (runs entire pipeline)
├── comprehensive_study.R          # Step 1: Model fitting
├── create_visualizations.R        # Step 2: Visualization
├── generate_report.R              # Step 3: Report generation
└── results/                       # Output directory (created automatically)
    ├── full_results.csv
    ├── method_comparison.csv
    ├── descriptive_stats.csv
    ├── all_results.rds
    ├── fitted_models.rds
    ├── study_summary.rds
    ├── ANALYTICAL_REPORT.md
    └── plots/
        ├── 01_aic_comparison.png
        ├── 02_bic_comparison.png
        ├── 03_rmse_comparison.png
        ├── 04_computation_time.png
        ├── 05_kurtosis_comparison.png
        ├── 06_skewness_comparison.png
        ├── 07_performance_heatmap.png
        ├── 08_method_differences.png
        ├── 09_best_model_diagnostics.png
        └── 10_summary_statistics.png
```

## Requirements

### R Packages

```r
# Core analysis
install.packages("EstemPMM")       # PMM2 implementation

# Visualization
install.packages("ggplot2")        # Plotting
install.packages("gridExtra")      # Multi-panel plots
install.packages("RColorBrewer")   # Color palettes
install.packages("tidyr")          # Data reshaping
install.packages("dplyr")          # Data manipulation

# Statistical tests
install.packages("moments")        # Skewness/kurtosis
install.packages("tseries")        # ADF test
```

### System Tools (Optional)

For HTML/PDF report generation:
```bash
# Install pandoc (for document conversion)
sudo apt install pandoc           # Ubuntu/Debian
brew install pandoc               # macOS

# For PDF generation (optional)
sudo apt install texlive-latex-base texlive-latex-extra
```

## Usage Examples

### Run Complete Study

```r
# Execute entire pipeline
source("analysis/run_full_study.R")
```

**Expected output:**
```
STEP 1: Running comprehensive model comparison...
  Fitting 12 models (6 specifications × 2 methods)
  ✓ Model fitting completed successfully

STEP 2: Creating visualizations...
  Generating 10 publication-quality plots
  ✓ Visualization creation completed successfully

STEP 3: Generating analytical report...
  Creating comprehensive markdown report
  ✓ Report generation completed successfully

STUDY COMPLETED SUCCESSFULLY
Total Time: 180 seconds (3.0 minutes)
```

### Run Individual Steps

```r
# Step 1 only: Fit models and collect data
source("analysis/comprehensive_study.R")

# Step 2 only: Create visualizations (requires Step 1 output)
source("analysis/create_visualizations.R")

# Step 3 only: Generate report (requires Step 1 output)
source("analysis/generate_report.R")
```

### Access Results in R

```r
# Load complete results
results <- readRDS("analysis/results/all_results.rds")
View(results)

# Load fitted models for further analysis
models <- readRDS("analysis/results/fitted_models.rds")

# Access best model
summary_info <- readRDS("analysis/results/study_summary.rds")
best_model_name <- summary_info$best_model_bic
```

### View and Export Report

```r
# View in terminal
system("cat analysis/results/ANALYTICAL_REPORT.md")

# Convert to HTML
system("pandoc analysis/results/ANALYTICAL_REPORT.md \
  -o analysis/results/ANALYTICAL_REPORT.html \
  -s --toc --toc-depth=3 \
  --metadata title='ARIMA Comparison Study'")

# Convert to PDF (requires LaTeX)
system("pandoc analysis/results/ANALYTICAL_REPORT.md \
  -o analysis/results/ANALYTICAL_REPORT.pdf \
  -s --toc --toc-depth=3")

# Open HTML in browser
browseURL("analysis/results/ANALYTICAL_REPORT.html")
```

## Expected Results

Based on WTI crude oil price characteristics, the study typically finds:

### Data Characteristics
- **Non-stationarity**: Original series requires differencing (d=1)
- **Volatility**: Moderate to high (CV ≈ 15-25%)
- **Distribution**: Non-normal with heavy tails (kurtosis > 3)
- **Asymmetry**: Often right-skewed during price shocks

### Model Performance
- **Best models**: Usually ARIMA(1,1,1) or ARIMA(0,1,1) by BIC
- **Complexity**: Simple models often preferred (parsimony principle)
- **Fit quality**: AIC values typically 7000-8000, BIC 7000-8100

### Method Comparison
- **Competitiveness**: PMM2 and CSS-ML often within 1-2% of each other
- **Win rates**: Typically 40-60% split depending on metric
- **Robustness**: PMM2 handles heavy tails slightly better
- **Speed**: PMM2 often 1-2x faster for simple models

## Interpretation Guidelines

### Information Criteria
- **AIC difference > 2**: Marginally significant
- **AIC difference > 10**: Strongly significant
- **BIC preferred**: For large samples (n > 500)

### Model Selection
1. Start with BIC rankings (parsimony)
2. Check top 3 models for stability
3. Validate residual diagnostics (Ljung-Box p > 0.05)
4. Assess Q-Q plot normality
5. Consider practical constraints (interpretability, computation)

### Method Choice
- **Use PMM2 when**:
  - Heavy-tailed residuals (|kurtosis| > 3)
  - Robustness to outliers required
  - Alternative validation needed

- **Use CSS-ML when**:
  - Approximately normal residuals
  - Standard reporting required
  - Integration with existing workflows

## Troubleshooting

### Common Issues

**1. Data file not found**
```r
# Check data path
file.exists("data/DCOILWTICO.csv")

# If in wrong directory, specify full path
data_path <- "/full/path/to/data/DCOILWTICO.csv"
```

**2. Package installation errors**
```r
# Update R to latest version
# Install packages one by one to identify problem
install.packages("EstemPMM")
install.packages("ggplot2")
# ... etc
```

**3. Memory issues (large datasets)**
```r
# Increase memory limit
memory.limit(size = 8000)  # Windows
# Or run on subset of data
```

**4. Convergence failures**
```r
# Check results for failed models
results[!results$Convergence, ]

# Try different starting values or increase iterations
arima_pmm2(x, order = c(1,1,1), max_iter = 100)
```

## Customization

### Test Different Models

Edit `comprehensive_study.R`:
```r
model_specs <- list(
  "ARIMA(3,1,1)" = c(3, 1, 1),  # Add more AR terms
  "ARIMA(1,1,3)" = c(1, 1, 3),  # Add more MA terms
  "ARIMA(2,0,2)" = c(2, 0, 2)   # No differencing
)
```

### Modify Visualizations

Edit `create_visualizations.R`:
```r
# Change color scheme
method_colors <- c("CSS-ML" = "blue", "PMM2" = "red")

# Adjust plot dimensions
ggsave(..., width = 12, height = 8, dpi = 600)  # Higher resolution

# Add custom plots
p_custom <- ggplot(...) + ...
```

### Customize Report

Edit `generate_report.R`:
```r
# Add custom sections
cat("## 9. CUSTOM ANALYSIS\n\n")
cat("Your additional analysis here...\n\n")

# Modify formatting
cat("**Bold text** or *italic text*\n")
```

## Citation

If you use this study in your research, please cite:

```bibtex
@software{arima_comparison_pmm2,
  title = {Comprehensive ARIMA Study: PMM2 vs Classical Methods},
  author = {EstemPMM Development Team},
  year = {2025},
  url = {https://github.com/SZabolotnii/EstemPMM}
}
```

## Support

For questions or issues:
1. Check this README
2. Review the analytical report
3. Examine script comments
4. Open an issue on GitHub

## License

This analysis framework is part of the EstemPMM package and follows the same license.
