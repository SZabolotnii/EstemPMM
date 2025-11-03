# PMM2-ARIMA Reproducibility Package

This directory bundles everything reviewers need to regenerate the empirical results that feed into `latex/PMM2_ARIMA.tex`. The scripts reproduce both the WTI case-study analysis and the Monte Carlo experiments from Section 3.

## Repository Layout

- `data/` – immutable inputs delivered with the package (`DCOILWTICO.csv` from FRED).
- `scripts/` – R entry points that orchestrate analysis, plotting, and reporting.
- `results/` – auto-generated artefacts (CSV tables, RDS stores, plots, Markdown report, Monte Carlo summaries).
- `latex/` – manuscript sources (`PMM2_ARIMA.tex`, `references.bib`).
- `DESCRIPTION` – declared R dependencies for quick inspection/tooling.

## Requirements

- R ≥ 4.3 (tested with 4.5.1).
- R packages: `EstemPMM` (version 0.1.1), `ggplot2`, `gridExtra`, `RColorBrewer`, `tseries`, `knitr`, `MASS`.
  - Install CRAN packages: `install.packages(c("ggplot2", "gridExtra", "RColorBrewer", "tseries", "knitr", "MASS"))`
  - Install `EstemPMM` version 0.1.1 from the included archive:
    `install.packages("EstemPMM2-lib/EstemPMM_0.1.1.tar.gz", repos = NULL, type = "source")`
  - Alternatively, install from GitHub: `remotes::install_github("SZabolotnii/EstemPMM")`
- Optional: `pandoc` if you plan to convert the Markdown report to HTML/PDF, and a LaTeX distribution (`TeXLive`, `TinyTeX`, …) to compile `latex/PMM2_ARIMA.tex`.

## Reproducing the WTI Experiment

From the root of this directory, run:

```bash
Rscript scripts/run_full_study.R
```

The command will:

1. Fit all ARIMA specifications with CSS-ML vs PMM2 (`scripts/comprehensive_study.R`).
2. Regenerate the 10 publication plots (`scripts/create_visualizations.R`).
3. Rebuild the narrative Markdown report (`scripts/generate_report.R`).

You can also execute any step individually:

```bash
# Quick console-only sanity check
Rscript scripts/arima_oil_quick_demo.R

# Full modelling pipeline
Rscript scripts/comprehensive_study.R

# Visuals and report (after the previous step)
Rscript scripts/create_visualizations.R
Rscript scripts/generate_report.R
```

## Output Map

After a successful run you should see:

- `results/full_results.csv` – values used in Table \ref{tab:wti_comprehensive_results} of the manuscript.
- `results/method_comparison.csv` – deltas reported in Table \ref{tab:wti_method_comparison}.
- `results/descriptive_stats.csv` – descriptive panel in Section 4.
- `results/plots/*.png` – Figures 1–10 (numbered to match captions).
- `results/ANALYTICAL_REPORT.md` – extended prose that aligns with the LaTeX discussion.
- `results/monte_carlo/*.csv` – simulation summaries mirroring the Monte Carlo tables (per-model metrics, relative efficiencies, residual cumulants).
- `results/monte_carlo/article_comparison.csv` – expanded table that juxtaposes simulation estimates with published article values and adds deviations/ratios.

The numerical values in the LaTeX tables were exported by doubling the AIC/BIC scores for formatting purposes; the CSV files contain the raw outputs from R. Reviewers can regenerate the manuscript by re-running the scripts and then compiling `latex/PMM2_ARIMA.tex`.

## Monte Carlo Simulations

The dedicated driver `scripts/run_monte_carlo.R` reconstructs the simulation study (default: 2000 replications per combination, seed = 12345). Expect the full run to take 20–40 minutes depending on hardware.

```bash
# Full Monte Carlo run
Rscript scripts/run_monte_carlo.R

# Customization (fewer replications for quick verification)
Rscript scripts/run_monte_carlo.R --reps=200 --seed=20250101

# Variation of assumptions
Rscript scripts/run_monte_carlo.R --standardize-innov=false --css-method=CSS-ML
```

Options `--standardize-innov` and `--css-method` allow checking sensitivity to innovation scale and type of classical estimate (available values: `CSS`, `CSS-ML`, `ML`). The `--m-est=false` flag disables additional robust estimates, which are computed by default for AR components.

Key artifacts:

- `monte_carlo_metrics.csv` – complete cross-section of metrics (bias, var, MSE, RE, VR) for each model, parameter, and distribution.
- `arima110_summary.csv`, `arima110_re_vs_sample_size.csv` – data for ARIMA(1,1,0) tables and RE dependence on N.
- `arima011_summary.csv`, `arima111_summary.csv`, `arima210_summary.csv` – condensed comparisons for other configurations.
- `arima110_residual_cumulants.csv` – average cumulants of PMM2 residuals (Table \ref{tab:residual_cumulants}).
- `article_comparison.csv` – comparison of simulations with published values and new columns with deviations.
- The `M-EST` column in summary files contains results of Huber M-estimates for AR components, allowing comparison of PMM2/CSS with robust estimates.

## Bootstrap Confidence Intervals (NEW)

After running Monte Carlo simulations, you can add bootstrap confidence intervals to the results:

```bash
# Add bootstrap CIs to Monte Carlo results
Rscript scripts/add_confidence_intervals.R

# Customization
Rscript scripts/add_confidence_intervals.R --bootstrap-reps=2000 --confidence=0.95
```

This creates `monte_carlo_metrics_with_ci.csv` containing standard errors and confidence intervals for bias, variance, MSE, and relative efficiency. Addresses reviewer request for uncertainty quantification in Monte Carlo tables.

## Out-of-Sample Validation (NEW)

Evaluate forecasting performance on held-out data:

```bash
# Fixed train/test split (80/20) + rolling window forecasts
Rscript scripts/wti_out_of_sample.R

# Customization
Rscript scripts/wti_out_of_sample.R --train-fraction=0.8 --window-size=100
```

Outputs:
- `wti_fixed_split_validation.csv` – RMSE/MAE for train/test split
- `wti_rolling_window_validation.csv` – rolling window forecast metrics

Addresses reviewer request for out-of-sample validation.

## Enhanced Diagnostics (NEW)

Generate comprehensive diagnostic plots and statistics with p-values:

```bash
# Full diagnostic suite for WTI case study
Rscript scripts/wti_diagnostics.R

# For different model specification
Rscript scripts/wti_diagnostics.R --order=1,1,1
```

Outputs:
- `wti_diagnostics_statistics.csv` – Ljung-Box, Jarque-Bera, Shapiro-Wilk, ARCH tests with p-values
- `wti_qq_plots.png` – Q-Q plots for residual normality assessment
- `wti_acf_pacf.png` – autocorrelation diagnostics
- `wti_residual_time_series.png` – residual plots over time
- `wti_residual_histograms.png` – residual distributions with normal overlay

Addresses reviewer request for diagnostic plots with p-values.

## Software Versions

- **R:** Version 4.5.1 (2025-06-13)
- **EstemPMM:** Version 0.1.1 (included in `EstemPMM2-lib/EstemPMM_0.1.1.tar.gz`)
- **Platform:** macOS Sequoia 15.6.1 (aarch64-apple-darwin24.4.0)
- **Random seed:** All scripts use `set.seed(12345)` for reproducibility
- **Session info:** Complete R environment details available in `sessionInfo.txt`

## Housekeeping

- `results/` ships with a `.gitignore`; delete or replace the placeholder when publishing generated artefacts.
- Use `DESCRIPTION` with `renv` or `pak` if you wish to create a fully locked environment before sharing.
