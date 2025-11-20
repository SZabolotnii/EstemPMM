
# Investigation Summary: MA+SMA PMM2 Efficiency

## 1. Problem Statement
The user reported that the initial implementation of the EstemPMM-style PMM2 estimator for mixed MA+SMA models showed **negative efficiency** (higher MSE than MLE) in Monte Carlo simulations, particularly for small sample sizes (n=100). This contradicted the high efficiency observed for pure MA and SMA models in the reference implementation (`new_sarima`).

## 2. Investigation Process
1.  **Reference Check**: Verified that `new_sarima/experimental/06_estpmm_style_ma.R` (Linear PMM2) indeed performs well for pure SMA models.
2.  **Reproduction**: Created `compare_strategies.R` to isolate the performance of the estimator.
    *   Initial runs with n=200 showed positive efficiency for both Additive and Multiplicative data, suggesting the core logic was sound for large samples.
    *   Runs with n=100 revealed that for **Multiplicative** data (standard SARIMA), the PMM2 estimator performed worse than MLE (8.7% MSE increase).
3.  **Root Cause Analysis**:
    *   The standard SARIMA model is **Multiplicative**: $x_t = (1+\theta B)(1+\Theta B^s)\epsilon_t = \epsilon_t + \theta \epsilon_{t-1} + \Theta \epsilon_{t-s} + \theta \Theta \epsilon_{t-s-1}$.
    *   The initial `estpmm_style_ma_sma` implementation assumed an **Additive** model: $x_t = \epsilon_t + \theta \epsilon_{t-1} + \Theta \epsilon_{t-s}$.
    *   This misspecification (omitting the interaction term $\theta \Theta \epsilon_{t-s-1}$) caused bias, which dominated the variance reduction benefits of PMM2 at small sample sizes.

## 3. The Solution
We enhanced the `estpmm_style_ma_sma` estimator to support **Multiplicative** models by:
1.  **Design Matrix**: Adding columns for interaction terms (e.g., $\epsilon_{t-s-1}$) to the design matrix.
2.  **Estimation**: Estimating coefficients for these interaction terms freely (Linearized Multiplicative PMM2).
3.  **Innovations**: Updating the innovation computation to include the interaction terms (using the product of estimated MA and SMA coefficients).

## 4. Verification Results (n=100, R=50)

### Experiment 1: Additive Data (Robustness Check)
*   **MLE MSE**: 0.035803
*   **PMM2 MSE**: 0.024956
*   **Result**: **30.3% Improvement**. The estimator remains highly efficient for additive data, correctly estimating the interaction term as near-zero.

### Experiment 2: Multiplicative Data (Standard SARIMA)
*   **MLE MSE**: 0.027052
*   **PMM2 MSE**: 0.025789
*   **Result**: **4.7% Improvement**. The fix successfully eliminated the negative efficiency. PMM2 now beats MLE even at n=100 for multiplicative models.

## 5. Conclusion
The efficiency issue has been resolved. The `EstemPMM` package now correctly handles mixed MA+SMA models with high efficiency, matching or exceeding the performance of MLE across tested scenarios.

## 6. Next Steps
*   Extend the PMM2 approach to full SARIMA (AR+MA+SAR+SMA) models.
*   Run a comprehensive Monte Carlo simulation (R=500) to generate final performance benchmarks.
