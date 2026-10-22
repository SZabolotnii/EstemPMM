# EstemPMM News

## Development version (0.5.0.9000)

Not on CRAN. CRAN 0.5.0 stays the reference for published work and for
the JSS manuscript until the next release.

### Standard errors after an OLS fall-back; `summary()` of PMM2 fits

`vcov()` for `PMM2fit` and `PMM3fit` (available since 0.3.2) returns
`g * m2 * (X'X)^{-1}`. When the iteration fails, the fit returns the OLS
start, and scaling by `g` would understate that estimate's variance.
`lm_pmm2()` and `lm_pmm3()` now record the fall-back in the attribute
`"fallback"`, and `vcov()` then returns the OLS covariance
`m2 * (X'X)^{-1}`. `summary()` of a `PMM2fit` prints standard errors,
z values and p-values from `vcov()`. Bootstrap inference still needs
`formula =` and `data =`.

`dev_scripts/check_pmm2_se.R` ran 1000 replications under chi-squared(3),
Gamma(2), exponential, lognormal and Gaussian errors. At n >= 200 the
95% interval covers 0.93-0.95. At n = 30-50 it covers 0.88-0.94, because
the plug-in standard error is 8-28% too small. At n = 30 with five
regressors, `lm_pmm2()` falls back to OLS in 14-28% of the replications
on skewed laws. Its MSE is still 1.2-1.6 times lower than OLS there.

### New: general degree-3 PMM, `lm_pmm3(symmetric = FALSE)`

`lm_pmm3()` fits the symmetric estimator, which assumes
`gamma3 = gamma5 = 0`. With `symmetric = FALSE` it fits the general degree-3
estimator instead. The polynomial is `psi(e) = sum_k h_k (e^k - alpha_k)`,
with `h` solving the full Kunchenko normal system built from the first six
central moments of the OLS residuals. Its factor is
`g3 = 1 / (alpha2 * b' F^{-1} b)`. That factor never exceeds `g2` and equals
the symmetric `g3` when the odd moments vanish. The Newton iteration is
safeguarded like PMM2: backtracking on the score norm, and a fall-back to
OLS with a warning. `PMM3fit` gains the slots `symmetric` (default `TRUE`)
and `poly_coef`. For the general fit, `vcov()` is `g3 * m2 * (X'X)^{-1}`.

- On 14 data sets it matches the solver of the PMM-SLS-BRIDGE revision
  (AISM-D-26-00101, `common.R::pmm_reg` with one outer step) within 1.5e-14
  (`tests/testthat/fixtures/reference_aism_pmm3.rds`,
  `dev_scripts/make_reference_aism.R`).
- In `dev_scripts/check_pmm3_general.R` (y = 2 + 1.5x + e, x ~ U(0, 5),
  1000 replications), the OLS-to-estimator slope MSE ratio at n = 200 is:

  | Errors | PMM2 | PMM3 symmetric | PMM3 general |
  |---|---|---|---|
  | chi-squared(3) | 2.05 | 1.30 | 3.16 |
  | Gamma(2) | 1.69 | 1.17 | 2.22 |
  | Exponential | 2.21 | 1.37 | 3.78 |
  | Lognormal(0, 0.5) | 1.73 | 1.25 | 2.24 |
  | Uniform | 0.99 | 2.82 | 2.83 |
  | Gaussian | 0.99 | 0.98 | 0.96 |

  The 95% interval from `vcov()` covers 0.93-0.96 at n >= 200 and
  0.88-0.94 at n = 50. At n = 50 the plug-in standard error underestimates
  the spread by 10-20%.
- The symmetric-PMM3 warning for `|gamma3| > 0.3` now names
  `symmetric = FALSE`.

### New: weak-moment PMM, `lm_wpmm2()` and `lm_wpmm3()`

Formula-interface versions of the windowed (weak-moment) PMM2 and PMM3
estimators from gmdhpmm (EstemPMM-GMDH, `R/weak_pmm.R`). The raw residual
moments are replaced by moments under a Gaussian window centred on the
residual median, with width `sigma_mult` times the MAD (default 2.5). The
fit is computed by IRLS with clipped, redescending weights. The new class
`WPMMfit` inherits from `PMMfit` and has `coef`, `residuals`, `fitted`,
`predict`, `nobs`, `vcov`, `confint`, `summary` and `show` methods.
`pmm_lm()` accepts `method = "wpmm2"` and `"wpmm3"`, but `method = "auto"`
never selects them.

- On ten KG-2 designs the estimates match gmdhpmm 0.2.0 within 2e-9; the gap
  comes from the 1e-8 ridge penalty gmdhpmm adds
  (`tests/testthat/fixtures/reference_gmdhpmm.rds`,
  `dev_scripts/make_reference_gmdhpmm.R`).
- `vcov()` is the M-estimator sandwich. Its bread keeps the window slope
  `w'(r)`, which does not vanish under symmetric errors. In
  `dev_scripts/check_wpmm_se.R` (y = 1 + 2x + e, 1000 replications,
  n = 100 and 400, errors t3, Laplace, Gaussian, centred exponential, 10%
  contamination and Cauchy) the mean sandwich standard error of the slope
  is within 5% of the empirical SD, and the 95% intervals cover 0.926-0.955.
- Slope MSE relative to OLS at n = 400:

  | Errors | WPMM2 | WPMM3 |
  |---|---|---|
  | t3 | 0.55 | 0.54 |
  | Laplace | 0.75 | 0.73 |
  | Centred exponential | 0.31 | 0.52 |
  | 10% contamination | 0.11 | 0.11 |
  | Gaussian | 1.03 | 1.04 |

  The Gaussian row is the cost of the weak estimators; at n = 100 it rises
  to 1.08 for WPMM2 and 1.12 for WPMM3.
- The documented limitations: under skewed errors the intercept estimates a
  location that depends on the window width, not the error mean, while the
  slopes are unaffected. `kappa^w` in WPMM3 is the coefficient of the
  windowed law, not the variance-minimising one.

### `pmm_dispatch()`: calibrated rule

The 0.5.0 thresholds (0.3 / 1.0 / -0.7 / 0.95) had never been calibrated.
A Monte Carlo study covered 20 error laws, n from 30 to 1000 and 5000
replications per cell (`dev_scripts/calibrate_dispatch.R`,
`dev_scripts/analyze_dispatch.R`; summary in `mc_results_dispatch/README.md`).
It measured the slope MSE of the dispatcher's choice against the best of
OLS, PMM2 and PMM3. The study found two systematic losses:

- **Small-sample platykurtic residuals.** Bimodal and uniform residuals at
  n = 30-50 often show a spurious skewness above 0.3. The arm "asymmetric AND
  platykurtic -> PMM2" then fired, at up to 2.7 times the best MSE.
- **Symmetric heavy tails.** With logistic, Laplace and Student-t errors and
  n >= 100, PMM3 has 0.62 to 0.94 of the OLS variance, but PMM3 was offered
  only for `gamma4 < -0.7`.

Retuning the existing thresholds does not fix both losses without losing
PMM2 on skewed laws at n <= 50: at n = 30 a skewed law cannot be told
reliably from a symmetric one. Two targeted changes do fix them, through new
arguments:

- `clear_asymmetry_z = 2`: skewed platykurtic residuals go to PMM2 only when
  `|gamma3| > 2 * sqrt(6/n)`, otherwise to PMM3;
- `g3_threshold = 0.95`, `n_leptokurtic = 100`: symmetric residuals with
  `g3 < 0.95` go to PMM3 once `n >= 100`, whatever the sign of `gamma4`.

On held-out replications the mean log efficiency loss falls from 0.061 to
0.028 and the worst cell from 2.67 to 1.23. 28 of 120 cells improve by more
than 3% and one gets worse (exponential errors, n = 30: 1.01 -> 1.04).
`clear_asymmetry_z = 0, g3_threshold = 0` reproduces the 0.5.0 decisions
exactly; this was checked on 20 000 random samples.

**Changed example.** `mpg ~ horsepower + I(horsepower^2)` on `auto_mpg`
(n = 392, gamma4 = 1.3) is now sent to PMM3 rather than OLS. A pairs
bootstrap shows PMM3 and OLS variances within 1-3% of each other there,
because the errors are heteroscedastic and the predicted gain assumes
i.i.d. errors. The vignette now says so.

**Not covered.** The calibration used linear regression with one Gaussian
regressor; the thresholds have not been checked on time-series models. The
symmetry threshold stays at a fixed `|gamma3| = 0.3`. Sample skewness of
heavy-tailed errors is noisy even at large n (0.335 on one Laplace sample of
2000), and such residuals then go to OLS rather than PMM3.

**OPEN_ISSUES D-1.** PMM2 also beats OLS on symmetric heavy-tailed errors:
MSE ratio 0.72-0.85 on t3 and 0.94-0.99 on Laplace. For Laplace and t5 the
gain shrinks with n, so it is a finite-sample effect and not a breach of
`g2 = 1`. t3 has no fourth moment, so the `g2` theorem does not apply to it.

### Safeguarded PMM2 Newton iteration

`pmm2_algorithm()`, the solver behind `lm_pmm2()` and the fixed-design
time-series fitters, took the full Newton step unconditionally and returned
the last iterate even when it had not converged. With strongly skewed
errors and few observations per parameter it could run away silently: on
300 replications per cell (`n = 30, p = 8` and `n = 40, p = 15`, exponential,
lognormal and Gamma(2) errors) the 0.5.0 coefficient MSE was 2.7 to 12 000
times that of OLS.

The iteration now

- keeps the full Newton step whenever it reduces the norm of the estimating
  function, and otherwise halves it (backtracking line search);
- returns the initial estimate (OLS, or CSS for time series) with
  `convergence = FALSE` and a warning when the line search stalls or the
  iterate is not finite;
- issues a `warning()` instead of a `verbose`-only message when the
  iteration limit is reached;
- warns when a converged root has a Jacobian that is not positive definite.
  The estimating equation can have several roots, and only the one reached
  from a consistent start is the PMM2 estimate; in small samples that root
  can fail the test too, so the estimate is reported, not replaced.

On the same cells the MSE ratio against OLS is now 0.82 to 0.97. With
`n = 100, p = 5` the results are unchanged except for one runaway
replication removed.

**Reproducibility.** A fixture of 27 fits made with the CRAN 0.5.0 build
(`tests/testthat/fixtures/reference_0.5.0.rds`, generated by
`dev_scripts/make_reference_050.R`) is replayed by the test suite. All 25
fits that 0.5.0 reported as converged, including every time-series and WTI
fit, are reproduced: 23 bit-for-bit, one to 9e-16. The remaining one
(lognormal errors, `n = 30, p = 12`) changes on purpose. In 0.5.0 it
"converged" to a second root where half of the Jacobian weights were
negative (coefficient error 2.6 against 0.95 for OLS). It now returns the
OLS start with `convergence = FALSE`.

### `pmm_dispatch()`: the `g3` guard is correct as it stands

It has been suggested that `.g3_factor()` should reject the band
`0 < 6 + 9*gamma4 + gamma6 <= gamma4^2`. With `gamma4` and `gamma6`
computed from sample moments this band cannot be reached, because
`6 + 9*gamma4 + gamma6 - gamma4^2 = m6/m2^3 - (m4/m2^2)^2 >= 0` by the
Cauchy-Schwarz inequality. The boundary is attained only by a symmetric
two-point law, and there `g3 = 0` is right: the PMM3 estimating equation
holds at every observation and the estimate is exact. The guard is
unchanged; a test now pins `0 <= g3 <= 1`.

## Version 0.5.0 (2026-08-09)

### Major: recursive (non-linearised) estimation of the MA component

Until 0.4.0 every PMM2 time-series path with a moving-average component
built the design matrix **once** from the CSS residuals
(`arma_build_design()`, `ma_build_design()`) and held it fixed while the
PMM2 Newton iteration ran. For the AR columns this is exact -- the
regressors are observed data and do not depend on the parameter. For the
MA columns it is a **one-step linearisation of the innovation
recursion**, and it costs a factor of exactly `1 - theta^2` in
asymptotic variance relative to CSS. The penalty does *not* shrink as
the sample grows: on ARIMA(0,1,1) it stayed near -30% of variance for
`T` = 100, 200, 500 and 1000 alike.

#### New: `ma_solver = c("linearized", "recursive")`

`ts_pmm2()`, `ma_pmm2()`, `arma_pmm2()` and `arima_pmm2()` gained an
`ma_solver` argument. With `"recursive"` the innovations
`eps_t(beta)` and the score regressors
`x_t(beta) = -d eps_t / d beta` are recomputed by the exact recursion at
**every** candidate parameter, so there is no linearisation. When
`mu3 = 0` the PMM2 estimating equation then reduces *identically* to the
CSS first-order condition, and the `1 - theta^2` penalty disappears
exactly rather than asymptotically.

Measured effect (ARIMA(0,1,1), `theta = -0.5`, `T = 500`, 2000
replications, relative efficiency against CSS; higher is better):

| innovations  | linearised | recursive | theory `(2+g4)/(2+g4-g3^2)` |
|--------------|-----------:|----------:|----------------------------:|
| Gaussian     |     0.7723 |    0.9863 |                      1.0000 |
| Gamma(2,1)   |     1.2807 |    1.6986 |                      1.6667 |
| lognormal    |     1.3289 |    1.8116 |                      1.6822 |
| chi-squared  |     1.3538 |    1.7557 |                      1.8000 |

Variance reduction on the MA branch rises from about 23% to about 42%,
i.e. to parity with the pure AR branch. The Gaussian penalty falls from
-23% to -1.4%.

**The default remains `"linearized"`** so that results from 0.4.x
reproduce bit-for-bit. `"recursive"` is the recommended setting whenever
`q >= 1`; the default is scheduled to change in the next major release.

#### Scope of the recursive branch

Supported: non-seasonal `ma`, `arma` and `arima` model types with
`q >= 1` and any `p >= 0`, with or without an intercept, and any
differencing order `d`. Verified against numerical differentiation for
ARMA(2,2) with an intercept.

Not supported (silently keep the linearised path): the seasonal
`sar_pmm2()` / `sma_pmm2()` / `sarma_pmm2()` / `sarima_pmm2()` family.
Note that `sarima_pmm2()` has an unrelated argument that is *also*
called `ma_method`, taking `"mle"` / `"pmm2"`; the two are not
interchangeable.

Pure AR models are unaffected: with `q = 0` the recursion is the
identity map and reproduces the fixed lag design bit-for-bit, so
`ma_solver` is ignored and the estimates are unchanged.

#### New: robust solver

The recursive estimating equation is solved by a safeguarded Newton
iteration with

- the analytic Jacobian obtained from a **second** recursion
  (`d x_t / d beta`, closed form), with a Gauss-Newton fallback;
- backtracking line search on the norm of the estimating function;
- projection onto the stationary/invertible region -- any step that
  leaves it is halved;
- a CSS fallback, so the solver never returns `NA`.

On the benchmark grid above (ARIMA(0,1,1), `theta = -0.5`, `T = 500`,
2000 replications x 4 laws) the solver returned a finite, invertible
estimate and reported convergence on **100.0%** of replications for all
four laws. A bracketing `uniroot()` on `[-0.97, 0.97]` -- the obvious
alternative -- fails on up to 6% of replications under skewed
innovations (1995/1882/1952/1873 successes out of 2000 for
Gaussian/Gamma/lognormal/chi-squared).

#### New: closed-form standard errors for MA/ARMA/ARIMA

For recursive fits the score regressors are `F_{t-1}`-measurable, so the
estimating function is a martingale difference sequence and the sandwich
covariance collapses to

```
Sigma = (Delta / a) * M^{-1},
  a = mu4 - mu2^2,  Delta = mu2 * a - mu3^2,  M = E[x_t x_t'],
  Delta / a = mu2 * (1 - gamma3^2 / (gamma4 + 2)) = mu2 * g2
```

-- the same expression already used for AR models, with `X` replaced by
the exact score regressors. `vcov()` and `confint()` therefore now work
for non-seasonal `ma` / `arma` / `arima` fits produced with
`ma_solver = "recursive"`; no numerical Jacobian and no HAC correction
are involved. Agreement with an empirical numerical sandwich is within
5%.

For linearised fits `vcov()` and `confint()` still refuse -- the frozen
CSS design is not the derivative of the criterion, so the formula does
not apply there -- and direct users to `ts_pmm2_inference()`.

#### Known limitation: intercept standard errors

Neyman orthogonality with respect to the plug-in moments holds if and
only if `mu3 * E[x_t] = 0`. Every **slope** regressor has `E[x_t] = 0`,
so estimating `mu2`, `mu3`, `mu4` from residuals does not affect the
limiting distribution of the AR/MA coefficients. The **intercept**
regressor has `E[x_t] = 1`, so under skewed innovations (`mu3 != 0`)
estimating `mu2` contributes at first order and the intercept's standard
error is understated. This affects the linearised branch as well and is
not new in 0.5.0, but it is now documented; `vcov()` covers the slope
parameters only and warns when an intercept was estimated under visible
skewness.

#### Other

- New S4 slot `TS2fit@ma_solver` records which treatment produced the
  fit. Objects created before 0.5.0 and every seasonal fit default to
  `"linearized"` via the class prototype.
- New internal functions: `arma_recursion()`, `pmm2_recursive_score()`,
  `pmm2_recursive_solve()`, `pmm2_recursive_design()`,
  `arma_admissible()`, `ma_pmm2_fit_recursive()`.
- Note for anyone comparing branches: on ARMA/ARIMA fits with `d = 0`
  and `include.mean = TRUE`, the linearised branch reports the fitted
  offset alone as `@intercept` while the recursive branch reports the
  full mean (CSS anchor plus correction). The recursive convention is
  the correct one.
- 365 new assertions in `tests/testthat/test-pmm2-recursive-ma.R`,
  including an exact-identity check against CSS under symmetric
  innovations and a regression test freezing the `1 - theta^2`
  behaviour of the linearised branch.

## Version 0.4.0 (2026-05-28)

### Major: revised class hierarchy and unified user interface

Substantial refactor in response to the JSS reviewer report (May 2026)
on the accompanying manuscript. The package methodology is unchanged;
the user interface and S4 class graph are reorganised to address the
reviewer's concerns about a flat class structure, scattered top-level
functions, and missing methods.

#### New unified API (primary going forward)

- **`pmm_lm(formula, data, method = c("auto", "pmm2", "pmm3"), ...)`** --
  primary entry point for non-Gaussian linear regression. `method = "auto"`
  consults `pmm_dispatch()` on the OLS residuals and forwards to the
  recommended polynomial order.
- **`pmm_ar()`**, **`pmm_ma()`**, **`pmm_arma()`**, **`pmm_arima()`**,
  **`pmm_sarima()`** -- unified time-series entry points with a
  `method = c("pmm2", "pmm3")` selector.
- The older `lm_pmm2()` / `lm_pmm3()` / `ar_pmm2()` / `arima_pmm3()` /
  ... functions remain available without runtime warnings, but are now
  documented as the lower-level fitters that the unified `pmm_*` family
  wraps. Runtime deprecation warnings are scheduled for 0.4.1; full
  removal for 0.5.0.

#### New virtual S4 class hierarchy

- **`PMMfit`** (virtual root) -- holds the slots common to every fit
  object: `coefficients`, `residuals`, `convergence`, `iterations`,
  `call`.
- **`BasePMM2`** (virtual; existing) -- adds PMM2-family moments
  `m2`, `m3`, `m4`. Now inherits from `PMMfit`.
- **`BasePMM3`** (virtual; new) -- adds PMM3-family moments and
  cumulant coefficients `m2`, `m4`, `m6`, `gamma4`, `gamma6`,
  `g_coefficient`, `kappa`.
- **`PMMtsfit`** (virtual; new) -- adds the time-series slots
  `model_type`, `intercept`, `original_series`, `order`.
- `TS2fit` and `TS3fit` now inherit from both their respective
  `BasePMM*` and from `PMMtsfit` via S4 multiple inheritance, with the
  diamond resolving at `PMMfit`.

#### New methods on existing classes

- **`show()`** methods for `PMM2fit`, `PMM3fit`, `TS2fit`, `TS3fit`
  replace the default S4 slot-dump display with an `lm()`-style header
  (call + coefficients + convergence/iteration footer).
- **`vcov()`** and **`confint()`** are now provided for `PMM3fit` and
  for AR-model `TS3fit` objects, mirroring the PMM2 implementations.
  The PMM3 asymptotic covariance formula `V = g_3 sigma^2 (X'X)^{-1}`
  with `g_3 = 1 - gamma_4^2 / (6 + 9*gamma_4 + gamma_6)` is derived in
  `notes/pmm3_vcov_derivation.md` and validated by the Monte Carlo
  script in the accompanying JSS replication bundle (max relative error
  < 3% at n = 500).
- New helper **`pmm3_variance_matrices(X, m2, m4, m6)`** returns both
  the OLS and PMM3 asymptotic covariance matrices.

#### S4/S3 cleanup

- The custom `setMethod("AIC", ...)` and `setMethod("BIC", ...)`
  methods for `PMM2fit`, `PMM3fit`, `TS2fit`, and `TS3fit` are
  **removed**. `AIC()` and `BIC()` now dispatch through the generic
  S3 mechanism using the existing `logLik()` methods, which are
  converted from `setMethod()` to S3 (`logLik.PMM2fit` etc.) so that
  `stats::AIC.default`'s internal `UseMethod("logLik")` finds them.
  Numerical values are identical to the removed custom methods.
- `logLik.TS2fit()` now filters non-finite residuals to match the AIC
  convention previously used by the removed `AIC.TS2fit` method.

#### New `PMMdispatch` S3 class

- `pmm_dispatch()` now returns an S3 object of class `PMMdispatch`
  with dedicated `print`, `summary`, and `format` methods instead of a
  bare list. Field access via `$` is preserved for back-compat.

#### Numerical / code-quality fixes

- `pmm2_inference()` and `ts_pmm2_inference()` now compute two-tailed
  p-values via `2 * pnorm(abs(t), lower.tail = FALSE)` instead of the
  numerically less stable `2 * (1 - pnorm(abs(t)))`. The result is
  identical for moderate t-values; the rewrite avoids loss of precision
  at extreme `|t|`.
- Both inference functions now centre the residuals before
  resampling so that any small non-zero residual mean introduced by
  the PMM iteration does not shift the simulated responses.
- Both inference functions gained a `ci_method` argument with three
  choices: `"normal"` (the new default; symmetric Wald interval
  using the bootstrap standard deviation, always centred on the
  point estimate), `"percentile"` (the previous default; Efron's
  empirical percentile interval), and `"basic"` (Davison & Hinkley
  1997 pivotal interval). The returned data frame now also carries
  a `bias` column equal to `mean(boot) - estimate` so users can see
  when a particular method exhibits finite-sample bias -- this is
  especially relevant for block bootstrap on AR coefficients near
  the stationarity boundary, where the percentile interval can
  drift away from the point estimate.
- Explicit `importFrom(methods, show)` so the new `show()` methods
  load cleanly from a fresh namespace.

---

## Version 0.3.1 (2026-04-06)

### CRAN Resubmission

- Added `CLAUDE.md` and built tarballs to `.Rbuildignore` to eliminate
  non-standard top-level file NOTEs.
- No code changes from 0.3.0; version bump required for CRAN resubmission.

---

## Version 0.3.0 (2026-03-19)

### New Feature: PMM3 for Symmetric Platykurtic Errors

PMM3 (S=3) extends the Polynomial Maximization Method to handle symmetric
error distributions with negative excess kurtosis (platykurtic), such as
uniform, beta-symmetric, and truncated normal errors.

#### New Functions

- **`lm_pmm3()`** - Linear regression estimation using PMM3 (S=3) with
  Newton-Raphson solver. Includes adaptive kappa mode, step-size limiting,
  and divergence guard.
- **`pmm_dispatch()`** - Automatic method selection (OLS / PMM2 / PMM3)
  based on residual cumulant analysis.
- **`compute_moments_pmm3()`** - Compute central moments m2, m4, m6 and
  derived quantities (gamma4, gamma6, g3, kappa) from residuals.
- **`pmm3_variance_factor()`** - Theoretical variance reduction factor
  g3 = 1 - gamma4^2 / (6 + 9*gamma4 + gamma6).
- **`pmm_gamma6()`** - Sixth-order cumulant coefficient.
- **`test_symmetry()`** - Test residual symmetry to guide PMM2 vs PMM3 choice.

#### PMM3 Time Series Functions

- **`ts_pmm3()`** - General PMM3 time series estimation (AR/MA/ARMA/ARIMA).
- **`ar_pmm3()`** - AR model estimation using PMM3.
- **`ma_pmm3()`** - MA model estimation using PMM3.
- **`arma_pmm3()`** - ARMA model estimation using PMM3.
- **`arima_pmm3()`** - ARIMA model estimation using PMM3.

#### New S4 Classes

- **`PMM3fit`** - Standalone class for linear regression (no inheritance
  from BasePMM2) with slots for m2, m4, m6, gamma4, gamma6, g_coefficient,
  and kappa. Full S4 methods: `coef()`, `residuals()`, `fitted()`,
  `predict()`, `summary()`, `plot()`, `AIC()`.
- **`TS3fit`** - Base class for PMM3 time series, with subclasses
  `ARPMM3`, `MAPMM3`, `ARMAPMM3`, `ARIMAPMM3`. Full S4 methods:
  `coef()`, `residuals()`, `fitted()`, `predict()`, `summary()`,
  `plot()`, `AIC()`.

#### Documentation

- New vignette: "PMM3: Linear Regression for Symmetric Platykurtic Errors"
- New vignette: "PMM3 for Time Series: AR, MA, ARMA, and ARIMA Models"
- Updated package-level documentation with PMM3 and method selection sections

---

## Version 0.2.0 (2025-11-20)

### Major Update: Unified PMM2 Architecture

This release represents a significant architectural improvement based on comprehensive research comparing different PMM2 implementation strategies.

#### New Features

- **Unified PMM2 Framework** - Universal PMM2 estimator supporting any nonlinear regression model
  - `pmm2_nonlinear_onestep()` - One-step global correction (default, recommended)
  - `pmm2_nonlinear_iterative()` - Full iterative Newton-Raphson procedure
  - Automatic numerical Jacobian computation via `numDeriv` package
  - Works with AR, MA, ARMA, ARIMA, SAR, SMA, SARIMA models
  
- **Three PMM2 Variants** - New `pmm2_variant` parameter in all time series functions:
  - `"unified_global"` (default) - One-step correction, fast and stable
  - `"unified_iterative"` - Full iterative procedure for maximum accuracy
  - `"linearized"` - Specialized linear approach for MA/SMA models (EstemPMM-style)
  
- **Enhanced Numerical Stability**
  - Optional numerical Jacobian when analytical derivatives unavailable
  - Improved convergence diagnostics
  - Regularization options for ill-conditioned systems

#### Research-Based Improvements

Based on Monte Carlo simulations (R=50, n=200) comparing three approaches:

| Approach | AR(1) | MA(1) | SARIMA | Status |
|----------|-------|-------|---------|--------|
| **Unified Iterative** | -2.9% MSE | -19.9% MSE | **-16.4% MSE** | ✅ **Best overall** |
| **Unified One-step** | -2.2% MSE | **-23.0% MSE** | -15.6% MSE | ✅ **Fastest** |
| **Linearized (MA)** | N/A | -21.6% MSE | N/A | ✅ **MA specialist** |
| Direct Nonlinear | N/A | ❌ **Failed** | ❌ **Failed** | ⛔ **Removed** |

**Key findings:**
- Unified approaches provide consistent 3-23% MSE improvement
- One-step (global) variant offers best speed/accuracy tradeoff
- Linearized approach optimal for pure MA/SMA models

#### Breaking Changes

- **Removed Direct Nonlinear PMM2** - Proved unstable in research (17× worse MSE)
- Previous default behavior preserved with `pmm2_variant = "unified_global"`

#### API Changes

```r
# Old way (still works, uses unified_global by default)
ar_pmm2(y, order = 2)

# New explicit variant selection
ar_pmm2(y, order = 2, pmm2_variant = "unified_iterative")
ma_pmm2(y, order = 1, pmm2_variant = "linearized")  # Best for MA
arima_pmm2(y, order = c(1,0,1), pmm2_variant = "unified_global")  # Default
```

### Documentation

- Updated all function documentation with `pmm2_variant` parameter
- Added comparison table of PMM2 variants to README
- New vignette examples demonstrating variant selection
- Research reports documenting Monte Carlo validation

### Dependencies

- Added `numDeriv` to Suggests for numerical Jacobian computation

### Bug Fixes

- Fixed convergence issues in mixed SARIMA models
- Improved moment estimation for small samples
- Enhanced error messages for degenerate cases

### Performance

- One-step variant: ~50% faster than iterative
- Numerical Jacobian: minimal overhead (<10%) when analytical unavailable
- Memory usage optimized for large time series (n > 1000)

---

## Version 0.1.4 (Development - Superseded by 0.2.0)

### New Features

- **EstemPMM-style PMM2 Estimator for MA/SMA Models** - Advanced parameter estimation for moving average components
  - New `ma_method` argument in `sarima_pmm2()` with options `"mle"` (default) and `"pmm2"`
  - `estpmm_style_ma()` - PMM2 estimator for pure MA(q) models using CSS residuals as fixed regressors
  - `estpmm_style_sma()` - PMM2 estimator for pure SMA(Q) models
  - **`estpmm_style_ma_sma()` - PMM2 estimator for mixed MA+SMA models** ⭐ NEW
  - Full support for MA(q)+SMA(Q) combinations in `sarima_pmm2()` with `ma_solver="pmm2"`
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

    *   Fixed residuals padding in `sar_pmm2`, `sarma_pmm2`, and `sarima_pmm2` to prevent length mismatch errors.
    *   Fixed S4 class definitions to ensure proper method dispatch.
    *   Corrected multiplicative SAR specification in tests.
    *   **Major Improvement**: Enhanced `estpmm_style_ma_sma` to support Multiplicative SARIMA models by including interaction terms in the design matrix. This resolves efficiency issues for mixed MA+SMA models at small sample sizes.
    *   **New Feature**: Extended PMM2 support to full SARIMA models (AR+MA+SAR+SMA) with multiplicative interactions, demonstrating improved efficiency over MLE.


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
- Comparison functions: `compare_with_ols()`, `compare_ts_methods()`, `compare_ar_methods()`, `compare_ma_solvers()`, `compare_arma_solvers()`, `compare_arima_solvers()`

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

### Roadmap

**1.0.0 (Stable API):**
- API stabilization and backward compatibility guarantee
- Seasonal PMM3 models (sar_pmm3, sarima_pmm3)
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
