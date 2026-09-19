# CRAN Submission Comments -- EstemPMM 0.5.0

## Release summary

This is an update to the CRAN release 0.4.0. It adds a second solver for the
moving-average component of the PMM2 time-series estimators and the inference
that goes with it.

- New argument `ma_solver = c("linearized", "recursive")` on `ts_pmm2()`,
  `ma_pmm2()`, `arma_pmm2()` and `arima_pmm2()`. Up to 0.4.0 the design matrix
  for the MA columns was built once from the first-stage residuals and held
  fixed, which is a one-step linearisation of the innovation recursion and
  costs a factor of `1 - theta^2` in asymptotic variance. The `"recursive"`
  solver recomputes the innovations and the score regressors by the exact
  recursion at every candidate parameter, and that penalty disappears.
- The default remains `"linearized"`, so results obtained with 0.4.0 are
  reproduced unchanged unless the user opts in.
- Closed-form sandwich `vcov()` and `confint()` methods for recursive MA fits.
- The argument formerly spelled `ma_method` in the development sources was
  renamed to `ma_solver` before release, to avoid a partial-matching collision
  with `sarima_pmm2()`. No released version exported `ma_method`, so this
  breaks no published interface.
- 365 new test assertions covering the recursive branch.

## Test environments

| Platform                   | R version | Result                     |
| -------------------------- | --------- | -------------------------- |
| macOS 27.0 (arm64)         | 4.5.3     | 0 ERROR, 0 WARNING, 1 NOTE |

Checked on 2026-09-19 against the tarball being submitted.

## R CMD check results

`R CMD check --as-cran` with `_R_CHECK_CRAN_INCOMING_=true` and
`_R_CHECK_CRAN_INCOMING_REMOTE_=true`:

```
Status: 1 NOTE
* checking CRAN incoming feasibility ... OK
```

The single NOTE is a local toolchain note, not a property of the package:
HTML manual validation and math rendering were skipped because the local
HTML Tidy is not recent enough and package `V8` is unavailable. No ERRORs and
no WARNINGs were reported, and the tests, examples, `--run-donttest` examples
and vignettes all pass.

## Reverse dependencies

None. This package has no reverse dependencies on CRAN.
