# CRAN Submission Comments -- EstemPMM 0.3.2

## Release summary

This candidate updates EstemPMM after the 0.3.1 CRAN release. The main changes
are:

- PMM3 documentation and examples now consistently describe the supported
  symmetric platykurtic regime.
- Seasonal SARIMA vignette examples were aligned with the current
  `sarima_pmm2()` API.
- S4 generic coverage was completed for `BIC`, `logLik`, and `nobs` methods
  where the package provides Gaussian approximate information criteria.

## Test environments

| Platform                    | R version | Result                     |
| --------------------------- | --------- | -------------------------- |
| macOS Tahoe 26.4.1 (arm64) | 4.5.3     | 0 ERROR, 0 WARNING, 2 NOTE |

## R CMD check results

```
Status: 0 ERROR, 0 WARNING, 2 NOTEs
```

The two local NOTEs are environment/toolchain notes:

1. `unable to verify current time` during future-timestamp checking.
2. HTML manual validation/math rendering skipped because local HTML Tidy is not
   recent enough and package `V8` is unavailable.

No ERRORs or WARNINGs were reported.

## Reverse dependencies

None. This package has no reverse dependencies on CRAN.
