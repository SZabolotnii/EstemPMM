# CRAN Submission Comments

## Submission Type
This is a resubmission addressing feedback from CRAN auto-check service.

## Changes from previous submission
* Fixed URL in README.md: changed GNU license URL from www.gnu.org to
  opensource.org to avoid timeout issues
* Confirmed LICENSE file is properly referenced in DESCRIPTION (GPL-3 + file LICENSE)
* Verified CRAN_SUBMISSION_CHECKLIST.md is correctly excluded via .Rbuildignore

## Test environments
* Local: macOS Sequoia 15.6.1, R 4.5.1
* Win-builder (r-devel and r-release)
* GitHub Actions (via R-CMD-check):
  - Ubuntu (latest), R release
  - Ubuntu (latest), R devel
  - Windows (latest), R release
  - macOS (latest), R release

## R CMD check results
0 errors | 0 warnings | 3 notes

### Notes:
1. "New submission" - This is expected for a new package
2. "unable to verify current time" - This is a system-specific note and does not affect package functionality
3. "Possibly misspelled words in DESCRIPTION: PMM" - PMM is an abbreviation for "Polynomial Maximization Method",
   which is the core statistical method implemented in this package. It is properly defined in the Description field
   and is used consistently throughout the package documentation and scientific literature.

## Downstream dependencies
There are currently no downstream dependencies for this package.

## Additional Notes
EstemPMM implements the Polynomial Maximization Method (PMM) for parameter
estimation in regression models with non-Gaussian errors. The package provides
an alternative to ordinary least squares that achieves lower variance estimates
when errors exhibit significant skewness.

The method is based on peer-reviewed scientific publications:
- Kunchenko & Lega (1992) - foundational reference
- Zabolotnii et al. (2018, 2022, 2023) - applications to regression and time series
