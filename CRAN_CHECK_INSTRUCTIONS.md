# CRAN Submission: Current Checklist (v0.1.3)

The translation work called out in earlier notes has been completed. All roxygen blocks, vignettes, and README files are now in English, so the remaining CRAN-prep focus is on **validation**, **documentation freshness**, and **submission hygiene**.

## 1. Align Monte Carlo Evidence
- Canonical SMA benchmark lives in `run_sma_monte_carlo.R` (n = 120, γ-innovations). Keep this script as the reference when CRAN asks about statistical validation.  
- The combined script `monte_carlo_seasonal_comparison.R` is broader (multiple n and scenarios). Document both in README/NEWS so reviewers understand the longer run time.

## 2. Regenerate Documentation
```r
library(devtools)
document()      # roxygen → man/
build_vignettes()
```
Make sure `inst/doc/` is cleaned before running `build_vignettes()` (the directory is auto-populated by `R CMD build`).

## 3. Build + Check
```bash
R CMD build .
R CMD check --as-cran EstemPMM_0.1.3.tar.gz
```
Expected outcome: 0 errors, 0 warnings, ≤2 notes ("new submission" + possible spelling note for "PMM").

## 4. GitHub Actions / External Checks
- Push to GitHub to trigger `.github/workflows/R-CMD-check.yaml` (Ubuntu + macOS + Windows).  
- Optional but recommended: `devtools::check_win_devel()` and `rhub::check_for_cran()`.

## 5. Update Submission Artifacts
- `README.md` / `README_uk.md`: include CRAN install snippet (`install.packages("EstemPMM")`) and mention how to rebuild documentation/tests.
- `NEWS.md`: ensure 0.1.3 section highlights documentation refresh + seasonal Monte Carlo evidence.
- `cran-comments.md`: describe the environments and provide a short bullet list of notes (see template in file).
- `CRAN_SUBMISSION_CHECKLIST.md`: mark the new version + date.

## 6. Pre-Submission Clean-up
- Delete `inst/doc/*`, `vignettes/*.html`, and any `.DS_Store` before building.  
- Ensure `.Rbuildignore` excludes `docs/`, `test_results/`, `.github/`, Monte Carlo artifacts, and large PDFs.  
- Confirm no files > 5 MB remain under version control (CRAN auto-check complains about bloated tarballs).

## 7. Suggested Release Order
1. `document()` + `build_vignettes()`  
2. `devtools::check()` (quick sanity)  
3. `R CMD build`  
4. `R CMD check --as-cran` on the tarball  
5. Update `cran-comments.md` with actual check output  
6. Tag release / create `EstemPMM_0.1.3.tar.gz` for upload  
7. Submit via https://cran.r-project.org/submit.html (attach tarball + `cran-comments.md`)

## 8. Useful Resources
- CRAN policy: https://cran.r-project.org/web/packages/policies.html  
- Submission checklist: https://cran.r-project.org/web/packages/submission_checklist.html  
- R Packages book: https://r-pkgs.org/

_Last updated: 2025-11-14_
