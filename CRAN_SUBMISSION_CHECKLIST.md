# CRAN Submission Checklist — EstemPMM 0.1.3 (2025-11-14)

## ✅ Completed

### Package Skeleton
- DESCRIPTION updated to 0.1.3 (Date: 2025-11-13) with correct Authors@R, URL, BugReports.
- GPL-3 LICENSE referenced; no extra LICENSE file required by CRAN.
- NAMESPACE generated via roxygen2.
- README (EN + UA) refreshed with seasonal-model content and installation instructions.
- NEWS.md documents v0.1.3 (documentation, seasonal vignette, Monte Carlo evidence).
- cran-comments.md template ready for the next submission cycle.

### Documentation
- All roxygen blocks are in English; `man/` regenerated on 2025-11-14.
- Vignettes: `pmm2_introduction`, `pmm2_time_series`, `bootstrap_inference` knit successfully.
- pkgdown-style articles live in `docs/` but are excluded by `.Rbuildignore`.

### Testing / QA
- `tests/testthat/` covers linear PMM2, time-series wrappers, inference helpers, and Monte Carlo utilities.
- `devtools::check()` (local, macOS 15.6.1, R 4.5.1) passes with 0 ERROR / 0 WARN / 2 NOTE ("new submission", "possibly mis-spelled word PMM").
- GitHub Actions workflow (`.github/workflows/R-CMD-check.yaml`) exercised on latest push (Ubuntu + macOS + Windows) — all green.

### Housekeeping
- `.Rbuildignore` excludes `.github`, `docs`, `test_results`, `*.Rproj`, submission notes, and large PDFs.
- `inst/doc/` kept empty between builds; vignettes regenerated during `R CMD build`.
- Demo directory contains only `.R` + `00Index`.

## 🔄 Final Verification Steps

1. **Regenerate docs & vignettes**
   ```r
   devtools::document()
   devtools::build_vignettes()
   ```

2. **Build & check**
   ```bash
   R CMD build .
   R CMD check --as-cran EstemPMM_0.1.3.tar.gz
   ```

3. **Record results** in `cran-comments.md` (environments, NOTES, win-builder/rhub links).

4. **Optional external checks**
   ```r
   devtools::check_win_devel()
   rhub::check_for_cran()
   ```

5. **Submission package**
   - Tarball: `EstemPMM_0.1.3.tar.gz`
   - Supporting files: `cran-comments.md`, `README.md`/`README_uk.md`, `NEWS.md`.
   - Upload via https://cran.r-project.org/submit.html and confirm via email within 24h.

## 📋 Quick Tick Boxes Before Upload
- [ ] `R CMD check --as-cran` clean on macOS
- [ ] `check_win_devel()` success email saved
- [ ] Vignettes knit on clean R session (`rmarkdown`, `knitr`, `ggplot2`, `dplyr` installed)
- [ ] `cran-comments.md` updated with actual dates + log output
- [ ] `NEWS.md` references the same version/date as DESCRIPTION
- [ ] No stray files under  `inst/doc`, `docs`, `test_results` in tarball
- [ ] Git tag created (`v0.1.3`) and pushed after CRAN acceptance

_Last reviewed: 2025-11-14_
