# CRAN Submission: Critical Issues and Instructions

## Critical Issue Found: Ukrainian Documentation

**CRAN REQUIREMENT**: All package documentation must be in English.

### Problem
All R source files currently contain Ukrainian text in:
- Roxygen documentation comments (`#'`)
- Function parameter descriptions
- Examples
- Internal comments in user-facing functions

### Files Affected
The following R files contain Ukrainian documentation:
- `R/pmm2_package.R` - **FIXED**
- `R/pmm2_main.R`
- `R/pmm2_classes.R`
- `R/pmm2_common.R`
- `R/pmm2_inference.R`
- `R/pmm2_monte_carlo.R`
- `R/pmm2_ts_design.R`
- `R/pmm2_ts_main.R`
- `R/pmm2_ts_methods.R`
- `R/pmm2_utils.R`

## Required Actions

### 1. Translate All Documentation to English

You need to translate all roxygen documentation (`#'` comments) in the R files listed above from Ukrainian to English.

**Important**:
- Keep the same structure and roxygen tags (@param, @return, @export, etc.)
- Translate parameter descriptions, function descriptions, details, and examples
- Internal code comments (non-roxygen `#` comments) can remain in Ukrainian if they're not user-facing, but it's better to translate them too

### 2. Regenerate Documentation

After translating the R source files, regenerate the man pages:

```r
# Install roxygen2 if needed
install.packages("roxygen2")

# Regenerate documentation
roxygen2::roxygenize()

# Or use devtools
library(devtools)
document()
```

### 3. Run R CMD check --as-cran

```bash
# Build the package
R CMD build .

# Check with CRAN standards
R CMD check --as-cran EstemPMM_0.1.0.tar.gz
```

### 4. Fix Any Warnings/Errors

Review the check output and fix any issues. Common CRAN notes that are acceptable:
- "New submission" NOTE
- "Non-standard directory 'docs'" (can be ignored if docs is in .Rbuildignore)

### 5. Alternative: Use devtools

```r
# Install devtools if needed
install.packages("devtools")

# Run CRAN checks
library(devtools)
check(cran = TRUE)
```

### 6. Pre-build cleanup (required)

- Update `DESCRIPTION` with the new version number before calling `R CMD build` (e.g., 0.1.1) so the resulting tarball matches the intended release.
- Remove generated artifacts from previous vignette builds: delete the contents of `inst/doc/` (they will be regenerated automatically).
- Delete `demo/README.md` or move its content elsewhere—`demo/` may only contain `.R`/`.Rout` files and an optional `00Index`.
- Verify that only portable file names remain (ASCII letters/digits plus `_` or `.`, no spaces) to avoid "invalid file names" warnings during `R CMD check`.

## Files Already Modified

The following files have been updated for CRAN compliance:
1. **cran-comments.md** - Created with submission information
2. **.Rbuildignore** - Updated to exclude non-package files (.github, docs, etc.)
3. **R/pmm2_package.R** - Translated from Ukrainian to English
4. **.github/workflows/R-CMD-check.yaml** - Added GitHub Actions workflow for automated CRAN checks

## GitHub Actions Workflow

A GitHub Actions workflow has been added that will automatically run R CMD check --as-cran on multiple platforms:
- Ubuntu (R release, R devel, R oldrel-1)
- Windows (R release)
- macOS (R release)

This will run automatically when you push to branches matching `claude/**` or to `main`/`master`.

## Translation Example

**Before (Ukrainian)**:
```r
#' Pidhaniaie liniinu model za dopomohoiu PMM2
#'
#' @param formula Formula R dlia modeli
#' @param data data.frame, shcho mistyt zminni u formuli
#' @return Ob'iekt S4 \code{PMM2fit}
```

**After (English)**:
```r
#' Fit linear model using PMM2
#'
#' @param formula R formula for the model
#' @param data data.frame containing variables in the formula
#' @return S4 object of class \code{PMM2fit}
```

## Next Steps

1. Translate all roxygen documentation in the R files listed above
2. Run `roxygen2::roxygenize()` to regenerate man pages
3. Run `R CMD check --as-cran EstemPMM_*.tar.gz` locally
4. Fix any warnings or errors
5. Commit and push changes
6. Review GitHub Actions check results
7. When all checks pass, proceed with CRAN submission

## Useful Resources

- CRAN Repository Policy: https://cran.r-project.org/web/packages/policies.html
- Writing R Extensions: https://cran.r-project.org/doc/manuals/r-release/R-exts.html
- R Packages book: https://r-pkgs.org/
