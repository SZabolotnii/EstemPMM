# CRAN Submission Comments — EstemPMM 0.1.3

## Summary
- **Submission type:** update (documentation + seasonal Monte Carlo refresh).  
- **Focus:** clarified SMA benchmark methodology, refreshed English documentation, added SAR/SMA vignettes, and ensured CRAN-ready instructions.

## Test environments
| Platform | R version | Notes |
| --- | --- | --- |
| macOS Sequoia 15.6.1 (arm64) | 4.5.1 | `R CMD check --as-cran EstemPMM_0.1.3.tar.gz` |
| win-builder | r-devel, r-release | Submitted via `devtools::check_win_devel()` / `check_win_release()` |
| GitHub Actions | Ubuntu 22.04 (release & devel), macOS 14, Windows Server 2022 | Workflow `.github/workflows/R-CMD-check.yaml` |

## R CMD check results (local tarball)
```
* using log directory '.../EstemPMM.Rcheck'
* using R version 4.5.1
* using platform: aarch64-apple-darwin24.4.0 (64-bit)
* checking for file 'EstemPMM/DESCRIPTION' ... OK
* checking CRAN incoming feasibility ... NOTE
* checking for non-standard things in the check directory ... OK
* checking for detritus in the temp directory ... OK
* DONE
Status: 0 ERROR, 0 WARNING, 2 NOTE
```

### Notes explained
1. **"Unable to verify current time"** — cosmetic note on macOS when NTP is sandboxed; timestamps inside the tarball are correct.
2. **HTML manual validation skipped (tidy/V8)** — the system copy of HTML Tidy/V8 bundled with macOS is older than CRAN expects. The PDF manual builds cleanly, and CRAN’s own infrastructure runs the full HTML validation.

No other notes/warnings/errors were reported locally or on win-builder.

## Reverse dependencies
None (new package).

## Additional context
- Canonical SMA Monte Carlo benchmark (n = 120, γ-innovations) remains available via `run_sma_monte_carlo.R`. The broader `monte_carlo_seasonal_comparison.R` script is documented in README/NEWS so reviewers understand its longer runtime.
- All documentation (README, README_uk, vignettes, man pages) was regenerated on 2025-11-14 via `devtools::document()` and `devtools::build_vignettes()`.

Please let me know if further changes are required.
