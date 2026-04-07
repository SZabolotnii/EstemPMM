# CRAN Submission Comments — EstemPMM 0.3.1

## Resubmission

This is a resubmission of EstemPMM 0.3.0, addressing the 2 NOTEs from
incoming pre-tests (2026-03-20). No code changes — only packaging fixes.

### NOTE 1: Possibly misspelled word "platykurtic" in DESCRIPTION

"Platykurtic" is a standard statistical term describing distributions with
negative excess kurtosis (e.g., uniform, beta-symmetric). It appears in
statistics textbooks (e.g., Westfall 2014, "Kurtosis as Peakedness") and
R documentation (e.g., `?moments::kurtosis`). `inst/WORDLIST` includes it.

### NOTE 2: Non-standard file/directory found at top level: 'figure'

The `figure/` directory (vignette build artefacts) was already excluded via
`.Rbuildignore` and is absent from the built tarball. Additionally, we now
exclude `CLAUDE.md` and `*.tar.gz` from `.Rbuildignore` to prevent any
future top-level NOTEs.

## Test environments

| Platform                       | R version | Result                    |
| ------------------------------ | --------- | ------------------------- |
| macOS Tahoe 26.3.1 (arm64)    | 4.5.2     | 0 ERROR, 0 WARNING, 0 NOTE |
| win-builder                    | R-devel   | 0 ERROR, 0 WARNING, 1 NOTE (spelling only) |
| Debian (CRAN incoming pretest) | R-devel   | 0 ERROR, 0 WARNING, 1 NOTE (spelling only) |

The remaining NOTE ("platykurtic") is a false positive — the word is a
valid statistical term, documented in `inst/WORDLIST`.

## R CMD check results

```
Status: 0 ERROR, 0 WARNING, 0 NOTE
```

## Reverse dependencies

None. This package has no reverse dependencies on CRAN.
