# LaTeX Documentation for PMM2-ARIMA Paper

This directory contains the LaTeX source files for the research paper:

**"Application of the Polynomial Maximization Method for ARIMA Parameter Estimation Under Non-Gaussian Innovations"**

## 📁 Directory Structure

```
LaTeX/
├── PMM2_ARIMA_Part1.tex    # Part 1: Introduction & Literature Review
├── references.bib           # BibTeX bibliography (40+ entries)
├── compile.sh              # Compilation script
├── README.md               # This file
└── templates.tex           # Template reference (not used directly)
```

## 📄 Current Status

### ✅ Completed Sections

**Part 1: Introduction and Literature Review** (`PMM2_ARIMA_Part1.tex`)
- Bilingual abstract (Ukrainian/English)
- Keywords (both languages)
- Section 1: Introduction
  - 1.1 Motivation and Problem Statement
  - 1.2 Limitations of Classical Methods
  - 1.3 Existing Approaches: Short Review
  - 1.4 Polynomial Maximization Method: Alternative Approach
  - 1.5 Research Gap and Contributions
  - 1.6 Paper Structure

### 🚧 To Be Completed

- **Part 2:** Methodology (Section 2)
- **Part 3:** Empirical Results (Section 3)
- **Part 4:** Discussion (Section 4)
- **Part 5:** Conclusion (Section 5)
- **Appendices:** Technical details, proofs, additional tables

## 🔧 Compilation Instructions

### Prerequisites

Ensure you have a LaTeX distribution installed:
- **Linux:** `sudo apt-get install texlive-full`
- **macOS:** Install MacTeX from https://www.tug.org/mactex/
- **Windows:** Install MiKTeX from https://miktex.org/

Required packages (should be included in full distributions):
- `babel` (Ukrainian and English support)
- `amsmath`, `amsthm`, `amsfonts`, `amssymb` (math support)
- `hyperref` (cross-references and links)
- `cite` (bibliography)

### Method 1: Using the Compilation Script (Recommended)

```bash
cd LaTeX
./compile.sh PMM2_ARIMA_Part1
```

The script will:
1. Run `pdflatex` (first pass)
2. Run `bibtex` to process bibliography
3. Run `pdflatex` (second pass)
4. Run `pdflatex` (final pass)
5. Optionally clean up auxiliary files
6. Optionally open the PDF

### Method 2: Manual Compilation

```bash
cd LaTeX
pdflatex PMM2_ARIMA_Part1.tex
bibtex PMM2_ARIMA_Part1
pdflatex PMM2_ARIMA_Part1.tex
pdflatex PMM2_ARIMA_Part1.tex
```

### Method 3: Using latexmk (if available)

```bash
cd LaTeX
latexmk -pdf PMM2_ARIMA_Part1.tex
```

## 📚 Bibliography

The `references.bib` file contains **40+ references** organized by category:

1. **Books and Monographs** (6 entries)
   - Box & Jenkins (2015), Hyndman & Athanasopoulos (2021)
   - Kunchenko works (1991, 2002, 2006)

2. **Financial Time Series** (5 entries)
   - Studies on fat tails, skewness, and heavy-tailed distributions

3. **Recent Studies (2024-2025)** (3 entries)
   - SARIMAX with skew-normal errors
   - ARFIMA-ANN hybrid models

4. **Methodological Papers** (20+ entries)
   - Robust methods (M-estimators, M-Whittle)
   - Quantile regression and LAD
   - Bayesian approaches
   - PMM applications and theory

All entries include DOI or arXiv links where available.

## 📊 Document Statistics

**Part 1 Statistics:**
- Pages: ~15-20 (estimated)
- Sections: 1 main section, 6 subsections
- Citations: 35+ unique references
- Equations: 1 major formula (relative efficiency)
- Lists: Multiple itemized and enumerated lists

## 🔍 Key Features

### Mathematical Content
- Properly formatted equations using `amsmath`
- Theorem environments defined but not yet used
- References to equations using `\ref` and `\label`

### Multilingual Support
- Full Ukrainian and English abstracts
- Proper UTF-8 encoding
- `babel` package for language switching

### Hyperlinks
- Colored links for citations (blue)
- URLs for online resources (cyan)
- Internal cross-references (blue)

## ⚠️ Known Issues

1. **Missing Sections:** Parts 2-5 not yet created
2. **Bibliography Warnings:** Some entries may generate warnings during first compilation (normal)
3. **Page Numbers:** Table of contents page numbering will be finalized once all sections are complete

## 🔨 Troubleshooting

### "Undefined control sequence" errors
- Ensure all required packages are installed
- Check for typos in LaTeX commands

### Bibliography not appearing
- Ensure `references.bib` is in the same directory
- Run the full compilation sequence (pdflatex → bibtex → pdflatex × 2)

### Ukrainian characters not displaying
- Ensure UTF-8 encoding
- Install Ukrainian language support: `sudo apt-get install texlive-lang-cyrillic`

### Citations showing as [?]
- Run `bibtex` and then `pdflatex` at least twice more
- Check that citation keys in `.tex` match those in `.bib`

## 📝 Citation Format

Citations use the `\cite{}` command:

```latex
% Single citation
\cite{box2015time}

% Multiple citations (comma-separated)
\cite{kunchenko2002polynomial,kunchenko2006stochastic}
```

## 🚀 Next Steps

1. **Part 2: Methodology**
   - Mathematical formulation of ARIMA models
   - PMM2 algorithm description
   - Asymptotic theory

2. **Part 3: Empirical Results**
   - Monte Carlo simulation design
   - Results tables and figures
   - Efficiency comparisons

3. **Part 4: Discussion**
   - Interpretation of results
   - Practical recommendations
   - Limitations

4. **Part 5: Conclusion**
   - Summary of contributions
   - Future research directions

## 📧 Contact

For questions or issues related to this LaTeX documentation:
- **Author:** Serhii Zabolotnii
- **Email:** zabolotniua@gmail.com
- **Institution:** Industrial Research Institute for Automation and Measurements PIAP, Warsaw, Poland

## 📜 License

This LaTeX source is part of the PMM2-ARIMA research project.

---

**Generated:** 2025-10-28
**Last Updated:** 2025-10-28
**Status:** Part 1 Complete, Parts 2-5 In Progress
