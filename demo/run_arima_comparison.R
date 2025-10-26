# Quick script to run ARIMA comparison on oil price data
# This generates the full HTML report

library(rmarkdown)

# Render the R Markdown report
output_file <- render(
  "arima_oil_comparison.Rmd",
  output_format = "html_document",
  output_dir = ".",
  quiet = FALSE
)

cat(sprintf("\nReport generated: %s\n", output_file))
cat("Open it in your web browser to view the results.\n")

# Also save as PDF if LaTeX is available
tryCatch({
  pdf_file <- render(
    "arima_oil_comparison.Rmd",
    output_format = "pdf_document",
    output_dir = ".",
    quiet = TRUE
  )
  cat(sprintf("PDF version: %s\n", pdf_file))
}, error = function(e) {
  cat("PDF generation skipped (LaTeX not available)\n")
})
