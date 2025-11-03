#!/usr/bin/env Rscript
# Script to convert PNG plots to PDF format for better LaTeX quality

library(magick)

# Input and output directories
input_dir <- "figures"
output_dir <- "figures"

# Get all PNG files
png_files <- list.files(input_dir, pattern = "\\.png$", full.names = TRUE)

cat("Converting", length(png_files), "PNG files to PDF...\n")

for (png_file in png_files) {
  # Create output filename
  pdf_file <- sub("\\.png$", ".pdf", png_file)

  cat("Converting:", basename(png_file), "->", basename(pdf_file), "\n")

  # Read PNG and write as PDF
  img <- image_read(png_file)
  image_write(img, path = pdf_file, format = "pdf")
}

cat("\nConversion complete! All files saved to:", output_dir, "\n")
