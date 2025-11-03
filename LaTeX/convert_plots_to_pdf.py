#!/usr/bin/env python3
"""Convert PNG plots to PDF format for better LaTeX quality"""

import os
from pathlib import Path

try:
    from PIL import Image
except ImportError:
    print("PIL not available, trying to install...")
    os.system("pip3 install pillow")
    from PIL import Image

# Directories
figures_dir = Path("figures")
png_files = sorted(figures_dir.glob("*.png"))

print(f"Converting {len(png_files)} PNG files to PDF...\n")

for png_file in png_files:
    pdf_file = png_file.with_suffix('.pdf')

    print(f"Converting: {png_file.name} -> {pdf_file.name}")

    # Open PNG and convert to RGB (PDF doesn't support RGBA)
    img = Image.open(png_file)
    if img.mode == 'RGBA':
        # Create white background
        rgb_img = Image.new('RGB', img.size, (255, 255, 255))
        rgb_img.paste(img, mask=img.split()[3])  # Use alpha channel as mask
        img = rgb_img
    elif img.mode != 'RGB':
        img = img.convert('RGB')

    # Save as PDF
    img.save(pdf_file, 'PDF', resolution=300.0)

print(f"\n✓ Conversion complete! All files saved to: {figures_dir}")
