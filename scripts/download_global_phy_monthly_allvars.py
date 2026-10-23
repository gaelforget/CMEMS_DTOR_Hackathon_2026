#!/usr/bin/env python3
"""
Download the native, complete global monthly files for:
GLOBAL_MULTIYEAR_PHY_001_030

Years requested: 2023–2026

This uses copernicusmarine.get(), not subset(), so each downloaded file
retains the full global grid, variables, and vertical levels provided
by Copernicus Marine.
"""

from pathlib import Path
import copernicusmarine

# Monthly global-physics reanalysis dataset
DATASET_ID = "cmems_mod_glo_phy_my_0.083deg_P1M-m"

# Directory in which the original Copernicus file structure is retained
OUTPUT_DIR = Path("/orcd/compute/gforget/001/GLORYS")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# The request below searches the dataset's internal yearly directories.
# Existing downloaded files are normally skipped automatically.
for year in range(2023, 2027):
    print(f"\nDownloading monthly global files for {year} ...")

    try:
        copernicusmarine.get(
            dataset_id=DATASET_ID,
            filter=f"*/{year}/*",
            output_directory=str(OUTPUT_DIR),
        )

    except Exception as exc:
        print(f"Download failed or no files were available for {year}:")
        print(exc)

print("\nDone.")
