"""
Convert BW25113 chemical-media growth curve Excel files to GUIbiont-compatible CSVs.

Source: figshare DOI 10.1038/s41597-025-05356-3
Input:  ../chemical-media-dataset/xlsx_raw/BW25113_Growth_Round0{1-7}.xlsx
Output: data/Round0{1-7}.csv  (Time_h, Curve00001, Curve00002, ...)

Each CSV follows the GUIbiont format:
  - First column: Time_h (numeric, hours)
  - Remaining columns: one growth curve per column (header = curve label)
  - No trailing empty rows
"""

import csv
import sys
from pathlib import Path

try:
    import openpyxl
except ImportError:
    sys.exit("openpyxl is required: pip install openpyxl")

XLSX_DIR = Path("/media/aivuk/64fce268-2613-4033-b39f-537ae2d28805/roms/pinheiroTech/chemical-media-dataset/xlsx_raw")
OUT_DIR  = Path(__file__).parent / "data"
ROUNDS   = range(1, 8)


def convert_round(round_num: int) -> None:
    src = XLSX_DIR / f"BW25113_Growth_Round{round_num:02d}.xlsx"
    dst = OUT_DIR  / f"Round{round_num:02d}.csv"

    print(f"  {src.name} → {dst.name}", end="", flush=True)

    wb = openpyxl.load_workbook(src, read_only=True, data_only=True)
    ws = wb.worksheets[0]

    rows = list(ws.iter_rows(values_only=True))
    wb.close()

    # Header: rename "Time (h)" → "Time_h" for safe column naming
    header = list(rows[0])
    header[0] = "Time_h"

    # Drop trailing rows where time is None (Round03 has ~900 empty rows)
    data_rows = [r for r in rows[1:] if r[0] is not None]

    dst.parent.mkdir(parents=True, exist_ok=True)
    with dst.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(data_rows)

    n_curves = len(header) - 1
    n_tp     = len(data_rows)
    t_max    = data_rows[-1][0]
    print(f"  ({n_curves} curves × {n_tp} timepoints, t_max={t_max}h)")


def main() -> None:
    print(f"Reading from: {XLSX_DIR}")
    print(f"Writing to:  {OUT_DIR}\n")

    if not XLSX_DIR.exists():
        sys.exit(f"Source directory not found: {XLSX_DIR}")

    for r in ROUNDS:
        convert_round(r)

    print("\nDone.")


if __name__ == "__main__":
    main()
