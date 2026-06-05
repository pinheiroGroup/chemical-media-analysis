"""
Convert BW25113 chemical-media growth curve Excel files to GUIbiont-compatible
experiment folders.

Source: figshare DOI 10.1038/s41597-025-05356-3
Input:  BW25113_Growth_Round0{1-7}.xlsx (location via --xlsx-dir, $DATA_SRC, or default)
Output: data/bw25113_round0{1-7}/data_channel_1.csv
        data/bw25113_round0{1-7}/annotation_clean.csv

Each output folder follows the GUIbiont Clean_data layout:
  - data_channel_1.csv: Time_h, Curve00001, Curve00002, ...
  - annotation_clean.csv: one row per curve (no header), `name,g,,,,,`
    where "g" marks the well as a real (non-blank, non-excluded) growth curve.
Drop these folders into $GUIBIONT/Clean_data/ to make the experiments visible
in the Batch Fit / Plot Growth tabs.
"""

import argparse
import csv
import os
import sys
from pathlib import Path

try:
    import openpyxl
except ImportError:
    sys.exit("openpyxl is required: pip install openpyxl")

DEFAULT_XLSX_DIR = Path(__file__).parent.parent / "chemical-media-dataset" / "xlsx_raw"
OUT_DIR = Path(__file__).parent / "data"
ROUNDS  = range(1, 8)


def exp_name(round_num: int) -> str:
    return f"bw25113_round{round_num:02d}"


# GUIbiont's batch fitter requires at least this many non-empty OD readings
# per curve (src/routes/fitting.jl: `length(valid_indices) < 10`). Curves
# below the threshold are dropped here so they don't show up as failures.
MIN_DATA_POINTS = 10

# Drop curves whose OD never changes by more than this — blanks, dead cultures,
# instrument-noise-only wells. No growth signal to fit, and they'd just
# contaminate per-condition means in the ML pipeline.
DEFAULT_MIN_AMPLITUDE = 0.05


def _is_numeric(v) -> bool:
    if v is None or v == "":
        return False
    try:
        float(v)
        return True
    except (TypeError, ValueError):
        return False


def convert_round(round_num: int, xlsx_dir: Path, min_amplitude: float = DEFAULT_MIN_AMPLITUDE) -> None:
    src     = xlsx_dir / f"BW25113_Growth_Round{round_num:02d}.xlsx"
    exp_dir = OUT_DIR / exp_name(round_num)
    data_dst = exp_dir / "data_channel_1.csv"
    ann_dst  = exp_dir / "annotation_clean.csv"

    print(f"  {src.name} → {exp_dir.name}/", end="", flush=True)

    wb = openpyxl.load_workbook(src, read_only=True, data_only=True)
    ws = wb.worksheets[0]
    rows = list(ws.iter_rows(values_only=True))
    wb.close()

    header = list(rows[0])
    header[0] = "Time_h"
    data_rows = [r for r in rows[1:] if r[0] is not None]

    # Per-column filtering. Two reasons to drop a curve column:
    #   1. Fewer than MIN_DATA_POINTS numeric readings (empty/very-short).
    #   2. Amplitude (max - min) below `min_amplitude` (flat / blank / dead).
    keep_idx = [0]  # always keep the time column
    n_empty  = 0
    n_flat   = 0
    for j in range(1, len(header)):
        vals = [float(r[j]) for r in data_rows if j < len(r) and _is_numeric(r[j])]
        if len(vals) < MIN_DATA_POINTS:
            n_empty += 1
            continue
        if min_amplitude > 0 and (max(vals) - min(vals)) < min_amplitude:
            n_flat += 1
            continue
        keep_idx.append(j)

    n_dropped = n_empty + n_flat

    filtered_header = [header[j] for j in keep_idx]
    filtered_rows   = [[r[j] if j < len(r) else "" for j in keep_idx] for r in data_rows]

    exp_dir.mkdir(parents=True, exist_ok=True)
    with data_dst.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(filtered_header)
        writer.writerows(filtered_rows)

    curve_names = filtered_header[1:]
    with ann_dst.open("w", newline="") as f:
        writer = csv.writer(f)
        for name in curve_names:
            writer.writerow([name, "g", "", "", "", "", ""])

    n_curves = len(curve_names)
    n_tp     = len(data_rows)
    t_max    = data_rows[-1][0]
    drop_parts = []
    if n_empty: drop_parts.append(f"{n_empty} empty")
    if n_flat:  drop_parts.append(f"{n_flat} flat")
    drop_note = f", dropped {' + '.join(drop_parts)}" if drop_parts else ""
    print(f"  ({n_curves} curves × {n_tp} timepoints, t_max={t_max}h{drop_note})")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument(
        "--xlsx-dir",
        type=Path,
        default=Path(os.environ.get("DATA_SRC", DEFAULT_XLSX_DIR)),
        help="Directory containing BW25113_Growth_Round0{1-7}.xlsx "
             "(default: $DATA_SRC, else ../chemical-media-dataset/xlsx_raw)",
    )
    parser.add_argument(
        "--min-amplitude",
        type=float,
        default=DEFAULT_MIN_AMPLITUDE,
        help=f"Drop curves with (max - min) OD below this. Catches blanks / dead cultures. "
             f"Default {DEFAULT_MIN_AMPLITUDE}; pass 0 to keep everything.",
    )
    args = parser.parse_args()
    xlsx_dir = args.xlsx_dir.expanduser().resolve()

    print(f"Reading from: {xlsx_dir}")
    print(f"Writing to:  {OUT_DIR}")
    print(f"Min amplitude: {args.min_amplitude} (0 = no filter)\n")

    if not xlsx_dir.exists():
        sys.exit(f"Source directory not found: {xlsx_dir}")

    for r in ROUNDS:
        convert_round(r, xlsx_dir, min_amplitude=args.min_amplitude)

    print("\nDone.")


if __name__ == "__main__":
    main()
