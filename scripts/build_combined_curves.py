"""
Interpolate all 7 per-round growth curve CSVs to a common time grid and
concatenate into a single file for upload to GUIbiont's Clustering tab.

Different rounds have different numbers of timepoints and slightly different
time axes.  This script linearly interpolates each curve to a shared 97-point
grid spanning 0–48 h so GUIbiont receives a uniform matrix.

Input:  data/Round01.csv … data/Round07.csv
Output: results/all_curves_combined.csv
        Columns: Time_h, <curve labels from all rounds>
"""

import sys
from pathlib import Path

try:
    import numpy as np
    import pandas as pd
except ImportError:
    sys.exit("numpy and pandas are required: pip install numpy pandas")

DATA_DIR    = Path(__file__).parent.parent / "data"
RESULTS_DIR = Path(__file__).parent.parent / "results"
ROUNDS      = range(1, 8)

T_START  = 0.0
T_END    = 48.0
N_POINTS = 97


def interpolate_round(df: pd.DataFrame, t_grid: np.ndarray) -> pd.DataFrame:
    t = df["Time_h"].values
    curves = df.drop(columns=["Time_h"])
    out = {}
    for col in curves.columns:
        y = pd.to_numeric(curves[col], errors="coerce").values
        # Only interpolate over the range where we have data
        valid = ~np.isnan(y)
        if valid.sum() < 2:
            out[col] = np.full(len(t_grid), np.nan)
        else:
            out[col] = np.interp(t_grid, t[valid], y[valid],
                                 left=np.nan, right=np.nan)
    return pd.DataFrame(out)


def main() -> None:
    t_grid = np.linspace(T_START, T_END, N_POINTS)

    all_curves = []
    total = 0
    for r in ROUNDS:
        path = DATA_DIR / f"Round{r:02d}.csv"
        if not path.exists():
            print(f"  WARNING: {path.name} not found — skipping round {r:02d}")
            continue
        df = pd.read_csv(path)
        interp = interpolate_round(df, t_grid)
        all_curves.append(interp)
        total += len(interp.columns)
        print(f"  Round{r:02d}: {len(interp.columns)} curves interpolated")

    if not all_curves:
        sys.exit("No round CSVs found in data/. Run preprocess.py first.")

    combined = pd.concat(all_curves, axis=1)
    combined.insert(0, "Time_h", t_grid)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    out = RESULTS_DIR / "all_curves_combined.csv"
    combined.to_csv(out, index=False)
    print(f"\nCombined: {total} curves × {N_POINTS} timepoints (0–{T_END} h)")
    print(f"Written to {out}")


if __name__ == "__main__":
    main()
