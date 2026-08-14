"""
Combine per-round batch fit CSVs into a single results/batch_fit_results.csv.

Input:  results/batch_fit_round01.csv ... results/batch_fit_round07.csv
        (downloaded from GUIbiont Batch Fit tab, one per round)
Output: results/batch_fit_results.csv
        Columns: round, label, gr, exit_lag_rate, N_max, shape, aicc, loss,
                 n_timepoints, converged
"""

import sys
from pathlib import Path

try:
    import pandas as pd
except ImportError:
    sys.exit("pandas is required: pip install pandas")

RESULTS_DIR = Path(__file__).parent.parent / "results"
ROUNDS = range(1, 8)


def main() -> None:
    frames = []
    for r in ROUNDS:
        path = RESULTS_DIR / f"batch_fit_round{r:02d}.csv"
        if not path.exists():
            print(f"  WARNING: {path.name} not found -- skipping round {r:02d}")
            continue
        df = pd.read_csv(path)
        df.insert(0, "round", f"Round{r:02d}")
        frames.append(df)
        n_conv = df["converged"].sum() if "converged" in df.columns else len(df)
        print(f"  Round{r:02d}: {len(df)} curves, {n_conv} converged")

    if not frames:
        sys.exit("No per-round CSV files found in results/. Run batch fitting in GUIbiont first.")

    combined = pd.concat(frames, ignore_index=True)

    out = RESULTS_DIR / "batch_fit_results.csv"
    combined.to_csv(out, index=False)
    n_conv = combined["converged"].sum() if "converged" in combined.columns else len(combined)
    print(f"\nCombined: {len(combined)} curves total, {n_conv} converged")
    print(f"Written to {out}")


if __name__ == "__main__":
    main()
