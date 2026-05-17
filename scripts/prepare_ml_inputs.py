"""
Join batch fit results with medium chemical composition and produce the two
CSV files needed by GUIbiont's ML Analysis tab.

Workflow (mirrors ml_analysis.jl):
  1. Load results/batch_fit_results.csv  (from GUIbiont Batch Fit, combined)
  2. Load BW25113_GrowthDataEvaluation.xlsx  → curve label → ConditionID
  3. Load BW25113_Medium composition.xlsx    → ConditionID → 44 compound concentrations
  4. Join and aggregate: mean gr, exit_lag_rate, N_max per condition
  5. Write GUIbiont-ready CSVs

Usage:
  python scripts/prepare_ml_inputs.py \\
      --fits      results/batch_fit_results.csv \\
      --composition  /path/to/BW25113_Medium\\ composition.xlsx \\
      --reference    /path/to/BW25113_GrowthDataEvaluation.xlsx

Output:
  results/guibiont_ml_inputs/fit_results.csv     (condition, gr, exit_lag_rate, N_max)
  results/guibiont_ml_inputs/feature_matrix.csv  (condition, <44 compounds>)
"""

import argparse
import sys
from pathlib import Path

try:
    import numpy as np
    import pandas as pd
except ImportError:
    sys.exit("numpy and pandas are required: pip install numpy pandas openpyxl")

OUT_DIR = Path(__file__).parent.parent / "results" / "guibiont_ml_inputs"


def load_fits(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "converged" in df.columns:
        df = df[df["converged"].astype(bool)]
    # GUIbiont's batch-fit CSV uses "well" for the curve name; the join below
    # expects "label" (matches the column in BW25113_GrowthDataEvaluation.xlsx).
    if "label" not in df.columns and "well" in df.columns:
        df = df.rename(columns={"well": "label"})
    print(f"  Fit results: {len(df)} converged curves")
    return df


def load_evaluation(path: Path) -> pd.DataFrame:
    df = pd.read_excel(path)
    df.columns = df.columns.str.strip()
    # Identify label and condition columns flexibly
    label_col = next((c for c in df.columns if "label" in c.lower()), None)
    cond_col  = next((c for c in df.columns if "condition" in c.lower()), None)
    if label_col is None or cond_col is None:
        sys.exit(f"Could not find Label/Condition columns in {path.name}. "
                 f"Columns found: {list(df.columns)}")
    df = df.rename(columns={label_col: "label", cond_col: "ConditionID"})
    df["label"]       = df["label"].astype(str).str.strip()
    df["ConditionID"] = df["ConditionID"].astype(str).str.strip()
    print(f"  Evaluation: {len(df)} curve-to-condition mappings")
    return df[["label", "ConditionID"]]


def load_composition(path: Path) -> tuple[pd.DataFrame, list[str]]:
    df = pd.read_excel(path)
    # Compound headers in the figshare workbook embed a literal newline before
    # the unit ("K2HPO4\n(mM)"). Collapse internal whitespace so the resulting
    # CSV stays single-line per row (and GUIbiont's column UI renders cleanly).
    df.columns = df.columns.str.replace(r"\s+", " ", regex=True).str.strip()
    cond_col = next((c for c in df.columns if "condition" in c.lower()), df.columns[0])
    df = df.rename(columns={cond_col: "ConditionID"})
    df["ConditionID"] = df["ConditionID"].astype(str).str.strip()
    # Compound columns: everything except ConditionID and the per-curve metadata
    # columns ("Label", "Assay ID") that the figshare workbook also carries.
    meta = {"ConditionID", "Label", "Assay ID"}
    compound_cols = [c for c in df.columns if c not in meta]
    for col in compound_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    print(f"  Composition: {len(df)} conditions × {len(compound_cols)} compounds")
    return df, compound_cols


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fits",        required=True, type=Path)
    parser.add_argument("--composition", required=True, type=Path)
    parser.add_argument("--reference",  required=True, type=Path)
    args = parser.parse_args()

    print("Loading data …")
    fits    = load_fits(args.fits)
    eval_df = load_evaluation(args.reference)
    comp_df, compound_cols = load_composition(args.composition)

    print("\nJoining …")
    joined = fits.merge(eval_df, on="label", how="left").dropna(subset=["ConditionID"])
    print(f"  After label→condition join: {len(joined)} curves")

    joined = joined.merge(comp_df, on="ConditionID", how="left").dropna(subset=[compound_cols[0]])
    print(f"  After composition join: {len(joined)} curves")

    print("\nAggregating per condition …")
    param_cols = [c for c in ["gr", "exit_lag_rate", "N_max"] if c in joined.columns]
    agg_params = joined.groupby("ConditionID")[param_cols].mean().reset_index()

    # Compound concentrations are constant within a condition — take first value
    agg_comp = joined.groupby("ConditionID")[compound_cols].first().reset_index()

    cond_params = agg_params.rename(columns={"ConditionID": "condition"})
    cond_comp   = agg_comp.rename(columns={"ConditionID": "condition"})
    print(f"  {len(cond_params)} unique conditions")

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    fit_out = OUT_DIR / "fit_results.csv"
    cond_params.to_csv(fit_out, index=False)
    print(f"\nWritten: {fit_out}")

    feat_out = OUT_DIR / "feature_matrix.csv"
    cond_comp.to_csv(feat_out, index=False)
    print(f"Written: {feat_out}")

    print("\nReady for GUIbiont ML Analysis tab:")
    print(f"  Fit results CSV  : {fit_out}")
    print(f"  Feature matrix   : {feat_out}")
    print(f"  Label column     : condition")
    print(f"  Parameters       : {', '.join(param_cols)}")


if __name__ == "__main__":
    main()
