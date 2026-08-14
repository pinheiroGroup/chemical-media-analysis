"""
Build the per-curve feature matrix consumed by `run_ml_via_guibiont.py`.

Every retained growth curve gets one row carrying the 44 compound
concentrations of the medium it was grown in:

  BW25113_GrowthDataEvaluation.xlsx  ->  curve label -> ConditionID
  BW25113_Medium composition.xlsx    ->  ConditionID -> 44 compound columns

The two tables are joined on ConditionID and the result is written per curve
label -- replicate curves of the same medium therefore repeat that medium's
concentration vector. `run_ml_via_guibiont.py` does its own aggregation by
medium afterwards, so the matrix handed to it must be per curve, not per
condition.

This is deliberately independent of any fit-result table: it is driven purely
by the evaluation and composition workbooks, so every curve present in the
evaluation mapping gets a row regardless of whether its fit converged. That
makes it broader than `prepare_ml_inputs.py`, which filters to converged aHPM
fits and aggregates to one row per condition for the manual ML Analysis tab.

Usage:
  python scripts/build_percurve_feature_matrix.py \\
      --composition "$DATA_SRC/BW25113_Medium composition.xlsx" \\
      --reference   "$DATA_SRC/BW25113_GrowthDataEvaluation.xlsx"

Output:
  results/guibiont_ml_inputs/feature_matrix.csv  (label, <44 compounds>)
"""

import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

try:
    import pandas as pd
except ImportError:
    sys.exit("pandas is required: pip install pandas openpyxl")

# Reuse the workbook parsing from prepare_ml_inputs so the two pathways can
# never drift apart on column detection or header normalisation.
from prepare_ml_inputs import load_composition, load_evaluation

HERE     = Path(__file__).resolve().parent
REPO     = HERE.parent
OUT_DIR  = REPO / "results" / "guibiont_ml_inputs"
DEFAULT_XLSX_DIR = Path(
    os.environ.get("DATA_SRC", REPO.parent / "chemical-media-dataset" / "xlsx_raw")
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument(
        "--composition",
        type=Path,
        default=DEFAULT_XLSX_DIR / "BW25113_Medium composition.xlsx",
        help="BW25113_Medium composition.xlsx "
             "(default: $DATA_SRC, else ../chemical-media-dataset/xlsx_raw)",
    )
    parser.add_argument(
        "--reference",
        type=Path,
        default=DEFAULT_XLSX_DIR / "BW25113_GrowthDataEvaluation.xlsx",
        help="BW25113_GrowthDataEvaluation.xlsx "
             "(default: $DATA_SRC, else ../chemical-media-dataset/xlsx_raw)",
    )
    args = parser.parse_args()

    for path in (args.composition, args.reference):
        if not path.exists():
            sys.exit(f"Input workbook not found: {path}")

    print("Loading data ...")
    eval_df = load_evaluation(args.reference)
    comp_df, compound_cols = load_composition(args.composition)

    # The composition workbook is itself stored one row per curve, so its
    # ConditionID column repeats. Collapse it to one row per medium before the
    # join, otherwise the merge below is many-to-many and explodes.
    #
    # Two conditions (Cond00577, Cond00578) carry more than one distinct value
    # for a single compound under the same ConditionID in the upstream
    # workbook. `first()` resolves those the same way `prepare_ml_inputs.py`
    # does, keeping both pathways consistent.
    n_rows_before = len(comp_df)
    comp_df = comp_df.groupby("ConditionID", as_index=False)[compound_cols].first()
    if len(comp_df) != n_rows_before:
        print(f"  Collapsed composition to {len(comp_df)} unique media")

    print("\nJoining label -> ConditionID -> composition ...")
    merged = eval_df.merge(comp_df, on="ConditionID", how="inner")
    unmatched = len(eval_df) - len(merged)
    if unmatched:
        print(f"  {unmatched} curve(s) had no matching ConditionID in the "
              f"composition workbook and were dropped")

    # One row per curve, deterministic order so re-running reproduces the file
    # byte for byte.
    out = (merged[["label"] + compound_cols]
           .sort_values("label", kind="stable")
           .reset_index(drop=True))

    # `load_composition` collapses the workbook's embedded newlines to spaces
    # ("Glucose (mM)"). The published result tables name compounds with
    # underscores ("Glucose_(mM)"), and these names flow straight through
    # run_ml_via_guibiont.py into the correlation and importance CSVs, so
    # convert to the published form here.
    out = out.rename(columns={c: c.replace(" ", "_") for c in compound_cols})

    dupes = int(out["label"].duplicated().sum())
    if dupes:
        print(f"  WARNING: {dupes} duplicate curve label(s) in the evaluation "
              f"workbook -- downstream joins may double-count them")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUT_DIR / "feature_matrix.csv"
    out.to_csv(out_path, index=False)

    print(f"\nWritten: {out_path}")
    print(f"  {len(out)} curves x {len(compound_cols)} compounds "
          f"({out['label'].nunique()} unique labels, "
          f"{merged['ConditionID'].nunique()} media)")


if __name__ == "__main__":
    main()
