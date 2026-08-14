"""
Build the per-curve feature matrix consumed by `run_ml_via_guibiont.py`.

Every growth curve retained during conversion gets one row carrying the 44
compound concentrations of the medium it was grown in:

  BW25113_GrowthDataEvaluation.xlsx  ->  curve label -> ConditionID
  BW25113_Medium composition.xlsx    ->  ConditionID -> 44 compound columns

The two tables are joined on ConditionID and the result is written per curve
label -- replicate curves of the same medium therefore repeat that medium's
concentration vector. `run_ml_via_guibiont.py` does its own aggregation by
medium afterwards, so the matrix handed to it must be per curve, not per
condition.

The conversion manifest supplies the retained label set. This keeps the ML
input aligned with the 13,400 curves that reached fitting, while the complete
13,608-record fit-status table remains available separately.

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

try:
    import pandas as pd
except ImportError:
    sys.exit("pandas is required: pip install pandas openpyxl")

HERE     = Path(__file__).resolve().parent
REPO     = HERE.parent
OUT_DIR  = REPO / "results" / "guibiont_ml_inputs"
MANIFEST = REPO / "data" / "conversion_manifest.csv"
DEFAULT_XLSX_DIR = Path(
    os.environ.get("DATA_SRC", REPO.parent / "chemical-media-dataset" / "xlsx_raw")
)


def load_evaluation(path: Path) -> pd.DataFrame:
    df = pd.read_excel(path)
    df.columns = df.columns.str.strip()
    label_col = next((c for c in df.columns if "label" in c.lower()), None)
    condition_col = next(
        (c for c in df.columns if "condition" in c.lower()), None)
    if label_col is None or condition_col is None:
        sys.exit(f"Could not find Label/Condition columns in {path.name}")
    df = df.rename(columns={label_col: "label",
                            condition_col: "ConditionID"})
    df["label"] = df["label"].astype(str).str.strip()
    df["ConditionID"] = df["ConditionID"].astype(str).str.strip()
    return df[["label", "ConditionID"]]


def load_composition(path: Path) -> tuple[pd.DataFrame, list[str]]:
    df = pd.read_excel(path)
    df.columns = df.columns.str.replace(r"\s+", " ", regex=True).str.strip()
    condition_col = next(
        (c for c in df.columns if "condition" in c.lower()), df.columns[0])
    df = df.rename(columns={condition_col: "ConditionID"})
    df["ConditionID"] = df["ConditionID"].astype(str).str.strip()
    metadata = {"ConditionID", "Label", "Assay ID"}
    compound_cols = [c for c in df.columns if c not in metadata]
    for col in compound_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df, compound_cols


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

    for path in (args.composition, args.reference, MANIFEST):
        if not path.exists():
            sys.exit(f"Required input not found: {path}")

    print("Loading data ...")
    eval_df = load_evaluation(args.reference)
    comp_df, compound_cols = load_composition(args.composition)
    manifest = pd.read_csv(MANIFEST)
    retained = manifest.loc[
        manifest["conversion_status"] == "retained", ["label"]]
    retained["label"] = retained["label"].astype(str).str.strip()
    eval_df = retained.merge(eval_df, on="label", how="left",
                             validate="one_to_one")
    if eval_df["ConditionID"].isna().any():
        labels = eval_df.loc[eval_df["ConditionID"].isna(), "label"].head()
        sys.exit(f"Missing ConditionID for labels: {labels.tolist()}")

    # The composition workbook is itself stored one row per curve, so its
    # ConditionID column repeats. Collapse it to one row per medium before the
    # join, otherwise the merge below is many-to-many and explodes.
    #
    # Two conditions (Cond00577, Cond00578) carry more than one distinct value
    # for a single compound under the same ConditionID in the upstream
    # workbook. The published analysis consistently keeps the first value.
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

    # One row per retained curve, in deterministic order so re-running
    # reproduces the file byte for byte.
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
    out.to_csv(out_path, index=False, lineterminator="\n")

    print(f"\nWritten: {out_path}")
    print(f"  {len(out)} curves x {len(compound_cols)} compounds "
          f"({out['label'].nunique()} unique labels, "
          f"{merged['ConditionID'].nunique()} media)")


if __name__ == "__main__":
    main()
