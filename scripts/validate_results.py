"""Validate every chemical-media count and headline result used in the paper."""

from pathlib import Path
import sys

import pandas as pd


REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data"
RESULTS = REPO / "results"
ML = RESULTS / "ml_results_modelfree"
INPUTS = RESULTS / "guibiont_ml_inputs"

ROUND_COUNTS = {
    1: (2628, 2628, 2602),
    2: (2640, 2640, 2635),
    3: (960, 960, 951),
    4: (2640, 2574, 2567),
    5: (2640, 2498, 2494),
    6: (1320, 1320, 1320),
    7: (780, 780, 780),
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> None:
    manifest = pd.read_csv(DATA / "conversion_manifest.csv")
    fits = pd.read_csv(RESULTS / "batch_fit_results_loglin.csv")
    features = pd.read_csv(INPUTS / "feature_matrix.csv")
    parameters_by_medium = pd.read_csv(INPUTS / "parameter_matrix_by_medium.csv")
    features_by_medium = pd.read_csv(INPUTS / "feature_matrix_by_medium.csv")
    correlations = pd.read_csv(ML / "correlations_loglin.csv")
    cv_summary = pd.read_csv(ML / "cv_r2_summary.csv")
    cv_folds = pd.read_csv(ML / "cv_r2_folds.csv")

    require(len(manifest) == 13608, "manifest must contain 13,608 records")
    require(len(fits) == 13608, "fit output must contain 13,608 records")
    require(not manifest[["round", "label"]].duplicated().any(),
            "manifest contains duplicate curve keys")
    require(not fits[["round", "label"]].duplicated().any(),
            "fit output contains duplicate curve keys")
    require(set(map(tuple, manifest[["round", "label"]].to_numpy())) ==
            set(map(tuple, fits[["round", "label"]].to_numpy())),
            "manifest and fit output have different curve keys")

    status_counts = fits["fit_status"].value_counts().to_dict()
    require(status_counts == {
        "fitted": 13349,
        "excluded_insufficient_data": 208,
        "no_positive_mu_max": 51,
    }, f"unexpected fit-status counts: {status_counts}")
    require(fits["ConditionID"].notna().all(),
            "every original record must have a medium identifier")

    for round_number, (original, retained, usable) in ROUND_COUNTS.items():
        rows = fits[fits["round"] == round_number]
        actual = (len(rows),
                  int((~rows["fit_status"].str.startswith("excluded_")).sum()),
                  int((rows["fit_status"] == "fitted").sum()))
        require(actual == (original, retained, usable),
                f"round {round_number}: expected {(original, retained, usable)}, "
                f"got {actual}")

    retained_labels = set(manifest.loc[
        manifest["conversion_status"] == "retained", "label"])
    require(len(features) == 13400, "feature matrix must contain 13,400 curves")
    require(features["label"].is_unique, "feature labels must be unique")
    require(set(features["label"]) == retained_labels,
            "feature matrix must contain exactly the retained curves")
    require(len(features.columns) == 45,
            "feature matrix must contain label plus 44 compounds")

    require(len(parameters_by_medium) == len(features_by_medium) == 1026,
            "both submitted matrices must contain 1,026 media")
    require(list(parameters_by_medium.columns) == ["label", "gr", "N_max_emp"],
            "parameter matrix must contain only the two manuscript targets")
    require(set(parameters_by_medium["label"]) ==
            set(features_by_medium["label"]),
            "submitted parameter and compound matrices have different media")
    require(len(features_by_medium.columns) == 45,
            "submitted compound matrix must contain 44 compounds")

    require(len(correlations) == 44, "correlation output must contain 44 compounds")
    expected_top = {"gr": ("Glucose_(mM)", 0.79),
                    "N_max_emp": ("GlutamicAcid/HCl_(mM)", 0.83)}
    for target, (compound, rounded_value) in expected_top.items():
        top = correlations.loc[correlations[target].idxmax()]
        require(top["compound"] == compound,
                f"unexpected top Spearman compound for {target}: {top['compound']}")
        require(round(float(top[target]), 2) == rounded_value,
                f"unexpected Spearman value for {target}: {top[target]}")

    expected_cv = {"gr": 0.70, "N_max_emp": 0.87}
    require(set(cv_summary["target"]) == set(expected_cv),
            "CV summary has unexpected targets")
    for target, rounded_value in expected_cv.items():
        row = cv_summary[cv_summary["target"] == target].iloc[0]
        require(round(float(row["cv_r2_mean"]), 2) == rounded_value,
                f"unexpected CV R-squared for {target}: {row['cv_r2_mean']}")
        require(int((cv_folds["target"] == target).sum()) == 5,
                f"{target} must have five CV folds")

    print("All chemical-media manuscript checks passed.")
    print("13,608 original; 13,400 retained; 13,349 usable; 1,026 media; 44 compounds.")


if __name__ == "__main__":
    try:
        main()
    except (AssertionError, FileNotFoundError, KeyError) as exc:
        print(f"VALIDATION FAILED: {exc}", file=sys.stderr)
        raise SystemExit(1)
