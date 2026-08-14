#!/usr/bin/env python3
"""Run GUIbiont's ML-downstream (impurity + permutation + CV R^2 + Spearman)
on the chemical-media screen through /api/ml-downstream and persist the
canonical CSVs.

This is the API-driven replacement for the earlier sklearn-via-Python
script (rerun_ml_with_modelfree.py): every number now comes from the same
RF that the GUIbiont UI displays -- Breiman defaults (sqrt-p feature
subsampling, 0.7 partial sampling, 100 trees, max depth 5) via Kinbiont's
DecisionTree.jl backend, not sklearn.

Inputs:
    results/batch_fit_results_loglin.csv  (model-free targets per curve)
    results/guibiont_ml_inputs/feature_matrix.csv  (44 compound concs per curve)

Outputs (results/ml_results_modelfree/):
    correlations_loglin.csv          # Spearman rho per (compound, target)
    feature_importance_{gr,lag_loglin,N_max_emp}.csv          # impurity
    perm_importance_{gr,lag_loglin,N_max_emp}.csv             # permutation
    cv_r2_summary.csv                # per-target CV R^2 mean +/- std
    cv_r2_folds.csv                  # five fold-level R^2 values per target
    ml_downstream_response.json      # complete, unmodified endpoint response

Run:
    /usr/bin/python run_ml_via_guibiont.py
"""
from __future__ import annotations

import json
import os
import sys
import urllib.request
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
RES  = HERE / "results"
OUT  = RES / "ml_results_modelfree"
OUT.mkdir(exist_ok=True)

LOGLIN_CSV = RES / "batch_fit_results_loglin.csv"
FEAT_CSV   = RES / "guibiont_ml_inputs" / "feature_matrix.csv"
# Curve label -> medium (Condition ID) mapping. Each medium has ~13 biological
# replicate curves; aggregating by medium before ML makes cross-validation hold
# out whole media (formulations) rather than individual replicate curves.
EVAL_XLSX  = Path(os.environ.get(
    "DATA_SRC", HERE.parent / "chemical-media-dataset" / "xlsx_raw")
) / "BW25113_GrowthDataEvaluation.xlsx"

API = os.environ.get("GUIBIONT_API", "http://localhost:8080")
TARGETS = ("gr", "N_max_emp")


def main() -> None:
    # --- Aggregate replicate curves by medium (Condition ID) before the ML call.
    # Random per-curve cross-validation leaks: every medium contributes ~13
    # replicate curves with an identical feature vector, so replicates of a
    # test-fold medium would also appear in the training fold. Averaging to one
    # row per medium removes this by construction and matches the manuscript's
    # "replicate runs were averaged beforehand" / "held-out formulations".
    ll = pd.read_csv(LOGLIN_CSV)
    ll["label"] = ll["label"].astype(str).str.strip()
    feat_pc = pd.read_csv(FEAT_CSV)
    feat_pc["label"] = feat_pc["label"].astype(str).str.strip()

    ev = pd.read_excel(EVAL_XLSX)
    lc = next(c for c in ev.columns if c.strip().lower() == "label")
    cc = next(c for c in ev.columns if "condition" in c.lower())
    ev = ev.rename(columns={lc: "label", cc: "ConditionID"})
    ev["label"] = ev["label"].astype(str).str.strip()
    ev["ConditionID"] = ev["ConditionID"].astype(str).str.strip()
    ev = ev[["label", "ConditionID"]]

    # Targets: mean over the converged replicate curves of each medium.
    llm = ll.merge(ev, on="label", how="inner")
    fit_df = (llm.groupby("ConditionID")[["gr_loglin", "lag_loglin", "N_max_emp"]]
                 .mean().reset_index()
                 .rename(columns={"ConditionID": "label", "gr_loglin": "gr"}))
    # Features: one row per medium (identical across its replicate curves).
    compound_cols = [c for c in feat_pc.columns if c != "label"]
    fmap = feat_pc.merge(ev, on="label", how="inner")
    feat = (fmap.groupby("ConditionID")[compound_cols].first().reset_index()
                .rename(columns={"ConditionID": "label"}))

    # Submit and publish the exact same set of media in both matrices. Three
    # conditions have no retained/fitted curves and therefore no usable target
    # values; the endpoint's inner join excluded them implicitly in older runs,
    # leaving a 1,029-row parameter file beside a 1,026-row feature file.
    common_labels = fit_df[["label"]].merge(feat[["label"]], on="label", how="inner")
    fit_df = common_labels.merge(fit_df, on="label", how="left")
    feat = common_labels.merge(feat, on="label", how="left")
    print(f"Aggregated {len(ll)} replicate curves -> {len(fit_df)} matched media "
          f"(features: {len(feat)} rows)")

    # Save the matrices actually submitted, not just the per-curve inputs.
    # Replicate averaging happens here, so these are the rows the random
    # forest and the cross-validation folds see; Supplementary Data S4 needs
    # exactly them to be reproducible.
    agg_dir = RES / "guibiont_ml_inputs"
    agg_dir.mkdir(parents=True, exist_ok=True)
    fit_df.to_csv(agg_dir / "parameter_matrix_by_medium.csv", index=False)
    feat.to_csv(agg_dir / "feature_matrix_by_medium.csv", index=False)

    payload = dict(
        fit_csv        = fit_df.to_csv(index=False),
        label_col      = "label",
        feature_matrix = feat.to_csv(index=False),
        params         = list(TARGETS),
    )
    req = urllib.request.Request(
        f"{API}/api/ml-downstream", method="POST",
        headers={"Content-Type": "application/json"},
        data=json.dumps(payload).encode())
    with urllib.request.urlopen(req, timeout=900) as r:
        body = json.loads(r.read())

    print(f"Got {body['n_wells']} rows joined; saving canonical CSVs.")

    # Preserve the complete endpoint response. In particular, cv_r2 contains
    # the individual fold values in addition to the displayed mean and SD.
    # Keeping this JSON prevents tabular exporters from silently discarding
    # fields returned by GUIbiont.
    response_path = OUT / "ml_downstream_response.json"
    response_path.write_text(json.dumps(body, indent=2, sort_keys=True) + "\n")
    print(f"  wrote {response_path}")

    # Spearman rho matrix: one row per compound, one column per target.
    corr_rows = []
    for row in body["correlations"]:
        out = {"compound": row["feature"]}
        for t in TARGETS:
            out[t] = row.get(t)
        corr_rows.append(out)
    pd.DataFrame(corr_rows).to_csv(OUT / "correlations_loglin.csv", index=False)
    print(f"  wrote {OUT / 'correlations_loglin.csv'}")

    # Importance + permutation per target.
    cv_rows = []
    cv_fold_rows = []
    for t in TARGETS:
        imp  = pd.DataFrame(body["importance"][t])
        imp  = imp.rename(columns={"feature": "compound"})
        imp.to_csv(OUT / f"feature_importance_{t}.csv", index=False)
        print(f"  wrote {OUT / f'feature_importance_{t}.csv'}  "
              f"top: {imp.iloc[0]['compound']} ({imp.iloc[0]['importance']:.3f})")

        perm = pd.DataFrame(body["permutation_importance"][t])
        perm = perm.rename(columns={"feature": "compound"})
        perm.to_csv(OUT / f"perm_importance_{t}.csv", index=False)
        print(f"  wrote {OUT / f'perm_importance_{t}.csv'}  "
              f"top: {perm.iloc[0]['compound']} "
              f"({perm.iloc[0]['permutation_importance']:.3f})")

        cv = body["cv_r2"][t]
        cv_rows.append(dict(
            dataset="chemical-media", medium="-", target=t,
            cv_r2_mean=cv["mean"], cv_r2_std=cv["std"], n=cv["n"],
        ))
        for fold_index, fold_r2 in enumerate(cv.get("folds", []), start=1):
            cv_fold_rows.append(dict(
                dataset="chemical-media", medium="-", target=t,
                fold=fold_index, cv_r2=fold_r2, n=cv["n"],
            ))

    pd.DataFrame(cv_rows).to_csv(OUT / "cv_r2_summary.csv", index=False)
    print(f"  wrote {OUT / 'cv_r2_summary.csv'}")
    pd.DataFrame(cv_fold_rows).to_csv(OUT / "cv_r2_folds.csv", index=False)
    print(f"  wrote {OUT / 'cv_r2_folds.csv'}")


if __name__ == "__main__":
    main()
