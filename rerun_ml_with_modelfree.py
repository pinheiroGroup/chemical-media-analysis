#!/usr/bin/env python3
"""Re-run GUIbiont's ML-downstream on the chemical-media screen using the
fully model-free target set: log-linear μ_max, Buchanan tangent-intercept
lag (lag_loglin), and empirical N_max (q95 of smoothed OD). All three
targets are produced by Kinbiont's log-linear estimator via the GUIbiont
/api/batch-fit-loglin endpoint — no parametric fit is consulted.

Inputs:
    results/guibiont_ml_inputs/feature_matrix.csv
    results/batch_fit_results_loglin.csv
        (after running run_loglin_via_guibiont.py against the upgraded
         GUIbiont server, this file carries lag_loglin and N_max_emp)

Outputs (kept separate from ml_results/ so the previous aHPM-based run
is preserved):
    results/ml_results_modelfree/correlations_loglin.csv
    results/ml_results_modelfree/feature_importance_gr.csv
    results/ml_results_modelfree/feature_importance_lag_loglin.csv
    results/ml_results_modelfree/feature_importance_N_max_emp.csv

Same Spearman+RF hyperparameters as GUIbiont-repo/src/ml_downstream.jl
(100 trees, max_depth=5, sample_fraction=0.7, seed=42).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor

HERE    = Path(__file__).resolve().parent
RESULTS = HERE / "results"
ML_OUT  = RESULTS / "ml_results_modelfree"
ML_OUT.mkdir(exist_ok=True)

FEAT_CSV   = RESULTS / "guibiont_ml_inputs" / "feature_matrix.csv"
LOGLIN_CSV = RESULTS / "batch_fit_results_loglin.csv"

RF_KW = dict(n_estimators=100, max_depth=5, max_samples=0.7,
             random_state=42, n_jobs=-1)


def spearman_table(X, y, feat_names, param_name):
    rows = []
    for j, n in enumerate(feat_names):
        xj = X[:, j]
        mask = np.isfinite(xj) & np.isfinite(y)
        rho, p = (np.nan, np.nan)
        if mask.sum() >= 3:
            rho, p = spearmanr(xj[mask], y[mask])
        rows.append({"compound": n, param_name: float(rho),
                     f"{param_name}_p": float(p)})
    return pd.DataFrame(rows)


def rf_importance_table(X, y, feat_names):
    mask = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    if mask.sum() < 10:
        return pd.DataFrame(columns=["compound", "importance"])
    rf = RandomForestRegressor(**RF_KW).fit(X[mask], y[mask])
    imp = rf.feature_importances_
    df = pd.DataFrame({"compound": feat_names, "importance": imp})
    return df.sort_values("importance", ascending=False).reset_index(drop=True)


def main():
    feat = pd.read_csv(FEAT_CSV)
    ll   = pd.read_csv(LOGLIN_CSV)

    required = {"label", "gr_loglin", "lag_loglin", "N_max_emp"}
    missing = required - set(ll.columns)
    if missing:
        sys.exit(
            f"ERROR: {LOGLIN_CSV.name} is missing {sorted(missing)}.\n"
            "Re-run run_loglin_via_guibiont.py against the upgraded "
            "GUIbiont server first.")

    fit = ll[["label", "gr_loglin", "lag_loglin", "N_max_emp"]].rename(
        columns={"gr_loglin": "gr"})
    print(f"Fit table: {len(fit)} curves (model-free μ + lag + N_max)")

    joined = feat.merge(fit, on="label", how="inner")
    print(f"Joined with feature matrix: {len(joined)} curves")

    feature_names = [c for c in feat.columns if c != "label"]
    X = joined[feature_names].to_numpy(float)

    results_corr = []
    for pname, col in [("gr",         "gr"),
                       ("lag_loglin", "lag_loglin"),
                       ("N_max_emp",  "N_max_emp")]:
        y = joined[col].to_numpy(float)
        sp = spearman_table(X, y, feature_names, pname)
        results_corr.append(sp.set_index("compound")[[pname]])
        imp = rf_importance_table(X, y, feature_names)
        out = ML_OUT / f"feature_importance_{pname}.csv"
        imp.to_csv(out, index=False)
        print(f"  {pname:>11}: top-5 importance → " + ", ".join(
            f"{r['compound']} ({r['importance']:.3f})"
            for _, r in imp.head(5).iterrows()))

    corr = pd.concat(results_corr, axis=1).reset_index()
    corr.to_csv(ML_OUT / "correlations_loglin.csv", index=False)
    print(f"\nWrote {ML_OUT / 'correlations_loglin.csv'} "
          f"({corr.shape[0]} compounds × 3 parameters)")


if __name__ == "__main__":
    main()
