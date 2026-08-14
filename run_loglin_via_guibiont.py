#!/usr/bin/env python3
"""Run /api/batch-fit-loglin on all seven chemical-media experiments.

The API sees the 13,400 curves retained during conversion.  The conversion
manifest is joined back afterwards so the published output accounts for all
13,608 original records, including the 208 excluded before fitting.

Outputs:
    results/batch_fit_results_loglin.csv   # one row per original record

Run:
    /usr/bin/python run_loglin_via_guibiont.py
"""
from __future__ import annotations

import json
import os
import time
import urllib.request
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
MANIFEST_CSV = HERE / "data" / "conversion_manifest.csv"
OUT_LOGLIN_CSV = RESULTS / "batch_fit_results_loglin.csv"
EVAL_XLSX = Path(os.environ.get(
    "DATA_SRC", HERE.parent / "chemical-media-dataset" / "xlsx_raw")
) / "BW25113_GrowthDataEvaluation.xlsx"

API = os.environ.get("GUIBIONT_API", "http://localhost:8080")

# Log-lin params -- same as the Keio run for cross-study consistency. The
# 97-timepoint chemical-media curves are even denser than the 200-point
# Keio mean curves, so the same window settings remain appropriate.
LOGLIN_PARAMS = {
    "blank_subtraction":       False,
    "type_of_smoothing":       "rolling_avg",
    "pt_avg":                  7,
    "pt_smoothing_derivative": 7,
    "pt_min_size_of_win":      7,
    "type_of_win":             "maximum",
    "threshold_of_exp":        0.9,
    "skip_flat_threshold":     0.0,
}


def _post(path, body):
    req = urllib.request.Request(
        f"{API}{path}", method="POST",
        headers={"Content-Type": "application/json"},
        data=json.dumps(body).encode())
    with urllib.request.urlopen(req) as r:
        return r.status, json.loads(r.read())


def _get(path):
    req = urllib.request.Request(f"{API}{path}")
    with urllib.request.urlopen(req) as r:
        return r.status, json.loads(r.read())


def submit_and_wait(experiment: str, poll_s: float = 4.0):
    payload = {"experiment": experiment, **LOGLIN_PARAMS}
    s, body = _post("/api/batch-fit-loglin", payload)
    if s != 200:
        raise RuntimeError(f"{experiment}: HTTP {s} {body}")
    job_id = body["job_id"]
    total  = body["total"]
    print(f"  {experiment}: job_id={job_id}, total={total}")
    start = time.time()
    while True:
        s, p = _get(f"/api/batch-fit/progress/{job_id}")
        if p.get("status") == "done":
            elapsed = time.time() - start
            print(f"  {experiment}: done in {elapsed:.1f}s -- "
                  f"summary={p['summary']}")
            return p
        if p.get("status") in {"error", "failed"}:
            raise RuntimeError(f"{experiment}: job failed: {p}")
        time.sleep(poll_s)


def collect_round(job, round_num: int):
    """One row per curve the round contained, successful or not.

    Curves the estimator skipped (flat-signal threshold) or failed on are
    written with empty descriptors and a fit_status saying why, so the output
    accounts for every retained record instead of silently shrinking to the
    curves that worked.
    """
    rows = []
    for r in job["results"]:
        # GUIbiont returns a normal (non-exception) result even when the
        # fit didn't converge to a meaningful positive mu_max -- e.g. a
        # flat/data-starved curve whose slope rejected below the epsilon+R^2
        # floor (analysis.jl:fit_well_loglin). Route those into the same
        # "no_positive_mu_max" bucket as the ones Kinbiont itself throws on,
        # instead of mislabeling them "fitted".
        converged = bool(r.get("loglin_converged", False))
        rows.append({
            "round":           round_num,
            "label":           r["well"],
            "fit_status":      "fitted" if converged else "no_positive_mu_max",
            "gr_loglin":       r.get("gr_loglin"),
            "gr_loglin_se":    r.get("gr_loglin_se"),
            "gr_max_sliding":  r.get("gr_max_sliding"),
            "t_exp_start":     r.get("t_exp_start_loglin"),
            "t_exp_end":       r.get("t_exp_end_loglin"),
            "doubling_time":   r.get("doubling_time_loglin"),
            "R_squared":       r.get("R_squared_loglin"),
            # Model-free companions returned by the same estimator. NaN when
            # the upgraded Kinbiont/GUIbiont is not running.
            "lag_loglin":      r.get("lag_loglin"),
            "N_max_emp":       r.get("N_max_emp"),
            "loglin_converged": bool(r.get("loglin_converged", False)),
        })

    empty = {k: "" for k in ("gr_loglin", "gr_loglin_se", "gr_max_sliding",
                             "t_exp_start", "t_exp_end", "doubling_time",
                             "R_squared", "lag_loglin", "N_max_emp")}
    summary = job.get("summary") or {}
    for well in summary.get("skipped", []) or job.get("skipped", []) or []:
        label = well if isinstance(well, str) else well.get("well", "")
        rows.append({"round": round_num, "label": label,
                     "fit_status": "skipped_flat", **empty,
                     "loglin_converged": False})
    for err in summary.get("errors", []) or job.get("errors", []) or []:
        text = err if isinstance(err, str) else json.dumps(err)
        label = text.split("'")[1] if "'" in text else text
        # Two distinct outcomes the paper reports separately: curves the
        # conversion step drops for having too few numeric OD readings, and
        # curves that were fitted but yielded no finite positive mu_max.
        if "insufficient data points" in text:
            status = "excluded_insufficient_data"
        elif "must be finite and positive" in text:
            status = "no_positive_mu_max"
        else:
            status = "failed"
        rows.append({"round": round_num, "label": label,
                     "fit_status": status, **empty,
                     "loglin_converged": False})
    return rows


def main():
    if not MANIFEST_CSV.exists():
        raise FileNotFoundError(
            f"Missing {MANIFEST_CSV}. Run preprocess.py before fitting.")
    if not EVAL_XLSX.exists():
        raise FileNotFoundError(
            f"Missing {EVAL_XLSX}. Set DATA_SRC to the raw-data directory.")

    rows = []
    print("Submitting batch-fit-loglin for each chemical-media round...")
    t0 = time.time()
    for n in range(1, 8):
        exp = f"bw25113_round{n:02d}"
        job = submit_and_wait(exp)
        rows.extend(collect_round(job, n))
    print(f"\nTotal wall time: {time.time()-t0:.1f}s")

    api_rows = pd.DataFrame(rows)
    manifest = pd.read_csv(MANIFEST_CSV)
    expected = manifest[manifest["conversion_status"] == "retained"]

    key = ["round", "label"]
    expected_keys = set(map(tuple, expected[key].itertuples(index=False,
                                                            name=None)))
    actual_keys = set(map(tuple, api_rows[key].itertuples(index=False,
                                                          name=None)))
    if expected_keys != actual_keys:
        missing = sorted(expected_keys - actual_keys)[:5]
        extra = sorted(actual_keys - expected_keys)[:5]
        raise RuntimeError(
            "API results do not match the conversion manifest: "
            f"missing={missing}, extra={extra}")

    # Attach the medium identifier promised in Supplementary Data S3.
    evaluation = pd.read_excel(EVAL_XLSX)
    evaluation.columns = evaluation.columns.str.strip()
    label_col = next(c for c in evaluation.columns if "label" in c.lower())
    condition_col = next(c for c in evaluation.columns
                         if "condition" in c.lower())
    evaluation = evaluation.rename(
        columns={label_col: "label", condition_col: "ConditionID"})
    evaluation["label"] = evaluation["label"].astype(str).str.strip()
    evaluation["ConditionID"] = (
        evaluation["ConditionID"].astype(str).str.strip())
    evaluation = evaluation[["label", "ConditionID"]].drop_duplicates()

    # Left-joining onto the complete manifest restores the 208 records that
    # never reached the API. Their descriptor fields remain empty and their
    # conversion status becomes the published fit status.
    df_ll = manifest.merge(api_rows, on=key, how="left", validate="one_to_one")
    excluded = df_ll["conversion_status"] != "retained"
    df_ll.loc[excluded, "fit_status"] = df_ll.loc[
        excluded, "conversion_status"]
    df_ll.loc[excluded, "loglin_converged"] = False
    df_ll = df_ll.merge(evaluation, on="label", how="left",
                        validate="many_to_one")
    if df_ll["ConditionID"].isna().any():
        labels = df_ll.loc[df_ll["ConditionID"].isna(), "label"].head().tolist()
        raise RuntimeError(f"Missing ConditionID for labels: {labels}")

    df_ll = (df_ll.drop(columns="conversion_status")
             .sort_values(key, kind="stable")
             .reset_index(drop=True))
    columns = ["round", "label", "ConditionID", "fit_status",
               "gr_loglin", "gr_loglin_se", "gr_max_sliding",
               "t_exp_start", "t_exp_end", "doubling_time", "R_squared",
               "lag_loglin", "N_max_emp", "loglin_converged"]
    df_ll = df_ll[columns]
    df_ll.to_csv(OUT_LOGLIN_CSV, index=False, lineterminator="\n")
    retained_count = int((~df_ll["fit_status"].str.startswith(
        "excluded_")).sum())
    converged_count = int(df_ll["loglin_converged"].sum())
    print(f"\nWrote {OUT_LOGLIN_CSV}  ({len(df_ll)} original records, "
          f"{retained_count} retained, {converged_count} usable fits)")

    print(f"\nmu_loglin quantiles: "
          f"p05={df_ll['gr_loglin'].quantile(0.05):.3f}  "
          f"p50={df_ll['gr_loglin'].median():.3f}  "
          f"p95={df_ll['gr_loglin'].quantile(0.95):.3f}  "
          f"p99={df_ll['gr_loglin'].quantile(0.99):.3f}  "
          f"max={df_ll['gr_loglin'].max():.3f}")

    # Surface the model-free companions so it's obvious whether the
    # upgraded Kinbiont was running (NaN otherwise).
    for col, label in [("lag_loglin", "lag_loglin (h)"),
                       ("N_max_emp", "N_max_emp (OD)")]:
        if col in df_ll.columns and df_ll[col].notna().any():
            print(f"{label} quantiles: "
                  f"p05={df_ll[col].quantile(0.05):.3f}  "
                  f"p50={df_ll[col].median():.3f}  "
                  f"p95={df_ll[col].quantile(0.95):.3f}")
        else:
            print(f"{label}: all NaN -- old GUIbiont/Kinbiont still running?")

if __name__ == "__main__":
    main()
