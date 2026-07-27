#!/usr/bin/env python3
"""Run /api/batch-fit-loglin on all chemical-media Round experiments and merge
the resulting μ_loglin into the existing aHPM batch_fit_results.csv (which
keeps the parametric λ and K).

Outputs:
    results/batch_fit_results_loglin.csv   # μ_loglin only, one row per curve
    results/batch_fit_results_merged.csv   # full merge with existing aHPM table

Run:
    /usr/bin/python run_loglin_via_guibiont.py
"""
from __future__ import annotations

import json
import os
import sys
import time
import urllib.request
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
EXISTING_FIT_CSV = RESULTS / "batch_fit_results.csv"
OUT_LOGLIN_CSV   = RESULTS / "batch_fit_results_loglin.csv"
OUT_MERGED_CSV   = RESULTS / "batch_fit_results_merged.csv"

API = os.environ.get("GUIBIONT_API", "http://localhost:9090")

# Log-lin params — same as the Keio run for cross-study consistency. The
# 97-timepoint chemical-media curves are even denser than the 200-point
# Keio mean curves, so the same window settings remain appropriate.
LOGLIN_PARAMS = {
    "pt_avg":                  7,
    "pt_smoothing_derivative": 7,
    "pt_min_size_of_win":      7,
    "threshold_of_exp":        0.9,
    "skip_flat_threshold":     0.0,
}


def _post(path, body):
    req = urllib.request.Request(
        f"{API}{path}", method="POST",
        headers={"Content-Type": "application/json"},
        data=json.dumps(body).encode())
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.status, json.loads(r.read())


def _get(path):
    req = urllib.request.Request(f"{API}{path}")
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.status, json.loads(r.read())


def submit_and_wait(experiment: str, poll_s: float = 4.0,
                    max_s: float = 1200.0):
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
            print(f"  {experiment}: done in {elapsed:.1f}s — "
                  f"summary={p['summary']}")
            return p
        if time.time() - start > max_s:
            raise TimeoutError(f"{experiment} did not finish in {max_s}s")
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
        rows.append({
            "round":           round_num,
            "label":           r["well"],
            "fit_status":      "fitted",
            "gr_loglin":       r.get("gr_loglin"),
            "gr_loglin_se":    r.get("gr_loglin_se"),
            "gr_max_sliding":  r.get("gr_max_sliding"),
            "t_exp_start":     r.get("t_exp_start_loglin"),
            "t_exp_end":       r.get("t_exp_end_loglin"),
            "doubling_time":   r.get("doubling_time_loglin"),
            "R_squared":       r.get("R_squared_loglin"),
            # Model-free companions added to fitting_one_well_Log_Lin so all
            # three RF targets share a single estimator definition. NaN when
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
    rows = []
    print("Submitting batch-fit-loglin for each chemical-media round…")
    t0 = time.time()
    for n in range(1, 8):
        exp = f"chem_round{n:02d}"
        job = submit_and_wait(exp)
        rows.extend(collect_round(job, n))
    print(f"\nTotal wall time: {time.time()-t0:.1f}s")

    df_ll = pd.DataFrame(rows)
    df_ll.to_csv(OUT_LOGLIN_CSV, index=False)
    print(f"\nWrote {OUT_LOGLIN_CSV}  ({len(df_ll)} rows, "
          f"{int(df_ll['loglin_converged'].sum())} converged)")

    print(f"\nμ_loglin quantiles: "
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
            print(f"{label}: all NaN — old GUIbiont/Kinbiont still running?")

    # ── Merge with existing aHPM batch fit ─────────────────────────────────
    if not EXISTING_FIT_CSV.exists():
        print(f"WARNING: no existing fit table at {EXISTING_FIT_CSV} — "
              "skipping merge")
        return
    df_ahpm = pd.read_csv(EXISTING_FIT_CSV)
    merged = df_ahpm.merge(df_ll, on=["round", "label"], how="left")
    merged.to_csv(OUT_MERGED_CSV, index=False)
    print(f"\nWrote {OUT_MERGED_CSV}  ({len(merged)} rows; "
          f"{merged['gr_loglin'].notna().sum()} have μ_loglin)")

    # Quick comparison aHPM gr vs log-lin gr (on matched rows)
    matched = merged.dropna(subset=["gr", "gr_loglin"])
    print(f"\nOn {len(matched)} matched curves:")
    print(f"  aHPM gr median = {matched['gr'].median():.3f},  "
          f"max = {matched['gr'].max():.3f}")
    print(f"  Log-lin gr_loglin median = {matched['gr_loglin'].median():.3f},  "
          f"max = {matched['gr_loglin'].max():.3f}")
    print(f"  Spearman ρ(aHPM, log-lin) = "
          f"{matched[['gr', 'gr_loglin']].corr('spearman').iloc[0,1]:.3f}")
    n_blown = int((matched["gr"] > 3).sum())
    print(f"  aHPM curves above physiological ceiling (gr > 3): "
          f"{n_blown}/{len(matched)}")


if __name__ == "__main__":
    main()
