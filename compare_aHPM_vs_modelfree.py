#!/usr/bin/env python3
"""Side-by-side comparison of RF feature importances and Spearman ρ between
the old run (log-lin μ + aHPM λ/K — saved in ml_results_aHPM_lag_K/) and
the new run (fully model-free targets — written by rerun_ml_with_modelfree.py
into ml_results_modelfree/).

Outputs:
    results/aHPM_vs_modelfree_comparison.csv   long-format table with both
                                               importances per (param, compound)
    results/aHPM_vs_modelfree_summary.md       short markdown report with
                                               top-N movers per parameter
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

HERE    = Path(__file__).resolve().parent
RESULTS = HERE / "results"
OLD_DIR = RESULTS / "ml_results_aHPM_lag_K"
NEW_DIR = RESULTS / "ml_results_modelfree"

OUT_CSV = RESULTS / "aHPM_vs_modelfree_comparison.csv"
OUT_MD  = RESULTS / "aHPM_vs_modelfree_summary.md"

# (label, old filename, old corr column, new filename, new corr column)
PARAMS = [
    ("gr",
     "feature_importance_gr.csv",        "gr",
     "feature_importance_gr.csv",        "gr"),
    ("lag",
     "feature_importance_exit_lag_rate.csv", "exit_lag_rate",
     "feature_importance_lag_loglin.csv",    "lag_loglin"),
    ("N_max",
     "feature_importance_N_max.csv",     "N_max",
     "feature_importance_N_max_emp.csv", "N_max_emp"),
]


def _md_table(df: pd.DataFrame) -> str:
    """Tabulate-free GitHub-flavored markdown table writer."""
    if df.empty:
        return "_(empty)_"
    cols = list(df.columns)
    def fmt(v):
        if pd.isna(v):
            return ""
        if isinstance(v, float):
            return f"{v:.4g}"
        return str(v)
    head = "| " + " | ".join(cols) + " |"
    sep  = "|" + "|".join(["---"] * len(cols)) + "|"
    body = "\n".join(
        "| " + " | ".join(fmt(v) for v in row) + " |"
        for row in df.itertuples(index=False, name=None)
    )
    return "\n".join([head, sep, body])


def _load_imp(path: Path, source: str) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=["compound", f"imp_{source}",
                                     f"rank_{source}"])
    df = pd.read_csv(path)
    df = df.rename(columns={"importance": f"imp_{source}"})
    df[f"rank_{source}"] = df[f"imp_{source}"].rank(ascending=False,
                                                    method="min").astype(int)
    return df[["compound", f"imp_{source}", f"rank_{source}"]]


def _load_corr(path: Path, col: str) -> pd.Series:
    if not path.exists() or col not in pd.read_csv(path, nrows=0).columns:
        return pd.Series(dtype=float, name=col)
    return pd.read_csv(path).set_index("compound")[col]


def main():
    old_corr = _load_corr(OLD_DIR / "correlations_loglin.csv", None) \
        if False else None  # placeholder; correlations handled per-param below
    long_rows = []

    md_lines = ["# aHPM vs model-free comparison", ""]
    md_lines.append("Old run: `ml_results_aHPM_lag_K/` "
                    "(log-lin μ paired with aHPM exit_lag_rate / N_max).")
    md_lines.append("New run: `ml_results_modelfree/` "
                    "(log-lin μ, Buchanan lag, q95 N_max — all Kinbiont).")
    md_lines.append("")

    for label, old_file, old_col, new_file, new_col in PARAMS:
        old_imp = _load_imp(OLD_DIR / old_file, "old")
        new_imp = _load_imp(NEW_DIR / new_file, "new")
        merged = old_imp.merge(new_imp, on="compound", how="outer")
        merged.insert(0, "parameter", label)
        # Bring in Spearman ρ from each side's correlations CSV.
        try:
            old_c = pd.read_csv(OLD_DIR / "correlations_loglin.csv")
            new_c = pd.read_csv(NEW_DIR / "correlations_loglin.csv")
            if old_col in old_c.columns:
                merged = merged.merge(
                    old_c[["compound", old_col]].rename(columns={old_col: "rho_old"}),
                    on="compound", how="left")
            else:
                merged["rho_old"] = float("nan")
            if new_col in new_c.columns:
                merged = merged.merge(
                    new_c[["compound", new_col]].rename(columns={new_col: "rho_new"}),
                    on="compound", how="left")
            else:
                merged["rho_new"] = float("nan")
        except FileNotFoundError:
            merged["rho_old"] = float("nan")
            merged["rho_new"] = float("nan")
        long_rows.append(merged)

        # Markdown summary: top movers by rank change.
        if "rank_old" in merged.columns and "rank_new" in merged.columns:
            m = merged.dropna(subset=["rank_old", "rank_new"]).copy()
            m["delta_rank"] = (m["rank_new"] - m["rank_old"]).astype(int)
            m["delta_imp"]  = m["imp_new"] - m["imp_old"]
            md_lines.append(f"## Parameter: `{label}`")
            md_lines.append("")
            md_lines.append("Top 5 old:")
            md_lines.append("")
            md_lines.append(m.sort_values("rank_old").head(5)[
                ["compound", "imp_old", "imp_new",
                 "rank_old", "rank_new", "delta_rank"]
            ].pipe(_md_table))
            md_lines.append("")
            md_lines.append("Top 5 new:")
            md_lines.append("")
            md_lines.append(m.sort_values("rank_new").head(5)[
                ["compound", "imp_old", "imp_new",
                 "rank_old", "rank_new", "delta_rank"]
            ].pipe(_md_table))
            md_lines.append("")
            md_lines.append("Biggest rank shifts (|Δrank| ≥ 5):")
            md_lines.append("")
            big = m[m["delta_rank"].abs() >= 5].sort_values("delta_rank")
            if len(big) == 0:
                md_lines.append("_None._")
            else:
                md_lines.append(big[
                    ["compound", "imp_old", "imp_new",
                     "rank_old", "rank_new", "delta_rank"]
                ].pipe(_md_table))
            md_lines.append("")

    long_df = pd.concat(long_rows, ignore_index=True)
    long_df.to_csv(OUT_CSV, index=False)
    OUT_MD.write_text("\n".join(md_lines))
    print(f"Wrote {OUT_CSV} ({len(long_df)} rows)")
    print(f"Wrote {OUT_MD}")


if __name__ == "__main__":
    main()
