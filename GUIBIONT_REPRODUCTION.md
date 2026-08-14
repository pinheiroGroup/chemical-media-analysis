# Reproducing the BW25113 chemical media screen with GUIbiont

This guide walks through replicating the batch fitting, clustering, and ML
analysis from the paper. Data preparation and the two result-producing steps
(log-linear fitting and the ML downstream analysis) run as Python scripts that
drive the GUIbiont HTTP API; clustering and the optional manual ML path run in
the browser.

The scripted path is the one that produced the numbers reported in the paper,
and it is the path that reproduces end to end. Where a GUIbiont tab offers the
same analysis interactively, that is noted as an equivalent alternative.

**Assumed paths** (adjust if your setup differs):

| Path | Contents |
|---|---|
| `data/` | Per-round CSVs produced by preprocessing |
| `results/` | Output directory |
| `$GUIBIONT` | Root of the GUIbiont server directory |
| `$DATA_SRC` | Directory containing the raw BW25113 xlsx files and `BW25113_Medium composition.xlsx` |

---

## Prerequisites

- GUIbiont server running (`julia --project=. --threads=auto web_server.jl` from the GUIbiont repo).
  This serves on **port 8080**, which is what the scripts in this repo target by
  default. To use a different port, start the server with `PORT=<n>` and point the
  scripts at it with `GUIBIONT_API=http://localhost:<n>`.
- Python ≥ 3.9 with `pandas`, `numpy`, `openpyxl`, `scipy`
- Install dependencies once: `pip install pandas numpy openpyxl scipy`

Two environment variables are honoured throughout:

| Variable | Meaning | Default |
|---|---|---|
| `DATA_SRC` | Directory holding the raw BW25113 xlsx files | `../chemical-media-dataset/xlsx_raw` |
| `GUIBIONT_API` | Base URL of the running GUIbiont server | `http://localhost:8080` |

---

## Step 1 — Convert raw xlsx to GUIbiont experiment folders

The raw data are seven BioTek plate-reader Excel files (`BW25113_Growth_Round01.xlsx` … `Round07.xlsx`).
The script used in this step converts each Excel file to a self-contained GUIbiont experiment folder: a
`data_channel_1.csv` (Time + one column per curve) and a companion
`annotation_clean.csv` (one row per curve marked as a real growth well).

Run from the repo root (point at the directory holding the raw xlsx files):

```bash
python preprocess.py --xlsx-dir "$DATA_SRC" --min-amplitude 0
# or rely on the default xlsx-dir: ../chemical-media-dataset/xlsx_raw
```

The manuscript analysis keeps all curves at this stage — `--min-amplitude 0`
disables the flat-curve filter described below, deferring exclusion to the
per-curve fit status recorded in Steps 4 and 5. Omitting the flag applies the
default `0.05` threshold instead, which drops more curves upstream and will
not reproduce the paper's counts.

Output (per round `NN` = 01 … 07):

```
data/bw25113_roundNN/
├── data_channel_1.csv   # Time_h, Curve00001, Curve00002, …
└── annotation_clean.csv # CurveNNNNN, g, , , , ,    (no header)
```

Two filters run during conversion (both reported in the per-round summary):

- **Empty columns** dropped (rounds 04 and 05 had 66 and 142 empty wells respectively).
- **Flat curves** dropped — any curve whose `max(OD) - min(OD)` is below
  `--min-amplitude` (default `0.05`). These are blanks, dead cultures, or
  instrument-noise-only wells; fitting them produces meaningless parameters
  that pollute the per-medium means in Step 8. Lower the threshold
  (e.g., `0.02`) to keep marginal slow-growers, or set to `0` to disable.

> **Why both files?** GUIbiont's `/api/experiment/{name}/info` requires both
> `data_channel_1.csv` *and* `annotation_clean.csv` to exist — otherwise the
> Batch Fit tab returns 404 and no wells render.

---

## Step 2 — Register experiments in GUIbiont

GUIbiont's Batch Fit tab reads from `$GUIBIONT/Clean_data/<experiment>/`.
The folders produced in Step 1 are already in the right shape — just copy them in:

```bash
cp -r data/bw25113_round0[1-7] "$GUIBIONT/Clean_data/"
```

Restart the GUIbiont server (or refresh the experiment list) so the new
experiments appear in the dropdown.

---

## Step 3 — Batch fitting (GUIbiont interface)

Repeat the following for each of the seven rounds (`bw25113_round01` … `bw25113_round07`).

1. Open the **Batch Fit** tab
2. From the experiment dropdown, select `bw25113_round01`
3. Under **Model selection**, leave **Single model** selected and pick `aHPM` from the model dropdown  
   (to fit several models and pick the best by AICc, use **Compare models — pick best by AICc** and tick the desired checkboxes instead)
4. *Optional* — tick **Blank subtraction** if the round contains blank wells you want subtracted (method defaults to *Point-by-point*)
5. **Optimizer** — switch the mode dropdown from **Single** to **Best of N**:
   - Deterministic: tick `LN_COBYLA` (default)
   - Stochastic: tick `BBO_adaptive_de_rand_1_bin_radiuslimited` (default)
   - **Runs per stochastic:** `3`
   
   This runs COBYLA once + BBO three times per curve and keeps the fit with the lowest RMSE.
   COBYLA gives a reproducible baseline; BBO's three restarts escape local minima where COBYLA gets stuck.
   On Curve02586 we measured ~3× lower error vs single-COBYLA, with most of the win coming from BBO finding a better basin.

   Leave max iterations at `20000`, tolerance at `1e-6`, and **Skip flat ≤** at `0.05`
   (the upstream filter in `preprocess.py` should already have removed flats; this is a safety net).
6. Under **Wells**, click **Select all**
7. Click **⚡ Run Batch Fit**  
   Runtime: best-of-N is ~4× a single fit. On 2,500 curves with `--threads=auto` on a modern laptop, expect ~15–30 min per round. The progress bar tracks completed wells.
8. When the results table appears, check the summary bar:  
   `✓ N fitted   ⊘ M skipped (flat)`  
   Expect ≥99% of non-flat curves to converge. Click **📥 Download CSV** → save as `results/batch_fit_round01.csv`.

Repeat for rounds 02–07, saving as `results/batch_fit_round02.csv` … `results/batch_fit_round07.csv`.

> **Why Best of N over a single optimizer?** Local optimizers (COBYLA, BOBYQA) are reproducible but
> can get trapped in basins around bad seeds. BBO's stochastic global search escapes those, but is
> noisy — one run might be excellent, another terrible. Running several attempts and keeping the
> best fit by RMSE gets the strengths of both. The new `loss` and `optimizer_used` columns in the
> CSV record which attempt won per curve.

> **Note:** The Batch Fit tab does not expose smoothing or stationary-phase trimming controls — fitting runs on the raw curves as supplied. The fitter does report a derived `stationary_phase_start` column in the output.

---

## Step 4 — Combine fit results (external script)

Merge the seven per-round CSVs into a single file, add the round label, and filter to converged curves:

```bash
python scripts/combine_fits.py
```

This script concatenates `results/batch_fit_round0{1-7}.csv`, adds a `round` column, and writes:

- `results/batch_fit_results.csv` — one row per fitted curve with columns `round, experiment, well, model, gr, exit_lag_rate, N_max, shape, stationary_phase_start, aic, loss, optimizer_used`

The `loss` column is the RMSE of the chosen fit against the raw OD over the growth window — useful for spotting per-curve outliers.  
The `optimizer_used` column records which optimizer (and which restart of BBO) produced the winning fit for each curve; useful for sanity-checking that Best-of-N is actually doing work.

---

## Step 5 — Log-linear fitting (external script, GUIbiont API)

The aHPM fit from Step 3 supplies the parametric λ and K. The growth rates the
paper reports come from GUIbiont's **log-linear** estimator, which is run over
every round through `/api/batch-fit-loglin`:

```bash
python run_loglin_via_guibiont.py
```

The script submits each of the seven `bw25113_roundNN` experiments registered in
Step 2, polls until the job finishes, and records one row per curve — including
curves that were skipped or failed to converge, so the output accounts for every
retained record rather than silently shrinking to the curves that worked.

Requires the GUIbiont server from the Prerequisites to be running, and
`results/batch_fit_results.csv` from Step 4 for the merge.

Output:

| File | Contents |
|---|---|
| `results/batch_fit_results_loglin.csv` | Per-curve log-linear `gr_loglin`, `lag_loglin`, `N_max_emp`, `R_squared`, `fit_status` (13,400 rows) |
| `results/batch_fit_results_merged.csv` | Step 4's aHPM table left-joined with the above |

Expected on the full screen: **13,400 curves**, of which **13,349** converge
(`loglin_converged`). Rows are written sorted by `round, label`, so re-running
against unchanged inputs reproduces the files byte for byte.

> **Note:** `lag_loglin` and `N_max_emp` are model-free companions added to
> `fitting_one_well_Log_Lin`. If they come back all-NaN, an older
> Kinbiont/GUIbiont is running — the script prints a warning to that effect.

---

## Step 6 — Clustering (GUIbiont interface)

Clustering is run on the combined set of all 13,608 raw growth curves.
First, build a single interpolated CSV across all rounds:

```bash
python scripts/build_combined_curves.py
```

This interpolates each round to a common 97-point time grid (0–48 h) and outputs
`results/all_curves_combined.csv` (Time_h + 13,608 curve columns).
The script uses boundary-constant extrapolation so the matrix is **NaN-free**
— GUIbiont's `/api/cluster-sweep` does not sanitise NaN and would otherwise
silently fail with an empty result.

Then in GUIbiont:

### 6a — Cluster sweep (find optimal k)

1. Open the **Clustering** tab → leave the **From File** toggle selected (default)
2. Upload `results/all_curves_combined.csv` via the file picker (or paste its absolute path into *Or enter server-side file path*)
3. Click **⚙ Advanced** to expand the options panel and set:
   - **Smoothing** → keep the default `LOWESS` with bandwidth fraction `0.05`
     *(do not pick `Rolling average` — it triggers a GUIbiont server bug:
     the sweep endpoint forgets to pass `smooth_pt_avg` and returns HTTP 500
     "Sweep failed: undefined")*
   - **Cluster method** → `K-means` (default)
   - Under **Non-growing pre-screen**, tick **Enable** and set τ (tolerance) to `1.5`
4. In the **Find best k** row, set the max-k input to `10` and click **🔍 Sweep k**
5. Inspect the **WCSS elbow plot** — for this dataset the elbow lands at **k = 4**
   (WCSS drops 102756 → 50439 between k=3 and k=4, then flattens; Davies-Bouldin
   also halves at that step and Calinski-Harabasz peaks there)

### 6b — Run clustering at optimal k

1. Set **k =** to the elbow value
2. Click **▶ Run Clustering**
3. Inspect the cluster grid — one cluster typically captures slow-growing / non-growing conditions (the sentinel cluster from the pre-screen)
4. Click **📥 Export all (CSV)** → save as `results/clustering/cluster_assignments.csv`

---

## Step 7 — Build the per-curve feature matrix (external script)

Give every curve the 44 compound concentrations of the medium it was grown in:

```bash
python scripts/build_percurve_feature_matrix.py \
    --composition "$DATA_SRC/BW25113_Medium composition.xlsx" \
    --reference   "$DATA_SRC/BW25113_GrowthDataEvaluation.xlsx"
# or rely on the $DATA_SRC-relative defaults
```

The script maps curve `label` → `ConditionID` via
`BW25113_GrowthDataEvaluation.xlsx`, then `ConditionID` → 44 compounds via
`BW25113_Medium composition.xlsx`, and writes **one row per curve**. Replicate
curves of the same medium therefore repeat that medium's concentration vector;
Step 8 does its own aggregation by medium.

Output:

| File | Contents |
|---|---|
| `results/guibiont_ml_inputs/feature_matrix.csv` | `label` + 44 compound concentrations, one row per curve (13,608 rows × 44 compounds, 1,029 media) |

Notes:

- This step is driven purely by the two workbooks and is **independent of fit
  results**, so every curve in the evaluation mapping gets a row — including the
  208 curves dropped upstream by the Step 1 data-point filter. Step 8 inner-joins
  against the log-linear table, so the surplus rows fall out there.
- The composition workbook is itself stored one row per curve, so its
  `ConditionID` column repeats; the script collapses it to one row per medium
  before joining. Two conditions (`Cond00577`, `Cond00578`) carry more than one
  distinct value for a single compound under the same `ConditionID` upstream —
  the first is kept, matching `prepare_ml_inputs.py`'s behaviour.
- Rows are sorted by `label`, so the output is byte-reproducible.

---

## Step 8 — ML downstream analysis (external script, GUIbiont API)

This is the step that produced the machine-learning numbers in the paper. It
posts the per-medium matrices to `/api/ml-downstream`, so every value comes from
the same random forest the GUIbiont UI displays (Kinbiont's DecisionTree.jl
backend, Breiman defaults: sqrt-p feature subsampling, 0.7 partial sampling, 100
trees, max depth 5) — not from sklearn.

```bash
python run_ml_via_guibiont.py
```

Inputs: `results/batch_fit_results_loglin.csv` (Step 5),
`results/guibiont_ml_inputs/feature_matrix.csv` (Step 7), and
`$DATA_SRC/BW25113_GrowthDataEvaluation.xlsx` for the curve → medium mapping.

Replicate curves are averaged to **one row per medium** before the ML call.
This matters: each medium contributes ~13 replicate curves with an identical
feature vector, so per-curve cross-validation would leak replicates of a
test-fold medium into the training fold. Aggregating first makes the folds hold
out whole formulations, as the manuscript describes.

Expected console output:

```
Aggregated 13400 replicate curves -> 1026 matched media (features: 1026 rows)
Got 1026 rows joined; saving canonical CSVs.
```

Reference values to check against:

| Quantity | Expected |
|---|---|
| Media in the final matrices | 1,026 |
| CV R² — `gr` | ≈ 0.70 |
| CV R² — `N_max_emp` | ≈ 0.87 |
| Top Spearman ρ for `gr` | `Glucose_(mM)`, ρ ≈ 0.79 |
| Top Spearman ρ for `N_max_emp` | `GlutamicAcid/HCl_(mM)`, ρ ≈ 0.83 |

Outputs land in `results/ml_results_modelfree/` — see the summary table below.
The complete, unmodified endpoint response is preserved in
`ml_downstream_response.json`, which carries the individual CV fold values
alongside the displayed mean and SD.

### Alternative — the ML Analysis tab (equivalent, manual)

The same analysis can be driven interactively. This path is aggregated **per
condition** rather than per curve, so it uses a different input pair, produced by
the legacy preparation script:

```bash
python scripts/prepare_ml_inputs.py \
    --fits results/batch_fit_results.csv \
    --composition "$DATA_SRC/BW25113_Medium composition.xlsx" \
    --reference   "$DATA_SRC/BW25113_GrowthDataEvaluation.xlsx"
```

That writes `results/guibiont_ml_inputs/fit_results_by_condition.csv` and
`feature_matrix_by_condition.csv` (one row per condition, filtered to converged
aHPM fits). Then:

1. Open the **🔬 ML Analysis** tab
2. **Step 1 — Batch-fit results CSV**: upload `fit_results_by_condition.csv`
3. When the label-column dropdown appears, select `condition`
4. **Step 2 — Feature matrix CSV**: upload `feature_matrix_by_condition.csv`
5. **Step 3 — Params for RF importance**: tick `gr`, `exit_lag_rate`, `N_max`
6. Click **🔬 Run ML Analysis**

The interface shows Spearman rank correlations, random-forest feature importance
(top 15), and partial dependence plots for the top 5 features per parameter.

> **Do not mix the two pathways.** The scripted path needs the **per-curve**
> `feature_matrix.csv` from Step 7; the manual path needs the **per-condition**
> `*_by_condition.csv` files. They are written under distinct names precisely so
> one cannot shadow the other. Feeding a per-condition matrix to
> `run_ml_via_guibiont.py` produces an empty join.

This manual path reports aHPM-based parameters (`gr`, `exit_lag_rate`, `N_max`),
not the log-linear targets, so its numbers are not directly comparable to the
table above.

---

## Summary of files produced

| File | How | Step |
|---|---|---|
| `data/bw25113_round0{1-7}/` | `preprocess.py` | 1 |
| `results/batch_fit_round0{1-7}.csv` | GUIbiont Batch Fit tab (×7) | 3 |
| `results/batch_fit_results.csv` | `scripts/combine_fits.py` | 4 |
| `results/batch_fit_results_loglin.csv` | `run_loglin_via_guibiont.py` | 5 |
| `results/batch_fit_results_merged.csv` | `run_loglin_via_guibiont.py` | 5 |
| `results/all_curves_combined.csv` | `scripts/build_combined_curves.py` | 6 |
| `results/clustering/cluster_assignments.csv` | GUIbiont Clustering tab | 6 |
| `results/guibiont_ml_inputs/feature_matrix.csv` | `scripts/build_percurve_feature_matrix.py` | 7 |
| `results/guibiont_ml_inputs/parameter_matrix_by_medium.csv` | `run_ml_via_guibiont.py` | 8 |
| `results/guibiont_ml_inputs/feature_matrix_by_medium.csv` | `run_ml_via_guibiont.py` | 8 |
| `results/ml_results_modelfree/correlations_loglin.csv` | `run_ml_via_guibiont.py` | 8 |
| `results/ml_results_modelfree/feature_importance_*.csv` | `run_ml_via_guibiont.py` | 8 |
| `results/ml_results_modelfree/perm_importance_*.csv` | `run_ml_via_guibiont.py` | 8 |
| `results/ml_results_modelfree/cv_r2_summary.csv` | `run_ml_via_guibiont.py` | 8 |
| `results/ml_results_modelfree/cv_r2_folds.csv` | `run_ml_via_guibiont.py` | 8 |
| `results/ml_results_modelfree/ml_downstream_response.json` | `run_ml_via_guibiont.py` | 8 |

Produced only by the optional manual ML path (Step 8, alternative):

| File | How |
|---|---|
| `results/guibiont_ml_inputs/fit_results_by_condition.csv` | `scripts/prepare_ml_inputs.py` |
| `results/guibiont_ml_inputs/feature_matrix_by_condition.csv` | `scripts/prepare_ml_inputs.py` |
