# Reproducing the BW25113 chemical media screen with GUIbiont

This guide walks through replicating the batch fitting, clustering, and ML
analysis from the paper using the GUIbiont web interface.
One Python script handles data preparation; one handles result aggregation;
everything else runs in the browser.

**Assumed paths** (adjust if your setup differs):

| Path | Contents |
|---|---|
| `data/` | Per-round CSVs produced by preprocessing |
| `results/` | Output directory |
| `$GUIBIONT` | Root of the GUIbiont server directory |
| `$DATA_SRC` | Directory containing the raw BW25113 xlsx files and `BW25113_Medium composition.xlsx` |

---

## Prerequisites

- GUIbiont server running (`julia --project=. --threads=auto web_server.jl` from the GUIbiont repo)
- Python ≥ 3.9 with `pandas`, `numpy`, `openpyxl`, `scipy`
- Install dependencies once: `pip install pandas numpy openpyxl scipy`

---

## Step 1 — Convert raw xlsx to GUIbiont experiment folders

The raw data are seven BioTek plate-reader Excel files (`BW25113_Growth_Round01.xlsx` … `Round07.xlsx`).
This script converts each to a self-contained GUIbiont experiment folder: a
`data_channel_1.csv` (Time + one column per curve) and a companion
`annotation_clean.csv` (one row per curve marked as a real growth well).

Run from the repo root (point at the directory holding the raw xlsx files):

```bash
python preprocess.py --xlsx-dir "$DATA_SRC"
# or rely on the default: ../chemical-media-dataset/xlsx_raw
```

Output (per round `NN` = 01 … 07):

```
data/bw25113_roundNN/
├── data_channel_1.csv   # Time_h, Curve00001, Curve00002, …
└── annotation_clean.csv # CurveNNNNN, g, , , , ,    (no header)
```

Total: 13,400 curves across 7 rounds (2628 / 2640 / 960 / 2574 / 2498 / 1320 / 780 per round).
Empty columns in the source xlsx are dropped automatically (rounds 04 and 05 had 66 and 142
empty wells respectively) — they would otherwise be reported as fit failures by GUIbiont.

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
5. Leave the optimizer at **BOBYQA (default)**, max iterations `20000`, tolerance `1e-6`
6. Under **Wells**, click **Select all**
7. Click **⚡ Run Batch Fit**  
   (runtime depends on hardware and Julia thread count — expect ≥ 95% convergence per round)
8. When the results table appears, click **📥 Download CSV** in the results bar → save as `results/batch_fit_round01.csv`

Repeat for rounds 02–07, saving as `results/batch_fit_round02.csv` … `results/batch_fit_round07.csv`.

> **Note:** The Batch Fit tab does not expose smoothing or stationary-phase trimming controls — fitting runs on the raw curves as supplied. The fitter does report a derived `stationary_phase_start` column in the output.

---

## Step 4 — Combine fit results (external script)

Merge the seven per-round CSVs into a single file, add the round label, and filter to converged curves:

```bash
python scripts/combine_fits.py
```

This script concatenates `results/batch_fit_round0{1-7}.csv`, adds a `round` column, and writes:

- `results/batch_fit_results.csv` — all 13,608 rows with columns `round, label, gr, exit_lag_rate, N_max, shape, aicc, loss, n_timepoints, converged`

---

## Step 5 — Clustering (GUIbiont interface)

Clustering is run on the combined set of all 13,400 raw growth curves.
First, build a single interpolated CSV across all rounds:

```bash
python scripts/build_combined_curves.py
```

This interpolates each round to a common 97-point time grid (0–48 h) and outputs
`results/all_curves_combined.csv` (Time_h + 13,400 curve columns).
The script uses boundary-constant extrapolation so the matrix is **NaN-free**
— GUIbiont's `/api/cluster-sweep` does not sanitise NaN and would otherwise
silently fail with an empty result.

Then in GUIbiont:

### 5a — Cluster sweep (find optimal k)

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

### 5b — Run clustering at optimal k

1. Set **k =** to the elbow value
2. Click **▶ Run Clustering**
3. Inspect the cluster grid — one cluster typically captures slow-growing / non-growing conditions (the sentinel cluster from the pre-screen)
4. Click **📥 Export all (CSV)** → save as `results/clustering/cluster_assignments.csv`

---

## Step 6 — Prepare ML feature matrix (external script)

Link each fitted curve to its medium's chemical composition (44 compounds).
The script reads `$DATA_SRC/BW25113_Medium composition.xlsx` and
`$DATA_SRC/BW25113_GrowthDataEvaluation.xlsx` (for reference K/r values),
aggregates fit parameters per condition (mean across replicates and rounds),
then joins with compound concentrations:

```bash
python scripts/prepare_ml_inputs.py \
    --fits results/batch_fit_results.csv \
    --composition "$DATA_SRC/BW25113_Medium composition.xlsx" \
    --reference  "$DATA_SRC/BW25113_GrowthDataEvaluation.xlsx"
```

Output:

| File | Contents |
|---|---|
| `results/guibiont_ml_inputs/fit_results.csv` | Per-condition mean `gr`, `exit_lag_rate`, `N_max` (1,026 rows) |
| `results/guibiont_ml_inputs/feature_matrix.csv` | 44 compound concentrations per condition (1,026 rows) |

---

## Step 7 — ML Analysis (GUIbiont interface)

1. Open the **🔬 ML Analysis** tab
2. **Step 1 — Batch-fit results CSV**: click **📂 Choose file** and upload `results/guibiont_ml_inputs/fit_results.csv`
3. When the label-column dropdown appears, select `condition`
4. **Step 2 — Feature matrix CSV**: click **📂 Choose file** and upload `results/guibiont_ml_inputs/feature_matrix.csv`
5. **Step 3 — Params for RF importance**: tick `gr`, `exit_lag_rate`, `N_max`
6. Click **🔬 Run ML Analysis**

The interface shows:
- **Spearman rank correlations** — bar chart per growth parameter (switch via the *Growth parameter* dropdown)
- **Random-forest feature importance (top 15)** — top compounds per growth parameter
- **Partial dependence plots (top 5 features per parameter)** — marginal effect of top compounds

Key findings to look for:
- `exit_lag_rate` is more strongly correlated with reference growth rate (ρ ≈ 0.77) than `gr` alone
- Cystine is the top predictor of `exit_lag_rate` (ρ ≈ 0.79)
- Glucose and iron are top predictors of `gr`
- Citrate shows an iron-dependent effect on lag exit

---

## Summary of files produced

| File | How |
|---|---|
| `data/Round0{1-7}.csv` | `preprocess.py` |
| `results/batch_fit_round0{1-7}.csv` | GUIbiont Batch Fit tab (×7) |
| `results/batch_fit_results.csv` | `scripts/combine_fits.py` |
| `results/all_curves_combined.csv` | `scripts/build_combined_curves.py` |
| `results/clustering/cluster_assignments.csv` | GUIbiont Clustering tab |
| `results/guibiont_ml_inputs/fit_results.csv` | `scripts/prepare_ml_inputs.py` |
| `results/guibiont_ml_inputs/feature_matrix.csv` | `scripts/prepare_ml_inputs.py` |
| ML figures and correlation tables | GUIbiont ML Analysis tab |
