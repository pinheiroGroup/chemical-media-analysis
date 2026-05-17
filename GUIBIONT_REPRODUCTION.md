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

## Step 1 — Convert raw xlsx to CSV

The raw data are seven BioTek plate-reader Excel files (`BW25113_Growth_Round01.xlsx` … `Round07.xlsx`).
This script converts each to a GUIbiont-compatible CSV with a `Time_h` column and one column per curve.

Run from the repo root:

```bash
python preprocess.py
```

Output: `data/Round01.csv` … `data/Round07.csv`  
Format: `Time_h, Curve00001, Curve00002, …` — one column per growth curve, rows = time points.  
Total: 13,608 curves across 7 rounds (276 / 324 / 110 / 309 / 471 / 231 / 146 per round).

---

## Step 2 — Register experiments in GUIbiont

GUIbiont's Batch Fit tab reads from `$GUIBIONT/Clean_data/<experiment>/`.
Create one experiment directory per round and copy the preprocessed CSV into each:

```bash
for i in 01 02 03 04 05 06 07; do
    mkdir -p "$GUIBIONT/Clean_data/bw25113_round${i}"
    cp "data/Round${i}.csv" "$GUIBIONT/Clean_data/bw25113_round${i}/Round${i}.csv"
done
```

---

## Step 3 — Batch fitting (GUIbiont interface)

Repeat the following for each of the seven rounds (`bw25113_round01` … `bw25113_round07`).

1. Open the **Batch Fit** tab
2. Select experiment `bw25113_round01`
3. Open **Advanced options**:
   - Smooth method: `rolling average`, window: `14`
   - Check **Cut stationary phase**
     - Percentile threshold: `0.05`
     - Smooth derivative window: `10`
     - Window size: `5`
4. **Models**: check `aHPM` only (deselect logistic, gompertz, baranyi_richards)
5. Click **Run**  
   (runtime depends on hardware and Julia thread count — expect ≥ 95% convergence per round)
6. Click **Download CSV** → save as `results/batch_fit_round01.csv`

Repeat for rounds 02–07, saving as `results/batch_fit_round02.csv` … `results/batch_fit_round07.csv`.

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

Clustering is run on the combined set of all 13,608 raw growth curves.
First, build a single interpolated CSV across all rounds:

```bash
python scripts/build_combined_curves.py
```

This interpolates each round to a common 97-point time grid (0–48 h) and outputs
`results/all_curves_combined.csv` (Time_h + 13,608 curve columns).

Then in GUIbiont:

### 5a — Cluster sweep (find optimal k)

1. Open the **Clustering** tab → click **File** mode
2. Upload `results/all_curves_combined.csv`
3. Open **Advanced options**:
   - Smooth method: `rolling average`, window: `14`
   - Check **Pre-screen constant curves**, tolerance: `1.5`
   - Cluster method: `kmeans`
4. Click **Sweep** (k = 2 to 10)
5. Look at the **WCSS elbow plot** — note the suggested k (expect 4–6)

### 5b — Run clustering at optimal k

1. Set **k** to the elbow value
2. Click **Run**
3. Inspect the cluster grid — one cluster typically captures slow-growing / non-growing conditions
4. Click **Export all (CSV)** → save as `results/clustering/cluster_assignments.csv`

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

1. Open the **ML Analysis** tab
2. **Fit results CSV**: upload `results/guibiont_ml_inputs/fit_results.csv`
3. **Label column**: `condition`
4. **Feature matrix CSV**: upload `results/guibiont_ml_inputs/feature_matrix.csv`
5. **Parameters to analyse**: select `gr`, `exit_lag_rate`, `N_max`
6. Click **Run**

The interface shows:
- **Spearman correlations** — ranked bar chart of compound–parameter associations
- **Random forest feature importance** — top compounds per growth parameter
- **Partial dependence plots** — marginal effect of top 5 compounds

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
