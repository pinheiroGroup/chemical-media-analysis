# Reproducing the BW25113 chemical-media analysis with GUIbiont

This repository contains the pipeline and outputs used for the chemical-media
results in the GUIbiont manuscript: conversion, log-linear fitting and
downstream analysis of medium composition.

## Prerequisites

- GUIbiont v1.1.1 with Kinbiont.jl v1.5.1, running on Julia 1.12.6.
- Python 3.9 or newer with `pandas`, `numpy` and `openpyxl`.
- The raw BW25113 workbooks from the published dataset.

The commands below assume that GUIbiont is in `$GUIBIONT`. `DATA_SRC` may
point to the directory containing the raw workbooks; without it, the scripts
use `../chemical-media-dataset/xlsx_raw`. The API URL defaults to
`http://localhost:8080` and can be overridden with `GUIBIONT_API`.

## 1. Convert the seven workbooks

```bash
python preprocess.py
```

The paper retains every curve having at least ten numeric OD measurements.
The default amplitude threshold is therefore zero. A non-zero
`--min-amplitude` is available for other analyses, but does not reproduce the
manuscript.

This writes seven GUIbiont experiment folders under `data/` and
`data/conversion_manifest.csv`, which records every original curve and its
conversion status. Expected totals are:

| Records | Count |
|---|---:|
| Original | 13,608 |
| Retained for analysis | 13,400 |
| Excluded for fewer than ten numeric measurements | 208 |

The exclusions are 66 records from round 4 and 142 from round 5. The generated
experiment folders committed here are byte-reproducible from the workbooks.

## 2. Register the experiments and start GUIbiont

Copy `data/bw25113_round01` through `data/bw25113_round07` into
`$GUIBIONT/Clean_data/`, then start the server from the GUIbiont repository:

```bash
julia --project=. --threads=auto web_server.jl
```

## 3. Run log-linear fitting

```bash
python run_loglin_via_guibiont.py
```

The script submits all seven retained experiment folders to
`/api/batch-fit-loglin` with the manuscript settings: rolling-average
smoothing, `pt_avg = 7`, `pt_smoothing_derivative = 7`,
`pt_min_size_of_win = 7`, `threshold_of_exp = 0.9`, no blank subtraction and
`skip_flat_threshold = 0`.

The API fits 13,400 curves. The script then joins those results back to the
complete conversion manifest and curve-to-medium mapping, producing
`results/batch_fit_results_loglin.csv` with one row for each of the 13,608
original records:

| Fit status | Count |
|---|---:|
| `fitted` | 13,349 |
| `no_positive_mu_max` | 51 |
| `excluded_insufficient_data` | 208 |

Rows are sorted by round and curve label. Each row includes its medium
identifier; excluded records have empty fit descriptors.

## 4. Build the retained per-curve compound matrix

```bash
python scripts/build_percurve_feature_matrix.py
```

The script maps each retained curve to its medium and the concentrations of 44
compounds. It writes
`results/guibiont_ml_inputs/feature_matrix.csv` with 13,400 rows from 1,026
media. The 208 excluded records, including all records from three media, are
not downstream-analysis inputs.

## 5. Run the downstream analysis

```bash
python run_ml_via_guibiont.py
```

Replicates are averaged by medium before the matrices are submitted to
`/api/ml-downstream`, ensuring that cross-validation holds out complete medium
formulations. The endpoint receives matching parameter and compound matrices
with 1,026 rows. The targets are maximum growth rate (`gr`) and empirical
saturation OD (`N_max_emp`).

Expected manuscript-level checks are:

| Quantity | Expected |
|---|---:|
| Fivefold CV R-squared, `gr` | 0.70 |
| Fivefold CV R-squared, `N_max_emp` | 0.87 |
| Highest Spearman association with `gr` | Glucose, 0.79 |
| Highest Spearman association with `N_max_emp` | Glutamic acid/HCl, 0.83 |

The complete endpoint response is retained alongside deterministic tabular
exports. Run `python scripts/validate_results.py` to check every count, label
set, target, fold count and headline value used by the manuscript.

## Canonical outputs

- `results/batch_fit_results_loglin.csv`
- `results/guibiont_ml_inputs/feature_matrix.csv`
- `results/guibiont_ml_inputs/parameter_matrix_by_medium.csv`
- `results/guibiont_ml_inputs/feature_matrix_by_medium.csv`
- `results/ml_results_modelfree/correlations_loglin.csv`
- `results/ml_results_modelfree/feature_importance_{gr,N_max_emp}.csv`
- `results/ml_results_modelfree/perm_importance_{gr,N_max_emp}.csv`
- `results/ml_results_modelfree/cv_r2_summary.csv`
- `results/ml_results_modelfree/cv_r2_folds.csv`
- `results/ml_results_modelfree/ml_downstream_response.json`

These are the only derived files used to build Supplementary Data S3 and S4
and to verify the numerical statements in `Guibiont.tex` and `SM.tex`.
