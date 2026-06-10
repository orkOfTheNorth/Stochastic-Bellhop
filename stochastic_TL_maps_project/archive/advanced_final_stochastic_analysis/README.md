# Stochastic Bellhop — Underwater Acoustics UQ Pipeline

Uncertainty quantification (UQ) of acoustic transmission loss (TL) in a shallow-water waveguide using the Bellhop Gaussian-beam ray tracer.  
Three uncertain physical inputs are propagated through Bellhop and analysed with four complementary methods:

| Input parameter | Nominal | Training range |
|---|---|---|
| Frequency | 10 000 Hz | ±1 % |
| Source depth z_S | 5 m | ±1 % |
| Temperature shift | 0 °C | ±1 °C |

Environment: depth D = 35 m, range up to 50 km, FOM threshold = 100 dB.

---

## תוכן עניינים — Table of Contents

1. [Quick start](#1-quick-start)
2. [Project structure](#2-project-structure)
3. [Monte Carlo](#3-monte-carlo)
4. [Delta Method](#4-delta-method)
5. [Decision Tree / Random Forest](#5-decision-tree--random-forest)
   - 5a. [Training & maps](#5a-training--maps-run_dtm)
   - 5b. [Explainability — PDPs & influence maps](#5b-explainability--pdps--influence-maps-run_dt_explainm)
   - 5c. [Feature importance (OOB permutation)](#5c-feature-importance-oob-permutation-run_dt_vardecompm)
   - 5d. [Out-of-distribution test](#5d-out-of-distribution-test-run_dt_oodm)
6. [Comparison](#6-comparison)
7. [Waveguide validation](#7-waveguide-validation)
8. [Utilities](#8-utilities)
9. [Run order](#9-run-order)
10. [Output files summary](#10-output-files-summary)

---

## 1. Quick start

```matlab
cd bellhop_09122025/2026-05-26
run_all          % runs everything in the correct order
```

Or run each step individually — see [Run order](#9-run-order).

---

## 2. Project structure

```
2026-05-26/
├── run_all.m                    — master script (runs all methods)
├── bellhopCached.m              — Bellhop wrapper with disk caching
├── viewer.m                     — interactive TL field viewer
│
├── Monte_Carlo/
│   ├── run_MC.m                 — Latin hypercube MC sampling
│   ├── show_MC.m                — interactive MC viewer (GUI)
│   ├── show_LN3_fits.m          — LN3 fit quality diagnostics
│   ├── TL_cache/                — cached Bellhop TL fields (single precision)
│   └── results/                 — MC statistics (.mat files)
│
├── Delta_Method/
│   ├── run_delta.m              — finite-difference Jacobians + Chebyshev bound
│   └── results/
│
├── Decision_Tree/
│   ├── run_DT.m                 — train Random Forest on MC TL cache
│   ├── run_DT_explain.m         — PDPs + spatial influence maps
│   ├── run_DT_vardecomp.m       — OOB permutation feature importance
│   ├── run_DT_OOD.m             — out-of-distribution generalization test
│   └── results/
│
├── Comparison/
│   ├── run_comparison.m         — 5-method shadow probability comparison
│   └── results/
│
└── Waveguide/
    ├── run_waveguide_comparison.m  — analytical normal modes vs Bellhop
    ├── waveguide_analytical.m      — normal mode TL computation
    └── results/
```

---

## 3. Monte Carlo

**Script:** `Monte_Carlo/run_MC.m`

Runs N = 100 Latin hypercube samples over (freq, z_S, temp) subsets.  
Seven parameter subsets are computed; the full 3-parameter subset is `TL_zS_freq_temp`.

| Output | Description |
|---|---|
| `TL_cache/TL_zS_freq_temp.mat` | Raw TL fields, single precision, Nz × Nr × N |
| `results/MC_zS_freq_temp.mat` | E[TL], Var[TL], P(TL>FOM) empirical, Chebyshev bound |

**Key methods computed per grid point:**
- `MC_EX` — sample mean E[TL]
- `MC_Var` — sample variance Var[TL]
- `MC_PrFOM` — empirical P(TL > FOM) = fraction of samples above threshold
- `Cheb_lb` — Chebyshev lower bound on P(TL > FOM): `1 − Var / (Var + (E−FOM)²)`
- LN3 fit — 3-parameter lognormal P(TL > FOM) via `ln3prob()`

**Diagnostics:** `show_LN3_fits.m` generates a 5×5 tiled histogram grid with LN3 PDF overlaid, 95 % model CI band, and a spatial KS-statistic map (fit quality).

---

## 4. Delta Method

**Script:** `Delta_Method/run_delta.m`

Finite-difference (Δ = 1 %) Jacobians of TL with respect to each uncertain input.  
The Chebyshev shadow-probability bound is computed from the linearised variance.

| Output | Description |
|---|---|
| `results/delta_jacobians.mat` | J_freq, J_zS, J_temp — dTL/dparam maps (dB/unit) |
| `results/delta_zS_freq_temp.mat` | Delta Chebyshev bound Cheb_lb_s |
| `figures/` | Jacobian maps, variance maps, shadow-probability map |

**Physical interpretation of Jacobians:**  
Each map shows how strongly TL changes at each (r, z) point for a unit perturbation of that parameter.  High-magnitude zones are where that parameter is most acoustically sensitive.

---

## 5. Decision Tree / Random Forest

### 5a. Training & maps — `run_DT.m`

Trains a 50-tree Random Forest (TreeBagger) to classify `TL > FOM` from 5 inputs:

```
X = [freq (Hz) | z_S (m) | temp shift (°C) | range (km) | depth (m)]
y = 1 if TL > FOM  (shadow zone),  0 otherwise
```

Training uses 80 % of the 100 LHS samples (split by sample index, not grid point, to prevent data leakage).

| Output | Description |
|---|---|
| `results/RF_results.mat` | RF object, metrics, RF_nominal_map, RF_mc_avg, grid axes |
| `figures/RF_metrics.png` | Confusion matrix, ROC curve, OOB feature importance |
| `figures/RF_maps.png` | 4-panel shadow-probability maps: RF nominal, RF MC-avg, MC empirical, Delta Chebyshev |
| `figures/RF_agreement.png` | RF vs MC density scatter + spatial difference map |

**Performance (corrected MC data):** Accuracy 99.4 %, AUC 0.9997, RMSE vs MC = 0.0066.

---

### 5b. Explainability — PDPs & influence maps — `run_DT_explain.m`

Two techniques to explain where and how each parameter drives the RF predictions:

**Partial Dependence Plots (PDPs)**  
Sweep each feature across a wide physical range (well beyond training bounds) while holding others fixed at nominal.  Average RF prediction shows the global marginal effect.  Grey band marks the training range.

**Spatial marginal influence maps**  
`influence(r,z) = P(shadow | param=HIGH) − P(shadow | param=LOW)`  
Red = high parameter value → more shadow there.  Blue = high value → less shadow.

| Output | Description |
|---|---|
| `figures/RF_PDP.png` | 5-panel PDP (one per feature) |
| `figures/RF_influence_maps.png` | 3-panel spatial influence maps (freq, z_S, temp) |
| `results/RF_explain.mat` | pdp_vals, infl_maps |

---

### 5c. Feature importance (OOB permutation) — `run_DT_vardecomp.m`

**Method — OOB Permutation Importance:**  
For each of the 5 input features, the feature's values are randomly shuffled across all out-of-bag (OOB) samples and the drop in OOB classification accuracy is measured.  Larger drop = the RF relies more on that feature.

Results (corrected MC, N = 100):

| Rank | Feature | OOB Importance | % of total |
|---|---|---|---|
| 1 | Range (km) | 571.6 | 77.7 % |
| 2 | Depth (m) | 122.2 | 16.6 % |
| 3 | Frequency (Hz) | 15.4 | 2.1 % |
| 4 | Source depth z_S | 14.3 | 1.9 % |
| 5 | Temp shift (°C) | 11.9 | 1.6 % |

> Range and Depth dominate because acoustic shadow zones have a fixed spatial structure — knowing *where* you are is the strongest single predictor.  
> Among the uncertain physical inputs alone: Frequency 37 %, z_S 34 %, Temp 28 %.

| Output | Description |
|---|---|
| `figures/RF_permutation_importance.png` | Sorted bar chart (raw + normalised) |
| `results/feature_importance_table.csv` | Machine-readable table (opens in Excel) |
| `results/feature_importance_table.txt` | Formatted human-readable report |

---

### 5d. Out-of-distribution test — `run_DT_OOD.m`

Tests RF generalisation on 5 physically valid but far-outside-training scenarios:

| Scenario | Parameters |
|---|---|
| Low frequency | f = 5 kHz |
| High frequency | f = 15 kHz |
| Deep source | z_S = 20 m |
| Cold water | T = −5 °C |
| Combined | f = 7 kHz, z_S = 12 m, T = +3 °C |

Figures generated: OOD shadow maps, accuracy table, ROC curves per scenario, shadow boundary overlays, accuracy vs training distance scatter.

**Key finding:** RF generalises well to temperature and source depth perturbations but degrades at different frequencies, because shadow-zone geometry is wavelength-dependent (interference fringes scale with λ).

---

## 6. Comparison

**Script:** `Comparison/run_comparison.m`

Side-by-side comparison of all 5 shadow-probability methods on the full 3-parameter uncertainty set:

1. Delta Chebyshev bound
2. MC Chebyshev bound
3. MC empirical P(TL > FOM)
4. MC LN3 fit
5. RF MC-averaged prediction

| Output | Description |
|---|---|
| `figures/compare_shadow_5methods_RF.png` | 1×5 panel shadow map |
| `figures/compare_shadow_methods.png` | Earlier 4-method comparison |
| `results/comparison_results.mat` | All method outputs on the same grid |

---

## 7. Waveguide validation

**Script:** `Waveguide/run_waveguide_comparison.m`

Validates Bellhop against an exact analytical solution: normal-mode sum in a rigid-bottom constant-SVP waveguide.  Both models share identical boundary conditions so any residual difference is purely numerical (Gaussian beams vs exact mode superposition).

At f = 10 kHz, D = 35 m, ≈ 467 propagating modes:

- **RMSE ≈ 10 dB** — Bellhop's Gaussian beams smooth out sharp modal interference fringes
- **Drift < 0.15 dB/km** — after Thorp volume absorption is matched in both models

This sets a floor: MC variance > 10 dB at some grid points may partially reflect Bellhop numerical noise rather than physical uncertainty.

| Output | Description |
|---|---|
| `figures/waveguide_TL_maps.png` | Side-by-side TL maps (analytical vs Bellhop) |
| `figures/waveguide_TL_comparison.png` | Range profiles at 4 depths |
| `figures/waveguide_difference.png` | Calibrated difference map |

---

## 8. Utilities

| File | Purpose |
|---|---|
| `bellhopCached.m` | Runs Bellhop; caches result to disk keyed by parameter hash. Avoids re-running identical simulations. |
| `viewer.m` | Interactive range-depth TL viewer with hover tooltips |
| `show_MC.m` | Interactive MC viewer — hover for statistics, click for histogram + LN3 fit |

---

## 9. Run order

```
1.  run_MC.m              — generate TL cache (slowest, ~hours for N=100)
2.  run_delta.m           — Delta method Jacobians (fast, ~minutes)
3.  run_DT.m              — train Random Forest on TL cache (~5 min)
4.  run_DT_explain.m      — PDPs + influence maps (~3 min)
5.  run_DT_vardecomp.m    — OOB permutation importance (<1 min)
6.  run_DT_OOD.m          — OOD test on retrained RF (~10 min)
7.  run_comparison.m      — 5-method comparison (~1 min)
8.  run_waveguide_comparison.m — analytical validation (~2 min)
9.  show_LN3_fits.m       — LN3 fit diagnostics (~3 min)
```

All scripts self-contained: `cd` to the script's folder is handled automatically via `fileparts(mfilename('fullpath'))`.

---

## 10. Output files summary

| Folder | Key figures | Key .mat files |
|---|---|---|
| `Monte_Carlo/` | (interactive viewers) | `TL_cache/TL_zS_freq_temp.mat`, `results/MC_zS_freq_temp.mat` |
| `Monte_Carlo/LN3_fits/` | `LN3_grid_fits.png`, `LN3_KS_map.png` | — |
| `Delta_Method/` | Jacobian maps, shadow map | `results/delta_jacobians.mat` |
| `Decision_Tree/` | `RF_metrics.png`, `RF_maps.png`, `RF_PDP.png`, `RF_influence_maps.png`, `RF_permutation_importance.png`, OOD figures | `results/RF_results.mat`, `results/RF_explain.mat`, `results/feature_importance_table.csv` |
| `Comparison/` | `compare_shadow_5methods_RF.png` | `results/comparison_results.mat` |
| `Waveguide/` | `waveguide_TL_maps.png`, `waveguide_difference.png` | `results/waveguide_results.mat` |
