# Decision Tree / Random Forest — Folder Guide

## What this folder does

Trains a **Random Forest classifier** on the TL samples already computed by
the Monte Carlo runs, then uses it to predict the shadow zone probability at
any combination of (freq, zS, temp_shift, range, depth) in milliseconds —
without running Bellhop again.

## Why a Random Forest?

- Learns the non-linear shadow-zone geometry from data (no analytical formula needed).
- Once trained, predicts P(TL>FOM) at a new parameter set in milliseconds.
- OOB (out-of-bag) importance reveals which input parameter drives shadow-zone
  variability the most.

## Inputs to the classifier (5 features)

| Feature | Unit | Training range |
|---------|------|----------------|
| Frequency | Hz | 9900–10100 |
| Source depth zS | m | 4.95–5.05 |
| Temperature shift | °C | −1 to +1 |
| Range | km | full grid |
| Depth | m | full grid |

## Key outputs

| File | What you see |
|------|-------------|
| `figures/RF_metrics.png` | Confusion matrix, ROC curve, OOB feature importance |
| `figures/RF_maps.png` | 4-panel: RF nominal map, RF MC-averaged, MC empirical, Delta Chebyshev |
| `figures/RF_agreement.png` | Density scatter (RF vs MC) + spatial difference map |
| `figures/RF_OOD_maps.png` | **OOD test**: Bellhop ground truth vs RF prediction at far-outside parameters |
| `figures/RF_OOD_table.png` | Accuracy summary across OOD scenarios |
| `results/RF_results.mat` | Trained RF model, metrics, maps |
| `results/RF_OOD.mat` | OOD test TL fields, RF probabilities, accuracy per scenario |

## Interpreting the RF vs MC comparison

The RF is an emulator, not a physics model. Within its training range (±1%
perturbations), it matches MC to R²≈1. Outside that range the RF extrapolates
using the closest training examples, which degrades as parameters diverge from
the training distribution — particularly for large frequency changes, since the
shadow-zone geometry depends on wavelength.

## Scripts

| Script | Purpose |
|--------|---------|
| `run_DT.m` | Train the RF and produce all in-distribution figures |
| `run_DT_OOD.m` | Out-of-distribution test on 5 physically valid far-outside scenarios |
