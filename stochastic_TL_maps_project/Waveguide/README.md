# Waveguide — Folder Guide

## What this folder does

Validates Bellhop against an **analytical normal-mode solution** for an ideal
waveguide (constant sound speed, rigid bottom).  This is a ground-truth
comparison: both models should agree on the same physics, so any residual
difference quantifies Bellhop's numerical approximation error.

## The analytical model

The waveguide model assumes:
- Pressure-release surface at z = 0 (free surface)
- Rigid (100% reflecting) bottom at z = D = 35 m
- Constant sound speed c = 1500 m/s
- Volume absorption via the Thorp formula (≈ 1.16 dB/km at 10 kHz)

Under these conditions, acoustic pressure is an exact sum of propagating modes:

    P(r,z) = Σₘ  [1/(k_xm · D)] · sin(k_zm · zS) · sin(k_zm · z) · H₀⁽¹⁾(k_xm · r)

where k_zm = (m + 0.5)π/D are the mode wavenumbers and H₀⁽¹⁾ is the
Hankel function approximated in the far field as √(2/πkr)·exp(j(kr−π/4)).

## Why "waveguide difference"?

The **difference map** (analytical − Bellhop) tells you:
- **Mean offset** (removed by calibration at r=5 km): just a reference-level
  convention difference, not a physics error.
- **Residual spatial pattern**: reflects Bellhop's Gaussian beam approximation
  vs the exact modal sum.  Typical RMSE ≈ 10 dB at 10 kHz, D=35 m.
- **Range drift**: if absorption is not matched between models, the difference
  grows linearly with range.  After adding the Thorp formula to the analytical
  model, residual drift is < 0.15 dB/km.

## Key outputs

| File | What you see |
|------|-------------|
| `figures/waveguide_TL_maps.png` | Side-by-side TL maps (analytical vs Bellhop), same colour scale |
| `figures/waveguide_TL_comparison.png` | Range profiles at 4 depths: analytical (blue) vs Bellhop (red) |
| `figures/waveguide_difference.png` | Difference map after level calibration; red/blue = Bellhop over/under-predicts TL |
| `results/waveguide_results.mat` | TL fields, grids, calibration offset |

## Comparison to stochastic results

Because the analytical model has no random parameters, it cannot replace the
MC or Delta methods directly.  Its role is:
1. Confirm that Bellhop is set up correctly (rigid bottom, same SVP, same geometry).
2. Provide intuition for how mode-interference patterns form at 10 kHz in 35 m water.
3. Serve as a fast substitute for sensitivity analysis ideas that can later be
   validated with the full stochastic Bellhop runs.

## Scripts

| Script | Purpose |
|--------|---------|
| `run_waveguide_comparison.m` | Runs both models, calibrates levels, produces all figures |
| `waveguide_analytical.m` | Normal-mode TL function (Nz×Nr output) |
