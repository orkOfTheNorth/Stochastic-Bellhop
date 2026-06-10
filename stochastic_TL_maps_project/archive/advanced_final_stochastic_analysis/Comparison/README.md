# Comparison — Folder Guide

## What this folder does

Loads the results from **Delta Method** and **Monte Carlo** (and optionally
the Random Forest) and directly compares them to reveal where the two
analytical approaches agree and where they diverge.

## Figures produced

### All 7 subsets

| File | What you see |
|------|-------------|
| `compare_EX_diff_all.png` | **E[TL] difference (MC − Delta)** across all 7 subsets. Red = MC predicts higher TL; blue = Delta predicts higher TL. Large values reveal where the Delta method's linear approximation breaks down (near shadow-zone boundaries with high curvature). |
| `compare_Var_diff_all.png` | **Var[TL] difference (MC − Delta)** across all 7 subsets. Positive (red) = MC variance is higher → Delta underestimates uncertainty due to nonlinearity. |
| `compare_shadow_all.png` | **Shadow-zone membership at 95%** encoded as a 4-colour map: grey=neither, blue=Delta only, red=MC only, green=both agree. Regions where they disagree are the most uncertain. |

### Full subset (zS + freq + temp) — detailed

| File | What you see |
|------|-------------|
| `compare_shadow_4methods_full.png` | 4-panel binary shadow map: Delta Chebyshev, MC Chebyshev, MC Empirical, MC LN3. All at the 95% threshold. Use this to judge how conservative each method is. |
| `compare_shadow_5methods_RF.png` | Same as above but with the **Random Forest** as a 5th panel (requires `run_DT.m` to have been run first). Shows how the trained RF shadow boundary compares to physics-based methods. |
| `compare_EX_Var_full_interactive.png` | Full-subset E[TL] and Var[TL] difference maps with an **interactive data cursor** — click anywhere in MATLAB to read exact values at that range/depth. |

## Key interpretation guide

| Observation | Meaning |
|-------------|---------|
| Large E[TL] difference near shadow boundary | Strong nonlinearity: Delta linearisation is inaccurate here |
| High Var[TL] difference for temp subsets | May reflect the SVP-depth bug in run_MC.m (fixed — re-run MC) |
| Delta Chebyshev much larger shadow zone than MC | Delta is conservative: overestimates variance in nonlinear regions |
| RF shadow zone matches MC empirical | RF has learned the shadow-zone geometry from data |
| RF shadow zone differs from MC | RF extrapolation error or smooth probability vs binary threshold |

## Scripts

| Script | Purpose |
|--------|---------|
| `run_comparison.m` | Generates all figures; requires Delta and MC results |
