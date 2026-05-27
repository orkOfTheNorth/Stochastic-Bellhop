# Delta Method — Folder Guide

## What this folder does

Estimates TL uncertainty using **first-order sensitivity (the Delta method)**:
rather than running hundreds of Bellhop simulations, it runs Bellhop at the
nominal parameter values plus small finite-difference perturbations, then
propagates input variance through the Jacobian to estimate output variance.

## The idea

For a smooth output y = f(x₁, x₂, x₃):

    Var[y] ≈ (∂y/∂x₁)² Var[x₁] + (∂y/∂x₂)² Var[x₂] + (∂y/∂x₃)² Var[x₃]

Each partial derivative is estimated by a finite difference:

    ∂TL/∂freq ≈  (TL(freq+Δ) − TL(freq−Δ)) / (2Δ)

This requires only 2×3 = 6 additional Bellhop runs (one pair per parameter)
instead of N=100 full runs.

## Key outputs

| File | What you see |
|------|-------------|
| `figures/delta_<subset>.png` | **Left**: Delta-estimated Var[TL]. **Right**: Chebyshev lower bound on P(TL>FOM). White line = FOM contour. |
| `results/delta_<subset>.mat` | `TL_nom`, `Jac_freq/zS/temp` (Jacobian maps), `Var_delta`, `Cheb_lb_s`. |
| `results/delta_jacobians.mat` | Combined Jacobians used by Comparison and Decision_Tree. |

## Chebyshev bound

The Cantelli inequality gives a conservative (guaranteed) lower bound on
shadow probability without assuming a distribution:

    P(TL > FOM) ≥ 1 − Var[TL] / (Var[TL] + (E[TL] − FOM)²)   when E[TL] > FOM

- **Solid white line** on figures: FOM=100 dB boundary in mean TL
- **Dashed white line**: 95% Chebyshev contour — region guaranteed (with ≥95% probability) to be in the shadow zone

## Scripts

| Script | Purpose |
|--------|---------|
| `run_delta.m` | Main Delta method runner — generates all Jacobians and uncertainty maps |
