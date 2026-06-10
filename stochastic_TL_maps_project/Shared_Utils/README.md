# Shared_Utils

Shared MATLAB utility functions used by all UQ method scripts.
Add this directory (and subdirectories) to the MATLAB path with `addpath(genpath('Shared_Utils'))`.

## File Index

| File | Purpose |
|------|---------|
| `loadConfig.m` | Read `config.json` from the project root and return a struct |
| `subsetDefs.m` | Canonical 7 uncertainty-parameter subsets (single, pairs, triple) — single source of truth |
| `computeVarianceBounds.m` | Compute per-parameter perturbation bounds and variances for a given distribution config entry |
| `lhsSample.m` | Generate N Latin Hypercube samples for [freq, zS, svp_level] given a distribution |
| `chebyshevBound.m` | Chebyshev lower bound on P(TL > FOM): `max(0, 1 - Var/(Var+d²))` |
| `computeJacobian.m` | Central finite-difference Jacobian `(TL_plus - TL_minus) / (2h)` [dB/unit] |
| `computeHessianDiag.m` | Diagonal Hessian `d²TL/dθ²` via second-order central finite difference |
| `makeSVPNoise.m` | Generate a perturbed Sound Velocity Profile using coupled surface-temperature + MLD model |
| `ln3fit.m` | Fit a 3-parameter lognormal to TL samples and return P(TL > FOM) with fitted parameters |
| `ln3prob.m` | Thin wrapper around `ln3fit` returning only the probability scalar |
| `deltaLN3prob.m` | Moment-matched LN3 P(shadow) from Delta method moments — no MC samples needed |
| `pce_basis.m` | Evaluate orthogonal PCE basis polynomials Φ₀…Φ_K (Legendre or Hermite) at nodes ξ |
| `gaussQuadNodes.m` | Gauss quadrature nodes and weights for PCE coefficient projection |
| `plotTLmap.m` | Standard TL heatmap (range × depth) on a given axes with FOM overlay |
| `shadowCategoryMap.m` | Plot categorical shadow-probability map with N threshold colour bands |
| `overlayBathymetry.m` | Draw seafloor line and grey fill on an existing axes |
| `saveFigPNG.m` | Save a figure as a 300-dpi PNG (no PDF) |
| `saveSVPGraph.m` | Save a sound-speed vs depth PNG for SVP perturbation validation |
| `redblue.m` | Red–white–blue diverging colormap |
| `gaussBlur2d.m` | Separable 2-D Gaussian blur with independent σ per axis |
| `lpf2d.m` | Ideal rectangular 2-D low-pass filter in the FFT domain |
| `progressBar.m` | ASCII progress bar with ETA printed to the command window |
