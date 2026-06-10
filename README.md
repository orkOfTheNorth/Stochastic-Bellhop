# Stochastic Bellhop — Underwater Acoustics UQ

Uncertainty quantification of transmission loss (TL) fields computed by the
[Bellhop](https://oalib-acoustics.org/AcousticsToolbox/) ray-tracing model.
Uncertain inputs are source frequency, source depth, and sound velocity profile.
Three analytical UQ methods are compared against a Monte Carlo ground truth across
multiple propagation scenarios and input distributions.

## UQ Methods

| Method | Key script | Summary |
|--------|-----------|---------|
| Monte Carlo | `run_MC.m` | N=50 LHS samples, parallel Bellhop runs, empirical stats |
| Delta (2nd-order Taylor) | `run_Delta.m` | Jacobian + diagonal Hessian via central FD; no extra Bellhop runs |
| PCE | `run_pce.m` | Polynomial Chaos fit (order 1–10) to cached MC samples |
| LN3 | `run_LN3.m` | 3-parameter lognormal fit to MC samples per pixel |

## How to Run

All scripts live in `stochastic_TL_maps_project/`. Run from that directory.

```matlab
% Re-run analysis on existing MC cache (PCE + LN3 + comparison figures):
matlab -batch "run_analysis_only"

% Full pipeline from scratch (runs Bellhop for MC, then analysis):
matlab -batch "run_full_pipeline"
```

Each step is cache-safe — if output files already exist the step is skipped.

## Key Outputs

```
Methods/MC/<scen>/<dist>/results/    — MC_freq.mat, MC_zS.mat, MC_svp.mat
Methods/Delta/<scen>/<dist>/results/ — delta_jacobians.mat, delta_<subset>.mat
Methods/PCE/<scen>/<dist>/results/   — pce_results.mat
Methods/LN3/<scen>/<dist>/figures/   — LN3_<subset>.png
Methods/*/figures/                   — shadow maps, variance maps, comparison PNGs
```

## Repository Layout

```
stochastic_TL_maps_project/   — main scripts and config.json
  Shared_Utils/               — shared MATLAB utilities (see Shared_Utils/README.md)
  Bellhop/                    — Bellhop binary and env-file builders
  Methods/                    — all output .mat and .png files (git-ignored large data)
  Cache/                      — Bellhop run cache (~16 GB, git-ignored)
  analytical/                 — standalone analytical validation scripts
research_sandbox/             — experimental / scratch scripts
```
