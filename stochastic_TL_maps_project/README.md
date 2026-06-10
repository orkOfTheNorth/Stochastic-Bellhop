# Stochastic Bellhop — TL Field UQ Project

Uncertainty quantification (UQ) of underwater acoustic transmission loss (TL) fields
computed by the Bellhop ray-tracing model. Three UQ methods are compared across five
propagation scenarios and five input-uncertainty distributions.

## UQ Methods

| Method | Script | Description |
|--------|--------|-------------|
| Monte Carlo (MC) | `run_MC.m` | N=50 Latin Hypercube samples, two-level caching |
| Delta / Taylor | `run_Delta.m` | First-order Taylor expansion via central-FD Jacobians |
| PCE | `run_pce.m` | Polynomial Chaos Expansion fitted to cached MC samples |

Also: `run_LN3.m` (3-parameter lognormal fits), `run_Comparison.m` (method cross-comparison).

## Scenarios

Defined in `config.json`. Each varies the bathymetry profile:

- **baseline** — flat 35 m depth, 50 km range
- **deep_water** — flat 2500 m depth, 50 km range
- **shallow_water** — flat 50 m depth, 5 km range
- **downslope** — 50 m → 550 m over 10 km
- **upslope** — 550 m → 50 m over 10 km

## Uncertain Parameters (Subsets)

Single: `freq`, `zS`, `svp`; pairs and triple combinations (7 subsets total).
Nominal values: freq=10 kHz, zS=5 m, SVP from Mackenzie formula.

## Input Distributions

`Uniform_1pct`, `Uniform_5pct`, `Uniform_10pct`, `Normal_5pct`, `Normal_10pct`
(percentage refers to the half-range or 2-sigma level relative to the nominal value).

## Directory Structure

```
config.json            — all nominal values, scenarios, distributions, MC settings
run_MC.m               — Monte Carlo runner
run_Delta.m            — Delta method runner
run_pce.m              — PCE runner
run_all.m              — sequential launcher for all methods
findings.tex           — LaTeX document summarising key results and figures
Shared_Utils/          — shared MATLAB utility functions (see Shared_Utils/README.md)
Bellhop/               — Bellhop binary and environment-file builders
Methods/               — outputs: MC/, Delta/, PCE/, Comparison/, LN3/ (see Methods/README.md)
Cache/                 — Bellhop run cache (~16 GB .mat files, git-ignored)
analytical/            — standalone analytical validation scripts
```
