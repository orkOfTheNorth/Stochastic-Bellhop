# Stochastic Bellhop — Underwater Acoustics UQ

Uncertainty quantification of transmission loss (TL) fields computed by the
[Bellhop](https://oalib-acoustics.org/AcousticsToolbox/) ray-tracing model.
Uncertain inputs are source frequency, source depth, and sound velocity profile.
Three analytical UQ methods are compared against a Monte Carlo ground truth across
multiple propagation scenarios and input distributions.

## UQ Methods

| Method | Key script | Summary |
|--------|-----------|---------|
| Monte Carlo | `run_MC.m` | N=300 base samples (Normal 10%); IS-recycled for Normal 1%/5% — no extra Bellhop runs |
| Delta (GH orders 1–10) | `run_Delta.m` | Gauss-Hermite quadrature sweep; L2 convergence vs MC reference |
| PCE | `run_pce.m` | Polynomial Chaos up to order 20; combined LOO cross-validation K\* selection |
| LN3 | `run_LN3.m` | 3-parameter lognormal fit; pixel-wise CDF RMSE maps |

## How to Run

All scripts live in `stochastic_TL_maps_project/`. The orchestrator runs all stages in sequence:

**Windows:**
```bat
run_server.bat
```

**Linux / macOS:**
```bash
./run_server.sh
```

Each stage is cache-safe — if output files already exist the step is skipped.
See [SERVER_GUIDE.md](SERVER_GUIDE.md) for full deployment instructions including Linux binary setup.

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
