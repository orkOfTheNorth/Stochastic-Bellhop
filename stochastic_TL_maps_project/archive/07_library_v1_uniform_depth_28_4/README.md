# Bellhop Function Library — v1 (uniform_depth_28_4)

Complete snapshot of the first structured Bellhop function library, originally at
`bellhop_09122025/uniform_depth_28_4/Code/Functions/`.

The analysis scripts that used this library are cross-listed in
[`02_delta_method_research/`](../02_delta_method_research/) and
[`03_monte_carlo_early/`](../03_monte_carlo_early/).

## Library functions

| File | Description |
|------|-------------|
| `simpleBellhopHazat.m` | **Main orchestrator.** Creates env/bathy files, runs `bellhop.exe`, reads output, returns TL matrix. Still in use as `early_stochastic_methods/code/Functions/simpleBellhopHazat.m`. |
| `svpMaker.m` | Generates summer/winter/constant SVP using Wilson's formula. Superseded by `Functions/makeSVPNoise.m` for UQ (which adds MLD perturbation). |
| `bathymetryMaker.m` | Constant-depth and linear-slope bathymetry builder. Still in use unchanged. |
| `createBellhopEnv.m` | Writes `bellhop.env` file. Still in use unchanged. |
| `createBellhopBty.m` | Writes `bellhop.bty` file. Still in use unchanged. |
| `bellhop.m` | Wrapper to execute `bellhop.exe`. Still in use unchanged. |
| `read_shd.m` | Reads Bellhop `.shd` output (dispatches to binary/mat/ascii). Still in use unchanged. |
| `read_shd_bin.m` | Low-level binary `.shd` reader. Still in use unchanged. |
| `plotshd.m` | Plots TL field from `.shd` file. Still in use unchanged. |
| `caxisrev.m` | Utility: set/reverse colorbar axis. |
| `colors.m` | Reads a PNG file as a colormap. |

## v1 vs current

The current library lives at `../early_stochastic_methods/code/Functions/`
and is functionally identical. The main differences introduced in `2026-05-26/`:

- `simpleBellhopHazat` now accepts `CustomSVP` and `CustomBathymetry` options
- `bellhopCached` wraps it with two-level caching
- `makeSVPNoise` replaces `svpMaker` for UQ perturbations (adds MLD coupling)
- `computeJacobian`, `ln3prob`, `shadowThresholdMaps`, etc. are new additions
