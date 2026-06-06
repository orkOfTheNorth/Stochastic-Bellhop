# Archive — Stochastic Bellhop Project

Legacy and experimental scripts from before the main `advanced_final_stochastic_analysis/` pipeline.
All active work lives in `../advanced_final_stochastic_analysis/`.

## Folder Structure

| Folder | Contents |
|--------|----------|
| [`01_early_experiments/`](01_early_experiments/) | First Bellhop wrapper runs; `.fig` output; original `.zip` snapshot |
| [`02_delta_method_research/`](02_delta_method_research/) | All delta method / derivative explorations: linearity tests, coupling checks, hybrid approaches |
| [`03_monte_carlo_early/`](03_monte_carlo_early/) | Early standalone MC over water depth (before stochastic_analysis pipeline) |
| [`04_wave_simulation_interactive/`](04_wave_simulation_interactive/) | Interactive educational GUIs: Lloyd mirror, waveguide modal decomposition |
| [`06_stubs_incomplete/`](06_stubs_incomplete/) | Empty / single-line stub files — kept for reference only |
| [`07_library_v1_uniform_depth_28_4/`](07_library_v1_uniform_depth_28_4/) | Complete snapshot of the v1 Bellhop function library (`simpleBellhopHazat`, `svpMaker`, etc.) and presentation scripts |

## Historical context

### v1 workflow (`uniform_depth_28_4/`, ~April 2024)
First systematic UQ comparison: delta method vs Monte Carlo over uniform water depth U[5, 100m].
Scripts: `derivative_method.m` + `monte_carlo_method.m` + `Presentation_Main.m`.
Revealed that the linearity assumption breaks down for large depth ranges.

### Root-level exploration (various dates)
- `depedence_map.m` — tested pairwise coupling between (freq, depth, c_bot) parameters
- `hybrid_approch.m` — combined linear (c_bot) + stochastic (freq, depth) UQ
- `run_sensitivity.m` — simple finite-difference dI/dH over bathymetry
- `coupling_check_second_order.m` — computed cross-derivatives ∂²TL/∂z∂c at 4 corners to quantify Taylor approximation error

### Current pipeline (`../advanced_final_stochastic_analysis/`)
Supersedes all of the above.  See `../advanced_final_stochastic_analysis/README.md` for the full description.
