# Delta Method Research

All scripts exploring the delta method (first-order Taylor / finite-difference UQ).
These represent the intellectual lineage leading to `2026-05-26/Delta_Method/run_delta.m`.

## Files

| File | Description | Status |
|------|-------------|--------|
| `derivative_method.m` | Core delta method: computes ∂TL/∂z via central FD, propagates water-depth variance U[5,100m] onto a fixed master grid. First clean implementation. | Superseded |
| `run_sensitivity.m` | Simpler sensitivity script: dI/dH over bathymetry, single-parameter. | Superseded |
| `depedence_map.m` | **Key insight script.** Tests pairwise coupling between (freq, depth, c_bot): how much does ΔTL(Δfreq + Δdepth) differ from ΔTL(Δfreq) + ΔTL(Δdepth)? Validates the linearity assumption. | Reference |
| `combined_change_coupling.m` | Checks if first-order Taylor predicts combined (depth + c_bot) perturbations accurately by comparing expansion vs full simulation. | Reference |
| `coupling_check_second_order.m` | Computes cross-derivative ∂²TL/∂z∂c at all four corners of the (z, c_bot) perturbation space. Quantifies the magnitude of nonlinearity / coupling error in the delta approximation. | Reference |
| `hybrid_approch.m` | Hybrid UQ: treats c_bot as linear (delta method) and (freq, depth) stochastically (MC). Combines variance contributions via root-sum-square. | Concept — not adopted |
| `uniform_depth.m` | Comparative study: Delta (local) vs MC (global) for U[5,100m] depth. Shows they agree for small variance but diverge for wide distribution. | Reference |
| `derivative_vs_monte_carlo.m` | Loads and plots pre-computed Delta and MC results side by side with difference maps. Earliest version of the comparison visualisation. | Superseded |
| `Presentation_Main.m` | Aggregated visualisation used for presentations: Delta vs MC at two nominal depths (52.5m, 100m) with interactive cursor. | Superseded |

## Key findings from this era

1. **Coupling is small for 1-10% perturbations** — `depedence_map.m` showed that cross-term coupling was < 5% of the total variance for typical parameter ranges, justifying the additive variance formula.
2. **Linearity breaks down for large depth ranges** — `uniform_depth.m` showed Delta underestimates variance when the depth uncertainty spans an entire waveguide (5→100m).
3. **Second-order corrections matter for c_bot** — `coupling_check_second_order.m` found ∂²TL/∂z∂c was non-negligible, motivating future Hessian terms.

## Successor
`2026-05-26/Delta_Method/run_delta.m` — uses `computeJacobian.m` (swappable derivative function), `makeSVPNoise` (physically coupled SVP), and produces 7-subset analysis with multi-threshold shadow maps.
