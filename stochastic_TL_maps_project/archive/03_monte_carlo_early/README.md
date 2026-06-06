# Monte Carlo — Early Implementation

## Files

| File | Description |
|------|-------------|
| `monte_carlo_method.m` | Runs 50-iteration Monte Carlo with Latin Hypercube Sampling over water depth U[5, 100m]. Interpolates all TL runs onto a fixed master grid, then computes E[TL] and Var[TL]. First implementation of LHS-MC in the project. |

## Design notes
- Fixed N=50 (same number chosen for the 2026-05-26 scenarios)
- Grid interpolation was necessary because each Bellhop depth produces a different z-grid
- No caching — re-runs Bellhop from scratch every time

## Successor
`2026-05-26/Monte_Carlo/run_MC.m` — adds two-level caching (individual runs + full TL cube), LN3 probability computation for all 7 subsets, and visual progress bar.
