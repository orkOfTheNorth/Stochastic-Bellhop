# Methods

Output directory for all UQ method results. Each sub-directory corresponds to one
method and is organised by scenario and distribution.

## Structure

```
Methods/
├── MC/          — Monte Carlo results (see MC/README.md)
├── Delta/       — Delta / Taylor first-order results (see Delta/README.md)
├── PCE/         — Polynomial Chaos Expansion results (see PCE/README.md)
├── LN3/         — 3-parameter lognormal fit summaries
└── Comparison/  — Cross-method comparison figures
```

Each method directory follows the same two-level layout:

```
<Method>/<scenario>/<distribution>/
    results/    — .mat files with computed statistics
    figures/    — .png figures (TL maps, Chebyshev bounds, etc.)
```

## Scenarios

`baseline`, `deep_water`, `shallow_water`, `downslope`, `upslope`

## Distributions

`Uniform_1pct`, `Uniform_5pct`, `Uniform_10pct`, `Normal_5pct`, `Normal_10pct`

## Notes

- `.mat` result files are excluded from version control (see `.gitignore`).
  Re-generate them by running the corresponding `run_*.m` script.
- Figures (`.png`) are committed and tracked.
