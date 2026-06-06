# Wave Simulation — Interactive GUIs

Educational / exploratory visualisations. Not part of the UQ pipeline.

## Files

| File | Description |
|------|-------------|
| `sea_wave_simulation_method_of_images.m` | Interactive 3-D GUI for acoustic pressure field using the **method of images** (Lloyd mirror). Features: point vs line source toggle, adjustable sound speed slider, real-time field update, far-field zoom. Runs standalone — no Bellhop needed. |
| `sea_wave_simulation_modal.m` | Interactive **modal decomposition** for shallow-water waveguide. Frequency slider controls the number of propagating modes. Visualises individual mode shapes and the total pressure field. |

## Purpose
Developed to build physical intuition for:
- How source depth affects the Lloyd mirror interference pattern (→ motivates why z_S is a sensitive UQ parameter)
- How frequency controls the number of propagating modes in shallow water (→ intuition for why frequency matters in the shadow zone analysis)

These are standalone demos — they do not depend on any project functions and can be run at any time.
