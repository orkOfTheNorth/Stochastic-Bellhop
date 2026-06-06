# Analytical Acoustics — Exact Solutions

Analytical (closed-form) TL field models. No ray tracing — direct physics.

## Files

| File | Description |
|------|-------------|
| `method_of_images_shallow_interactive.m` | Lloyd's Mirror, shallow scale (10m×10m). Interactive: source type, sound speed, zoom. |
| `method_of_images_deep_ocean_interactive.m` | Lloyd's Mirror, deep-ocean scale (10km×2500m). Interactive: sound speed slider. |
| `waveguide_modal_interactive.m` | Normal-mode waveguide solution (D=2500m, 10km). Interactive: frequency slider 10–500 Hz. |
| `save_frequency_maps.m` | **Run this** to generate all PNG maps at multiple frequencies (no UI, batch export). |
| `figures/` | Output PNGs (created by `save_frequency_maps.m`). |

## Physics

### Method of Images (Lloyd's Mirror)

Pressure field for a point source near a pressure-release surface:

```
P(x,z) = exp(i k R1)/R1  −  exp(i k R2)/R2
```

where `R1` = distance from real source, `R2` = distance from image source (mirror reflection above surface). Surface reflection coefficient = -1.

TL = −20 log₁₀|P|

### Waveguide Normal Modes

Isovelocity waveguide (Dirichlet surface, Neumann bottom):

```
p(x,z) = Σ_m  [sin(k_zm z_s) · sin(k_zm z) · exp(−i k_xm x)] / (i k_xm D)
```

- `k_zm = (m + 0.5)π/D`  — vertical wavenumber
- `k_xm = √(k² − k_zm²)` — horizontal wavenumber (real = propagating, imaginary = evanescent)
- Propagating modes: f > (m+0.5)c/(2D)

## Frequency Maps

`save_frequency_maps.m` generates maps at:

| Simulation | Frequencies | Rationale |
|------------|-------------|-----------|
| Shallow MoI | 100, 300, 1000, 3000, 10000 Hz | λ spans depth×15 down to depth/70 |
| Deep Ocean MoI | 15, 50, 200, 1000, 10000 Hz | λ from 100m to 0.15m |
| Waveguide Modal | 15, 30, 75, 150, 300 Hz | 1→2→5→10→25 propagating modes |

## Relation to Bellhop Waveguide

`bellhop_09122025/Waveguide/` contains the same waveguide scenario solved numerically by Bellhop (ray tracing). Compare with `waveguide_modal_interactive.m` to see where rays and modes agree.
