# Monte Carlo — Folder Guide

## What this folder does

Runs **N=100 Monte Carlo samples** of the Bellhop acoustic model across
7 parameter subsets and summarises how transmission-loss (TL) uncertainty
depends on which physical parameters are perturbed.

## Parameter subsets

| Subset | What is randomised |
|--------|-------------------|
| `zS` | Source depth (±1% around 5 m) |
| `freq` | Frequency (±1% around 10 kHz) |
| `temp` | Near-surface temperature shift (−1 to +1 °C) |
| `zS_freq` | Both zS and freq |
| `zS_temp` | Both zS and temp |
| `freq_temp` | Both freq and temp |
| `zS_freq_temp` | All three (full uncertainty) |

## Key outputs

| File | What you see |
|------|-------------|
| `figures/MC_<subset>_summary.png` | **Left panel**: mean TL field (E[TL]). White solid line = FOM=100 dB boundary. **Right panel**: TL variance map (Var[TL]). High variance = region where parameter uncertainty matters most. |
| `results/MC_<subset>.mat` | `MC_EX` (mean TL), `MC_Var` (variance), `MC_PrFOM` (empirical P(TL>FOM)), `Cheb_lb` (Chebyshev lower bound on shadow probability). |
| `TL_cache/TL_<subset>.mat` | Full TL cube (Nz × Nr × N, single precision). Loaded instead of re-running if bounds are unchanged. |

## How to interpret the variance map

- **High variance near shadow-zone boundaries**: small parameter changes
  shift the acoustic interference pattern, making these regions sensitive.
- **Near-zero variance** in clearly detected or clearly undetected regions:
  the model output is robust there regardless of parameter uncertainty.
- **Temperature subsets** were affected by a bug (SVP depth=5000 m instead
  of 35 m) that is now fixed in `run_MC.m`. Re-run with `force_rerun=true`
  to regenerate corrected results.

## Scripts

| Script | Purpose |
|--------|---------|
| `run_MC.m` | Main MC runner — run this to generate all results |
| `show_MC.m` | Interactive viewer for exploring individual subsets |
| `show_LN3_fits.m` | Diagnostic: checks whether 3-parameter lognormal fits the TL distribution at each grid point |
