# Bellhop Environment

Bellhop acoustic propagation executable and MATLAB wrapper.

## Contents

| File | Description |
|------|-------------|
| `bellhop.exe` | Bellhop ray-tracing propagation program executable |
| `bellhop.m` | MATLAB wrapper — finds and runs `bellhop.exe` via `which()` |
| `bellhop.bty`, `bellhop.env`, etc. | Working run files from Bellhop test runs (auto-generated) |

## Usage

Add this folder to MATLAB's path so that `bellhop.m` and `bellhop.exe` are found together:

```matlab
addpath(fullfile(base_dir, '../bellhop_env'));
```

The `bellhopCached.m` function in `advanced_final_stochastic_analysis/` already handles this
via the `early_stochastic_methods/code/Functions/` addpath, which also contains a copy of the exe.
This folder is the authoritative location for the Bellhop executable.

## Version

AT&T Bellhop — ocean acoustic propagation using Gaussian beam ray tracing.
