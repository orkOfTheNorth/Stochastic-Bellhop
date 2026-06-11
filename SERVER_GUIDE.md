# Running the Stochastic-Bellhop Pipeline on a Windows Server

This guide covers everything needed to clone the repo and run the full UQ pipeline
on a Windows server over SSH. No prior knowledge of the codebase is required.

---

## Requirements

| Requirement | Notes |
|---|---|
| **MATLAB R2021a or later** | Must be installed and on the system PATH |
| **Parallel Computing Toolbox** | Optional but recommended — without it, `parfor` loops run sequentially (slower) |
| **Git** | To clone the repo |
| *GPU (NVIDIA, Compute Capability ≥ 3.5)* | Optional — CPU mode works without one |

LaTeX / pdflatex is **not** required to run the pipeline.

---

## Step 1 — Clone the repository

Open a terminal (CMD or PowerShell) and run:

```bat
git clone https://github.com/orkOfTheNorth/Stochastic-Bellhop.git
cd Stochastic-Bellhop
```

---

## Step 2 — Verify the environment

Run the setup check to confirm MATLAB can find everything it needs:

```bat
run_server.bat check_setup
```

Expected output (some items may be FAIL — see notes below):

```
[OK]   MATLAB version >= R2021a
[OK]   config.json present
[OK]   config.json parseable
[OK]   Shared_Utils on MATLAB path
[OK]   bellhop.exe on MATLAB path
[OK]   bellhopcuda.exe on MATLAB path
[OK]   Cache/ directory writable
[OK]   Methods/ directory writable
```

**Acceptable FAILs on a CPU-only server:**
- `Parallel Computing Toolbox` — parfor will run one sample at a time (slower but correct)
- `GPU available` — CPU bellhop.exe will be used automatically

If `bellhop.exe on MATLAB path` fails, the clone is missing files — check that you cloned
the full `PCE` branch: `git checkout PCE && git pull`.

---

## Step 3 — Run the pipeline

### Option A: Full pipeline (Monte Carlo + all analyses)

This runs everything from scratch. Takes several hours depending on hardware.

```bat
run_server.bat run_full_pipeline
```

### Option B: Analysis only (uses existing MC results from Cache/)

If Monte Carlo results are already in `stochastic_TL_maps_project/Cache/`, skip
the expensive MC step and run only PCE, LN3, and Comparison:

```bat
run_server.bat run_analysis_only
```

### Option C: Run individual stages

```bat
run_server.bat run_MC          # Monte Carlo sampling (slowest step)
run_server.bat run_pce         # Polynomial Chaos Expansion
run_server.bat run_LN3         # LN3 fit
run_server.bat run_Comparison  # Generate comparison figures
run_server.bat run_Delta       # Delta method
```

All modes log their output to `stochastic_TL_maps_project\run_<MODE>.log`.

---

## Step 4 — Check results

After the pipeline completes, results are in:

```
stochastic_TL_maps_project/Methods/Delta/        ← Delta method outputs
stochastic_TL_maps_project/Methods/PCE/          ← PCE outputs
stochastic_TL_maps_project/Methods/LN3/          ← LN3 outputs
stochastic_TL_maps_project/Methods/Comparison/   ← Cross-method comparison figures
```

Figures are saved as `.fig` (MATLAB) files. To view them on the server, open MATLAB
interactively and run `openfig('path/to/figure.fig')`.

---

## Step 5 — Fill findings.tex values (optional)

After a full pipeline run, you can extract the computed statistics for the
LaTeX report:

```bat
run_server.bat fill_findings
```

This prints copy-paste-ready LaTeX snippets to the console. If LaTeX is installed,
compile `stochastic_TL_maps_project/findings.tex` manually with:

```bat
cd stochastic_TL_maps_project
pdflatex findings.tex
pdflatex findings.tex
```

(Two passes are needed for cross-references.)

---

## GPU acceleration (automatic)

If the server has a compatible NVIDIA GPU (Compute Capability ≥ 3.5), the CUDA
build (`bellhopcuda.exe`) is selected automatically — no configuration needed.
The `check_setup` step shows which binary will be used.

---

## Resuming an interrupted run

All stages support resumption. If a run is interrupted, re-run the same command —
already-computed samples are detected via the `Cache/` directory and skipped.

---

## Troubleshooting

| Symptom | Fix |
|---|---|
| `matlab` not found | Add MATLAB bin folder to PATH: `set PATH=%PATH%;C:\Program Files\MATLAB\R202Xa\bin` |
| `bellhop.exe not found` | Verify you are on the PCE branch: `git checkout PCE` |
| Out of memory during MC | Reduce `N` in `config.json` → `"N": 500` (default is 1000) |
| Figures not saving | Confirm `Methods/` directory exists and is writable |
