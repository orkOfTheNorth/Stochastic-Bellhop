# Running the Stochastic-Bellhop Pipeline on a Server

This guide is for a partner picking up the repo cold — no prior codebase knowledge required.  
Branch: **`12.6-READY-FOR-SERVER-RUN`**

---

## Requirements

| Requirement | Windows | Linux |
|---|---|---|
| **MATLAB R2021a+** | Must be on PATH (`matlab` launches from CMD) | Must be on PATH |
| **Bellhop binary** | `bellhop.exe` (committed to repo) | `bellhop` — see [Step 0 (Linux)](#step-0-linux-only--get-the-bellhop-binary) |
| **Git** | Yes | Yes |
| **Python 3.9+** | Only for Streamlit UI | Only for Streamlit UI |
| **NVIDIA GPU ≥ CC 3.5** | Optional — auto-detected | Optional — auto-detected |
| Parallel Computing Toolbox | Optional — `parfor` falls back to sequential | Same |

LaTeX / pdflatex is **not** required to run the MATLAB pipeline.

---

## Step 1 — Clone and switch branch

```bash
git clone https://github.com/orkOfTheNorth/Stochastic-Bellhop.git
cd Stochastic-Bellhop
git checkout 12.6-READY-FOR-SERVER-RUN
```

---

## Step 0 (Linux only) — Get the Bellhop binary

The repo ships `bellhop.exe` (Windows). On Linux you need a native build.

### Option A — Download pre-built binary (easiest)

The Acoustics Toolbox project distributes Linux binaries:

```bash
# Download the AT toolbox (adjust URL to latest release)
wget http://oalib.hlsresearch.com/AcousticsToolbox/at.zip
unzip at.zip -d at_toolbox
# Copy the Linux bellhop binary into the repo Bellhop/ folder
cp at_toolbox/bin/bellhop stochastic_TL_maps_project/Bellhop/bellhop
chmod +x stochastic_TL_maps_project/Bellhop/bellhop
```

### Option B — Compile from source (most compatible)

Requires `gfortran` (install via `sudo apt install gfortran` or `sudo yum install gcc-gfortran`):

```bash
git clone https://github.com/A-New-BellHope/bellhop.git bellhop_src
cd bellhop_src
make
cp bellhop ../stochastic_TL_maps_project/Bellhop/bellhop
cd ..
```

### Option C — GPU / CUDA build (optional)

If the server has a CUDA-capable NVIDIA GPU, also place `bellhopcuda` (or `bellhop_cuda`) in
`stochastic_TL_maps_project/Bellhop/`. The MATLAB wrapper auto-detects it.

### Verify

```bash
stochastic_TL_maps_project/Bellhop/bellhop --help 2>&1 | head -3
# Should print the BELLHOP header line, not "command not found"
```

---

## Step 2 — Verify the environment

**Windows:**
```bat
run_server.bat check_setup
```

**Linux:**
```bash
chmod +x run_server.sh   # only needed once
./run_server.sh check_setup
```

Expected output:
```
[OK]   MATLAB version >= R2021a
[OK]   config.json present and parseable
[OK]   Shared_Utils/ on MATLAB path
[OK]   bellhop on MATLAB path
[OK]   Cache/ directory writable
[OK]   Methods/ directory writable
```

**Acceptable FAILs:**
- `Parallel Computing Toolbox` — parfor runs sequentially; correct but slower
- `GPU available` — CPU bellhop is used automatically

---

## Step 3 — Run the pipeline

### Full v2 pipeline (recommended)

**Windows:**
```bat
run_server.bat
```

**Linux:**
```bash
./run_server.sh
```

Both default to `run_all_v2`, which runs all 7 stages:

| # | Stage | Script | What it does |
|---|---|---|---|
| 1 | MC | `run_MC` | 300 Normal-10% samples; IS reweighting for 1% and 5% |
| 2 | Delta | `run_Delta` | GH quadrature orders 1–10; L2 convergence figures |
| 3 | PCE | `run_pce` | PCE up to order 20; combined LOO cross-validation; K\* |
| 4 | LN3 | `run_LN3` | LN3 histogram diagnostics + pixel-wise CDF RMSE maps |
| 5 | Comparison | `run_Comparison` | Side-by-side detection-probability maps |
| 6 | MC convergence | `run_MC_convergence` | EX/Var vs N scatter and derivative plots |
| 7 | Fill findings | `fill_findings` | Prints LaTeX-ready values for findings.tex |

**Runtime estimate:** 4–8 hours on a modern workstation (Step 1 MC dominates).  
**Log:** `stochastic_TL_maps_project/run_run_all_v2.log`

---

### Run individual stages

**Windows / Linux (just swap `run_server.bat` ↔ `./run_server.sh`):**

```bash
./run_server.sh run_MC               # Monte Carlo — slowest (~3-5 h)
./run_server.sh run_Delta            # Delta GH sweep  (~30 min)
./run_server.sh run_pce              # PCE             (~20 min)
./run_server.sh run_LN3              # LN3 RMSE maps   (~15 min)
./run_server.sh run_Comparison       # Comparison maps (~5 min)
./run_server.sh run_MC_convergence   # Convergence diagnostics (~10 min)
./run_server.sh fill_findings        # Print LaTeX values (seconds)
```

**All stages skip already-completed work.** If interrupted, re-run the same command and
it picks up where it left off via the `Cache/` and `Methods/` skip-if-done guards.

---

## Step 4 — Check results

Figures (`.png`) and data (`.mat`) land in:

```
stochastic_TL_maps_project/Methods/
  MC/              ← EX, Var, P_detect maps per scenario/distribution/param
  Delta/           ← GH-order sweep, L2 convergence, scatter figures
  PCE/             ← LOO-CV curves, variance scatter, K* selection
  LN3/             ← Histogram diagnostics, CDF RMSE heatmaps
  Comparison/      ← Side-by-side detection-probability maps
  MC_convergence/  ← EX/Var vs N scatter, mean curves, derivative plots
```

---

## Step 5 — Compile the report (optional)

```bash
cd stochastic_TL_maps_project
pdflatex findings.tex
pdflatex findings.tex   # second pass for cross-references
```

Output: `findings.pdf` (24 pages).

---

## Step 6 — Launch the Streamlit UI (optional)

Install Python dependencies once:

```bash
pip install -r stochastic_TL_maps_project/ui/requirements.txt
# On Linux you may need pip3 if pip points to Python 2:
# pip3 install -r stochastic_TL_maps_project/ui/requirements.txt
```

Launch:

```bash
cd stochastic_TL_maps_project
streamlit run ui/app.py
# On a headless server, add: --server.headless true --server.port 8501
```

**SSH tunnel** (run on your local machine):

```bash
ssh -L 8501:localhost:8501 user@server-address
```

Then open `http://localhost:8501` in your local browser.

**What the UI shows:**
- Sidebar: scenario, parameter, distribution, FOM (dB), detection threshold, overlay method
- Main panel: EX(TL) heatmap with green detection zones (P(detect) ≥ threshold)
- Q-TIP inspector: drag Range/Depth sliders → EX, Var, P(detect) for MC / Delta / PCE side by side + MC histogram

---

## GPU acceleration

The MATLAB wrapper (`Bellhop/bellhop.m`) auto-selects the CUDA build when:
- A CUDA binary is present (`bellhopcuda.exe` on Windows, `bellhopcuda` or `bellhop_cuda` on Linux)
- MATLAB's `gpuDevice()` returns compute capability ≥ 3.5

No configuration needed. `check_setup` shows which binary will be used.

---

## Pipeline design (for orientation)

**Importance sampling (IS):** Runs 300 Bellhop evaluations at `Normal_10pct` once, then
reweights analytically for `Normal_5pct` and `Normal_1pct` — no extra Bellhop calls.

**Detection probability:** All maps show P(TL < FOM) = P(detect).
Green zones = likely detection. FOM = 100 dB for all scenarios except `shallow_water` (20 dB).

**Scenarios:** `deep_water`, `shallow_water`, `upslope`, `downslope`.

**Distributions:** `Normal_1pct`, `Normal_5pct`, `Normal_10pct`.

---

## Troubleshooting

| Symptom | Fix |
|---|---|
| `matlab: command not found` | `export PATH=$PATH:/usr/local/MATLAB/R2024a/bin` (adjust year) |
| `bellhop: command not found` inside MATLAB | See Step 0 — Linux binary not placed in `Bellhop/` |
| `bellhop.exe` error on Linux | Expected — the MATLAB wrapper ignores `.exe` on Linux automatically |
| Permission denied on `run_server.sh` | `chmod +x run_server.sh` |
| Out of memory during MC | Reduce `N` in `config.json` → `"N": 100` (default 300) |
| Figures not saving | `Methods/` is created automatically; check disk space with `df -h` |
| Streamlit `ModuleNotFoundError` | `pip3 install -r stochastic_TL_maps_project/ui/requirements.txt` |
| PCE skipped unexpectedly | MC results must exist first: `Methods/MC/<scen>/<dist>/results/MC_<param>.mat` |
| `run_server.sh: line N: $'\r': command not found` | Windows line endings — fix with: `sed -i 's/\r//' run_server.sh` |
