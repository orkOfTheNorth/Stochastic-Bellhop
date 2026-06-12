# Running the Stochastic-Bellhop Pipeline on a Windows Server

This guide is for a partner picking up the repo cold. No prior codebase knowledge needed.  
Branch: **`12.6-READY-FOR-SERVER-RUN`**

---

## Requirements

| Requirement | Notes |
|---|---|
| **MATLAB R2021a or later** | Must be on the system PATH (`matlab` launches from CMD) |
| **Parallel Computing Toolbox** | Optional — without it `parfor` runs sequentially (slower, still correct) |
| **Git** | To clone the repo |
| **Python 3.9+** | Only needed for the Streamlit UI (not for the MATLAB pipeline) |
| *GPU — NVIDIA Compute Capability ≥ 3.5* | Optional — CPU mode works automatically |

LaTeX / pdflatex is **not** required to run the MATLAB pipeline.

---

## Step 1 — Clone the repository and switch branch

```bat
git clone https://github.com/orkOfTheNorth/Stochastic-Bellhop.git
cd Stochastic-Bellhop
git checkout 12.6-READY-FOR-SERVER-RUN
```

---

## Step 2 — Verify the environment

```bat
run_server.bat check_setup
```

Expected output:

```
[OK]   MATLAB version >= R2021a
[OK]   config.json present and parseable
[OK]   Shared_Utils/ on MATLAB path
[OK]   bellhop.exe on MATLAB path
[OK]   Cache/ directory writable
[OK]   Methods/ directory writable
```

**Acceptable FAILs on a CPU-only machine:**
- `Parallel Computing Toolbox` — parfor runs sequentially; no accuracy impact
- `GPU available` — `bellhop.exe` (CPU) is used automatically

---

## Step 3 — Run the pipeline

### Recommended: full v2 pipeline

```bat
run_server.bat run_all_v2
```

This is the **default** — calling `run_server.bat` with no argument does the same thing.

It runs all seven stages in sequence:

| Stage | Script | What it does |
|---|---|---|
| 1 | `run_MC` | 300 Normal-10% samples per scenario/param; IS reweighting for 1% and 5% |
| 2 | `run_Delta` | Gauss–Hermite quadrature orders 1–10; L2 convergence figures |
| 3 | `run_pce` | PCE up to order 20; combined LOO cross-validation; K\* selection |
| 4 | `run_LN3` | LN3 histogram diagnostics + pixel-wise CDF RMSE maps |
| 5 | `run_Comparison` | Side-by-side detection-probability maps for all methods |
| 6 | `run_MC_convergence` | EX/Var convergence scatter and derivative plots |
| 7 | `fill_findings` | Extracts computed statistics and prints LaTeX-ready values |

Total runtime estimate on a modern workstation: **4–8 hours** (dominated by Step 1 MC).

Log: `stochastic_TL_maps_project\run_run_all_v2.log`

---

### Run individual stages

Any stage can be run in isolation (e.g. to resume after a failure or re-run one step):

```bat
run_server.bat run_MC               # Monte Carlo — slowest step (~3-5 h)
run_server.bat run_Delta            # Delta GH sweep  (~30 min)
run_server.bat run_pce              # PCE             (~20 min)
run_server.bat run_LN3              # LN3 RMSE maps   (~15 min)
run_server.bat run_Comparison       # Comparison maps (~5 min)
run_server.bat run_MC_convergence   # Convergence diagnostics (~10 min)
run_server.bat fill_findings        # Print LaTeX values (seconds)
```

**All stages skip already-completed work** — if interrupted, re-run the same command and it
picks up where it left off via the `Cache/` and `Methods/` skip-if-done guards.

---

## Step 4 — Check results

After the pipeline completes, figures (`.png`) and data (`.mat`) are in:

```
stochastic_TL_maps_project/Methods/
  MC/              ← EX, Var, P_detect maps per scenario/distribution/param
  Delta/           ← GH-order sweep, L2 convergence, scatter figures
  PCE/             ← LOO-CV curves, variance scatter, K* selection
  LN3/             ← Histogram diagnostics, CDF RMSE heatmaps
  Comparison/      ← Side-by-side detection-probability maps
  MC_convergence/  ← EX/Var vs N scatter, mean curves, derivative plots
```

Figures are saved as `.png` and open in any image viewer.

---

## Step 5 — Compile the report (optional)

If MiKTeX / pdflatex is installed:

```bat
cd stochastic_TL_maps_project
pdflatex findings.tex
pdflatex findings.tex
```

Two passes are needed for cross-references. Output: `findings.pdf`.

---

## Step 6 — Launch the Streamlit UI (optional)

The interactive map viewer requires Python. Install dependencies once:

```bat
pip install -r stochastic_TL_maps_project\ui\requirements.txt
```

Then launch:

```bat
cd stochastic_TL_maps_project
streamlit run ui\app.py
```

Opens at `http://localhost:8501`.  
On a headless server, SSH-tunnel the port first:

```bat
ssh -L 8501:localhost:8501 user@server
```

Then open `http://localhost:8501` in your local browser.

**What the UI shows:**
- Sidebar: scenario, parameter, distribution, FOM (dB), detection threshold, overlay method
- Main panel: EX(TL) heatmap with colour-coded detection zones (green = P(detect) ≥ threshold)
- Q-TIP inspector: drag the Range/Depth sliders to any pixel — shows EX, Var, P(detect) for
  MC, Delta 1st/2nd/LN3, and PCE side by side, plus the MC TL histogram at that pixel

---

## GPU acceleration (automatic)

If the server has a compatible NVIDIA GPU, `bellhopcuda.exe` is selected automatically — no
configuration needed. `check_setup` shows which binary will be used.

---

## Resuming an interrupted run

Re-run the same command. Every stage checks for existing output before recomputing:
- **MC**: skips a param/scenario if `Cache/<scen>/<dist>/TL_<param>.mat` already exists
- **Delta / PCE / LN3 / Comparison**: skips if the corresponding figure or `.mat` result exists
- **run_all_v2**: skips any stage whose completion marker is present

---

## Pipeline design notes (for orientation)

**Importance sampling (IS):** The pipeline runs 300 Bellhop evaluations for `Normal_10pct`
(the widest distribution), then reweights them analytically to obtain `Normal_5pct` and
`Normal_1pct` statistics — no extra Bellhop calls. IS weights are
`w_i ∝ exp(x_i²/2 · (1/σ_large² − 1/σ_target²))`.

**Detection probability:** All maps show `P(TL < FOM) = P(detect)` — green zones are where
detection is likely, not where signal is lost. FOM defaults to 100 dB for all scenarios
except `shallow_water` (FOM = 20 dB).

**Scenarios:** `deep_water`, `shallow_water`, `upslope`, `downslope` — no baseline scenario.

**Distributions:** `Normal_1pct`, `Normal_5pct`, `Normal_10pct` — no Uniform distributions.

---

## Troubleshooting

| Symptom | Fix |
|---|---|
| `matlab` not found | Add MATLAB bin to PATH: `set PATH=%PATH%;C:\Program Files\MATLAB\R202Xa\bin` |
| `bellhop.exe not found` | Confirm you are on the correct branch: `git checkout 12.6-READY-FOR-SERVER-RUN` |
| Out of memory during MC | Reduce `N` in `config.json` → `"N": 100` (default is 300) |
| Figures not saving | Confirm `Methods/` directory is writable; it is created automatically if absent |
| Streamlit `ModuleNotFoundError` | Run `pip install -r stochastic_TL_maps_project\ui\requirements.txt` |
| PCE skipped unexpectedly | Check that MC results exist: `Methods/MC/<scen>/<dist>/results/MC_<param>.mat` |
