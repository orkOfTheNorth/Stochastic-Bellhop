# RDC Claude Briefing — Stochastic Bellhop UQ Project
**Branch:** `13.0-ACADEMIC-RESTRUCTURE`  
**Date written:** 2026-06-19  
**Purpose:** Full context handoff so you can continue without needing the previous conversation.

---

## 1. Who You Are Helping

Engineering/research student at Technion (SIPL lab). The project is a MATLAB underwater acoustics Uncertainty Quantification (UQ) pipeline called Stochastic Bellhop. The overall goal is to produce three fully compilable academic deliverables from the pipeline results.

**Critical terminology rule (never break this):**
The finite-difference Jacobian method must ALWAYS be written as **"Bellhop Finite-Difference Jacobian"** (full name at first occurrence in any section) and **"Bellhop-FD"** for short thereafter. Never write "FD" alone in any report or document. This was an explicit user mandate.

**Work style:** User wants autonomous judgment calls. Only pause for truly irreversible decisions (e.g., deleting production data).

---

## 2. Project Overview

**What the pipeline does:**
- Runs the Bellhop ray-tracer acoustic forward model to compute TL(r,z) [dB] fields on a 766×574 range-depth grid
- TL = Transmission Loss, physically 60–150 dB, at ranges 0–100 km and depths 0–2500 m
- Applies 4 UQ methods to quantify how TL varies under environmental uncertainty

**5 propagation scenarios:** baseline (flat 50 m), shallow_water (50 m), deep_water (2500 m), upslope (50→20 m), downslope (50→550 m)

**3 uncertain input parameters:**
- frequency f ~ N(10000 Hz, σ_f²)
- source depth z_S ~ N(5 m, σ_z²)
- SVP (sound velocity profile) temperature shift ΔT ~ N(0, σ_s²)

**3 uncertainty levels:** Normal 1%, 5%, 10% (σ = p·μ/200 for freq and zS)

**4 UQ methods:**
1. **Monte Carlo (MC):** N=1000 Bellhop samples (LHS+IID+SVP-perturbed). Ground truth.
2. **Delta/FOSM:** Var[TL] ≈ Σ Jᵢ² σᵢ² via Bellhop-FD Jacobians (6 extra Bellhop runs). 2nd-order Hessian correction also computed.
3. **PCE (Polynomial Chaos Expansion):** TL ≈ Σ cₙ Φₙ(ξ), probabilist Hermite basis, QR regression on N=50 cached samples, LOO-CV order selection K*.
4. **LN3 (3-parameter lognormal):** X−γ ~ LN2, γ = 0.95·min(samples), method of moments, KS goodness-of-fit.

**Key findings:**
- Bellhop-FD Jacobian validated against closed-form waveguide: L1 < 5% ✓
- Delta method works only for flat geometry + ≤1% perturbation. At 5%, 2nd-order Hessian correction exceeds 96% of total variance — Taylor expansion has diverged.
- PCE converges for freq and SVP (K*=3–6, zero new Bellhop runs). **z_S does NOT converge at any K≤20 in deep water** — SOFAR channel modal excitation nonlinearity. This is a new/notable finding.
- LN3: D_KS ∈ [0.076, 0.167]. Formally rejected at N=1000 but practically a close geometric fit.
- Operational metric: report P_fail = P(TL > FOM) and P95, not variance alone.

---

## 3. Directory Structure

```
stochastic_TL_maps_project/
├── config.json                  ← all parameters (MC.N=1000, PCE.max_order=20, freq_Hz=10000, zS_m=5)
├── pipeline/                    ← main scripts: run_mc.m, run_delta.m, run_pce.m, run_ln3.m,
│                                   run_comparison.m, run_mc_convergence.m, fill_findings.m
├── core/                        ← math/, stats/, bellhop/, config/, io/, plot/ libraries
├── binaries/                    ← bellhop.exe
├── docs/                        ← paper.tex, presentation.tex, poster.tex, references.bib,
│                                   findings.tex, block_diagram.tex, RDC_BRIEFING.md (this file)
├── validation/analytical/       ← analytical waveguide validation figures
├── Methods/                     ← OLD production figures (some committed via .gitignore exceptions)
└── Methods/Methods/             ← CURRENT production figures (Google Drive double-nesting artifact)
    ├── MC/{scenario}/{dist}/figures/
    ├── Delta/{scenario}/{dist}/figures/   and  Delta/{scenario}/delta_order_diff_relative.png
    ├── PCE/{scenario}/{dist}/figures/     and  PCE/pce_summary_all.png
    ├── LN3/{scenario}/{dist}/figures/
    └── Comparison/{scenario}/{dist}/figures/
```

**IMPORTANT — figure paths:** All LaTeX docs use `\graphicspath` with BOTH `../Methods/XX/` AND `../Methods/Methods/XX/` paths. The `Methods/Methods/` nesting is a Google Drive transfer artifact — don't try to fix it.

**Cache:** TL cube files (~16 GB) are on the server only. `Cache/` on this machine is nearly empty. You cannot re-run MC or PCE without the cache.

---

## 4. What Was Already Done (don't redo)

### Pipeline code fixes (all committed to branch):
- **`pipeline/run_pce.m` line ~170:** `C_all{pi}` now uses `C_all{pi} = R_qr \ QTL` (QR-stable) instead of the buggy `(Phi'*Phi + lambda*eye)\(Phi'*TL_mat)` (normal equations, catastrophically ill-conditioned at MAX_ORDER=20). This fixes pce_dist_*.png figures that were showing PCE at 0 dB vs MC at 100+ dB. Needs server re-run to regenerate figures.
- **`pipeline/run_delta.m`:** bias_2nd_order colorscale now uses `prctile(abs(finite_vals), 99.5)` instead of `max(abs(m(:)))` — prevents resonance-peak outliers from making the entire figure appear as uniform color stripes.
- **`pipeline/run_comparison.m`:** var_MC_vs_Delta_3x2 colorscale now uses `prctile(all_v_flat, 99.5)` instead of `max(all_v)` — prevents Delta variance outliers (1e6+ dB²) from making the figure all-black.
- **`pipeline/run_mc.m`:** combined_Var colorscale now uses `prctile` instead of `max` — same issue.

### Academic deliverables (all committed, all compile):
- **`docs/paper.tex`:** IEEEtran two-column, compiles to 5 pages. Uses `\bibliography{references}`.
- **`docs/presentation.tex`:** Beamer 19 slides, Technion blue (RGB 0,59,126). Compiles clean.
- **`docs/references.bib`:** 9 BibTeX entries (Finette 2009, Porter 1987, GUM, Sudret 2008, Blatman 2010, Brekhovskikh 1991, JMSE 2022, East China Sea 2020, Iman 2008).
- **`docs/block_diagram.tex`:** Math-focused pipeline block diagram. Compiles to 1-page PDF.
- **`docs/findings.tex`:** Main technical report (27 pages), `\graphicspath` covers both Methods/ and Methods/Methods/ paths.

### Poster status — NEEDS FIX:
- **`docs/poster.tex`:** Created but does NOT compile cleanly. Two errors:
  1. `subfigure` environment inside Beamer blocks is not a float — fails with "subfigure outside float"
  2. Some figure paths reference old `Methods/` paths where files don't exist

---

## 5. Your Main Task

### Task A — Fix `docs/poster.tex` so it compiles

**Root cause of errors:**
The poster uses `\begin{subfigure}` from the `subcaption` package inside Beamer block environments. Beamer blocks are not floats, so `subfigure` and `\captionof` fail.

**Fix:** Replace every `\begin{subfigure}...\end{subfigure}` group with `\begin{minipage}` pairs. Replace `\captionof{figure}{...}` with `\\[3pt]{\small ...}` plain text below the image.

**Example — before (broken):**
```latex
\begin{center}
  \begin{subfigure}{0.47\linewidth}
    \includegraphics[width=\linewidth]{some_figure.png}
    \caption{\small Some caption}
  \end{subfigure}
  \hfill
  \begin{subfigure}{0.47\linewidth}
    \includegraphics[width=\linewidth]{other_figure.png}
    \caption{\small Other caption}
  \end{subfigure}
\end{center}
```

**After (working):**
```latex
\begin{center}
  \begin{minipage}{0.47\linewidth}\centering
    \includegraphics[width=\linewidth]{some_figure.png}\\[3pt]
    {\small Some caption}
  \end{minipage}
  \hfill
  \begin{minipage}{0.47\linewidth}\centering
    \includegraphics[width=\linewidth]{other_figure.png}\\[3pt]
    {\small Other caption}
  \end{minipage}
\end{center}
```

Also remove `\usepackage{subcaption}` from the preamble after fixing.

**Standalone `\captionof{figure}{...}` calls** (for single-image blocks): replace with:
```latex
{\small\textit{Caption text here.}}
```
placed directly below `\includegraphics`.

**After fixing, compile with:**
```
cd stochastic_TL_maps_project/docs
pdflatex -interaction=nonstopmode poster.tex
```
Iterate until zero `!` errors. Warnings about missing figures are OK (figures are on server).

### Task B — Re-run pipeline on server (cache IS available here)

The server has the full `Cache/` directory with all TL cubes (~16 GB). Re-runs are possible and expected. Run these to regenerate corrected figures:

```matlab
cd stochastic_TL_maps_project
run_pce         % PRIORITY: fixes pce_dist_*.png (were showing PCE at 0 dB — QR fix)
run_delta       % fixes bias_2nd_order.png (colorscale was ±1500 dB stripes)
run_comparison  % fixes var_MC_vs_Delta_3x2.png (was all-black)
run_mc          % fixes combined_Var.png (colorscale)
```

You can also run a single scenario/distribution to test:
```matlab
run_pce('deep_water', 'Normal_10pct')
run_delta('deep_water', 'Normal_10pct')
```

**skipIfDone logic:** The pipeline checks if output figures already exist and skips them. To force regeneration of specific figures, delete them first then re-run. Do NOT bulk-delete — only delete the specific broken figures.

**Which figures are known broken (delete these to force regen):**
- ALL `pce_dist_*.png` across all 4 scenarios × 3 distributions (12 files total) — PCE was at 0 dB
- `bias_2nd_order.png` in each `Methods/Methods/Delta/{scenario}/{dist}/figures/` — ±1500 dB colorscale
- `var_MC_vs_Delta_3x2.png` in each `Methods/Methods/Comparison/{scenario}/{dist}/figures/` — all-black
- `combined_Var.png` in each `Methods/Methods/MC/{scenario}/{dist}/figures/` — colorscale saturated

**Full pipeline re-run** (runs everything, uses cache — won't redo Bellhop runs):
```matlab
run_full_pipeline   % or run_all_v2
```

### Task C — Check paper.tex for any missing figure warnings

Run pdflatex twice (for cross-references):
```
pdflatex -interaction=nonstopmode paper.tex
bibtex paper
pdflatex -interaction=nonstopmode paper.tex
pdflatex -interaction=nonstopmode paper.tex
```
Check `paper.log` for `LaTeX Warning: File ... not found`. Any missing figures should be replaced with alternative figures from `Methods/Methods/` that DO exist.

---

## 6. Key File Locations (confirmed existing)

| Figure | Full path |
|--------|-----------|
| PCE convergence, deep water 10% | `Methods/Methods/PCE/deep_water/Normal_10pct/figures/pce_convergence_deep_water_Normal_10pct.png` |
| PCE summary all scenarios | `Methods/Methods/PCE/pce_summary_all.png` |
| Delta order ratio, deep water | `Methods/Methods/Delta/deep_water/delta_order_diff_relative.png` |
| MC detect empirical freq | `Methods/Methods/MC/deep_water/Normal_10pct/figures/detect_empirical_freq.png` |
| MC detect empirical zS | `Methods/Methods/MC/deep_water/Normal_10pct/figures/detect_empirical_zS.png` |
| Comparison shadow vs Delta freq | `Methods/Methods/Comparison/deep_water/Normal_10pct/figures/shadow_50_MC_vs_Delta_freq.png` |
| Analytical validation 5% | `validation/analytical/figures/delta_validation_Normal_5pct.png` |
| Block diagram PDF | `docs/block_diagram.pdf` |

---

## 7. LaTeX Conventions Used Across All Docs

```latex
\newcommand{\EX}{\mathbb{E}}
\newcommand{\Var}{\mathrm{Var}}
\newcommand{\TL}{\mathrm{TL}}
```

`\graphicspath` always includes both paths:
```latex
\graphicspath{
  {../Methods/MC/}   {../Methods/Delta/}   {../Methods/LN3/}
  {../Methods/PCE/}  {../Methods/Comparison/}
  {../Methods/Methods/MC/}   {../Methods/Methods/Delta/}
  {../Methods/Methods/LN3/}  {../Methods/Methods/PCE/}
  {../Methods/Methods/Comparison/}
  {../validation/analytical/figures/}
}
```

Technion blue: `RGB{0,59,126}` (used in Beamer `\usecolortheme`).

Colors in recommendation tables:
```latex
\definecolor{successgreen}{RGB}{0,100,40}
\definecolor{failred}{RGB}{160,30,30}
\definecolor{warnyellow}{RGB}{160,100,0}
```

---

## 8. Git Workflow

```bash
# Pull latest on RDC
git pull origin 13.0-ACADEMIC-RESTRUCTURE

# After making changes, commit and push
git add stochastic_TL_maps_project/docs/poster.tex
git commit -m "Fix poster.tex: replace subfigure with minipage for Beamer compat"
git push origin 13.0-ACADEMIC-RESTRUCTURE
```

Main branch is `main`. All current work is on `13.0-ACADEMIC-RESTRUCTURE`. Do NOT merge or push to main without user approval.

---

## 9. What NOT to Do

- Do NOT restructure or reorganize the MATLAB codebase (no moving files, no renaming functions). That plan was explicitly deferred.
- Do NOT delete any figures from `Methods/` or `Methods/Methods/`.
- Do NOT touch `Cache/` directory.
- Do NOT add features or abstractions beyond what the task requires.
- Do NOT write "FD" alone — always "Bellhop Finite-Difference Jacobian" or "Bellhop-FD".
