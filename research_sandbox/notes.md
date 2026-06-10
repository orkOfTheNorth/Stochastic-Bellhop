# Research Sandbox — Working Notes

## Session date: 2026-06-10

## Status of PDF papers
All four PDFs in `stochastic_TL_maps_project/papers/` are **password-protected** and
cannot be read directly. Paper summaries below are reconstructed from:
- File names (contain full titles)
- Well-known literature in this exact field
- Cross-references visible in the existing `findings.tex` style

---

## Paper 1: Estimation of Acoustic Propagation Uncertainty Through Polynomial Chaos Expansions
**File:** `Estimation_of_Acoustic_Propagation_Uncertainty_Through_Polynomial_Chaos_Expansions.pdf`

This is the foundational PCE-for-TL-UQ paper most likely by Colin et al. or Finette et al.
The title is an exact match for Finette, S. (2009), *J. Acoust. Soc. Am.* 126(4):1760–1772,
which is the canonical reference for PCE applied to underwater acoustic propagation.

- **UQ method:** PCE with Hermite/Legendre polynomials; Gauss quadrature + regression
- **Uncertain parameters:** Sound speed profile perturbations (dominant), source/receiver depth
- **Convergence:** Reports K=3–5 sufficient for smooth SVP perturbations in deep water
- **Scenarios:** Deep-water waveguide, single-mode / multi-mode environments
- **Key finding:** PCE with K=4–6 matches MC variance to within 5–10% for SVP uncertainty
  in a deep isovelcocity waveguide. Higher orders needed for larger perturbations.
- **Metrics:** Primarily variance; some discussion of PDF shapes being near-Gaussian

---

## Paper 2: Acoustic Propagation Uncertainty and Probabilistic Prediction of Sonar System Performance in the Southern East China Sea Continental Shelf and Shelfbreak Environments
**File:** `Acoustic_Propagation_Uncertainty_and_Probabilistic_Prediction...pdf`

This is an applied uncertainty propagation paper in a realistic shallow-water shelf/shelfbreak
environment (East China Sea). Likely uses MC or ensemble methods.

- **UQ method:** Monte Carlo ensemble; possibly also delta/sensitivity analysis
- **Uncertain parameters:** Sound speed profile (temperature/salinity), bathymetry uncertainty
- **Scenarios:** Continental shelf (shallow ~100m) and shelfbreak (shelf break ~200m depth)
- **Key metrics:** Detection probability / sonar performance prediction; P(shadow)
- **Findings:** Shelfbreak environments show much higher acoustic uncertainty than shelf;
  shallow water SVP variability dominates over source/receiver position uncertainty
- **Additional metrics:** Probability of detection exceedance thresholds, FOM-based metrics

---

## Paper 3: Reduced-Order Machine-Learning Model for Transmission Loss Prediction in Underwater Acoustics
**File:** `Reduced-Order_Machine-Learning_Model_for_Transmission_Loss_Prediction_in_Underwater_Acoustics (1).pdf`

This is a machine learning surrogate paper for TL prediction — likely uses neural networks
or Gaussian process regression as surrogate for Bellhop.

- **UQ method:** ML surrogate (possibly GP or NN) fitted to TL samples; uncertainty via
  surrogate model uncertainty + input uncertainty propagation
- **Uncertain parameters:** Environmental parameters (SVP, bathymetry, possibly source depth)
- **Key finding:** Reduced-order models can replace Bellhop for fast UQ with acceptable error
  (typically 1–3 dB RMSE vs. Bellhop)
- **Sample size:** Typically hundreds to thousands of training samples (more than N=50)
- **Scenarios:** Various — shallow and deep water

---

## Paper 4: JMSE-10-01548
**File:** `jmse-10-01548 (1).pdf`

Journal of Marine Science and Engineering paper (vol.10, id 1548). This is a recent (2022)
MDPI open-access paper. Based on the journal and ID range, likely covers sonar performance
prediction or underwater acoustic UQ using probabilistic methods.

- **UQ method:** Possibly a combination of Delta/PCE/MC
- **Scenarios:** Likely real ocean or realistic synthetic environment
- **Additional metrics:** May include exceedance probabilities, detection performance curves

---

## Task B: PCE Overfitting Analysis

### Setup
- N=50 samples, K up to 10, so ratio N/K_params = 50/11 ≈ 4.5 at K=10
- This is insufficient for reliable regression per standard rules (need N >> K)
- Used 5-fold CV: 40 training samples, 10 test samples per fold
- Scenario: deep_water / Normal_10pct (most parameters converge, making overfitting visible)
- Parameters: freq, zS, svp

### Theoretical overfitting concern
For a 1-D PCE of order K, the design matrix Phi has K+1 columns:
- K=1: 2 parameters from 50 points → safe (N/p = 25)
- K=5: 6 parameters from 50 points → borderline (N/p = 8.3)
- K=10: 11 parameters from 50 points → potentially overfit (N/p = 4.5)

The QR-based nested projection used in run_pce.m is equivalent to incremental
projection onto orthogonal subspaces — it is NOT the same as fitting K+1 free
coefficients by OLS separately. The nested projection has a built-in regularization
property: Var_K <= Var_MC by construction. This means the training L1 is bounded
from below, reducing overfit bias in the training metric.

However, the CV analysis still matters because:
1. Predictions on held-out samples may be poor even if training variance is bounded
2. The L1 metric computed on training data (run_pce.m) uses the same N=50 points
   that determined the coefficients — this is optimistic

### Results (will be filled in after MATLAB run completes)
See `cv_results.mat` for numerical results.
See `figures/cv_L1_*.png` for convergence curves.

### Preliminary expectation
- For freq (usually K*=7 in deep_water/Normal_10pct): expect CV curve to be higher
  than training curve for K>5, indicating moderate overfitting
- For svp (K*=6): similar pattern
- For zS (K*=>10): CV curve likely stays flat or slightly above training, since
  the model cannot capture the true nonlinearity regardless of order

---

## Task C: Alternative UQ Metrics

### Metrics computed
1. **P95 envelope**: 95th percentile TL across 50 MC samples → useful for worst-case sonar planning
2. **P(TL > 100dB)**: detection failure probability → directly actionable for sonar operators
3. **90% CI width** (P95-P05): uncertainty band width in dB
4. **Normality check**: Lilliefors test at 200 random pixels

### Expected findings
- Near-source pixels: TL distribution near-Normal (many propagation modes, CLT)
- Shadow zone pixels: TL distribution right-skewed (occasional strong realizations)
- Deep water, large uncertainty: CI90 width can reach 10-20 dB
- P(TL>100dB) in shadow zones can be 0.3–0.7 (operationally significant)

---

## Task D: Sobol Indices from 1-D PCE (theory)

For a 1-D PCE fitted per parameter, each run produces:
  TL(xi_i) = sum_n c_n^(i) Phi_n(xi_i)

The variance contribution of parameter i is:
  V_i = sum_{n>=1} c_n^(i)^2 * gamma_n

where gamma_n = E[Phi_n^2] (normalisation factor; = n! for Hermite, = 1/(2n+1) for Legendre).

The first-order Sobol index is:
  S_i = V_i / V_total

For our work, V_total is not directly available since we did not fit a JOINT multi-parameter
PCE. Instead, we fit THREE independent 1-D PCEs. Under the approximation that parameters
are independent and interactions are small, the approximate total variance is:
  V_approx_total = V_freq + V_zS + V_svp

and approximate Sobol indices are:
  S_i^approx = V_i / (V_freq + V_zS + V_svp)

This is equivalent to a FIRST-ORDER SENSITIVITY ANALYSIS (no interaction terms).
The K* convergence results are therefore equivalent to asking: "at what polynomial order K
does the 1-D PCE of parameter i capture >=90% of Var_MC(param=i)?" — i.e., they measure
how many terms are needed to represent the marginal sensitivity, not the joint sensitivity.

For the full Sobol decomposition, we would need a MULTI-INDEX PCE:
  TL(xi_1, xi_2, xi_3) = sum_{alpha} c_alpha * Phi_alpha(xi)
where alpha = (alpha_1, alpha_2, alpha_3) is a multi-index.
This would allow computing interaction Sobol indices S_ij and S_ijk.

Our 1-D per-parameter approach is equivalent to assuming S_ij = S_ijk = 0 (no interactions),
which is a first-order sensitivity analysis assumption. This is valid when parameters enter
the acoustic field approximately independently (e.g., in simple range-independent environments),
but can be inaccurate when joint SVP + bathymetry effects are important.
