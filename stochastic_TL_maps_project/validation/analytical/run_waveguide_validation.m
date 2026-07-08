%% run_waveguide_validation.m
% Pekeris מוליך גלים (Pekeris Waveguide) — normal-mode analytical solution
% validated against Bellhop for an isovelocity shallow-water channel.
%
% ════════════════════════════════════════════════════════════════════════════
% JENSEN REFERENCES (Computational Ocean Acoustics, 2nd ed., Springer 2011)
%   §2.4.5   pp. 119–131  — Pekeris waveguide, dispersion relation, modes
%   §5.2     pp. 338–340  — Modal expansion, normalisation, pressure field
%   §5.9     pp. 386–389  — Perturbation theory for modal attenuation
% ════════════════════════════════════════════════════════════════════════════
%
% Physical setup — Jensen Fig. 2.25 (p. 119):
%   z = 0      pressure-release surface  (p = 0)
%   0 < z < D  water column, speed c1, density ρ1  (isovelocity)
%   z = D      ocean floor (interface)
%   z > D      semi-infinite fluid bottom, speed c2 > c1, density ρ2
%
% This project's bottom parameters come from config.json "nominal.geo"
% = [c_bottom  rho_ratio  alpha_dB_per_lambda]; createBellhopEnv.m writes
% geo(1) VERBATIM as the bottom half-space sound speed of the .env file,
% so here c2 = geo(1) = 1483.5 m/s, rho2/rho1 = geo(2) = 1.63,
% bottom attenuation alpha = geo(3) = 0.07 dB/lambda.
%
% ════════════════════════════════════════════════════════════════════════════
% DERIVATION ROAD-MAP  (maps each code block to Jensen equations)
%
% Step 1 — Helmholtz → Hankel transform (Jensen §2.4.5, eq. 2.176, p. 119)
%   Fourier-transform in time: ∂/∂t → -iω.  Result: Helmholtz equation
%       ∇²p + k²p = δ(r-rS),   k = ω/c = 2πf/c.
%   Hankel-transform in range r → kr:  depth equation per layer is
%       ψ'' + kz²ψ = 0,   kz² = k² - kr²   (harmonic oscillator in z).
%   → code: k1, k2
%
% Step 2 — Trapped-mode band (Jensen §2.4.5, eq. 2.187, p. 121)
%   Only modes with  k2 < kr < k1  are trapped (water propagating,
%   bottom evanescent).  Requires c2 > c1.
%   Maximum vertical wavenumber for a trapped mode:
%       kz_max = sqrt(k1² - k2²)
%   Number of trapped modes at frequency f:
%       M = floor(kz_max · D / π + 1/2)     (from eq. 2.190, p. 124)
%   → code: kzmax, M
%
% Step 3 — Dispersion relation (Jensen §2.4.5, eq. 2.185, p. 121)
%   Apply boundary conditions at z=0 and z=D to the 3-amplitude system
%   (Jensen eq. 2.182, p. 120).  The modal wavenumbers are the zeros of
%   the determinant:
%       det = 2i [ρ1·kz2·sin(kz1·D) + i·ρ2·kz1·cos(kz1·D)] = 0
%   For a TRAPPED mode kz2 = i·γ with γ = sqrt(kr²-k2²) > 0, so:
%       ρ1·γ·sin(kz1·D) + ρ2·kz1·cos(kz1·D) = 0
%   Dividing by ρ1 gives the singularity-free form used in the code:
%       g(kz) = γ·sin(kz·D) + (ρ2/ρ1)·kz·cos(kz·D) = 0    ← eq. 2.185
%   Equivalent to Jensen eq. 2.186:  tan(kz1·D) = -(ρ2/ρ1)·kz1/γ
%   The m-th root lies in  kz·D ∈ ((m-½)π, m·π), exactly one per interval.
%   → code: g(kz), fzero loop
%
% Step 4 — Mode normalisation (Jensen §5.2, eq. 5.6, p. 339)
%   Modes must satisfy the density-weighted orthonormality:
%       ∫₀^∞  Ψm²(z)/ρ(z) dz = 1
%   Water part (0<z<D, ρ1=1):
%       Nw = ∫₀^D sin²(kz·z) dz = D/2 - sin(2kz·D)/(4kz)
%   Bottom evanescent tail (z>D, ρ2):
%       Nb = sin²(kz·D)/(2γρ2)   [from ∫_D^∞ sin²(kzD)e^{-2γ(z-D)}/ρ2 dz]
%   Total norm: Nm = Nw + Nb
%   Normalised mode shape: Ψm(z) = sin(kz·z) / sqrt(Nm)
%   → code: Nw, Nb, Nm, phi_z, phi_s
%
% Step 5 — Modal attenuation (Jensen §5.9, eqs. 5.175–5.180, pp. 386–388)
%   First-order perturbation theory: add imaginary part α(z) to wavenumber.
%   For material absorption:  δk² = 2i·α(z)·ω/c(z)  (Jensen eq. 5.174).
%   Perturbed eigenvalue gives modal attenuation Im(krm) (Jensen eq. 5.175):
%       αm = (1/krm) ∫₀^∞ α(z)·ω/c(z) · Ψm²(z)/ρ(z) dz
%   Bottom contribution (Jensen eq. 5.180, evanescent tail):
%       δαm = Ψm²(D)·αb·ω / (2·krm·γ·cb·ρb)     ← dominant term here
%   Thorp volume attenuation in water added the same way.
%   → code: alpha2_Np, alphaW_Np, Im_krm, krm_c
%
% Step 6 — Modal pressure field (Jensen §5.2, eqs. 5.13–5.14, p. 340)
%   Exact form (eq. 5.13):
%       p(r,z) = [i / (4ρ(zS))] · Σm Ψm(zS)·Ψm(z)·H0^(1)(krm·r)
%   Far-field asymptotic (Hankel: H0^(1)(x) ≈ sqrt(2/πx)·e^{i(x-π/4)}):
%       p(r,z) ≈ [i / (ρ(zS)·sqrt(8πr))] · e^{-iπ/4}
%                · Σm  Ψm(zS)·Ψm(z)·e^{ikrm·r} / sqrt(krm)   (eq. 5.14)
%   → code: S = phi_z * (amp .* exp(1i*krm_c*r))
%
% Step 7 — Transmission Loss (Jensen §5.2, eqs. 5.15–5.17, p. 340)
%   Jensen ref pressure: p0(r=1) = 1/4  (3-D point source, eq. 5.16)
%   Bellhop/our convention: p0(r=1) = 1/(4π)  (unit-pressure source at 1m)
%   → factor 4π vs Jensen's 4 means the code prefactor becomes sqrt(2π/r):
%       TL = -20·log10( sqrt(2π/r) · |S| )     with ρ(zS)=1, S as above
%   Jensen eq. 5.17 gives sqrt(2/r)/ρ(zS); the π comes from the Bellhop
%   source normalisation (p0 = 1/(4π) instead of 1/4).
%   → code: TL_an = -20*log10( sqrt(2*pi./r) .* abs(S) )
%
% ════════════════════════════════════════════════════════════════════════════
% WHY DOES THE STOCHASTIC TL DISTRIBUTION FIT LN3 (3-parameter lognormal)?
% (for a reader who knows probability but has never studied ocean acoustics)
%
% The acoustic pressure at receiver (r,z) is a SUM OF COMPLEX PHASORS
% (Jensen eq. 5.14, derived above):
%
%   p(r,z; θ) = C(r) · Σ_{m=1}^{M}  A_m(θ) · exp(i · krm(θ) · r)
%
% where θ is a physical parameter (SVP shape, source depth zS, frequency f)
% that is UNCERTAIN and drawn from a Normal distribution:  θ ~ N(μ_θ, σ_θ²).
% Each mode m contributes:
%   - a PHASE     φm(θ) = krm(θ)·r   (wavenumber × range, grows with r)
%   - an AMPLITUDE A_m(θ) = Ψm(zS;θ)·Ψm(z;θ)/√krm(θ)  (mode excitation)
%
% ── STEP 1: Log-amplitude is approximately Normal (delta method) ──────────
%
%   Define  U(θ) = ln|p(r,z;θ)|.  Taylor-expand around the nominal θ₀:
%
%       U(θ) ≈ U(θ₀) + (θ-θ₀) · dU/dθ |_{θ₀}
%               ↑                  ↑
%            mean            slope × perturbation
%
%   Since θ-θ₀ ~ N(0, σ_θ²), and the expression above is LINEAR in (θ-θ₀):
%
%       U  ~  N( U(θ₀),  (dU/dθ)²·σ_θ² )       [delta method]
%
%   Therefore:  |p| = e^U  is LOGNORMAL.
%   And:        TL = -20/ln(10) · U  is NORMAL in the small-σ limit.
%
% ── STEP 2: BUT TL has a hard lower bound ────────────────────────────────
%
%   Regardless of how θ changes, sound must travel at least the geometric
%   spreading distance r.  The MINIMUM TL is set by the dominant mode alone:
%
%       γ  ≈  TL_geometric(r) = 20·log10(r) + 10·log10(2π/r)_correction
%          ≈  20·log10(r)   [dB]
%
%   This is a physical lower bound: TL ≥ γ always.
%   The EXCESS attenuation  X = TL - γ ≥ 0  represents interference loss.
%
% ── STEP 3: The excess X = TL - γ is LOGNORMAL ───────────────────────────
%
%   Why?  The interference term (from the double sum in |p|²) has a
%   MULTIPLICATIVE structure.  Write the intensity:
%
%   I = |p|² = |C|²/r  ·  [A₁² + A₂² + ... + 2A₁A₂cos(Δφ₁₂) + ...]
%                                           ↑ interference cross-terms
%
%   At DESTRUCTIVE INTERFERENCE (Δφ = π):  I → 0 → TL → ∞.  This creates
%   a very long RIGHT TAIL in the TL distribution.
%
%   Taking logarithm:  ln(X) = ln(TL-γ) ≈ sum over pairs of ln|interference|
%   → by CLT:  ln(X) is approximately Normal
%   ⟹  X = TL - γ  ~  Lognormal(μL, σL)   [LN2]
%   ⟹  TL  ~  LN3(γ, μL, σL)               [3-parameter lognormal]
%
% ── WHAT THE THREE LN3 PARAMETERS MEAN PHYSICALLY ────────────────────────
%
%   γ  (location/shift)  = geometric spreading floor = lower bound on TL
%   μL (log-mean)        = mean ln(TL-γ), grows with range and σ_θ
%   σL (log-std)         = spread of TL around the mean; driven by:
%                            - number of modes M (more modes = more nulls)
%                            - uncertainty level σ_θ (larger σ → more variation)
%                            - range r (phases randomise more at long range)
%
% ── LIMITING CASES ────────────────────────────────────────────────────────
%
%   σ_θ → 0  (nearly deterministic):     LN3 → 2-param lognormal → Normal(TL)
%   σ_θ → ∞, M → ∞ (many random modes): |p| → Rayleigh → TL → Gumbel
%   Practical 1–10% parameter spread:    LN3 fits best (see KS table in report)
%
% CONCLUSION: LN3 is not a heuristic fit.  It is theoretically motivated by
%   (a) the modal-sum structure of the acoustic field (Jensen §5.2),
%   (b) the delta method applied to log-amplitude (linear stats),
%   (c) the multiplicative interference structure creating a right-skewed,
%       bounded-below distribution for the excess TL above γ.
% ════════════════════════════════════════════════════════════════════════════
%
% Outputs:
%   validation/analytical/figures/waveguide_validation_50modes.png
%   validation/analytical/results/waveguide_validation_50modes.mat

ROOT = fileparts(fileparts(fileparts(mfilename('fullpath'))));
try; cd(ROOT); catch; end
addpath(genpath(fullfile(ROOT,'core')));
addpath(fullfile(ROOT,'binaries'));

cfg = jsondecode(fileread(fullfile(ROOT,'config.json')));

%% ── Parameters (shared with Bellhop run) ─────────────────────────────────
scen = cfg.scenarios;
if ~iscell(scen), scen = num2cell(scen); end
sc = scen{cellfun(@(s) strcmp(s.name,'shallow_water'), scen)};
if isfield(sc,'FOM_dB') && ~isempty(sc.FOM_dB), FOM = sc.FOM_dB; else, FOM = cfg.nominal.FOM_dB; end

freq = cfg.nominal.freq_Hz;      % 10000 Hz (nominal)
zS   = cfg.nominal.zS_m;         % 5 m source depth (nominal)
geo  = cfg.nominal.geo(:).';     % [c2  rho2/rho1  alpha (dB/lambda)]
D    = sc.maxDepth_m;            % 50 m  — water depth
maxR = sc.maxR_m;                % 5000 m

c1   = 1479.4;                   % isovelocity water speed → exactly 50 modes
c2   = geo(1);                   % 1483.5 m/s  bottom half-space
rho  = geo(2);                   % rho2/rho1 = 1.63
a2dl = geo(3);                   % bottom attenuation [dB/wavelength]

% ── Jensen §2.4.5 eq. (2.176) p.119: medium wavenumbers k = ω/c ──────────
k1 = 2*pi*freq/c1;   % water wavenumber
k2 = 2*pi*freq/c2;   % bottom wavenumber
assert(c2 > c1, 'Pekeris trapped modes require c2 > c1  (Jensen §2.4.5 p.121).');

% ── Jensen §2.4.5 eq. (2.187) p.121 + eq. (2.190) p.124 ─────────────────
% Trapped modes exist only for  k2 < kr < k1.
% Maximum vertical wavenumber for any trapped mode:
kzmax = sqrt(k1^2 - k2^2);
% Number of trapped modes at this frequency (Jensen eq. 2.190):
M = floor(kzmax*D/pi + 0.5);
fprintf('Pekeris מוליך גלים (waveguide): f=%g Hz, D=%g m, c1=%g, c2=%g, rho2/rho1=%g -> %d trapped modes\n', ...
        freq, D, c1, c2, rho, M);

%% ── Dispersion relation — Jensen §2.4.5 eq. (2.185) p.121 ───────────────
% det = 2i[ρ1·kz2·sin(kz1·D) + i·ρ2·kz1·cos(kz1·D)] = 0
% For trapped mode kz2 = i·γ, γ = sqrt(kr²-k2²) > 0.
% Substituting and dividing by ρ1:
%   g(kz) = γ·sin(kz·D) + (ρ2/ρ1)·kz·cos(kz·D) = 0
% Equivalent to tan(kz·D) = -(ρ2/ρ1)·kz/γ  (negative: root in ((m-½)π, m·π))
g = @(kz) sqrt(max(kzmax^2 - kz.^2, 0)) .* sin(kz*D) + rho .* kz .* cos(kz*D);

kz = zeros(M,1);
for m = 1:M
    lo = ((m-0.5)*pi/D) * (1 + 1e-12) + 1e-12;
    hi = min(m*pi/D, kzmax) * (1 - 1e-12);
    assert(g(lo)*g(hi) < 0, 'No sign change for mode %d — bracket error.', m);
    kz(m) = fzero(g, [lo hi]);
end
krm = sqrt(k1^2 - kz.^2);   % horizontal wavenumbers (real part), Jensen §2.4.5
gam = sqrt(krm.^2 - k2^2);  % bottom decay constants γm > 0

%% ── Mode normalisation — Jensen §5.2 eq. (5.6) p.339 ────────────────────
% Require:  ∫₀^∞ Ψm²(z)/ρ(z) dz = 1
%   Water integral (0<z<D, ρ1=1):   Nw = ∫₀^D sin²(kz·z) dz
%   Bottom tail  (z>D, ρ2=rho):     Nb = sin²(kz·D)·∫_D^∞ e^{-2γ(z-D)}/ρ2 dz
Nw = (D/2 - sin(2*kz*D)./(4*kz));          % analytic water integral
Nb = sin(kz*D).^2 ./ (2*gam*rho);          % analytic bottom-tail integral
Nm = Nw + Nb;                               % total norm  (Nw + Nb = Jensen eq.5.6)

%% ── Modal attenuation — Jensen §5.9 eqs. (5.175)–(5.180) pp.386–388 ─────
% First-order perturbation: Im(krm) = (1/krm) ∫ α(z)·ω/c(z)·Ψm²(z)/ρ(z) dz
% Bottom contribution (Jensen eq. 5.180, evanescent tail Ψm(z>D)):
%   δαm = Ψm²(D)·αb·ω / (2·krm·γ·cb·ρb)   [analytic tail integral]
alpha2_Np = a2dl * freq / c2 / 8.6858896;          % dB/lambda → Np/m in bottom
% Thorp volume attenuation in water (CVWT option in Bellhop .env):
f_kHz     = freq/1000;
thorp_dBkyd = 40*f_kHz^2/(4100+f_kHz^2) + 0.1*f_kHz^2/(1+f_kHz^2);
alphaW_Np = thorp_dBkyd / 914.4 / 8.6858896;       % Np/m in water
% Combine: Jensen eq. (5.177) for bottom + (5.176) for water column
Im_krm = (k2*alpha2_Np ./ krm) .* (Nb ./ Nm) ...  % bottom fraction (Jensen 5.180)
       + (k1*alphaW_Np ./ krm) .* (Nw ./ Nm);     % water fraction  (Jensen 5.176)
krm_c  = krm + 1i*Im_krm;                          % complex kr (Jensen §2.4.5.5)

%% ── Bellhop reference run (same physics: isovelocity c1, same geo) ───────
svp_type  = sprintf('const_%.6g', c1);
sim_pars  = {freq, maxR, zS, 0, 0, svp_type, sc.bathy_type, geo, FOM};
cache_dir = fullfile(ROOT, 'Cache', char(sc.name), 'bellhop_raw');
[TL_bh, r_grid, z_grid] = bellhopCached(sim_pars, cache_dir);
r = r_grid(:).';   z = z_grid(:);

%% ── Analytical mode sum — Jensen §5.2 eqs. (5.13)–(5.14) p.340 ──────────
% Normalised mode shapes at all receiver depths and at source depth:
phi_z  = sin(z * kz.') ./ sqrt(Nm.');   % Nz × M  — Ψm(z)  (Jensen eq.5.6)
phi_s  = sin(kz * zS) ./ sqrt(Nm);     % M  × 1  — Ψm(zS) (Jensen eq.5.6)
% Modal amplitude factor  Ψm(zS)/√krm  (Jensen eq. 5.14 numerator):
amp    = phi_s ./ sqrt(krm_c);         % M × 1
r_safe = max(r, 1);                    % guard r = 0 column
% S = Σm  Ψm(zS)·Ψm(z)·e^{ikrm·r}/√krm   (Jensen eq. 5.14 sum)
S      = phi_z * (amp .* exp(1i * krm_c * r_safe));  % Nz × Nr

%% ── Transmission Loss — Jensen §5.2 eqs. (5.15)–(5.17) p.340 ─────────────
% Jensen eq. 5.17 (ρ(zS)=1):  TL = -20·log10( sqrt(2/r)·|S| )
% Bellhop convention p0 = 1/(4π) instead of Jensen's 1/4 → factor π:
%       TL = -20·log10( sqrt(2π/r) · |S| )
TL_an  = -20*log10( sqrt(2*pi ./ r_safe) .* abs(S) + eps );
TL_an(:, r < 1) = NaN;   % far-field asymptotic invalid at r ~ 0

%% ── Error metrics (mask interference nulls + leaky-mode near field) ──────
r_min = 500;
valid = (ones(numel(z),1) * (r >= r_min)) & isfinite(TL_an) & isfinite(TL_bh) ...
        & (TL_bh < 110) & (TL_an < 110);
d         = TL_an - TL_bh;
mean_abs  = mean(abs(d(valid)));
med_abs   = median(abs(d(valid)));
max_abs   = max(abs(d(valid)));
bias      = mean(d(valid));
fprintf(['Pekeris מוליך גלים validation (%d modes, r >= %g m, nulls masked):\n' ...
         '  mean |dTL| = %.2f dB | median = %.2f dB | max = %.2f dB | bias = %+.2f dB\n'], ...
        M, r_min, mean_abs, med_abs, max_abs, bias);

%% ── Figure ───────────────────────────────────────────────────────────────
fig_dir = fullfile(ROOT,'validation','analytical','figures');
res_dir = fullfile(ROOT,'validation','analytical','results');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
if ~exist(res_dir,'dir'), mkdir(res_dir); end

r_km = r/1000;
tl_lo = prctile(TL_bh(valid), 2);  tl_hi = prctile(TL_bh(valid), 98);
cm_tl = colors('bar.png');

fig = figure('Visible','off','Position',[10 10 1500 850]);

ax1 = subplot(2,2,1);
imagesc(ax1, r_km, z, TL_an, [tl_lo tl_hi]);
set(ax1,'YDir','reverse'); colormap(ax1, cm_tl); cb = colorbar(ax1); cb.Label.String = 'TL (dB)';
xlabel(ax1,'Range (km)'); ylabel(ax1,'Depth (m)');
title(ax1, sprintf('Analytical normal modes — %d trapped modes (Jensen §5.2 eq.5.14)', M));

ax2 = subplot(2,2,2);
imagesc(ax2, r_km, z, TL_bh, [tl_lo tl_hi]);
set(ax2,'YDir','reverse'); colormap(ax2, cm_tl); cb = colorbar(ax2); cb.Label.String = 'TL (dB)';
xlabel(ax2,'Range (km)'); ylabel(ax2,'Depth (m)');
title(ax2, 'Bellhop (isovelocity SVP)');

ax3 = subplot(2,2,3);
dmap = d;  dmap(~valid) = NaN;
imagesc(ax3, r_km, z, dmap, [-10 10]);
set(ax3,'YDir','reverse'); colormap(ax3, redblue(256)); cb = colorbar(ax3); cb.Label.String = '\DeltaTL (dB)';
xlabel(ax3,'Range (km)'); ylabel(ax3,'Depth (m)');
title(ax3, sprintf('Analytical - Bellhop | mean |\\Delta| = %.2f dB', mean_abs), 'Interpreter','tex');

ax4 = subplot(2,2,4);
[~, iz] = min(abs(z - D/2));
plot(ax4, r_km, TL_bh(iz,:), '-',  'Color',[0.85 0.33 0.10], 'LineWidth',1.0); hold(ax4,'on');
plot(ax4, r_km, TL_an(iz,:), '--', 'Color',[0.00 0.45 0.74], 'LineWidth',1.0);
set(ax4,'YDir','reverse'); grid(ax4,'on');
xlabel(ax4,'Range (km)'); ylabel(ax4,'TL (dB)');
legend(ax4, {'Bellhop','Analytical modes (Jensen §5.2)'}, 'Location','southwest');
title(ax4, sprintf('TL slice at z = %.0f m', z(iz)));
xlim(ax4, [0 max(r_km)]); ylim(ax4, [tl_lo-5 tl_hi+15]);

sgtitle(sprintf(['Pekeris מוליך גלים (Jensen §2.4.5–§5.2) — %d normal modes\n' ...
                 'f = %g kHz, D = %g m, z_S = %g m, c_1 = %.1f m/s, c_2 = %.1f m/s, ' ...
                 '\\rho_2/\\rho_1 = %.2f, \\alpha_2 = %.2f dB/\\lambda'], ...
                M, freq/1000, D, zS, c1, c2, rho, a2dl), 'Interpreter','tex');

saveFigPNG(fig, fullfile(fig_dir,'waveguide_validation_50modes'));
close(fig);
fprintf('Saved: %s\n', fullfile(fig_dir,'waveguide_validation_50modes.png'));

%% ── Save results ─────────────────────────────────────────────────────────
save(fullfile(res_dir,'waveguide_validation_50modes.mat'), ...
     'TL_an','TL_bh','r_grid','z_grid','kz','krm','Im_krm','Nm','M', ...
     'freq','zS','D','c1','c2','rho','a2dl','r_min', ...
     'mean_abs','med_abs','max_abs','bias','-v7.3');
fprintf('run_waveguide_validation DONE\n');
