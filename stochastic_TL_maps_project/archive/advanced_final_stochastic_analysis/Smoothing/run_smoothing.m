%% run_smoothing.m  –  Smoothing the Delta Chebyshev Probability Map
%
% Two families of smoothing for the Delta-method shadow-probability map:
%
%   1. IDEAL 2-D LOW-PASS FILTER (FFT-based)
%      Keeps the lowest CUTOFF fraction of spatial frequencies in each axis.
%      cutoff = 0.90 → keeps 90% → removes only top 10% (mild, user request)
%      Also shows 0.70 and 0.50 for comparison.
%      Artefact: Gibbs ringing at sharp edges (esp. the shadow-zone boundary).
%
%   2. GAUSSIAN BLUR (convolution)
%      Isotropic in physical metres: σ_m applied to both depth and range
%      after converting to grid pixels.  No ringing, purely lowpass.
%      Grid spacing: Δz ≈ 0.06 m/px, Δr ≈ 65.4 m/px
%      So the same σ_m looks VERY different per axis in pixel space →
%      σ_px_z = σ_m / Δz  (many pixels in depth for a given σ_m)
%      σ_px_r = σ_m / Δr  (very few pixels in range for the same σ_m)
%
% OUTPUT
%   Fig 1: LPF comparison  (original + 3 cutoffs side by side)
%   Fig 2: Gaussian overview  (original + 4 σ values, 1×5 panel)
%   Fig 3: Gaussian detail grid  (σ values × 2 maps: Cheb_lb and Var[TL])
%
% Data: loads full-subset (z_S+Freq+SVP) from Delta_Method/results.
%       Falls back to smaller subsets if the full one is missing.

clear; close all; clc;
try, cd(fileparts(mfilename('fullpath'))); catch; end

base_dir = fullfile('..');
addpath(genpath(fullfile(base_dir, 'Functions')));

if ~exist('figures','dir'), mkdir('figures'); end

%% ── LOAD DELTA RESULTS ───────────────────────────────────────────────────
jac_file = fullfile(base_dir,'Delta_Method','results','delta_jacobians.mat');
dat_file = fullfile(base_dir,'Delta_Method','results','delta_zS_freq_svp.mat');

if ~isfile(jac_file) || ~isfile(dat_file)
    error('Run Delta_Method/run_delta.m first. Expected:\n  %s\n  %s', jac_file, dat_file);
end

J   = load(jac_file, 'TL_nom','r_km','z_m');
D   = load(dat_file, 'Cheb_lb_s','Var_TL','r_km','z_m');

TL_nom   = J.TL_nom;
Cheb_lb  = D.Cheb_lb_s;
Var_TL   = D.Var_TL;
r_km     = D.r_km;
z_m      = D.z_m;

% Grid spacing in physical units
dz = z_m(2) - z_m(1);           % m per depth pixel
dr = (r_km(2) - r_km(1)) * 1000; % m per range pixel
[Nz, Nr] = size(Cheb_lb);
FOM = 100;

fprintf('Grid: %d × %d  (Δz=%.3fm/px  Δr=%.1fm/px)\n', Nz, Nr, dz, dr);

THRESHOLDS = [0.70 0.80 0.90 0.95];

%% ══════════════════════════════════════════════════════════════════════════
%  PART 1 – IDEAL LOW-PASS FILTER (FFT)
%% ══════════════════════════════════════════════════════════════════════════
% cutoff_frac = fraction of frequency bandwidth kept per dimension.
% Mask: |f_norm| ≤ cutoff_frac/2  (f_norm ∈ [-0.5, 0.5])
% cutoff = 0.9 → removes frequencies with |f_norm| > 0.45

cutoffs = [0.90, 0.70, 0.50];
lpf_labels = {'LPF 90%  (keeps 90%)','LPF 70%  (keeps 70%)','LPF 50%  (keeps 50%)'};

% Pre-compute LPF results
lpf_cheb = cell(1, numel(cutoffs));
lpf_var  = cell(1, numel(cutoffs));
for k = 1:numel(cutoffs)
    lpf_cheb{k} = lpf2d(Cheb_lb, cutoffs(k));
    lpf_var{k}  = lpf2d(Var_TL,  cutoffs(k));
end

% ── Fig 1a: LPF on Cheb probability map ───────────────────────────────────
fig1a = figure('Name','LPF – Chebyshev Shadow Map','Position',[30 30 1800 800]);

% Panel 1: original
ax = subplot(1, 4, 1);
shadowThresholdMaps(ax, r_km, z_m, Cheb_lb, 'Original (no filter)', THRESHOLDS);

for k = 1:3
    ax = subplot(1, 4, k+1);
    shadowThresholdMaps(ax, r_km, z_m, lpf_cheb{k}, lpf_labels{k}, THRESHOLDS);
end
sgtitle(sprintf(['Delta Chebyshev P(TL>%ddB) | Ideal LPF Comparison\n' ...
    'Contours: 70%% (grey) / 80%% (orange) / 90%% (red) / 95%% (blue)'], FOM), ...
    'FontSize', 11);
saveFig(fig1a, fullfile('figures','smooth_lpf_cheb'));

% ── Fig 1b: LPF on Var[TL] ────────────────────────────────────────────────
fig1b = figure('Name','LPF – Var[TL]','Position',[60 60 1800 800]);

subplot(1,4,1);
pcolor(r_km, z_m, Var_TL); shading interp; set(gca,'YDir','reverse');
colormap(gca,hot); colorbar; title('Var[TL] Original'); xlabel('Range (km)'); ylabel('Depth (m)');

for k = 1:3
    subplot(1,4,k+1);
    pcolor(r_km, z_m, lpf_var{k}); shading interp; set(gca,'YDir','reverse');
    colormap(gca,hot); colorbar;
    title(sprintf('Var[TL] | %s',lpf_labels{k}),'FontSize',9);
    xlabel('Range (km)'); ylabel('Depth (m)');
end
sgtitle('Delta Var[TL] | Ideal LPF Comparison');
saveFig(fig1b, fullfile('figures','smooth_lpf_var'));

% ── Fig 1c: LPF frequency content visualisation ───────────────────────────
fig1c = figure('Name','LPF – Frequency Spectrum','Position',[90 90 1500 600]);
F_orig = fftshift(fft2(Cheb_lb));
logmag = log1p(abs(F_orig));

subplot(1,4,1);
imagesc([-0.5 0.5],[-0.5 0.5],logmag); colormap(gca,hot); colorbar; axis xy;
title('|FFT| original (log scale)'); xlabel('Normalised f_r'); ylabel('Normalised f_z');

for k = 1:3
    subplot(1,4,k+1);
    F_filt = fftshift(fft2(lpf_cheb{k}));
    imagesc([-0.5 0.5],[-0.5 0.5],log1p(abs(F_filt))); colormap(gca,hot); colorbar; axis xy;
    title(sprintf('|FFT| %s',lpf_labels{k}),'FontSize',9);
    xlabel('Normalised f_r'); ylabel('Normalised f_z');
end
sgtitle('Frequency-Domain Content Before and After LPF');
saveFig(fig1c, fullfile('figures','smooth_lpf_spectrum'));

fprintf('LPF figures saved.\n');

%% ══════════════════════════════════════════════════════════════════════════
%  PART 2 – GAUSSIAN BLUR
%% ══════════════════════════════════════════════════════════════════════════
% σ specified in physical METRES → different number of pixels in z vs r
% because Δz ≈ 0.06 m/px but Δr ≈ 65.4 m/px.
%
% Shown sigma values (in metres):
%   0.30 m  → σ_z≈5px  σ_r≈0.005px  (sub-meter: smooths depth fringes only)
%   0.60 m  → σ_z≈10px σ_r≈0.009px
%   2.0  m  → σ_z≈33px σ_r≈0.031px  (scale of source depth)
%   5.0  m  → σ_z≈83px σ_r≈0.077px  (full water-column scale)
%
% For significant range smoothing, the σ must be >> 65.4 m.
% Additional set with σ specified per-axis independently:
%   σ_z=2m, σ_r=500m  → σ_px_z≈33, σ_px_r≈7.6  (physically motivated 2-D Gaussian)
%   σ_z=5m, σ_r=1000m → σ_px_z≈83, σ_px_r≈15.3

% Set A: isotropic σ in physical metres
sigma_iso_m = [0.30, 0.60, 2.0, 5.0];   % metres

% Set B: anisotropic σ – physically motivated (separate z and r smoothing)
sigma_aniso = [2, 500; 5, 1000];  % [σ_z_m, σ_r_m]

% ── Pre-compute all Gaussian results ──────────────────────────────────────
gauss_iso_cheb = cell(1, numel(sigma_iso_m));
gauss_iso_var  = cell(1, numel(sigma_iso_m));
for k = 1:numel(sigma_iso_m)
    sig_m = sigma_iso_m(k);
    sz = sig_m / dz;   % σ in depth pixels
    sr = sig_m / dr;   % σ in range pixels
    gauss_iso_cheb{k} = gaussBlur2d(Cheb_lb, sz, sr);
    gauss_iso_var{k}  = gaussBlur2d(Var_TL,  sz, sr);
    fprintf('Gaussian ISO σ=%.2fm  → σ_z=%.1fpx  σ_r=%.3fpx\n', sig_m, sz, sr);
end

gauss_aniso_cheb = cell(1, size(sigma_aniso,1));
for k = 1:size(sigma_aniso,1)
    sz = sigma_aniso(k,1) / dz;
    sr = sigma_aniso(k,2) / dr;
    gauss_aniso_cheb{k} = gaussBlur2d(Cheb_lb, sz, sr);
    fprintf('Gaussian ANISO σ_z=%.0fm σ_r=%.0fm → σ_z=%.0fpx σ_r=%.1fpx\n', ...
            sigma_aniso(k,1), sigma_aniso(k,2), sz, sr);
end

% ── Fig 2: Isotropic Gaussian on Cheb map (1×5 panel) ────────────────────
fig2 = figure('Name','Gaussian ISO – Cheb Shadow','Position',[30 100 2000 500]);

ax = subplot(1, numel(sigma_iso_m)+1, 1);
shadowThresholdMaps(ax, r_km, z_m, Cheb_lb, 'Original', THRESHOLDS);

for k = 1:numel(sigma_iso_m)
    ax = subplot(1, numel(sigma_iso_m)+1, k+1);
    sig_m = sigma_iso_m(k);
    shadowThresholdMaps(ax, r_km, z_m, gauss_iso_cheb{k}, ...
        sprintf('Gauss ISO σ=%.2fm\n(σ_z=%.0fpx, σ_r=%.2fpx)', ...
                sig_m, sig_m/dz, sig_m/dr), THRESHOLDS);
end
sgtitle(sprintf(['Chebyshev P(shadow) | Isotropic Gaussian Blur (σ in physical metres)\n' ...
    'Note: Δz=%.3fm/px  Δr=%.1fm/px → same σ_m affects axes very differently'], dz, dr), ...
    'FontSize', 10);
saveFig(fig2, fullfile('figures','smooth_gauss_iso_cheb'));

% ── Fig 3: Isotropic Gaussian on Var[TL] ─────────────────────────────────
fig3 = figure('Name','Gaussian ISO – Var','Position',[60 100 2000 500]);
subplot(1, numel(sigma_iso_m)+1, 1);
pcolor(r_km, z_m, Var_TL); shading interp; set(gca,'YDir','reverse');
colormap(gca,hot); colorbar; title('Var[TL] Original');
xlabel('Range (km)'); ylabel('Depth (m)');
for k = 1:numel(sigma_iso_m)
    subplot(1, numel(sigma_iso_m)+1, k+1);
    pcolor(r_km, z_m, gauss_iso_var{k}); shading interp; set(gca,'YDir','reverse');
    colormap(gca,hot); colorbar;
    title(sprintf('σ=%.2fm  (σ_z=%.0fpx)',sigma_iso_m(k),sigma_iso_m(k)/dz),'FontSize',9);
    xlabel('Range (km)'); ylabel('Depth (m)');
end
sgtitle('Var[TL] | Isotropic Gaussian Blur (in physical metres)');
saveFig(fig3, fullfile('figures','smooth_gauss_iso_var'));

% ── Fig 4: Anisotropic Gaussian (physically motivated) ────────────────────
% Separate σ for depth and range → most physically meaningful smoothing
fig4 = figure('Name','Gaussian ANISO – Cheb Shadow','Position',[90 100 1600 700]);

all_cheb = {Cheb_lb, gauss_aniso_cheb{:}};
all_titles = {'Original (no filter)'};
for k = 1:size(sigma_aniso,1)
    all_titles{end+1} = sprintf('Aniso Gauss\nσ_z=%.0fm (%.0fpx)  σ_r=%.0fm (%.0fpx)', ...
        sigma_aniso(k,1), sigma_aniso(k,1)/dz, ...
        sigma_aniso(k,2), sigma_aniso(k,2)/dr);
end

for k = 1:numel(all_cheb)
    ax = subplot(1, numel(all_cheb), k);
    shadowThresholdMaps(ax, r_km, z_m, all_cheb{k}, all_titles{k}, THRESHOLDS);
end
sgtitle(sprintf(['Chebyshev P(shadow) | Anisotropic Gaussian (separate σ_z and σ_r)\n' ...
    'Physically motivated: larger range smoothing than depth smoothing']), 'FontSize', 10);
saveFig(fig4, fullfile('figures','smooth_gauss_aniso_cheb'));

% ── Fig 5: Direct comparison – LPF 90% vs best Gaussian ──────────────────
fig5 = figure('Name','Best Comparison – LPF vs Gauss','Position',[100 50 1900 800]);

panels = {Cheb_lb, lpf_cheb{1}, gauss_iso_cheb{2}, gauss_aniso_cheb{1}, gauss_aniso_cheb{2}};
ttls = {'Original', ...
        sprintf('LPF 90%%\n(removes top 10%% freqs)'), ...
        sprintf('Gauss ISO σ=%.2fm\n(σ_z=%.0fpx σ_r=%.2fpx)', ...
                sigma_iso_m(2), sigma_iso_m(2)/dz, sigma_iso_m(2)/dr), ...
        sprintf('Gauss Aniso\nσ_z=%.0fm σ_r=%.0fm',sigma_aniso(1,1),sigma_aniso(1,2)), ...
        sprintf('Gauss Aniso\nσ_z=%.0fm σ_r=%.0fm',sigma_aniso(2,1),sigma_aniso(2,2))};

for k = 1:numel(panels)
    ax = subplot(1, numel(panels), k);
    shadowThresholdMaps(ax, r_km, z_m, panels{k}, ttls{k}, THRESHOLDS);
end
sgtitle(sprintf('Shadow Probability Smoothing Comparison | Full 3-Param | FOM=%ddB', FOM), ...
    'FontSize', 11);
saveFig(fig5, fullfile('figures','smooth_comparison_best'));

fprintf('\n=== Smoothing complete. Figures in Smoothing/figures/ ===\n');
fprintf('Files:\n');
fprintf('  smooth_lpf_cheb.png         — LPF at 90/70/50%% cutoff on P(shadow)\n');
fprintf('  smooth_lpf_var.png          — LPF on Var[TL]\n');
fprintf('  smooth_lpf_spectrum.png     — Frequency-domain effect of LPF\n');
fprintf('  smooth_gauss_iso_cheb.png   — Isotropic Gaussian (physical metres) on P(shadow)\n');
fprintf('  smooth_gauss_iso_var.png    — Isotropic Gaussian on Var[TL]\n');
fprintf('  smooth_gauss_aniso_cheb.png — Anisotropic Gaussian (separate z/r σ)\n');
fprintf('  smooth_comparison_best.png  — Side-by-side best options\n');

%% ── LOCAL FILTER FUNCTIONS ───────────────────────────────────────────────
function out = lpf2d(img, cutoff_frac)
% Ideal rectangular 2-D LPF in the FFT domain.
% cutoff_frac: fraction of full frequency bandwidth to KEEP per axis.
%   0.90 → keep |f_norm| ≤ 0.45  (remove top 10%)
%   0.50 → keep |f_norm| ≤ 0.25  (remove top 50%)
    [Nz, Nr] = size(img);
    F = fftshift(fft2(double(img)));

    fz = (-Nz/2 : Nz/2-1) / Nz;
    fr = (-Nr/2 : Nr/2-1) / Nr;
    [FR, FZ] = meshgrid(fr, fz);

    % Rectangular pass-band: keep both axes within ±(cutoff/2)
    mask = (abs(FZ) <= cutoff_frac/2) & (abs(FR) <= cutoff_frac/2);
    out  = real(ifft2(ifftshift(F .* mask)));
    out  = max(0, min(1, out));   % clamp to valid probability range [0,1]
end

function out = gaussBlur2d(img, sigma_z_px, sigma_r_px)
% 2-D Gaussian blur with potentially different σ in z (depth) and r (range).
% sigma_z_px, sigma_r_px: standard deviations in PIXELS.
% Uses separable convolution for efficiency.
    img = double(img);

    % Depth (row) direction
    if sigma_z_px >= 0.5
        rz  = ceil(3 * sigma_z_px);
        kz  = exp(-((-rz:rz).^2) / (2 * sigma_z_px^2));
        kz  = kz / sum(kz);
        img = conv2(kz(:), 1, img, 'same');
    end

    % Range (column) direction
    if sigma_r_px >= 0.5
        rr  = ceil(3 * sigma_r_px);
        kr  = exp(-((-rr:rr).^2) / (2 * sigma_r_px^2));
        kr  = kr / sum(kr);
        img = conv2(1, kr(:)', img, 'same');
    end

    out = max(0, min(1, img));
end
