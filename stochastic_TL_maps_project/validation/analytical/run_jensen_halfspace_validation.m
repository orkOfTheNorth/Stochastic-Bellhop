%% run_jensen_halfspace_validation.m
% Compares Bellhop against Jensen's EXACT method-of-images half-space
% solution (jensenHalfspaceImage.m, Jensen et al. 2011 Sec. 2.3.4 Eq. 2.76).
%
% Bellhop is made to simulate a true half-space (no waveguide, no bottom
% reflection at all) by setting the bottom geoacoustic parameters to
% exactly match the water column: geoAcoustics = [1, 1, 0] means
% c_bottom/c_water = 1, rho_bottom/rho_water = 1, attenuation = 0 — an
% acoustically invisible bottom, i.e. Bellhop just sees an unbounded
% homogeneous fluid, same as Jensen's idealization.
%
% Isovelocity SVP (svpType = 'const_1500') on both sides — this is a
% homogeneous-medium problem, no depth-varying sound speed.
%
% Grid matches the deep_water scenario size (maxR=50 km, maxDepth=2500 m).
%
% Output: validation/analytical/figures/jensen_halfspace_validation.png
%         (3 panels: Jensen TL, Bellhop TL, difference map)

ROOT = fileparts(fileparts(fileparts(mfilename('fullpath'))));
cd(ROOT);
addpath(genpath(fullfile(ROOT,'core')));
addpath(fullfile(ROOT,'binaries'));
addpath(fullfile(ROOT,'validation','analytical'));

%% Parameters (matching deep_water scenario + nominal source)
cfg   = loadConfig();
freq  = cfg.nominal.freq_Hz;      % 10000 Hz
zS    = cfg.nominal.zS_m;         % 5 m
c     = 1500;                     % isovelocity sound speed [m/s]
maxR  = 50000;                    % m, matches deep_water
FOM   = cfg.nominal.FOM_dB;

sim_pars = {freq, maxR, zS, 0, 0, 'const_1500', 'const_2500', [1 1 0], FOM};

fprintf('=== Jensen half-space (image method) vs Bellhop ===\n');
fprintf('freq=%g Hz, zS=%g m, c=%g m/s, geo=[1 1 0] (impedance-matched bottom)\n', freq, zS, c);

cache_dir = fullfile(ROOT, 'Cache', 'validation_halfspace');
[TL_bellhop, r_grid, z_grid] = bellhopCached(sim_pars, cache_dir);

%% Jensen exact solution on the same grid
[R, Z] = meshgrid(r_grid, z_grid);
p_jensen  = jensenHalfspaceImage(freq, zS, R, Z, c);
TL_jensen = pressureToTL(p_jensen, 1);   % p_ref = 1 at R=1m, matches Bellhop's own TL convention

%% Difference
% Exclude the r=0 column (R1=0 there is a true singularity of the
% point-source solution, TL -> +Inf, not a modelling discrepancy) AND any
% other non-finite pixels explicitly — 'omitnan' does NOT drop Inf values.
diff_map = TL_bellhop - TL_jensen;
valid_r  = r_grid > 0;
d_valid  = diff_map(:, valid_r);
d_valid  = d_valid(isfinite(d_valid));
fprintf('mean|diff| = %.3f dB, median|diff| = %.3f dB, max|diff| = %.3f dB  (r>0, finite only)\n', ...
    mean(abs(d_valid)), median(abs(d_valid)), max(abs(d_valid)));

% Sanity check: is the gap explained by Bellhop's built-in Thorp volume
% attenuation? (core/bellhop/createBellhopEnv.m hardcodes the SSP options
% string 'CVWT' — the trailing 'T' turns Thorp attenuation ON; Jensen's
% eq. 2.76 has NO volume attenuation term, by design — see jensenHalfspaceImage.m.)
% Thorp's formula (f in kHz, alpha in dB/km):
f_kHz = freq/1000;
alpha_thorp_dBkm = 0.11*f_kHz^2/(1+f_kHz^2) + 44*f_kHz^2/(4100+f_kHz^2) ...
                   + 2.75e-4*f_kHz^2 + 0.003;
far_field_gap = TL_bellhop(:, end-5:end-1) - TL_jensen(:, end-5:end-1);
far_field_extra_dB = mean(far_field_gap(isfinite(far_field_gap)));
r_far_km = mean(r_grid(end-5:end-1))/1000;
fprintf('Thorp predicted extra loss at %.1f km: %.1f dB.  Observed extra loss (Bellhop-Jensen) there: %.1f dB.\n', ...
    r_far_km, alpha_thorp_dBkm*r_far_km, far_field_extra_dB);

%% Figure
fig = figure('Position',[50 50 1500 420], 'Visible','off');
clim_tl = [50 100];

ax1 = subplot(1,3,1);
imagesc(ax1, r_grid/1000, z_grid, TL_jensen, clim_tl);
set(ax1,'YDir','reverse'); colormap(ax1,'jet'); colorbar(ax1);
xlabel(ax1,'Range (km)'); ylabel(ax1,'Depth (m)');
title(ax1,'Jensen (image method, exact)');

ax2 = subplot(1,3,2);
imagesc(ax2, r_grid/1000, z_grid, TL_bellhop, clim_tl);
set(ax2,'YDir','reverse'); colormap(ax2,'jet'); colorbar(ax2);
xlabel(ax2,'Range (km)');
title(ax2,'Bellhop (impedance-matched bottom)');

ax3 = subplot(1,3,3);
mx = max(abs(d_valid(:)), [], 'omitnan');
imagesc(ax3, r_grid/1000, z_grid, diff_map, [-mx mx]);
set(ax3,'YDir','reverse'); colormap(ax3, redblue(256)); colorbar(ax3);
xlabel(ax3,'Range (km)');
title(ax3, sprintf('Bellhop - Jensen (median|diff|=%.1f dB, r>0) — range-growing gap = Thorp attenuation', ...
    median(abs(d_valid(:)),'omitnan')));

sgtitle('Half-space validation: Bellhop vs Jensen (2011) Sec. 2.3.4, Eq. 2.76 (image method only, no volume attenuation)', 'Interpreter','none');

out_dir = fullfile(ROOT, 'validation', 'analytical', 'figures');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
out_png = fullfile(out_dir, 'jensen_halfspace_validation.png');
print(fig, out_png, '-dpng', '-r150');
close(fig);
fprintf('Saved: %s\n', out_png);
