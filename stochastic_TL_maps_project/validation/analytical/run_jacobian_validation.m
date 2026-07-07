%% run_jacobian_validation.m
% Validates the Delta method's core assumption: is the FINITE-DIFFERENCE
% Jacobian (core/math/computeJacobian.m — exactly what pipeline/run_delta.m
% uses) actually a good approximation of the TRUE (analytical) derivative
% dTL/dtheta?
%
% Since we have a closed-form solution for the half-space problem
% (validation/analytical/jensenHalfspaceImage.m, Jensen 2011 Eq. 2.76), we
% can differentiate it ANALYTICALLY in freq and zS, and compare against
% Bellhop's own finite-difference Jacobian (same central-difference recipe
% and step sizes as run_delta.m: cfg.jacobian_steps).
%
% Setup: deep-water-sized grid (maxR=50 km, maxDepth=2500 m), CONSTANT
% (isovelocity) sound speed — no SVP perturbation here, this isolates the
% differentiation question from the SVP-model question — and a fully
% penetrable (impedance-matched) bottom, same as run_jensen_halfspace_validation.m.
%
% Output: validation/analytical/figures/jacobian_validation.png
%         (2 rows: dTL/dfreq and dTL/dzS; 3 cols: analytical, numerical (FD), diff)

ROOT = fileparts(fileparts(fileparts(mfilename('fullpath'))));
cd(ROOT);
addpath(genpath(fullfile(ROOT,'core')));
addpath(fullfile(ROOT,'binaries'));
addpath(fullfile(ROOT,'validation','analytical'));

cfg   = loadConfig();
freq0 = cfg.nominal.freq_Hz;
zS0   = cfg.nominal.zS_m;
c     = 1500;
maxR  = 50000;
FOM   = cfg.nominal.FOM_dB;
h_freq = cfg.jacobian_steps.freq_Hz;   % 50 Hz
h_zS   = cfg.jacobian_steps.zS_m;      % 0.025 m

cache_dir = fullfile(ROOT, 'Cache', 'validation_halfspace');

fprintf('=== Jacobian validation: numerical (FD) vs analytical derivative ===\n');
fprintf('freq0=%g Hz, zS0=%g m, c=%g m/s, h_freq=%g Hz, h_zS=%g m\n', freq0, zS0, c, h_freq, h_zS);

%% ── Numerical (finite-difference) Jacobian, EXACT run_delta.m recipe ──────
sim_fp = {freq0+h_freq, maxR, zS0, 0, 0, 'const_1500', 'const_2500', [1 1 0], FOM};
sim_fm = {freq0-h_freq, maxR, zS0, 0, 0, 'const_1500', 'const_2500', [1 1 0], FOM};
sim_zp = {freq0, maxR, zS0+h_zS, 0, 0, 'const_1500', 'const_2500', [1 1 0], FOM};
sim_zm = {freq0, maxR, zS0-h_zS, 0, 0, 'const_1500', 'const_2500', [1 1 0], FOM};

[TL_fp, r_grid, z_grid] = bellhopCached(sim_fp, cache_dir);
TL_fm = bellhopCached(sim_fm, cache_dir);
TL_zp = bellhopCached(sim_zp, cache_dir);
TL_zm = bellhopCached(sim_zm, cache_dir);

J_freq_num = computeJacobian(TL_fp, TL_fm, h_freq);   % dB/Hz
J_zS_num   = computeJacobian(TL_zp, TL_zm, h_zS);      % dB/m

%% ── Analytical Jacobian (exact derivative of jensenHalfspaceImage) ─────────
[R, Z] = meshgrid(r_grid, z_grid);
k0 = 2*pi*freq0/c;
R1 = sqrt(R.^2 + (Z-zS0).^2);
R2 = sqrt(R.^2 + (Z+zS0).^2);
p  = exp(1i*k0*R1)./R1 - exp(1i*k0*R2)./R2;

% dp/dfreq :  d/df[exp(ikR)/R] = i*(2*pi/c)*exp(ikR)   (the 1/R cancels the
% R from d(kR)/df = (2*pi/c)*R)
dp_dfreq = 1i*(2*pi/c) * (exp(1i*k0*R1) - exp(1i*k0*R2));

% dp/dzS :  via dR1/dzS = (zS0-Z)/R1 ,  dR2/dzS = (Z+zS0)/R2
dR1_dzS = (zS0 - Z) ./ R1;
dR2_dzS = (Z + zS0) ./ R2;
dp_dzS = exp(1i*k0*R1).*dR1_dzS.*(1i*k0./R1 - 1./R1.^2) ...
       - exp(1i*k0*R2).*dR2_dzS.*(1i*k0./R2 - 1./R2.^2);

% TL = -20*log10(|p|)  =>  dTL/dparam = -(20/ln(10)) * Re( (1/p) * dp/dparam )
J_freq_an = -(20/log(10)) * real(dp_dfreq ./ p);   % dB/Hz
J_zS_an   = -(20/log(10)) * real(dp_dzS   ./ p);   % dB/m

%% ── Difference maps + stats (exclude r=0 singular column) ─────────────────
valid_r = r_grid > 0;

diff_freq = J_freq_num - J_freq_an;
d_freq    = diff_freq(:, valid_r);
d_freq    = d_freq(isfinite(d_freq));
fprintf('dTL/dfreq [dB/Hz]: mean|diff|=%.3g, median|diff|=%.3g\n', mean(abs(d_freq)), median(abs(d_freq)));

diff_zS = J_zS_num - J_zS_an;
d_zS    = diff_zS(:, valid_r);
d_zS    = d_zS(isfinite(d_zS));
fprintf('dTL/dzS   [dB/m]:  mean|diff|=%.3g, median|diff|=%.3g\n', mean(abs(d_zS)), median(abs(d_zS)));

% Relative error (normalized by the spread of the analytical Jacobian itself,
% restricted to the same valid finite r>0 pixels used above)
Ja_freq_valid = J_freq_an(:, valid_r); Ja_freq_valid = Ja_freq_valid(isfinite(Ja_freq_valid));
Ja_zS_valid   = J_zS_an(:, valid_r);   Ja_zS_valid   = Ja_zS_valid(isfinite(Ja_zS_valid));
rel_freq = 100 * median(abs(d_freq)) / median(abs(Ja_freq_valid));
rel_zS   = 100 * median(abs(d_zS))   / median(abs(Ja_zS_valid));
fprintf('Relative median error: dTL/dfreq = %.1f%%,  dTL/dzS = %.1f%%\n', rel_freq, rel_zS);

%% ── Figure ──────────────────────────────────────────────────────────────
fig = figure('Position',[50 50 1500 800], 'Visible','off');

rows = {J_freq_an, J_freq_num, diff_freq; J_zS_an, J_zS_num, diff_zS};
row_lbl = {'dTL/dfreq (dB/Hz)', 'dTL/dzS (dB/m)'};
col_lbl = {'Analytical (exact)', 'Numerical (FD, run\_delta.m recipe)', 'Numerical - Analytical'};

for row = 1:2
    % Robust color scale (99th percentile of |analytical|, not raw max) —
    % a handful of pixels sit exactly on interference nulls where the
    % analytical derivative is singular (|p|->0 there); a max-based scale
    % is dominated by those few pixels and washes out everything else.
    an_finite = rows{row,1}(isfinite(rows{row,1}));
    cscale = prctile(abs(an_finite), 99);
    for col = 1:3
        ax = subplot(2,3,(row-1)*3+col);
        data = rows{row,col};
        imagesc(ax, r_grid/1000, z_grid, data, [-cscale cscale]);
        colormap(ax, redblue(256));
        set(ax,'YDir','reverse'); colorbar(ax);
        xlabel(ax,'Range (km)');
        if col==1, ylabel(ax, row_lbl{row}); end
        if row==1, title(ax, col_lbl{col}, 'Interpreter','none'); end
    end
end
sgtitle('Delta-method sanity check: numerical (FD) vs analytical Jacobian, half-space, isovelocity', 'Interpreter','none');

out_dir = fullfile(ROOT, 'validation', 'analytical', 'figures');
if ~exist(out_dir,'dir'), mkdir(out_dir); end
out_png = fullfile(out_dir, 'jacobian_validation.png');
print(fig, out_png, '-dpng', '-r150');
close(fig);
fprintf('Saved: %s\n', out_png);
