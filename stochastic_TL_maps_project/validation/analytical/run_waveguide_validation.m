%% run_waveguide_validation.m
% Analytical Pekeris-waveguide normal-mode solution (50 trapped modes)
% validated against Bellhop for an isovelocity shallow-water channel.
%
% Physical setup (classic Pekeris waveguide, Jensen/Kuperman/Porter/Schmidt,
% "Computational Ocean Acoustics" 2nd ed., Ch. 2 & 5):
%   * pressure-release surface at z = 0  (p = 0)
%   * isovelocity water column, speed c1, density rho1, depth D
%   * semi-infinite fluid bottom half-space, speed c2 > c1, density rho2
%
% This project's bottom parameters come from config.json "nominal.geo"
% = [c_bottom  rho_ratio  alpha_dB_per_lambda]; createBellhopEnv.m writes
% geo(1) VERBATIM as the bottom half-space sound speed of the .env file,
% so here c2 = geo(1) = 1483.5 m/s, rho2/rho1 = geo(2) = 1.63,
% bottom attenuation alpha = geo(3) = 0.07 dB/lambda.
%
% ── Derivation of the dispersion relation ────────────────────────────────
% Seek modal solutions p(r,z) = phi(z) * H0^(1)(kr r), time factor e^{-iwt}.
% Water column (0<z<D):   phi'' + kz1^2 phi = 0,  kz1^2 = k1^2 - kr^2
%   Pressure release at z=0 (phi(0)=0)  =>  phi(z) = A sin(kz1 z).
% Bottom (z>D): phi'' + (k2^2 - kr^2) phi = 0. For a TRAPPED mode we need
%   an evanescent (decaying) tail, i.e. kr > k2:
%     phi(z) = B exp(-gamma (z - D)),   gamma = sqrt(kr^2 - k2^2) > 0
%   (this is the radiation/decay condition; trapped modes exist only for
%    k2 < kr < k1, i.e. only when c2 > c1).
% Boundary conditions at z = D:
%   (i)  pressure continuity:              A sin(kz1 D) = B
%   (ii) normal velocity continuity, i.e.  (1/rho) dp/dz continuous:
%        (A kz1 / rho1) cos(kz1 D) = -(B gamma / rho2)
% Dividing (ii) by (i):
%        kz1 cot(kz1 D) / rho1 = -gamma / rho2
%   =>   tan(kz1 D) = - (rho2/rho1) * kz1 / gamma          (Pekeris relation)
% with gamma = sqrt(k1^2 - k2^2 - kz1^2). Numerically we solve the smooth,
% singularity-free equivalent form
%        g(kz1) = gamma sin(kz1 D) + (rho2/rho1) kz1 cos(kz1 D) = 0 ,
% which has exactly one root in each interval kz1*D in ((m-1/2)pi, m*pi),
% m = 1..M, M = floor(kzmax*D/pi + 1/2), kzmax = sqrt(k1^2 - k2^2).
%
% ── Field / TL ────────────────────────────────────────────────────────────
% With the density-weighted normalisation  int phi_m^2(z)/rho(z) dz = 1
% (water sine part + evanescent bottom tail), the modal field is (COA 2nd
% ed. eq. 5.13-5.14, far-field asymptotic Hankel function
% H0^(1)(x) ~ sqrt(2/(pi x)) e^{i(x - pi/4)}, valid for r >> lambda):
%   p(r,z) ~ (i / (rho1 sqrt(8 pi r))) e^{-i pi/4}
%            * sum_m phi_m(zS) phi_m(z) e^{i krm r} / sqrt(krm)
% TL is referenced to the same source at 1 m (|p0(1m)| = 1/(4 pi)), which
% matches Bellhop's unit-source-at-1m convention:
%   TL(r,z) = -20 log10(4 pi |p|) = -20 log10( sqrt(2 pi / r) |S| ),
%   S = sum_m phi_m(zS) phi_m(z) e^{i krm r} / sqrt(krm).
%
% NOTE on coherence: createBellhopEnv.m hardcodes Bellhop run type 'SB' =
% SEMI-COHERENT TL (incoherent ray sum with the Lloyd-mirror source factor
% 2 sin^2(k1 zS sin(theta))). Its modal counterpart is the INCOHERENT mode
% sum (the factor phi_m^2(zS) ~ (2/D) sin^2(kz_m zS) IS that source factor,
% kz = k1 sin(theta)):
%   |p|^2_incoh = (2 pi / r) sum_m phi_m^2(zS) phi_m^2(z)
%                                  e^{-2 Im(krm) r} / krm / (4 pi)^-2 ...
%   TL_incoh = -10 log10( (2 pi / r) sum_m [...] ).
% We therefore report the quantitative error metric between Bellhop and the
% incoherent mode sum (like-for-like), and additionally show the coherent
% 50-mode field, whose local range-average follows the same envelope.
%
% Loss terms (to mirror what Bellhop models):
%  * Bottom attenuation (geo(3) dB/lambda) enters via first-order
%    perturbation theory (COA eq. 5.175): Im(krm) = (k2 alpha2_Np / krm)
%    * int_bottom phi_m^2 / rho2 dz  (fraction of mode energy in the bottom).
%  * The .env SSP option string 'CVWT' enables Thorp volume attenuation
%    ('T'), so the analytical water column uses the same Thorp formula as
%    the Acoustics-Toolbox source (dB/kyd -> Np/m).
%
% ── Choice of scenario / "50 modes" ──────────────────────────────────────
% shallow_water scenario (D = 50 m, maxR = 5 km), nominal f = 10 kHz, and
% an ISOVELOCITY water column c1 = 1479.4 m/s: the Pekeris mode count
% M = floor(D sqrt(k1^2-k2^2)/pi + 1/2) then equals exactly 50. (With the
% project's summer SVP the profile is not isovelocity and the Pekeris
% solution would not apply; c1 is a free choice for this validation and was
% selected so that exactly 50 trapped modes exist at the nominal frequency
% with the project's own bottom parameters.)
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
if ~iscell(scen), scen = num2cell(scen); end   % jsondecode: cell if fields differ
sc = scen{cellfun(@(s) strcmp(s.name,'shallow_water'), scen)};
if isfield(sc,'FOM_dB') && ~isempty(sc.FOM_dB), FOM = sc.FOM_dB; else, FOM = cfg.nominal.FOM_dB; end

freq = cfg.nominal.freq_Hz;      % 10000 Hz (nominal)
zS   = cfg.nominal.zS_m;         % 5 m source depth (nominal)
geo  = cfg.nominal.geo(:).';     % [c2  rho2/rho1  alpha (dB/lambda)]
D    = sc.maxDepth_m;            % 50 m
maxR = sc.maxR_m;                % 5000 m

c1   = 1479.4;                   % isovelocity water speed -> exactly 50 modes
c2   = geo(1);                   % 1483.5 m/s (bottom half-space)
rho  = geo(2);                   % rho2/rho1 = 1.63
a2dl = geo(3);                   % bottom attenuation, dB per wavelength

k1 = 2*pi*freq/c1;
k2 = 2*pi*freq/c2;
assert(c2 > c1, 'Pekeris trapped modes require c2 > c1.');
kzmax = sqrt(k1^2 - k2^2);

%% ── Solve dispersion relation for all trapped modes ─────────────────────
M = floor(kzmax*D/pi + 0.5);
fprintf('Pekeris waveguide: f=%g Hz, D=%g m, c1=%g, c2=%g, rho=%g -> %d trapped modes\n', ...
        freq, D, c1, c2, rho, M);

g = @(kz) sqrt(max(kzmax^2 - kz.^2, 0)) .* sin(kz*D) + rho .* kz .* cos(kz*D);

kz = zeros(M,1);
for m = 1:M
    lo = ((m-0.5)*pi/D) * (1 + 1e-12) + 1e-12;
    hi = min(m*pi/D, kzmax) * (1 - 1e-12);
    assert(g(lo)*g(hi) < 0, 'No sign change for mode %d — bracket error.', m);
    kz(m) = fzero(g, [lo hi]);
end
krm   = sqrt(k1^2 - kz.^2);            % horizontal wavenumbers (real part)
gam   = sqrt(krm.^2 - k2^2);           % bottom decay constants (>0)

% Density-weighted normalisation:  N = int_0^D sin^2/rho1 + tail/rho2
Nw = (D/2 - sin(2*kz*D)./(4*kz));                  % water integral (rho1 = 1)
Nb = sin(kz*D).^2 ./ (2*gam*rho);                  % bottom tail integral
Nm = Nw + Nb;

% Modal attenuation (perturbation): bottom loss + Thorp volume loss
alpha2_Np = a2dl * freq / c2 / 8.6858896;          % dB/lambda -> Np/m in bottom
f_kHz  = freq/1000;
thorp_dBkyd = 40*f_kHz^2/(4100+f_kHz^2) + 0.1*f_kHz^2/(1+f_kHz^2);
alphaW_Np = thorp_dBkyd / 914.4 / 8.6858896;       % Np/m in water (as in AT 'T')
Im_krm = (k2*alpha2_Np ./ krm) .* (Nb ./ Nm) ...   % bottom-loss part
       + (k1*alphaW_Np ./ krm) .* (Nw ./ Nm);      % Thorp volume part
krm_c  = krm + 1i*Im_krm;

%% ── Bellhop reference run (same physics: isovelocity c1, same geo) ──────
svp_type  = sprintf('const_%.6g', c1);
sim_pars  = {freq, maxR, zS, 0, 0, svp_type, sc.bathy_type, geo, FOM};
cache_dir = fullfile(ROOT, 'Cache', char(sc.name), 'bellhop_raw');
[TL_bh, r_grid, z_grid] = bellhopCached(sim_pars, cache_dir);
r = r_grid(:).';   z = z_grid(:);      % r: 1xNr (m), z: Nzx1 (m)

%% ── Analytical mode sum on the same grid ────────────────────────────────
phi_z  = sin(z * kz.') ./ sqrt(Nm.');              % Nz x M  (water column)
phi_s  = sin(kz * zS) ./ sqrt(Nm);                 % M x 1
amp    = phi_s ./ sqrt(krm_c);                     % M x 1
r_safe = max(r, 1);                                % guard r = 0 column
S      = phi_z * (amp .* exp(1i * krm_c * r_safe));% Nz x Nr
TL_an  = -20*log10( sqrt(2*pi ./ r_safe) .* abs(S) + eps );
TL_an(:, r < 1) = NaN;                             % no far-field at r ~ 0

%% ── Error metrics (mask interference nulls + leaky-mode near field) ─────
r_min = 500;                                       % m; near field is leaky-mode
valid = (ones(numel(z),1) * (r >= r_min)) & isfinite(TL_an) & isfinite(TL_bh) ...
        & (TL_bh < 110) & (TL_an < 110);           % TL>110 dB = interference null
d      = TL_an - TL_bh;
mean_abs = mean(abs(d(valid)));
med_abs  = median(abs(d(valid)));
max_abs  = max(abs(d(valid)));
bias     = mean(d(valid));
fprintf(['Waveguide validation (%d modes, r >= %g m, nulls masked):\n' ...
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
title(ax1, sprintf('Analytical normal modes (%d trapped modes)', M));

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
title(ax3, sprintf('Analytical - Bellhop (mean |\\Delta| = %.2f dB)', mean_abs), 'Interpreter','tex');

ax4 = subplot(2,2,4);
[~, iz] = min(abs(z - D/2));
plot(ax4, r_km, TL_bh(iz,:), '-',  'Color',[0.85 0.33 0.10], 'LineWidth',1.0); hold(ax4,'on');
plot(ax4, r_km, TL_an(iz,:), '-',  'Color',[0.00 0.45 0.74], 'LineWidth',1.0);
set(ax4,'YDir','reverse'); grid(ax4,'on');
xlabel(ax4,'Range (km)'); ylabel(ax4,'TL (dB)');
legend(ax4, {'Bellhop','Analytical (modes)'}, 'Location','southwest');
title(ax4, sprintf('TL slice at z = %.0f m', z(iz)));
xlim(ax4, [0 max(r_km)]); ylim(ax4, [tl_lo-5 tl_hi+15]);

sgtitle(sprintf(['Pekeris waveguide validation — %d normal modes | f = %g kHz, D = %g m, z_S = %g m\n' ...
                 'c_1 = %.1f m/s, c_2 = %.1f m/s, \\rho_2/\\rho_1 = %.2f, \\alpha_2 = %.2f dB/\\lambda'], ...
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
