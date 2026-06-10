%% run_delta.m  –  Multivariate Delta Method UQ  (2026-05-26 | Baseline scenario)
%
% MATHEMATICAL MODEL
% ──────────────────
% TL(θ) ≈ TL(θ₀) + J'(θ − θ₀)       (first-order Taylor expansion)
%
%   E[TL]   ≈  TL(θ₀)
%   Var[TL] ≈  Σᵢ∈S  Jᵢ² σᵢ²          (independent inputs)
%
% Jacobians via central FD  →  see Functions/computeJacobian.m
% SVP noise  →  Functions/makeSVPNoise.m  (surface temp + MLD coupling)
%
% Perturbation parameters (±1%):
%   freq      ~ U[9900, 10100] Hz        σ² ≈ 3333 Hz²
%   z_S       ~ U[4.95, 5.05] m          σ² ≈ 8.33e-4 m²
%   svp_noise ~ U[-0.25, +0.25] °C       σ² ≈ 0.00521 °C²
%     (shifts surface temp AND MLD; see makeSVPNoise.m)
%
% Outputs: results/delta_<subset>.mat  +  figures/delta_*.png

clear; close all; clc; warning('off');
try, cd(fileparts(mfilename('fullpath'))); catch; end

base_dir = fullfile('..');
addpath(genpath(fullfile(base_dir, 'Functions')));
addpath(genpath(fullfile(base_dir, '../../early_stochastic_methods/code/Functions')));
addpath(base_dir);   % for bellhopCached.m

if ~exist('results','dir'), mkdir('results'); end
if ~exist('figures','dir'), mkdir('figures'); end
cache_dir = fullfile(base_dir, 'cache');

%% ── PARAMETERS ──────────────────────────────────────────────────────────────
freq0 = 10000;   zS0 = 5;   maxR = 50000;   FOM = 100;
geo   = [0.989*1500, 1.63, 0.07];

freq_bnd = [freq0*0.99, freq0*1.01];
zS_bnd   = [zS0*0.99,  zS0*1.01];
svp_bnd  = [-0.25, 0.25];   % °C  (1% of 25°C surface temperature nominal)

var_freq = diff(freq_bnd)^2 / 12;
var_zS   = diff(zS_bnd)^2  / 12;
var_svp  = diff(svp_bnd)^2 / 12;

h_freq = 50;     % Hz
h_zS   = 0.025;  % m
h_svp  = 0.125;  % °C   (half the ±0.25 bound)

fprintf('=== DELTA METHOD  |  Baseline (const_35, 50km) ===\n');
fprintf('σ²_freq=%.2f Hz²   σ²_zS=%.6f m²   σ²_svp=%.6f °C²\n\n', ...
        var_freq, var_zS, var_svp);

set(0,'DefaultFigureVisible','off');

%% ── NOMINAL RUN ──────────────────────────────────────────────────────────────
fprintf('[1/7] Nominal ...\n');
sim_nom = {freq0, maxR, zS0, 0, 0, "summer", "const_35", geo, FOM};
[TL_nom, r_grid, z_grid] = bellhopCached(sim_nom, cache_dir);
r_km = r_grid / 1000;   z_m = z_grid;
max_depth = max(z_grid);

%% ── FREQUENCY PERTURBATIONS ──────────────────────────────────────────────────
fprintf('[2/7] freq+h ...\n');
sim = sim_nom; sim{1} = freq0 + h_freq;
[TL_fUp,~,~] = bellhopCached(sim, cache_dir);

fprintf('[3/7] freq-h ...\n');
sim{1} = freq0 - h_freq;
[TL_fDn,~,~] = bellhopCached(sim, cache_dir);

%% ── SOURCE-DEPTH PERTURBATIONS ───────────────────────────────────────────────
fprintf('[4/7] zS+h ...\n');
sim = sim_nom; sim{3} = zS0 + h_zS;
[TL_zUp,~,~] = bellhopCached(sim, cache_dir);

fprintf('[5/7] zS-h ...\n');
sim{3} = zS0 - h_zS;
[TL_zDn,~,~] = bellhopCached(sim, cache_dir);

%% ── SVP NOISE PERTURBATIONS ──────────────────────────────────────────────────
fprintf('[6/7] svp+h ...\n');
sim = sim_nom; sim{6} = "custom";
[TL_sUp,~,~] = bellhopCached(sim, cache_dir, 'CustomSVP', makeSVPNoise(+h_svp, max_depth));

fprintf('[7/7] svp-h ...\n');
[TL_sDn,~,~] = bellhopCached(sim, cache_dir, 'CustomSVP', makeSVPNoise(-h_svp, max_depth));

%% ── SVP PERTURBATION INSPECTION ─────────────────────────────────────────────
fig_svp = figure('Name','Delta SVP Perturbations','Visible','off','Position',[50 50 520 640]);
svp_nom  = makeSVPNoise(0,       max_depth);
svp_plus = makeSVPNoise(+h_svp,  max_depth);
svp_neg  = makeSVPNoise(-h_svp,  max_depth);
hold on;
plot(svp_neg(:,2),  svp_neg(:,1),  'b-', 'LineWidth',1.5, 'DisplayName',sprintf('-%.3f°C',h_svp));
plot(svp_nom(:,2),  svp_nom(:,1),  'k-', 'LineWidth',2.5, 'DisplayName','Nominal (0°C)');
plot(svp_plus(:,2), svp_plus(:,1), 'r-', 'LineWidth',1.5, 'DisplayName',sprintf('+%.3f°C',h_svp));
hold off;
set(gca,'YDir','reverse'); grid on; legend('Location','southeast');
xlabel('Sound Speed (m/s)'); ylabel('Depth (m)');
title(sprintf('SVP Perturbations for FD Jacobian  h_{svp}=%.3f deg C',h_svp));
saveFig(fig_svp, fullfile('figures','svp_perturbations'));
close(fig_svp);

set(0,'DefaultFigureVisible','on');
close all;

%% ── JACOBIANS ────────────────────────────────────────────────────────────────
J_freq = computeJacobian(TL_fUp, TL_fDn, h_freq);   % dB/Hz
J_zS   = computeJacobian(TL_zUp, TL_zDn, h_zS);    % dB/m
J_svp  = computeJacobian(TL_sUp, TL_sDn, h_svp);   % dB/°C

V_freq = J_freq.^2 * var_freq;
V_zS   = J_zS.^2   * var_zS;
V_svp  = J_svp.^2  * var_svp;

save(fullfile('results','delta_jacobians.mat'), ...
     'TL_nom','J_freq','J_zS','J_svp','V_freq','V_zS','V_svp', ...
     'r_grid','z_grid','r_km','z_m','FOM', ...
     'freq0','zS0','var_freq','var_zS','var_svp', ...
     'freq_bnd','zS_bnd','svp_bnd');
fprintf('Jacobians saved.\n');

%% ── 7 SUBSETS ────────────────────────────────────────────────────────────────
% active = [freq_active, zS_active, svp_active]
subset_names  = {'zS','freq','svp','zS_freq','zS_svp','freq_svp','zS_freq_svp'};
subset_labels = {'z_S only','Freq only','SVP only', ...
                 'z_S+Freq','z_S+SVP','Freq+SVP','z_S+Freq+SVP'};
active = {[0 1 0],[1 0 0],[0 0 1],[1 1 0],[0 1 1],[1 0 1],[1 1 1]};
Vpool  = {V_freq, V_zS, V_svp};

Var_delta = cell(1,7);
Cheb_lb   = cell(1,7);

for s = 1:7
    fl = active{s};
    Vs = fl(1)*Vpool{1} + fl(2)*Vpool{2} + fl(3)*Vpool{3};
    Var_delta{s} = Vs;

    d  = TL_nom - FOM;
    lb = zeros(size(TL_nom));
    mask = d > 0;
    lb(mask) = max(0,  1 - Vs(mask) ./ (Vs(mask) + d(mask).^2));
    Cheb_lb{s} = lb;

    TL_expected   = TL_nom;
    Var_TL        = Vs;
    Cheb_shadow95 = lb >= 0.95;
    Cheb_lb_s     = lb;
    fname = fullfile('results', sprintf('delta_%s.mat', subset_names{s}));
    save(fname, 'TL_expected','Var_TL','Cheb_lb_s','Cheb_shadow95', ...
         'r_grid','z_grid','r_km','z_m','FOM');
    fprintf('  Saved %s\n', fname);
end

%% ── FIGURES ──────────────────────────────────────────────────────────────────
TL_contour = movmean(movmean(TL_nom, 20, 2), 20, 1);
cmap_tl = jet;   cmap_var = hot;
thresholds = [0.70 0.80 0.90 0.95];

% Fig 1: Single-variable (2×3): E[TL] row, Var row
fig1 = figure('Name','Delta – Single-Variable','Position',[30 30 1500 720]);
for k = 1:3
    subplot(2,3,k);
    pcolor(r_km, z_m, TL_nom); shading interp; set(gca,'YDir','reverse');
    colormap(gca,cmap_tl); colorbar; clim([50 150]);
    hold on; contour(r_km,z_m,TL_contour,[FOM FOM],'w-','LineWidth',0.8); hold off;
    xlabel('Range (km)'); ylabel('Depth (m)');
    title(sprintf('Delta E[TL] | %s',subset_labels{k}));

    subplot(2,3,3+k);
    pcolor(r_km, z_m, Var_delta{k}); shading interp; set(gca,'YDir','reverse');
    colormap(gca,cmap_var); colorbar;
    xlabel('Range (km)'); ylabel('Depth (m)');
    title(sprintf('Delta Var[TL] | %s',subset_labels{k}));
end
sgtitle(sprintf('Delta Method – Single-Variable | f=%dHz zS=%.1fm FOM=%ddB',freq0,zS0,FOM));
saveFig(fig1, fullfile('figures','delta_single_variable'));

% Fig 2: Pair-variable
fig2 = figure('Name','Delta – Two-Variable','Position',[60 30 1500 720]);
for k = 1:3
    s = 3+k;
    subplot(2,3,k);
    pcolor(r_km, z_m, TL_nom); shading interp; set(gca,'YDir','reverse');
    colormap(gca,cmap_tl); colorbar; clim([50 150]);
    hold on; contour(r_km,z_m,TL_contour,[FOM FOM],'w-','LineWidth',0.8); hold off;
    xlabel('Range (km)'); ylabel('Depth (m)');
    title(sprintf('Delta E[TL] | %s',subset_labels{s}));

    subplot(2,3,3+k);
    pcolor(r_km, z_m, Var_delta{s}); shading interp; set(gca,'YDir','reverse');
    colormap(gca,cmap_var); colorbar;
    xlabel('Range (km)'); ylabel('Depth (m)');
    title(sprintf('Delta Var[TL] | %s',subset_labels{s}));
end
sgtitle(sprintf('Delta Method – Two-Variable | f=%dHz zS=%.1fm FOM=%ddB',freq0,zS0,FOM));
saveFig(fig2, fullfile('figures','delta_pair_variable'));

% Fig 3: Full 3-param
fig3 = figure('Name','Delta – Full 3-Param','Position',[90 30 1100 480]);
subplot(1,2,1);
pcolor(r_km, z_m, TL_nom); shading interp; set(gca,'YDir','reverse');
colormap(gca,cmap_tl); colorbar; clim([50 150]);
hold on; contour(r_km,z_m,TL_contour,[FOM FOM],'w-','LineWidth',0.8); hold off;
xlabel('Range (km)'); ylabel('Depth (m)');
title(sprintf('Delta E[TL] | All 3 Params | f=%dHz zS=%.1fm',freq0,zS0));
subplot(1,2,2);
pcolor(r_km, z_m, Var_delta{7}); shading interp; set(gca,'YDir','reverse');
colormap(gca,cmap_var); colorbar;
xlabel('Range (km)'); ylabel('Depth (m)');
title(sprintf('Var[TL] | σ²f=%.0f σ²z=%.5f σ²s=%.5f',var_freq,var_zS,var_svp));
sgtitle(sprintf('Delta – Full 3-Param | f=%dHz zS=%.1fm FOM=%ddB',freq0,zS0,FOM));
saveFig(fig3, fullfile('figures','delta_full_3param'));

% Fig 4: All-7 Mean TL (E[TL] = TL_nom for all subsets in Delta method)
N_s = 7;  fw = max(1400, 200*N_s);
fig4 = figure('Name','Delta – Mean TL all subsets','Visible','off','Position',[50 50 fw 380]);
for s = 1:N_s
    ax = subplot(1,N_s,s);
    pcolor(r_km,z_m,TL_nom); shading interp; set(ax,'YDir','reverse');
    colormap(ax,cmap_tl); colorbar; clim([50 150]);
    hold on; contour(r_km,z_m,TL_contour,[FOM FOM],'w-','LineWidth',0.8); hold off;
    xlabel('Range (km)'); if s==1, ylabel('Depth (m)'); end
    title(subset_labels{s},'FontSize',8);
end
sgtitle('Delta E[TL]  (= TL_{nom} for all subsets)');
saveFig(fig4, fullfile('figures','delta_all_subsets_mean_TL'));  close(fig4);

% Fig 5: All-7 Variance (1×7 row)
fig5 = figure('Name','Delta – Var all subsets','Visible','off','Position',[50 50 fw 380]);
for s = 1:N_s
    ax = subplot(1,N_s,s);
    pcolor(r_km,z_m,Var_delta{s}); shading interp; set(ax,'YDir','reverse');
    colormap(ax,cmap_var); colorbar;
    xlabel('Range (km)'); if s==1, ylabel('Depth (m)'); end
    title(subset_labels{s},'FontSize',8);
end
sgtitle('Delta Var[TL] — All 7 Subsets');
saveFig(fig5, fullfile('figures','delta_all_subsets_variance'));  close(fig5);

% Fig 6: All-7 Chebyshev P(shadow) (1×7 row)
fig6 = figure('Name','Delta – Cheb shadow all subsets','Visible','off','Position',[50 50 fw 380]);
for s = 1:N_s
    ax = subplot(1,N_s,s);
    shadowThresholdMaps(ax,r_km,z_m,Cheb_lb{s},subset_labels{s},thresholds);
end
sgtitle(sprintf('Delta Chebyshev P(shadow)  FOM=%ddB  Thresholds: 70/80/90/95%%',FOM));
saveFig(fig6, fullfile('figures','delta_all_subsets_cheb'));  close(fig6);

fprintf('\n=== Delta complete. Results in results/  Figures in figures/ ===\n');
