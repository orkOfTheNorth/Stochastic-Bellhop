%% run_delta.m  –  Multivariate Delta Method UQ (2026-05-26)
%
% MATHEMATICAL MODEL
% ──────────────────
% TL(θ) ≈ TL(θ₀) + J'(θ − θ₀)      (first-order Taylor expansion)
%
%   E[TL]   ≈  TL(θ₀)                 (nominal run)
%   Var[TL] ≈  J' Σ_θ J  =  Σᵢ∈S  (∂TL/∂θᵢ)² σᵢ²    (independent inputs)
%
% where J = [∂TL/∂freq, ∂TL/∂zS, ∂TL/∂temp]'  is the Jacobian at θ₀,
% computed via central finite differences:
%   ∂TL/∂θᵢ ≈ [TL(θ₀+hᵢeᵢ) – TL(θ₀−hᵢeᵢ)] / (2hᵢ)
%
% Input distributions (independent uniform U[a,b]):
%   freq  ~ U[9900, 10100] Hz    σ²_freq = (200)²/12  ≈ 3333   Hz²
%   z_S   ~ U[4.95,  5.05] m    σ²_zS   = (0.10)²/12 ≈ 8.3e−4 m²
%   temp  ~ U[−1,    +1]  °C    σ²_temp = (2)²/12    ≈ 0.333  °C²
%
% Chebyshev shadow zone (one-sided Cantelli inequality):
%   When E[TL] > FOM:
%     P(TL > FOM) ≥ max(0,  1 − Var / (Var + (E[TL]−FOM)²))
%   95% guaranteed shadow condition:   E[TL] − FOM  ≥  σ · √19 ≈ 4.36 σ
%
% Outputs: results/delta_<subset>.mat  and  figures/delta_*.png

clear; close all; clc; warning('off');

% Set CWD to script folder so Bellhop writes its files here
try, cd(fileparts(mfilename('fullpath'))); catch; end
addpath(genpath('../../uniform source depth analasis 13.5/code/Functions'));
if ~exist('results','dir'), mkdir('results'); end
if ~exist('figures','dir'), mkdir('figures'); end

%% ── PARAMETERS ──────────────────────────────────────────────────────────────
freq0 = 10000;   % Hz  (nominal)
zS0   = 5;       % m   (nominal source depth)
maxR  = 50000;   % m
FOM   = 100;     % dB  (figure of merit – shadow zone threshold)
geo   = [0.989*1500, 1.63, 0.07];

% Uniform bounds: ±1% for freq and zS, ±1°C for surface temperature
freq_bnd = [freq0*0.99, freq0*1.01];   % [9900, 10100] Hz
zS_bnd   = [zS0*0.99,  zS0*1.01];     % [4.95, 5.05] m
temp_bnd = [-1, 1];                    % °C  (surface temp shift)

% Input variances: σ² = (b−a)²/12  for U[a,b]
var_freq = diff(freq_bnd)^2 / 12;    % ≈ 3333 Hz²
var_zS   = diff(zS_bnd)^2  / 12;    % ≈ 8.33e-4 m²
var_temp = diff(temp_bnd)^2 / 12;   % ≈ 0.333 °C²

% Finite-difference step sizes (≈ 0.5% – well inside the ±1% range)
h_freq = 50;     % Hz
h_zS   = 0.025;  % m
h_temp = 0.20;   % °C

fprintf('=== DELTA METHOD  |  freq±1%%  zS±1%%  Temp±1°C ===\n');
fprintf('σ²_freq = %.2f Hz²     σ²_zS = %.6f m²     σ²_temp = %.4f °C²\n\n', ...
        var_freq, var_zS, var_temp);

%% ── SUPPRESS BELLHOP PLOT WINDOWS DURING RUNS ───────────────────────────────
set(0,'DefaultFigureVisible','off');

%% ── NOMINAL (1 of 7) ─────────────────────────────────────────────────────────
fprintf('[1/7] Nominal ...\n');
sim_nom = {freq0, maxR, zS0, 0, 0, "summer", "const_35", geo, FOM};
[TL_nom, r_grid, z_grid] = simpleBellhopHazat(sim_nom);

% Grid helpers (r_grid is a row vector in m; z_grid is a column vector in m)
r_km = r_grid / 1000;
z_m  = z_grid;

%% ── FREQUENCY PERTURBATIONS (2–3 of 7) ──────────────────────────────────────
fprintf('[2/7] freq + h ...\n');
sim = sim_nom; sim{1} = freq0 + h_freq;
[TL_fUp,~,~] = simpleBellhopHazat(sim);

fprintf('[3/7] freq − h ...\n');
sim{1} = freq0 - h_freq;
[TL_fDn,~,~] = simpleBellhopHazat(sim);

%% ── SOURCE-DEPTH PERTURBATIONS (4–5 of 7) ───────────────────────────────────
fprintf('[4/7] zS + h ...\n');
sim = sim_nom; sim{3} = zS0 + h_zS;
[TL_zUp,~,~] = simpleBellhopHazat(sim);

fprintf('[5/7] zS − h ...\n');
sim{3} = zS0 - h_zS;
[TL_zDn,~,~] = simpleBellhopHazat(sim);

%% ── TEMPERATURE PERTURBATIONS (6–7 of 7) ────────────────────────────────────
fprintf('[6/7] temp + h ...\n');
sim = sim_nom; sim{6} = "custom";
[TL_tUp,~,~] = simpleBellhopHazat(sim, 'CustomSVP', makeSVP(h_temp));

fprintf('[7/7] temp − h ...\n');
[TL_tDn,~,~] = simpleBellhopHazat(sim, 'CustomSVP', makeSVP(-h_temp));

set(0,'DefaultFigureVisible','on');
close all;

%% ── JACOBIAN MAPS ────────────────────────────────────────────────────────────
J_freq = (TL_fUp - TL_fDn) / (2*h_freq);   % dB / Hz
J_zS   = (TL_zUp - TL_zDn) / (2*h_zS);    % dB / m
J_temp = (TL_tUp - TL_tDn) / (2*h_temp);  % dB / °C

% Per-parameter variance contribution maps
V_freq = J_freq.^2 * var_freq;
V_zS   = J_zS.^2   * var_zS;
V_temp = J_temp.^2 * var_temp;

save(fullfile('results','delta_jacobians.mat'), ...
     'TL_nom','J_freq','J_zS','J_temp','V_freq','V_zS','V_temp', ...
     'r_grid','z_grid','r_km','z_m','FOM', ...
     'freq0','zS0','var_freq','var_zS','var_temp');
fprintf('Jacobians saved.\n');

%% ── 7 SUBSETS: compute Var[TL] and Chebyshev maps ───────────────────────────
%  active flags = [freq_active, zS_active, temp_active]
subset_names  = {'zS','freq','temp','zS_freq','zS_temp','freq_temp','zS_freq_temp'};
subset_labels = {'z_S only','Freq only','Temp only', ...
                 'z_S + Freq','z_S + Temp','Freq + Temp','z_S + Freq + Temp'};
active = {[0 1 0],[1 0 0],[0 0 1],[1 1 0],[0 1 1],[1 0 1],[1 1 1]};
Vpool  = {V_freq, V_zS, V_temp};   % indexed by [freq, zS, temp]

Var_delta = cell(1,7);
Cheb_lb   = cell(1,7);  % Cantelli lower bound on P(TL > FOM)

for s = 1:7
    fl = active{s};
    Vs = fl(1)*Vpool{1} + fl(2)*Vpool{2} + fl(3)*Vpool{3};
    Var_delta{s} = Vs;

    % Cantelli: P(TL > FOM) >= 1 − Var/(Var+(μ−FOM)²)  when μ > FOM
    d = TL_nom - FOM;          % positive → expected to be in shadow
    lb = zeros(size(TL_nom));
    mask = d > 0;
    lb(mask) = max(0,  1 - Vs(mask) ./ (Vs(mask) + d(mask).^2));
    Cheb_lb{s} = lb;

    % Save
    TL_expected   = TL_nom;
    Var_TL        = Vs;
    Cheb_shadow95 = lb >= 0.95;
    fname = fullfile('results', sprintf('delta_%s.mat', subset_names{s}));
    save(fname, 'TL_expected','Var_TL','Cheb_lb','Cheb_shadow95', ...
         'r_grid','z_grid','r_km','z_m','FOM');
    fprintf('  Saved %s\n', fname);
end

%% ── FIGURES ──────────────────────────────────────────────────────────────────
cmap_tl  = jet;
cmap_var = hot;
cmap_ch  = parula;

% ── Fig 1: Single-variable (2×3): Row1=Var, Row2=Chebyshev ───────────────────
fig1 = figure('Name','Delta – Single-Variable Analysis','Position',[30 30 1500 720]);
for k = 1:3
    s = k;   % subsets 1=zS, 2=freq, 3=temp
    lbl = subset_labels{s};

    subplot(2,3,k);
    pcolor(r_km, z_m, Var_delta{s}); shading interp; set(gca,'YDir','reverse');
    colormap(gca,cmap_var); colorbar;
    xlabel('Range (km)'); ylabel('Depth (m)');
    title(sprintf('\\DeltaMethod Var[TL] | %s | f=%dHz zS=%.1fm Temp=0°C', lbl,freq0,zS0));

    subplot(2,3,3+k);
    pcolor(r_km, z_m, Cheb_lb{s}*100); shading interp; set(gca,'YDir','reverse');
    colormap(gca,cmap_ch); colorbar; clim([0 100]);
    hold on;
    contour(r_km, z_m, Cheb_lb{s}, [0.95 0.95], 'w--','LineWidth',2);
    hold off;
    xlabel('Range (km)'); ylabel('Depth (m)');
    title(sprintf('Cantelli P(TL>%ddB) %% | %s | white=95%% boundary',FOM,lbl));
end
sgtitle(sprintf('Delta Method – Single-Variable UQ  |  f=%dHz, zS=%.1fm, FOM=%ddB',freq0,zS0,FOM));
saveas(fig1, fullfile('figures','delta_single_variable.png'));

% ── Fig 2: Pair-variable (2×3): Row1=Var, Row2=Chebyshev ─────────────────────
fig2 = figure('Name','Delta – Two-Variable Analysis','Position',[60 30 1500 720]);
for k = 1:3
    s = 3+k;   % subsets 4=zS+freq, 5=zS+temp, 6=freq+temp
    lbl = subset_labels{s};

    subplot(2,3,k);
    pcolor(r_km, z_m, Var_delta{s}); shading interp; set(gca,'YDir','reverse');
    colormap(gca,cmap_var); colorbar;
    xlabel('Range (km)'); ylabel('Depth (m)');
    title(sprintf('\\DeltaMethod Var[TL] | %s | f=%dHz zS=%.1fm', lbl,freq0,zS0));

    subplot(2,3,3+k);
    pcolor(r_km, z_m, Cheb_lb{s}*100); shading interp; set(gca,'YDir','reverse');
    colormap(gca,cmap_ch); colorbar; clim([0 100]);
    hold on;
    contour(r_km, z_m, Cheb_lb{s}, [0.95 0.95], 'w--','LineWidth',2);
    hold off;
    xlabel('Range (km)'); ylabel('Depth (m)');
    title(sprintf('Cantelli P(TL>%ddB) %% | %s | white=95%% boundary',FOM,lbl));
end
sgtitle(sprintf('Delta Method – Two-Variable UQ  |  f=%dHz, zS=%.1fm, FOM=%ddB',freq0,zS0,FOM));
saveas(fig2, fullfile('figures','delta_pair_variable.png'));

% ── Fig 3: Full 3-variable (1×3) ─────────────────────────────────────────────
s = 7;
fig3 = figure('Name','Delta – Full 3-Variable Analysis','Position',[90 30 1500 420]);

subplot(1,3,1);
pcolor(r_km, z_m, TL_nom); shading interp; set(gca,'YDir','reverse');
colormap(gca,cmap_tl); colorbar; clim([50 150]);
xlabel('Range (km)'); ylabel('Depth (m)');
title(sprintf('E[TL] = Nominal TL | f=%dHz, zS=%.1fm, Temp=0°C',freq0,zS0));

subplot(1,3,2);
pcolor(r_km, z_m, Var_delta{s}); shading interp; set(gca,'YDir','reverse');
colormap(gca,cmap_var); colorbar;
xlabel('Range (km)'); ylabel('Depth (m)');
title(sprintf('Var[TL] | All 3 Params | σ²_f=%.0fHz² σ²_z=%.4fm² σ²_T=%.3f°C²', ...
              var_freq,var_zS,var_temp));

subplot(1,3,3);
pcolor(r_km, z_m, Cheb_lb{s}*100); shading interp; set(gca,'YDir','reverse');
colormap(gca,cmap_ch); colorbar; clim([0 100]);
hold on;
contour(r_km, z_m, Cheb_lb{s}, [0.95 0.95], 'w--','LineWidth',2);
hold off;
xlabel('Range (km)'); ylabel('Depth (m)');
title(sprintf('Cantelli P(TL>%ddB) %% | All 3 Params | white=95%% boundary',FOM));

sgtitle(sprintf('Delta Method – Full 3-Parameter UQ  |  f=%dHz, zS=%.1fm, FOM=%ddB',freq0,zS0,FOM));
saveas(fig3, fullfile('figures','delta_full_3param.png'));

% ── Fig 4: All-7 Variance summary (2×4 layout) ───────────────────────────────
fig4 = figure('Name','Delta – All Subsets Variance','Position',[120 50 1700 900]);
for s = 1:7
    subplot(2,4,s);
    pcolor(r_km, z_m, Var_delta{s}); shading interp; set(gca,'YDir','reverse');
    colormap(gca,cmap_var); colorbar;
    xlabel('Range (km)'); ylabel('Depth (m)');
    title(sprintf('Var[TL] | %s', subset_labels{s}),'FontSize',9);
end
subplot(2,4,8); axis off;
sgtitle(sprintf('Delta Method – Variance Maps for All 7 Subsets | f=%dHz, zS=%.1fm, FOM=%ddB',freq0,zS0,FOM));
saveas(fig4, fullfile('figures','delta_all_subsets_variance.png'));

% ── Fig 5: All-7 Chebyshev shadow summary (2×4 layout) ───────────────────────
fig5 = figure('Name','Delta – All Subsets Chebyshev Shadow','Position',[150 50 1700 900]);
for s = 1:7
    subplot(2,4,s);
    pcolor(r_km, z_m, Cheb_lb{s}*100); shading interp; set(gca,'YDir','reverse');
    colormap(gca,cmap_ch); colorbar; clim([0 100]);
    hold on;
    contour(r_km, z_m, Cheb_lb{s}, [0.95 0.95], 'w--','LineWidth',1.5);
    hold off;
    xlabel('Range (km)'); ylabel('Depth (m)');
    title(sprintf('P(TL>%ddB) %% | %s',FOM,subset_labels{s}),'FontSize',9);
end
subplot(2,4,8); axis off;
sgtitle(sprintf(['Delta Method – Cantelli Lower Bound P(TL>%ddB) [%%] | All 7 Subsets\n' ...
                 'white dashed = 95%% guarantee boundary | f=%dHz, zS=%.1fm'],FOM,freq0,zS0));
saveas(fig5, fullfile('figures','delta_all_subsets_chebyshev.png'));

fprintf('\n=== Delta Method complete. Results in results/  Figures in figures/ ===\n');

%% ── LOCAL FUNCTION ───────────────────────────────────────────────────────────
function svp = makeSVP(temp_shift)
% Build a custom SVP (Wilson formula) with a uniform surface-temperature shift.
    depths  = linspace(0, 5000, 1000);
    base_td = [0 30 180 400 5000; 25 25 17 13.6 13.6]';  % [depth, temp]
    t  = interp1(base_td(:,1), base_td(:,2) + temp_shift, depths);
    S  = 37;  % salinity [ppt]
    sv = 1499.2 + 4.6*t - 0.055*t.^2 + 0.00029*t.^3 + ...
         (1.34 - 0.01*t).*(S-35) + 0.016*depths;
    svp = [depths.' sv.'];
end
