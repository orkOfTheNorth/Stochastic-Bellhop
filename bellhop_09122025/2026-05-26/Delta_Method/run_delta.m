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
addpath('..');   % for bellhopCached.m
if ~exist('results','dir'), mkdir('results'); end
if ~exist('figures','dir'), mkdir('figures'); end
cache_dir = fullfile('..','cache');   % shared across all methods

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
[TL_nom, r_grid, z_grid] = bellhopCached(sim_nom, cache_dir);

% Grid helpers (r_grid is a row vector in m; z_grid is a column vector in m)
r_km = r_grid / 1000;
z_m  = z_grid;

%% ── FREQUENCY PERTURBATIONS (2–3 of 7) ──────────────────────────────────────
fprintf('[2/7] freq + h ...\n');
sim = sim_nom; sim{1} = freq0 + h_freq;
[TL_fUp,~,~] = bellhopCached(sim, cache_dir);

fprintf('[3/7] freq − h ...\n');
sim{1} = freq0 - h_freq;
[TL_fDn,~,~] = bellhopCached(sim, cache_dir);

%% ── SOURCE-DEPTH PERTURBATIONS (4–5 of 7) ───────────────────────────────────
fprintf('[4/7] zS + h ...\n');
sim = sim_nom; sim{3} = zS0 + h_zS;
[TL_zUp,~,~] = bellhopCached(sim, cache_dir);

fprintf('[5/7] zS − h ...\n');
sim{3} = zS0 - h_zS;
[TL_zDn,~,~] = bellhopCached(sim, cache_dir);

%% ── TEMPERATURE PERTURBATIONS (6–7 of 7) ────────────────────────────────────
fprintf('[6/7] temp + h ...\n');
sim = sim_nom; sim{6} = "custom";
[TL_tUp,~,~] = bellhopCached(sim, cache_dir, 'CustomSVP', makeSVP(h_temp));

fprintf('[7/7] temp − h ...\n');
[TL_tDn,~,~] = bellhopCached(sim, cache_dir, 'CustomSVP', makeSVP(-h_temp));

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

    % Save (Cheb_lb saved as numeric matrix, not cell)
    TL_expected   = TL_nom;
    Var_TL        = Vs;
    Cheb_shadow95 = lb >= 0.95;
    Cheb_lb_s     = lb;   % numeric matrix for this subset only
    fname = fullfile('results', sprintf('delta_%s.mat', subset_names{s}));
    save(fname, 'TL_expected','Var_TL','Cheb_lb_s','Cheb_shadow95', ...
         'r_grid','z_grid','r_km','z_m','FOM');
    fprintf('  Saved %s\n', fname);
end

%% ── FIGURES ──────────────────────────────────────────────────────────────────
% Rules:
%  • EX maps  → one thin FOM boundary contour overlaid (smoothed to avoid fringe noise)
%  • Var maps → NO overlay lines (variance coloring speaks for itself)
%  • Shadow-category maps (saved only, no pop-up):
%      White  = non-shadow (Cheb_lb < 0.95)
%      Blue   = shadow zone at 95% Chebyshev confidence

% Smooth TL_nom for the FOM contour line only (no toolbox needed)
TL_contour = movmean(movmean(TL_nom, 20, 2), 20, 1);

cmap_tl  = jet;
cmap_var = hot;

% Shadow-category colormap: [white ; steel-blue]
cmap_shd = [1 1 1; 0.18 0.46 0.71];

% ── Fig 1: Single-variable (2×3): Row1=E[TL]+FOM line, Row2=Var[TL] ─────────
fig1 = figure('Name','Delta – Single-Variable Analysis','Position',[30 30 1500 720]);
for k = 1:3
    s   = k;
    lbl = subset_labels{s};

    subplot(2,3,k);
    pcolor(r_km, z_m, TL_nom); shading interp; set(gca,'YDir','reverse');
    colormap(gca,cmap_tl); colorbar; clim([50 150]);
    hold on;
    contour(r_km, z_m, TL_contour, [FOM FOM], 'w-', 'LineWidth',0.8);
    hold off;
    xlabel('Range (km)'); ylabel('Depth (m)');
    title(sprintf('Delta E[TL] | %s | f=%dHz zS=%.1fm  (white=FOM)',lbl,freq0,zS0));

    subplot(2,3,3+k);
    pcolor(r_km, z_m, Var_delta{s}); shading interp; set(gca,'YDir','reverse');
    colormap(gca,cmap_var); colorbar;
    xlabel('Range (km)'); ylabel('Depth (m)');
    title(sprintf('Delta Var[TL] | %s',lbl));
end
sgtitle(sprintf('Delta Method – Single-Variable UQ | f=%dHz, zS=%.1fm, FOM=%ddB', ...
    freq0,zS0,FOM));
saveFig(fig1, fullfile('figures','delta_single_variable'));

% ── Fig 2: Pair-variable (2×3): Row1=E[TL]+FOM line, Row2=Var[TL] ───────────
fig2 = figure('Name','Delta – Two-Variable Analysis','Position',[60 30 1500 720]);
for k = 1:3
    s   = 3+k;
    lbl = subset_labels{s};

    subplot(2,3,k);
    pcolor(r_km, z_m, TL_nom); shading interp; set(gca,'YDir','reverse');
    colormap(gca,cmap_tl); colorbar; clim([50 150]);
    hold on;
    contour(r_km, z_m, TL_contour, [FOM FOM], 'w-', 'LineWidth',0.8);
    hold off;
    xlabel('Range (km)'); ylabel('Depth (m)');
    title(sprintf('Delta E[TL] | %s | f=%dHz zS=%.1fm  (white=FOM)',lbl,freq0,zS0));

    subplot(2,3,3+k);
    pcolor(r_km, z_m, Var_delta{s}); shading interp; set(gca,'YDir','reverse');
    colormap(gca,cmap_var); colorbar;
    xlabel('Range (km)'); ylabel('Depth (m)');
    title(sprintf('Delta Var[TL] | %s',lbl));
end
sgtitle(sprintf('Delta Method – Two-Variable UQ | f=%dHz, zS=%.1fm, FOM=%ddB', ...
    freq0,zS0,FOM));
saveFig(fig2, fullfile('figures','delta_pair_variable'));

% ── Fig 3: Full 3-variable (1×2): E[TL]+FOM line, Var[TL] ───────────────────
s = 7;
fig3 = figure('Name','Delta – Full 3-Variable Analysis','Position',[90 30 1100 480]);

subplot(1,2,1);
pcolor(r_km, z_m, TL_nom); shading interp; set(gca,'YDir','reverse');
colormap(gca,cmap_tl); colorbar; clim([50 150]);
hold on;
contour(r_km, z_m, TL_contour, [FOM FOM], 'w-', 'LineWidth',0.8);
hold off;
xlabel('Range (km)'); ylabel('Depth (m)');
title(sprintf('Delta E[TL] | All 3 Params | f=%dHz, zS=%.1fm  (white=FOM)',freq0,zS0));

subplot(1,2,2);
pcolor(r_km, z_m, Var_delta{s}); shading interp; set(gca,'YDir','reverse');
colormap(gca,cmap_var); colorbar;
xlabel('Range (km)'); ylabel('Depth (m)');
title(sprintf('Delta Var[TL] | All 3 Params | σ²_f=%.0f σ²_z=%.5f σ²_T=%.3f', ...
              var_freq,var_zS,var_temp));

sgtitle(sprintf('Delta Method – Full 3-Parameter UQ | f=%dHz, zS=%.1fm, FOM=%ddB', ...
    freq0,zS0,FOM));
saveFig(fig3, fullfile('figures','delta_full_3param'));

% ── Fig 4: All-7 Variance summary (2×4), no overlay lines ────────────────────
fig4 = figure('Name','Delta – All Subsets Variance','Position',[120 50 1700 900]);
for s = 1:7
    subplot(2,4,s);
    pcolor(r_km, z_m, Var_delta{s}); shading interp; set(gca,'YDir','reverse');
    colormap(gca,cmap_var); colorbar;
    xlabel('Range (km)'); ylabel('Depth (m)');
    title(sprintf('Var[TL] | %s',subset_labels{s}),'FontSize',9);
end
subplot(2,4,8); axis off;
sgtitle(sprintf('Delta Method – Var[TL] All 7 Subsets | f=%dHz, zS=%.1fm, FOM=%ddB', ...
    freq0,zS0,FOM));
saveFig(fig4, fullfile('figures','delta_all_subsets_variance'));

% ── Fig 5: Shadow-category maps (2×4, saved only) ────────────────────────────
% White = non-shadow  |  Blue = Chebyshev 95% shadow
fig5 = figure('Name','Delta – Shadow Category Maps','Position',[150 50 1700 900], ...
              'Visible','off');   % save without popping up
for s = 1:7
    subplot(2,4,s);
    shd = double(Cheb_lb{s} >= 0.95);   % 0=safe, 1=shadow
    imagesc(r_km, z_m, shd); set(gca,'YDir','reverse');
    colormap(gca, cmap_shd); clim([0 1]);
    cb = colorbar; cb.Ticks = [0.25 0.75]; cb.TickLabels = {'Non-shadow','Shadow'};
    xlabel('Range (km)'); ylabel('Depth (m)');
    title(sprintf('Chebyshev Shadow (95%%) | %s',subset_labels{s}),'FontSize',9);
end
subplot(2,4,8); axis off;
sgtitle(sprintf(['Delta – Chebyshev 95%% Shadow Zone | All 7 Subsets\n' ...
    'White=safe  Blue=P(TL>%ddB)≥95%% | f=%dHz, zS=%.1fm'],FOM,freq0,zS0));
saveFig(fig5, fullfile('figures','delta_all_subsets_shadow_category'));
close(fig5);

fprintf('\n=== Delta Method complete. Results in results/  Figures in figures/ ===\n');

%% ── LOCAL FUNCTIONS ──────────────────────────────────────────────────────────
function saveFig(fig, base_path)
% Save as 300-dpi PNG + 300-dpi raster-in-PDF (fast for pcolor/surf figures).
    print(fig, base_path, '-dpng', '-r300');
    try
        exportgraphics(fig, [base_path '.pdf'], 'ContentType','image', 'Resolution',300);
    catch
        % skip PDF on older MATLAB versions
    end
end

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
