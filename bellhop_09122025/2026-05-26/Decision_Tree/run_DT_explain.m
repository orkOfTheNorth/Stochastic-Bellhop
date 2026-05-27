%% run_DT_explain.m  —  RF Explainability Analysis  (2026-05-26)
%
% Explains WHERE and HOW MUCH each physical parameter drives the RF shadow
% predictions using two complementary techniques:
%
%   1. Partial Dependence Plots (PDPs)
%      For each input feature, sweep its value across a wide physical range
%      while holding others fixed at nominal.  The RF's average predicted
%      P(TL>FOM) as a function of that feature reveals its global effect.
%      Grey band shows the training range (±1%); lines outside show how the
%      RF extrapolates.
%
%   2. Spatial marginal influence maps
%      At every grid point (r, z), compute:
%         influence(r,z) = P(shadow | param=HIGH) − P(shadow | param=LOW)
%      Positive (red)  = higher parameter value → more shadow there.
%      Negative (blue) = higher parameter value → less shadow there.
%      This is the data-driven equivalent of the Delta method Jacobian map.
%
% REQUIRES
%   results/RF_results.mat  (run run_DT.m first)
%
% OUTPUT
%   figures/RF_PDP.png            — PDPs for all 5 input features
%   figures/RF_influence_maps.png — spatial influence maps (freq, zS, temp)
%   results/RF_explain.mat

clear; close all; clc;
try, cd(fileparts(mfilename('fullpath'))); catch; end

if ~exist('results','dir'), mkdir('results'); end
if ~exist('figures','dir'), mkdir('figures'); end

fprintf('=== RF Explainability Analysis ===\n\n');

%% ── Load RF and reconstruct grid ─────────────────────────────────────────
rf_file = fullfile('results','RF_results.mat');
if ~isfile(rf_file)
    error('RF model not found — run run_DT.m first.');
end
res      = load(rf_file, 'RF', 'r_km_map', 'z_m_map');
RF       = res.RF;
r_km_map = res.r_km_map(:)';   % 1 × nri  (km)
z_m_map  = res.z_m_map(:)';    % 1 × nzi  (m)
nri      = numel(r_km_map);
nzi      = numel(z_m_map);
n_pts    = nzi * nri;

[RI, ZI] = meshgrid(1:nri, 1:nzi);
r_vec = r_km_map(RI(:))';   % n_pts × 1
z_vec = z_m_map(ZI(:))';    % n_pts × 1

freq0 = 10000;   zS0 = 5;   temp0 = 0;
FOM   = 100;

fprintf('RF: %d trees  |  Grid: %d depth × %d range = %d points\n', ...
    RF.NumTrees, nzi, nri, n_pts);

%% ── Build nominal feature matrix (all spatial pts, nominal params) ────────
X_nom = [repmat(freq0, n_pts, 1), ...
         repmat(zS0,   n_pts, 1), ...
         repmat(temp0, n_pts, 1), ...
         r_vec, z_vec];

%% ══════════════════════════════════════════════════════════════════════════
%% FIGURE 1 — Partial Dependence Plots
%  For each feature: fix all other features at nominal, sweep the target
%  feature across a wide range, and record the average predicted P(shadow).
%% ══════════════════════════════════════════════════════════════════════════
fprintf('\nComputing PDPs (sweeping each feature over wide range)...\n');

% {label, feature index, sweep values, training lo, training hi}
sweeps = { ...
    'Frequency (Hz)',       1, linspace(4000,  16000, 100)', 9900,  10100; ...
    'Source depth zS (m)',  2, linspace(1,     33,    100)', 4.95,  5.05;  ...
    'Temp shift (C)',       3, linspace(-6,    6,     100)', -1,    1;     ...
    'Range (km)',           4, linspace(0.2,   50,    100)', 0.1,   50;    ...
    'Depth (m)',            5, linspace(0.5,   34,    100)', 0,     35;    ...
};
nFeat = size(sweeps,1);
pdp_vals = cell(nFeat,1);

% Background: 800 random spatial grid points at nominal physical params
rng(7);
idx_bg = randperm(n_pts, min(800, n_pts));
X_bg   = [repmat(freq0, numel(idx_bg), 1), ...
          repmat(zS0,   numel(idx_bg), 1), ...
          repmat(temp0, numel(idx_bg), 1), ...
          r_vec(idx_bg), z_vec(idx_bg)];

for k = 1:nFeat
    feat_k = sweeps{k,2};
    vals_k = sweeps{k,3};
    pdp_k  = zeros(numel(vals_k), 1);
    for vi = 1:numel(vals_k)
        X_v = X_bg;
        X_v(:, feat_k) = vals_k(vi);
        [~, sc] = predict(RF, X_v);
        pdp_k(vi) = mean(sc(:,2));
    end
    pdp_vals{k} = pdp_k;
    fprintf('  [%d/5] %s done\n', k, sweeps{k,1});
end

fig1 = figure('Name','RF PDP','Position',[30 30 1500 440]);
for k = 1:nFeat
    ax = subplot(1, nFeat, k);
    hold on;

    % Grey fill: training range
    tr_lo = sweeps{k,4};  tr_hi = sweeps{k,5};
    fill([tr_lo tr_hi tr_hi tr_lo], [0 0 1 1], ...
        [0.82 0.82 0.82], 'FaceAlpha', 0.45, 'EdgeColor','none');

    plot(sweeps{k,3}, pdp_vals{k}, 'b-', 'LineWidth', 2);
    yline(0.5, 'k--', 'LineWidth', 0.8);
    hold off;

    xlabel(sweeps{k,1}, 'FontSize', 8);
    if k==1, ylabel('E[P(TL>FOM)]', 'FontSize', 9); end
    title(sweeps{k,1}, 'FontSize', 8, 'Interpreter','none');
    ylim([0 1]); grid on;
    text(0.05, 0.93, 'Training range', 'Units','norm', ...
        'FontSize', 6, 'Color',[0.4 0.4 0.4]);
end
sgtitle({'Partial Dependence Plots — RF average predicted P(TL>FOM) vs each input feature', ...
    'Grey = training range  |  Blue line = marginal effect  |  Dashed = 50% threshold'}, ...
    'FontSize', 10);
saveFig(fig1, fullfile('figures','RF_PDP'));
close(fig1);
fprintf('Saved RF_PDP\n');

%% ══════════════════════════════════════════════════════════════════════════
%% FIGURE 2 — Spatial Marginal Influence Maps
%  For each physical parameter:  influence(r,z) = P(shadow|high) − P(shadow|low)
%  This shows WHERE in the field each parameter matters most.
%  Analogous to the Delta method Jacobian but for shadow probability.
%% ══════════════════════════════════════════════════════════════════════════
fprintf('\nComputing spatial influence maps...\n');

% {label, feature index, low value, high value}
phys = { ...
    'Freq:  5kHz vs 15kHz',   1,  5000,  15000; ...
    'z_S:   2m   vs 25m',     2,  2,     25;    ...
    'Temp: -5C   vs +5C',     3, -5,     5;     ...
};
nP = size(phys,1);
infl_maps = zeros(nzi, nri, nP);

for p = 1:nP
    feat_p = phys{p,2};
    val_lo = phys{p,3};
    val_hi = phys{p,4};
    X_lo = X_nom;  X_lo(:, feat_p) = val_lo;
    X_hi = X_nom;  X_hi(:, feat_p) = val_hi;
    [~, sc_lo] = predict(RF, X_lo);
    [~, sc_hi] = predict(RF, X_hi);
    infl_maps(:,:,p) = reshape(sc_hi(:,2) - sc_lo(:,2), nzi, nri);
    fprintf('  [%d/3] %s done\n', p, phys{p,1});
end

fig2 = figure('Name','RF Influence Maps','Position',[40 40 1550 500]);
cmax_all = max(abs(infl_maps(:)));
for p = 1:nP
    ax = subplot(1, nP, p);
    pcolor(r_km_map, z_m_map, infl_maps(:,:,p)); shading interp;
    set(ax,'YDir','reverse');
    colormap(ax, redblue(256));
    clim([-cmax_all  cmax_all]);
    cb = colorbar; cb.Label.String = 'dP(shadow)'; cb.FontSize = 8;
    xlabel('Range (km)', 'FontSize', 8);
    if p==1, ylabel('Depth (m)', 'FontSize', 8); end
    title({phys{p,1}, '[high] - [low]'}, 'FontSize', 8, 'Interpreter','none');
end
sgtitle({'RF Spatial Influence Maps — dP(shadow) when switching each parameter from low to high', ...
    'Red = high param causes MORE shadow   Blue = high param causes LESS shadow'}, 'FontSize', 10);
saveFig(fig2, fullfile('figures','RF_influence_maps'));
close(fig2);
fprintf('Saved RF_influence_maps\n');

%% ── Save ─────────────────────────────────────────────────────────────────
save(fullfile('results','RF_explain.mat'), ...
     'pdp_vals','infl_maps','sweeps','phys','r_km_map','z_m_map', '-v7.3');

fprintf('\n=== Explainability complete. Figures in figures/ ===\n');

%% ── Local functions ───────────────────────────────────────────────────────
function cmap = redblue(n)
    n2=ceil(n/2); n1=n-n2;
    cmap=[linspace(0,1,n1)',linspace(0,1,n1)',ones(n1,1);
          ones(n2,1),linspace(1,0,n2)',linspace(1,0,n2)'];
end

function saveFig(fig, base_path)
    print(fig, base_path, '-dpng', '-r200');
    try
        exportgraphics(fig,[base_path '.pdf'],'ContentType','image','Resolution',200);
    catch
    end
end

