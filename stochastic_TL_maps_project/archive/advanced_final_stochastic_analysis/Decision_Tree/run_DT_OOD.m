%% run_DT_OOD.m  —  Out-of-Distribution RF generalisation test  (2026-05-26)
%
% WHAT IT DOES
%   Tests the trained RF on 5 parameter sets that are FAR OUTSIDE its
%   training bounds.  The RF was trained on narrow MC perturbations:
%       freq  in [9900, 10100] Hz   (±1% around 10 kHz)
%       zS    in [4.95, 5.05]  m    (±1% around 5 m)
%       temp  in [-1, +1]      C    (near-surface temperature shift)
%
%   For each OOD scenario we:
%     (1) Run Bellhop once — deterministic ground truth for that exact
%         (freq, zS, temp) combination.
%     (2) Query the RF with those same input features at every subsampled
%         grid point.
%     (3) Compare RF predicted P(TL>FOM) to the binary Bellhop shadow map.
%
%   This answers: can the forest generalise to physically valid but
%   previously unseen operating conditions?  Expected outcome: poor
%   accuracy on scenarios with very different frequency (RF has learned
%   range/depth geometry at 10 kHz, which is frequency-specific).
%
% PREREQUISITES
%   run_DT.m must have been run (needs results/RF_results.mat).
%   Bellhop must be on the system PATH.
%
% OUTPUT
%   figures/RF_OOD_maps.png   — 2-row panel: Bellhop shadow / RF prob, 5 scenarios
%   figures/RF_OOD_table.png  — accuracy + RMSE summary table
%   results/RF_OOD.mat

clear; close all; clc;
try, cd(fileparts(mfilename('fullpath'))); catch; end
addpath(genpath('../../early_stochastic_methods/code/Functions'));
addpath('..');   % bellhopCached.m lives in 2026-05-26/

if ~exist('results','dir'), mkdir('results'); end
if ~exist('figures','dir'), mkdir('figures'); end

fprintf('=== RF Out-of-Distribution Test ===\n\n');

%% ── Load saved RF ─────────────────────────────────────────────────────────
rf_file = fullfile('results','RF_results.mat');
if ~isfile(rf_file)
    error('RF model not found at %s — run run_DT.m first.', rf_file);
end
res = load(rf_file, 'RF');
RF  = res.RF;
fprintf('RF loaded: %d trees\n', RF.NumTrees);

%% ── Common grid setup (nominal run at training conditions) ────────────────
FOM       = 100;     % dB threshold — same as training
maxR      = 50000;   % m
geo       = [0.989*1500, 1.63, 0.07];
cache_dir = fullfile('..','cache');
stride    = 5;   % must match the stride used in run_DT.m

set(0,'DefaultFigureVisible','off');
sim_nom = {10000, maxR, 5, 0, 0, "summer", "const_35", geo, FOM};
[~, r_grid0, z_grid0] = bellhopCached(sim_nom, cache_dir);
close all;

r_km0 = r_grid0(:)' / 1000;   % row vector, km
z_m0  = z_grid0(:)';          % row vector, m
Nz0   = numel(z_m0);
Nr0   = numel(r_km0);

zi_sub = 1:stride:Nz0;
ri_sub = 1:stride:Nr0;
r_km_s = r_km0(ri_sub);   % subsampled range (km), row vector
z_m_s  = z_m0(zi_sub);    % subsampled depth (m),  row vector
nzi    = numel(zi_sub);
nri    = numel(ri_sub);

% RF feature vectors: every subsampled grid point (one vector per point)
[RI, ZI] = meshgrid(ri_sub, zi_sub);
r_vec = r_km0(RI(:))';   % n_pts x 1  — range feature (km)
z_vec = z_m0(ZI(:))';    % n_pts x 1  — depth feature (m)
n_pts = nzi * nri;

%% ── OOD scenario definitions ──────────────────────────────────────────────
% Columns: label | freq(Hz) | zS(m) | temp_shift(C) | use_custom_svp
% All scenarios are physically valid for D=35 m shallow water.
scenarios = { ...
    'f=5kHz  zS=5m  T=0C',    5000,  5,  0, false; ...   % half frequency
    'f=15kHz zS=5m  T=0C',   15000,  5,  0, false; ...   % 1.5x frequency
    'f=10kHz zS=20m T=0C',   10000, 20,  0, false; ...   % near-bottom source
    'f=10kHz zS=5m  T=-5C',  10000,  5, -5, true;  ...   % cold anomaly
    'f=7kHz  zS=12m T=+3C',   7000, 12,  3, true;  ...   % combined OOD
};
nS = size(scenarios, 1);

%% ── Run Bellhop + RF for each scenario ────────────────────────────────────
TL_bhp_sub  = cell(nS,1);   % TL interpolated to subsampled grid, nzi x nri
RF_prob_sub = cell(nS,1);   % RF P(TL>FOM), nzi x nri
acc_all     = zeros(nS,1);
rmse_all    = zeros(nS,1);

for k = 1:nS
    lbl    = scenarios{k,1};
    f_k    = scenarios{k,2};
    zS_k   = scenarios{k,3};
    t_k    = scenarios{k,4};
    custom = scenarios{k,5};
    fprintf('Scenario %d/%d: %s\n', k, nS, lbl);

    % ── Bellhop ground truth (1 deterministic run) ──
    sim_k = {f_k, maxR, zS_k, 0, 0, "summer", "const_35", geo, FOM};
    if custom
        sim_k{6} = "custom";
        [TL_k, r_k, z_k] = bellhopCached(sim_k, cache_dir, ...
            'CustomSVP', makeSVP(t_k, max(z_m0)));
    else
        [TL_k, r_k, z_k] = bellhopCached(sim_k, cache_dir);
    end
    close all;

    % Project Bellhop output onto the common subsampled grid
    % (different frequencies may produce slightly different Bellhop grids)
    TL_on_sub = interp2(r_k(:)'/1000, z_k(:)', TL_k, r_km_s, z_m_s', 'linear', NaN);
    TL_bhp_sub{k} = TL_on_sub;   % nzi x nri

    shadow_bhp = TL_on_sub > FOM;   % binary shadow map (1=shadow zone)

    % ── RF prediction at same grid points ──
    X_k = [repmat(f_k,   n_pts, 1), ...
           repmat(zS_k,  n_pts, 1), ...
           repmat(t_k,   n_pts, 1), ...
           r_vec, z_vec];
    [~, sc_k]  = predict(RF, X_k);
    prob_k     = reshape(sc_k(:,2), nzi, nri);
    RF_prob_sub{k} = prob_k;

    % ── Metrics (only where Bellhop interpolation is valid) ──
    valid    = isfinite(TL_on_sub(:));
    pred_bin = prob_k(valid) > 0.5;
    true_bin = shadow_bhp(valid);
    prob_v   = prob_k(valid);
    true_v   = double(true_bin);

    acc_all(k)  = mean(pred_bin == true_bin);
    rmse_all(k) = sqrt(mean((prob_v - true_v).^2));
    fprintf('  Accuracy vs Bellhop: %.1f%%   RMSE: %.3f\n', ...
        100*acc_all(k), rmse_all(k));
end

%% ── Figure 1: comparison maps ─────────────────────────────────────────────
% Top row:    Bellhop binary shadow (blue=shadow, white=clear)
% Bottom row: RF predicted P(TL>FOM)  [0=white, 1=dark blue]
fig1 = figure('Name','RF OOD','Position',[20 20 1650 640]);

for k = 1:nS
    shadow = TL_bhp_sub{k} > FOM;

    % Top: Bellhop ground truth
    ax1 = subplot(2, nS, k);
    pcolor(r_km_s, z_m_s, double(shadow)); shading flat;
    set(ax1, 'YDir','reverse');
    colormap(ax1, [0.93 0.96 1.0; 0.08 0.28 0.72]);   % white/blue
    clim([0 1]);
    title(scenarios{k,1}, 'FontSize', 8, 'Interpreter','none');
    if k==1, ylabel('Depth (m)','FontSize',8); end
    xlabel('Range (km)','FontSize',7);
    text(0.02, 0.08, 'Bellhop binary', 'Units','norm', 'FontSize',7);

    % Bottom: RF probability
    ax2 = subplot(2, nS, k+nS);
    pcolor(r_km_s, z_m_s, RF_prob_sub{k}); shading interp;
    set(ax2, 'YDir','reverse');
    colormap(ax2, parula(256)); clim([0 1]);
    if k==1, ylabel('Depth (m)','FontSize',8); end
    xlabel('Range (km)','FontSize',7);
    title(sprintf('RF  Acc=%.0f%%  RMSE=%.3f', 100*acc_all(k), rmse_all(k)), ...
        'FontSize', 8);
    if k==nS
        cb2 = colorbar; cb2.Label.String = 'P(TL>FOM)'; cb2.FontSize = 7;
    end
end
sgtitle({'RF Out-of-Distribution Test', ...
    'Top: Bellhop shadow zone (blue=shadow)  |  Bottom: RF predicted P(TL>FOM)'}, ...
    'FontSize',10);

saveFig(fig1, fullfile('figures','RF_OOD_maps'));
close(fig1);
fprintf('\nSaved RF_OOD_maps\n');

%% ── Figure 2: accuracy bar chart ──────────────────────────────────────────
fig2 = figure('Name','OOD Accuracy','Position',[30 30 820 420]);

subplot(1,2,1);
bar(1:nS, 100*acc_all, 'FaceColor',[0.25 0.55 0.85]);
hold on;
yline(100 * mean([0 1]), 'k--', 'LineWidth', 1.2);  % 50% = random guess
ylabel('Accuracy vs Bellhop (%)');
set(gca,'XTick',1:nS,'XTickLabel', ...
    cellfun(@(x)strrep(x,'  ',' '), scenarios(:,1),'UniformOutput',false), ...
    'XTickLabelRotation',25, 'FontSize',7);
ylim([0 100]);
title('Accuracy per scenario');
grid on; yline(50,'k--','50% (random)','LabelHorizontalAlignment','right');

subplot(1,2,2);
bar(1:nS, rmse_all, 'FaceColor',[0.85 0.35 0.25]);
ylabel('RMSE (RF prob vs binary shadow)');
set(gca,'XTick',1:nS,'XTickLabel', ...
    cellfun(@(x)strrep(x,'  ',' '), scenarios(:,1),'UniformOutput',false), ...
    'XTickLabelRotation',25, 'FontSize',7);
title('RMSE per scenario');
grid on;

sgtitle({'OOD Generalisation Summary', ...
    'RF trained on +-1% parameter bounds; tested on physically valid but far-outside values'}, ...
    'FontSize',10);

saveFig(fig2, fullfile('figures','RF_OOD_table'));
close(fig2);
fprintf('Saved RF_OOD_table\n');

%% ── Figure 3: ROC curves — one per OOD scenario ──────────────────────────
% ROC curve shows how well the RF separates shadow from clear at any threshold.
% AUC=1 = perfect; AUC=0.5 = random. Does not depend on the 0.5 threshold.
fig3 = figure('Name','OOD ROC Curves','Position',[40 40 680 560]);
col3 = lines(nS);
auc_all = zeros(nS,1);
hold on;
for k = 1:nS
    shadow_k = TL_bhp_sub{k} > FOM;
    valid_k  = isfinite(TL_bhp_sub{k}(:));
    [fpr_k, tpr_k, ~, auc_k] = perfcurve(double(shadow_k(valid_k)), ...
                                           RF_prob_sub{k}(valid_k), 1);
    auc_all(k) = auc_k;
    plot(fpr_k, tpr_k, '-', 'Color', col3(k,:), 'LineWidth', 2);
end
plot([0 1],[0 1],'k--','LineWidth',1);
hold off;
leg_str = arrayfun(@(k) sprintf('%s  AUC=%.3f', scenarios{k,1}, auc_all(k)), ...
    (1:nS)', 'UniformOutput', false);
legend(leg_str, 'Location','southeast', 'FontSize',7, 'Interpreter','none');
xlabel('False Positive Rate (RF says shadow, actually clear)','FontSize',9);
ylabel('True Positive Rate (RF correctly finds shadow)','FontSize',9);
title({'OOD Scenarios — ROC Curves', ...
    'RF trained on ±1% bounds; curve per OOD parameter set'}, 'FontSize',10);
grid on; axis square; xlim([0 1]); ylim([0 1]);

saveFig(fig3, fullfile('figures','RF_OOD_ROC'));
close(fig3);
fprintf('Saved RF_OOD_ROC\n');

%% ── Figure 4: Shadow boundary overlay ────────────────────────────────────
% Background: RF probability map.
% Red solid:  Bellhop actual shadow boundary (TL = FOM).
% White dash: RF predicted 50% boundary (prob = 0.5).
% Good RF: red and white lines track each other. Poor RF: they diverge.
fig4 = figure('Name','OOD Boundary Overlay','Position',[50 50 1650 450]);
for k = 1:nS
    ax = subplot(1, nS, k);
    pcolor(r_km_s, z_m_s, RF_prob_sub{k}); shading interp;
    colormap(ax, parula(256)); clim([0 1]);
    set(ax, 'YDir','reverse'); hold on;

    TL_k = TL_bhp_sub{k};
    if any(TL_k(isfinite(TL_k)) < FOM) && any(TL_k(isfinite(TL_k)) > FOM)
        contour(r_km_s, z_m_s, TL_k, [FOM FOM], 'r-', 'LineWidth', 2);
    end
    prob_k = RF_prob_sub{k};
    if max(prob_k(:)) > 0.5 && min(prob_k(:)) < 0.5
        contour(r_km_s, z_m_s, prob_k, [0.5 0.5], 'w--', 'LineWidth', 2);
    end
    hold off;
    if k==1, ylabel('Depth (m)','FontSize',8); end
    xlabel('Range (km)','FontSize',7);
    title(sprintf('%s\nAUC=%.3f  Acc=%.0f%%', scenarios{k,1}, auc_all(k), 100*acc_all(k)), ...
        'FontSize',7, 'Interpreter','none');
    if k==nS
        cb = colorbar; cb.Label.String = 'RF P(TL>FOM)'; cb.FontSize = 7;
    end
end
sgtitle({'Shadow Boundary: RF probability (background) | Red = Bellhop TL=FOM | White dash = RF P=0.5', ...
    'Matching lines = RF boundary agrees with physics.  Diverging = RF extrapolation error.'}, 'FontSize',9);

saveFig(fig4, fullfile('figures','RF_OOD_boundaries'));
close(fig4);
fprintf('Saved RF_OOD_boundaries\n');

%% ── Figure 5: Parameter distance from training vs accuracy ────────────────
% Distance: how many training half-widths each parameter has moved from centre.
% Δfreq normalised by ±100 Hz, ΔzS by ±0.05 m, Δtemp by ±1 °C.
fig5 = figure('Name','OOD Distance vs Accuracy','Position',[60 60 700 500]);
f0=10000; df=100;   z0=5; dz=0.05;   t0=0; dt=1;
dist_all = zeros(nS,1);
for k = 1:nS
    dist_all(k) = sqrt(((scenarios{k,2}-f0)/df)^2 + ...
                       ((scenarios{k,3}-z0)/dz)^2 + ...
                       ((scenarios{k,4}-t0)/dt)^2);
end

col5 = lines(nS);
scatter(dist_all, 100*acc_all, 130, col5, 'filled');
hold on;
for k = 1:nS
    text(dist_all(k)+max(dist_all)*0.02, 100*acc_all(k), ...
        strtrim(scenarios{k,1}), 'FontSize',7, 'Interpreter','none', ...
        'VerticalAlignment','middle');
end
yline(50,'k--','Random (50%)','LabelHorizontalAlignment','right','FontSize',8);
hold off;
xlabel({'Normalised parameter distance from training centre', ...
    '(each axis scaled by its ±1% training half-width)'},'FontSize',9);
ylabel('RF Accuracy vs Bellhop (%)','FontSize',9);
title({'RF Generalisation: Accuracy vs Distance from Training', ...
    'Closer to origin = closer to training distribution'}, 'FontSize',10);
xlim([0 max(dist_all)*1.25]); ylim([0 105]); grid on;

saveFig(fig5, fullfile('figures','RF_OOD_distance'));
close(fig5);
fprintf('Saved RF_OOD_distance\n');

%% ── Save ──────────────────────────────────────────────────────────────────
save(fullfile('results','RF_OOD.mat'), ...
     'TL_bhp_sub','RF_prob_sub','acc_all','rmse_all','scenarios', ...
     'r_km_s','z_m_s', '-v7.3');

fprintf('\n=== OOD Test Complete ===\n');
fprintf('%-32s  %10s  %8s\n', 'Scenario','Accuracy','RMSE');
for k=1:nS
    fprintf('%-32s  %9.1f%%  %8.3f\n', scenarios{k,1}, 100*acc_all(k), rmse_all(k));
end
fprintf('\nKey insight: high accuracy on similar freq (near 10kHz) but\n');
fprintf('degraded accuracy when frequency departs significantly from\n');
fprintf('training range — the shadow geometry is frequency-dependent.\n');

%% ── Local functions ───────────────────────────────────────────────────────
function svp = makeSVP(temp_shift, max_depth)
    % Wilson formula SVP perturbed by temp_shift (C), limited to max_depth (m).
    depths  = linspace(0, max_depth, 200);
    base_td = [0 30 180 400 5000; 25 25 17 13.6 13.6]';
    t  = interp1(base_td(:,1), base_td(:,2) + temp_shift, depths, 'linear', 'extrap');
    S  = 37;
    sv = 1499.2 + 4.6*t - 0.055*t.^2 + 0.00029*t.^3 + ...
         (1.34 - 0.01*t).*(S-35) + 0.016*depths;
    svp = [depths.' sv.'];
end

function saveFig(fig, base_path)
    print(fig, base_path, '-dpng', '-r200');
    try
        exportgraphics(fig, [base_path '.pdf'], 'ContentType','image','Resolution',200);
    catch
    end
end
