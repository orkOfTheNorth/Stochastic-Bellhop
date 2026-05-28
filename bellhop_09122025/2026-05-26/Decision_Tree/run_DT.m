%% run_DT.m  ──  Random Forest shadow predictor  (2026-05-26)
%
% WHAT IT DOES
%   Trains a Random Forest (50 decision trees) to predict whether TL > FOM
%   at any grid point, given 5 inputs:
%       freq (Hz), src-depth zS (m), temp-shift (°C), range (km), depth (m)
%
%   The forest is trained on Bellhop data already saved in TL_cache
%   so NO new Bellhop runs are needed.
%
% WHY A RANDOM FOREST?
%   - Handles the complex, non-linear acoustic shadow zones naturally.
%   - Once trained it predicts in milliseconds (much faster than Bellhop).
%   - Gives feature importance: which input drives TL exceedance the most?
%
% OUTPUTS (saved to Decision_Tree/)
%   figures/RF_metrics.png   – confusion matrix, ROC curve, feature importance
%   figures/RF_maps.png      – 4-panel shadow prob map: RF vs MC vs Delta
%   figures/RF_agreement.png – spatial agreement between RF prediction and MC
%   results/RF_results.mat   – metrics, maps (forest NOT saved — too large)

clear; close all; clc;
try, cd(fileparts(mfilename('fullpath'))); catch; end

if ~exist('results','dir'), mkdir('results'); end
if ~exist('figures','dir'), mkdir('figures'); end

fprintf('=== Random Forest Shadow Predictor ===\n\n');

%% ══════════════════════════════════════════════════════════════════════════
%% STEP 1 ── LOAD SAVED TL SAMPLES FROM TL_CACHE
%  TL_save is (Nz × Nr × N), single precision, TL in dB.
%  The cache also stores the parameter bounds used when running the MC,
%  which we need to reconstruct the sample parameters.
%% ══════════════════════════════════════════════════════════════════════════
fprintf('Step 1: Loading TL cache...\n');
tc   = load(fullfile('..','Monte_Carlo','TL_cache','TL_zS_freq_temp.mat'));
TL_all   = tc.TL_save;          % single, Nz × Nr × N
N        = tc.N_saved;
freq_bnd = tc.freq_bnd;
zS_bnd   = tc.zS_bnd;
temp_bnd = tc.temp_bnd;
[Nz, Nr, ~] = size(TL_all);
fprintf('  Grid: %d × %d  |  N = %d MC samples\n', Nz, Nr, N);

% Load grid axes and MC empirical shadow probabilities for comparison
mc_res  = load(fullfile('..','Monte_Carlo','results','MC_zS_freq_temp.mat'), ...
               'r_km','z_m','MC_PrFOM');
r_km    = mc_res.r_km(:)';   % 1 × Nr  (row vector)
z_m     = mc_res.z_m(:)';    % 1 × Nz  (row vector)
MC_PrFOM = mc_res.MC_PrFOM;  % Nz × Nr  (empirical P(TL>FOM) from MC)

FOM  = 100;   % figure of merit threshold (dB)
freq0 = 10000;
zS0   = 5;

%% ══════════════════════════════════════════════════════════════════════════
%% STEP 2 ── RECONSTRUCT WHICH (freq, zS, temp) EACH SAMPLE USED
%  The MC used rng(42) + lhsdesign(N,3).  Re-running that here gives
%  the exact same random samples (same seed, same N) so we recover
%  the input parameters without storing them separately.
%% ══════════════════════════════════════════════════════════════════════════
fprintf('Step 2: Reconstructing MC sample parameters (rng=42)...\n');
rng(42);
lhs       = lhsdesign(N, 3);
samp_freq = freq_bnd(1) + lhs(:,1) * diff(freq_bnd);   % N×1 Hz
samp_zS   = zS_bnd(1)   + lhs(:,2) * diff(zS_bnd);    % N×1 m
samp_temp = temp_bnd(1)  + lhs(:,3) * diff(temp_bnd); % N×1 °C shift

%% ══════════════════════════════════════════════════════════════════════════
%% STEP 3 ── BUILD FEATURE MATRIX AND BINARY LABELS
%  Subsample the spatial grid every 5 points in both dimensions so the
%  training set fits comfortably in memory.
%
%  Each row represents one (MC sample, grid point) combination:
%    X = [freq_Hz | zS_m | temp_shift_C | range_km | depth_m]
%    y = 1 if TL > FOM (shadow zone), 0 otherwise
%% ══════════════════════════════════════════════════════════════════════════
fprintf('Step 3: Building feature matrix...\n');
stride  = 5;                          % spatial subsampling step
zi_sub  = 1:stride:Nz;               % subsampled depth indices
ri_sub  = 1:stride:Nr;               % subsampled range indices
nzi     = numel(zi_sub);
nri     = numel(ri_sub);
n_pts   = nzi * nri;                 % grid points per sample
n_rows  = N * n_pts;                 % total rows in feature matrix

% Build a fixed grid of (range, depth) coordinates (same for every sample)
[RI, ZI] = meshgrid(ri_sub, zi_sub);   % nzi × nri
r_vec    = r_km(RI(:))';               % n_pts × 1
z_vec    = z_m(ZI(:))';                % n_pts × 1

% Allocate feature matrix and label vector
X = zeros(n_rows, 5);          % [freq, zS, temp, r_km, z_m]
y = false(n_rows, 1);          % 1 = shadow (TL > FOM)

for i = 1:N
    rows = (i-1)*n_pts + (1:n_pts);

    % 5 input features for this sample at every subsampled grid point
    X(rows,:) = [repmat(samp_freq(i), n_pts, 1), ...
                 repmat(samp_zS(i),   n_pts, 1), ...
                 repmat(samp_temp(i), n_pts, 1), ...
                 r_vec, z_vec];

    % Label: is TL above FOM at each grid point?
    tl_sub   = TL_all(zi_sub, ri_sub, i);   % nzi × nri slice
    y(rows)  = tl_sub(:) > FOM;
end

feature_names = {'Freq (Hz)','Src depth (m)','Temp shift (°C)','Range (km)','Depth (m)'};
fprintf('  Feature matrix: %d rows × %d features\n', n_rows, size(X,2));
fprintf('  Shadow fraction: %.1f%% of all labels\n', 100*mean(y));

%% ══════════════════════════════════════════════════════════════════════════
%% STEP 4 ── TRAIN / TEST SPLIT BY SAMPLE INDEX
%  We split by MC sample (not by grid point) to avoid data leakage:
%  a single Bellhop run should not appear in both train and test sets.
%  80% of samples → train,  20% of samples → test.
%% ══════════════════════════════════════════════════════════════════════════
n_train    = floor(0.8 * N);
n_test     = N - n_train;
fprintf('\nStep 4: %d training samples (%d rows), %d test samples (%d rows)\n', ...
    n_train, n_train*n_pts, n_test, n_test*n_pts);

train_rows = 1 : n_train * n_pts;
test_rows  = (n_train*n_pts + 1) : n_rows;

X_train = X(train_rows, :);
y_train = y(train_rows);
X_test  = X(test_rows,  :);
y_test  = y(test_rows);

%% ══════════════════════════════════════════════════════════════════════════
%% STEP 5 ── TRAIN RANDOM FOREST (pruned for generalisation)
%  TreeBagger creates 200 decision trees by bagging (bootstrap sampling).
%
%  PRUNING / TRIMMING rationale:
%    With only N=100 LHS samples (80 training), there are just 80 unique
%    (freq, zS, temp) parameter combinations.  A deep tree (MinLeafSize=10)
%    can memorise which of those 80 combos produce shadow at each spatial
%    point rather than learning a smooth, generalisable function.
%
%    Fixes applied:
%      MinLeafSize  10 → 200  Forces leaves to average ≥200 training rows,
%                             preventing the RF from splitting down to single
%                             parameter-combination neighbourhoods.
%      NumTrees     50 → 200  More trees compensate for the higher per-tree
%                             bias; the ensemble variance stays low.
%      NumPredictorsToSample  kept at 3 (= round(sqrt(5)))
%% ══════════════════════════════════════════════════════════════════════════
fprintf('\nStep 5: Training Random Forest (200 trees, MinLeafSize=200)...\n');
tic;
RF = TreeBagger(200, X_train, double(y_train), ...
    'Method',                 'classification', ...
    'OOBPredictorImportance', 'on', ...
    'MinLeafSize',            200,  ...
    'NumPredictorsToSample',  3);
t_train = toc;
fprintf('  Training done in %.1f s\n', t_train);

%% ══════════════════════════════════════════════════════════════════════════
%% STEP 6 ── EVALUATE ON TEST SET
%  predict() returns: class labels (cell array) + class probabilities.
%  Class 1 (index 2 in scores) = shadow zone probability.
%% ══════════════════════════════════════════════════════════════════════════
fprintf('\nStep 6: Evaluating on test set...\n');
[lbl_cell, scores_test] = predict(RF, X_test);
pred_class = cellfun(@str2double, lbl_cell);   % cell '0'/'1' → numeric 0/1
pred_proba = scores_test(:,2);                 % P(TL > FOM) per row

accuracy     = mean(pred_class == double(y_test));
[~,~,~,auc] = perfcurve(double(y_test), pred_proba, 1);

tp = sum(pred_class==1 & y_test==1);
fp = sum(pred_class==1 & y_test==0);
tn = sum(pred_class==0 & y_test==0);
fn = sum(pred_class==0 & y_test==1);

fprintf('  Accuracy : %.1f%%\n', 100*accuracy);
fprintf('  AUC      : %.4f\n',   auc);
fprintf('  Confusion: TP=%d  FP=%d  TN=%d  FN=%d\n', tp, fp, tn, fn);

%% ══════════════════════════════════════════════════════════════════════════
%% STEP 7 ── BUILD SHADOW PROBABILITY MAPS
%
%  Map (a) – RF at NOMINAL params (freq0, zS0, temp_shift=0):
%    What does the forest predict at every grid point if all inputs
%    are at their nominal (central) values?
%
%  Map (b) – RF AVERAGED OVER MC UNCERTAINTY:
%    For each grid point, run all N LHS parameter sets through the forest
%    and average the predicted probabilities.  This is the RF equivalent
%    of MC_PrFOM — the shadow probability marginalized over uncertainty.
%% ══════════════════════════════════════════════════════════════════════════
fprintf('\nStep 7: Computing prediction maps...\n');

% Subsampled map grid (same stride as training)
[RI_m, ZI_m] = meshgrid(ri_sub, zi_sub);
r_vec_m      = r_km(RI_m(:))';        % n_pts × 1
z_vec_m      = z_m(ZI_m(:))';         % n_pts × 1
n_map        = n_pts;

% (a) Nominal prediction map
X_nom = [repmat(freq0, n_map, 1), repmat(zS0, n_map, 1), ...
         zeros(n_map,1), r_vec_m, z_vec_m];
[~, sc_nom]     = predict(RF, X_nom);
RF_nominal_map  = reshape(sc_nom(:,2), nzi, nri);

% (b) MC-averaged map: stack all N×n_map predictions, then average
fprintf('  Predicting over N=%d samples × %d grid pts = %d total...\n', ...
    N, n_map, N*n_map);
X_mc = zeros(N * n_map, 5);
for i = 1:N
    rows = (i-1)*n_map + (1:n_map);
    X_mc(rows,:) = [repmat(samp_freq(i), n_map, 1), ...
                    repmat(samp_zS(i),   n_map, 1), ...
                    repmat(samp_temp(i), n_map, 1), ...
                    r_vec_m, z_vec_m];
end
[~, sc_mc]   = predict(RF, X_mc);
probs_mc     = reshape(sc_mc(:,2), n_map, N);   % rows=grid pts, cols=samples
RF_mc_avg    = reshape(mean(probs_mc, 2), nzi, nri);   % average over samples

fprintf('  Maps complete.\n');

%% ══════════════════════════════════════════════════════════════════════════
%% STEP 8 ── LOAD DELTA CHEBYSHEV MAP FOR COMPARISON
%% ══════════════════════════════════════════════════════════════════════════
D        = load(fullfile('..','Delta_Method','results','delta_zS_freq_temp.mat'), 'Cheb_lb_s');
Cheb_lb  = D.Cheb_lb_s;   % Nz × Nr

% Subsample all comparison maps to the same grid as the RF maps
MC_sub   = MC_PrFOM(zi_sub, ri_sub);   % MC empirical shadow probability
Cheb_sub = Cheb_lb(zi_sub, ri_sub);    % Delta Chebyshev lower bound

r_km_map = r_km(ri_sub);
z_m_map  = z_m(zi_sub);

%% ══════════════════════════════════════════════════════════════════════════
%% FIGURE 1 ── METRICS: confusion matrix, ROC curve, feature importance
%% ══════════════════════════════════════════════════════════════════════════
fig1 = figure('Name','RF Metrics','Position',[30 30 1350 440]);

% --- Confusion matrix ---
ax1 = subplot(1,3,1);
cm  = [tn fp; fn tp];
imagesc(cm); colormap(ax1, flipud(hot(64))); colorbar;
set(gca, 'XTick',1:2, 'XTickLabel',{'True: Clear','True: Shadow'}, ...
         'YTick',1:2, 'YTickLabel',{'Pred: Clear','Pred: Shadow'}, ...
         'FontSize', 9);
for r=1:2, for c=1:2
    text(c, r, num2str(cm(r,c),'%d'), 'HorizontalAlignment','center', ...
         'FontSize',11, 'Color','w', 'FontWeight','bold');
end; end
title(sprintf('Confusion Matrix\nAccuracy = %.1f%%  |  AUC = %.3f', 100*accuracy, auc));

% --- ROC curve ---
subplot(1,3,2);
[fpr, tpr] = perfcurve(double(y_test), pred_proba, 1);
plot(fpr, tpr, 'b-', 'LineWidth',2); hold on;
plot([0 1],[0 1],'k--','LineWidth',1);
xlabel('False Positive Rate'); ylabel('True Positive Rate');
title(sprintf('ROC Curve  (AUC = %.3f)', auc));
legend(sprintf('RF  AUC=%.3f', auc), 'Random', 'Location','southeast');
grid on; axis square;

% --- Feature importance ---
subplot(1,3,3);
imp = RF.OOBPermutedPredictorDeltaError;
barh(imp, 'FaceColor',[0.25 0.55 0.85]);
set(gca,'YTick',1:5,'YTickLabel',feature_names,'FontSize',9);
xlabel('OOB Permutation Importance');
title('Feature Importance');
grid on;

sgtitle(sprintf('Random Forest | %d trees | MinLeafSize=%d | %d train | %d test', ...
    RF.NumTrees, 200, n_train, n_test), 'FontSize',11);
saveFig(fig1, fullfile('figures','RF_metrics'));
close(fig1);
fprintf('Saved RF_metrics\n');

%% ══════════════════════════════════════════════════════════════════════════
%% FIGURE 2 ── SHADOW PROBABILITY MAPS: 4-panel comparison
%% ══════════════════════════════════════════════════════════════════════════
fig2 = figure('Name','Shadow Probability Maps','Position',[40 40 1600 500]);
panel_maps   = {RF_nominal_map,  RF_mc_avg,      MC_sub,       Cheb_sub};
panel_titles = {'RF — Nominal params', 'RF — MC-averaged (uncertainty)', ...
                'MC Empirical P(TL>FOM)', 'Delta Chebyshev lower bound'};

for p = 1:4
    ax = subplot(1,4,p);
    pcolor(r_km_map, z_m_map, panel_maps{p}); shading interp;
    set(ax,'YDir','reverse');
    colormap(ax, parula(256)); cb = colorbar;
    clim([0 1]);
    if p==4, ylabel(cb,'Shadow probability','FontSize',8); end
    xlabel('Range (km)','FontSize',8);
    if p==1, ylabel('Depth (m)','FontSize',8); end
    title(panel_titles{p},'FontSize',8);
end
sgtitle(sprintf('Shadow probability maps | FOM=%d dB | N=%d samples | stride=%d pts', ...
    FOM, N, stride), 'FontSize',10);
saveFig(fig2, fullfile('figures','RF_maps'));
close(fig2);
fprintf('Saved RF_maps\n');

%% ══════════════════════════════════════════════════════════════════════════
%% FIGURE 3 ── AGREEMENT: RF vs MC (scatter + spatial difference map)
%% ══════════════════════════════════════════════════════════════════════════
diff_map  = RF_mc_avg - MC_sub;           % positive = RF predicts more shadow
rmse_val  = sqrt(mean(diff_map(:).^2));
r2_val    = corr(MC_sub(:), RF_mc_avg(:))^2;
bias_val  = mean(diff_map(:));

fig3 = figure('Name','RF vs MC Agreement','Position',[50 50 1200 500]);

% Density scatter plot — 2D histogram so every point is visible
subplot(1,2,1);
histogram2(MC_sub(:), RF_mc_avg(:), [60 60], ...
    'DisplayStyle','tile','ShowEmptyBins','off');
colormap(gca, flipud(hot(256))); cb = colorbar;
ylabel(cb,'Point count');
hold on;
plot([0 1],[0 1],'w-','LineWidth',2);   % white 1:1 line visible on dark bg
xlabel('MC Empirical P(TL>FOM)'); ylabel('RF MC-averaged P(TL>FOM)');
title(sprintf('RF vs MC density  |  RMSE=%.3f   R²=%.3f', rmse_val, r2_val));
axis equal; xlim([0 1]); ylim([0 1]); grid on;

% Spatial difference map
subplot(1,2,2);
pcolor(r_km_map, z_m_map, diff_map); shading interp;
set(gca,'YDir','reverse');
colormap(gca, redblue(256)); cb = colorbar;
cb_lim = max(0.01, max(abs(diff_map(:))));
clim([-cb_lim  cb_lim]);
ylabel(cb,'RF − MC (probability)','FontSize',9);
xlabel('Range (km)'); ylabel('Depth (m)');
title(sprintf('Difference map (RF − MC)   bias=%.3f   RMSE=%.3f', bias_val, rmse_val));

sgtitle('Random Forest vs Monte Carlo Agreement','FontSize',11);
saveFig(fig3, fullfile('figures','RF_agreement'));
close(fig3);
fprintf('Saved RF_agreement\n');

%% ══════════════════════════════════════════════════════════════════════════
%% SAVE METRICS AND MAPS
%% ══════════════════════════════════════════════════════════════════════════
metrics = struct( ...
    'accuracy',       accuracy, ...
    'AUC',            auc, ...
    'RMSE_vs_MC',     rmse_val, ...
    'R2_vs_MC',       r2_val, ...
    'bias_vs_MC',     bias_val, ...
    'TP', tp, 'FP', fp, 'TN', tn, 'FN', fn, ...
    'n_train_samples', n_train, ...
    'n_test_samples',  n_test, ...
    'feature_names',  {feature_names});

save(fullfile('results','RF_results.mat'), ...
     'RF', 'metrics', 'RF_nominal_map', 'RF_mc_avg', 'r_km_map', 'z_m_map', '-v7.3');

fprintf('\n=== Summary ===\n');
fprintf('  Accuracy    : %.1f%%\n', 100*accuracy);
fprintf('  AUC         : %.4f\n',   auc);
fprintf('  RMSE vs MC  : %.4f\n',   rmse_val);
fprintf('  R²   vs MC  : %.4f\n',   r2_val);
fprintf('\nAll figures saved to figures/\n');

%% ══════════════════════════════════════════════════════════════════════════
%% LOCAL FUNCTIONS
%% ══════════════════════════════════════════════════════════════════════════
function saveFig(fig, base_path)
    print(fig, base_path, '-dpng', '-r200');
    try
        exportgraphics(fig, [base_path '.pdf'], 'ContentType','image','Resolution',200);
    catch
    end
end

function cmap = redblue(n)
    % Diverging colormap: blue → white → red  (white at centre = zero diff)
    n2   = ceil(n/2);
    n1   = n - n2;
    cmap = [ linspace(0,1,n1)', linspace(0,1,n1)', ones(n1,1); ...
             ones(n2,1), linspace(1,0,n2)', linspace(1,0,n2)' ];
end
