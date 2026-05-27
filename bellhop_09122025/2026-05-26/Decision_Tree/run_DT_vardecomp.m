%% run_DT_vardecomp.m  —  RF Permutation Feature Importance  (2026-05-26)
%
% METHOD: OOB Permutation Importance
% ────────────────────────────────────────────────────────────────────────
%   During RF training each tree is evaluated on its out-of-bag (OOB)
%   samples — the ~37% of training rows NOT used to build that tree.
%   This gives an unbiased accuracy estimate without a separate test set.
%
%   For each feature k (one at a time):
%     1. Take the OOB predictions from all trees
%     2. Randomly SHUFFLE the values of feature k across all OOB samples
%        (this breaks the real relationship between k and the shadow label)
%     3. Re-predict with the shuffled feature; measure accuracy drop
%     4. OOB_importance_k = baseline_accuracy − shuffled_accuracy
%
%   Larger drop = the RF relied on this feature to classify correctly.
%   Near-zero drop = the RF can do fine without this feature.
%
%   MATLAB computes this automatically via OOBPredictorImportance='on'
%   during TreeBagger training and stores it in RF.OOBPermutedPredictorDeltaError.
%
% NOTE ON RANGE AND DEPTH
%   Range (km) and Depth (m) are spatial COORDINATES, not uncertain
%   parameters.  Their high importance simply means "knowing where you are
%   in the field is the strongest predictor of whether you are in a shadow
%   zone" — which is physically obvious.  The physically interesting
%   question (answered separately in the variance decomposition) is which
%   of the UNCERTAIN inputs drives the most uncertainty.
%
% REQUIRES:  results/RF_results.mat  (from run_DT.m)
% OUTPUT:
%   figures/RF_permutation_importance.png
%   results/feature_importance_table.csv
%   results/feature_importance_table.txt

clear; close all; clc;
try, cd(fileparts(mfilename('fullpath'))); catch; end
if ~exist('results','dir'), mkdir('results'); end
if ~exist('figures','dir'), mkdir('figures'); end

fprintf('=== RF Permutation Feature Importance ===\n\n');

%% ── Load RF ───────────────────────────────────────────────────────────────
res = load(fullfile('results','RF_results.mat'), 'RF', 'metrics');
RF  = res.RF;

feat_names = {'Frequency (Hz)', 'Source depth z_S (m)', 'Temp shift (°C)', ...
              'Range (km)', 'Depth (m)'};
n_feat = numel(feat_names);

oob_imp  = RF.OOBPermutedPredictorDeltaError(:);   % 5×1  accuracy drop
oob_pct  = 100 * oob_imp / sum(oob_imp);           % normalised to 100%

% Sort descending for display
[oob_sorted, sort_idx] = sort(oob_imp, 'descend');
pct_sorted   = oob_pct(sort_idx);
names_sorted = feat_names(sort_idx);

fprintf('  %-26s  %12s  %8s\n', 'Feature', 'OOB Imp.', '% total');
fprintf('  %s\n', repmat('-',1,52));
for j = 1:n_feat
    fprintf('  %-26s  %12.5f  %7.1f%%\n', names_sorted{j}, oob_sorted(j), pct_sorted(j));
end

%% ── Colours: physical params (blue family) vs spatial coords (grey) ───────
% Features 1-3 are physical, 4-5 are spatial coordinates
is_physical = [true true true false false];   % original order
clr_phys    = [0.20 0.45 0.85;   % freq     — blue
               0.90 0.50 0.10;   % zS       — orange
               0.15 0.65 0.30];  % temp     — green
clr_spatial = [0.60 0.60 0.60;   % range    — grey
               0.75 0.75 0.75];  % depth    — light grey
all_colors  = [clr_phys; clr_spatial];   % 5×3 in original feature order

% Reorder colours to match sorted order
sorted_colors = all_colors(sort_idx, :);

%% ── FIGURE — Permutation Importance ──────────────────────────────────────
fig = figure('Name','RF Permutation Importance','Position',[50 50 950 540]);

% Left panel: sorted bar chart (raw OOB importance)
ax1 = subplot(1,2,1);
b1  = barh(n_feat:-1:1, oob_sorted, 'FaceColor','flat');
b1.CData = flipud(sorted_colors);
hold on;
% Percentage labels at end of each bar
for j = 1:n_feat
    text(oob_sorted(j) + max(oob_sorted)*0.015, n_feat-j+1, ...
         sprintf('%.1f%%', pct_sorted(j)), ...
         'VerticalAlignment','middle','FontSize',9,'FontWeight','bold', ...
         'Color', sorted_colors(j,:));
end
set(ax1,'YTick',1:n_feat,'YTickLabel',fliplr(names_sorted),'FontSize',9);
xlabel('OOB accuracy drop when shuffled','FontSize',9);
title({'\bfOOB Permutation Importance', 'sorted by importance'}, ...
    'FontSize',9,'Interpreter','tex');
grid on; axis tight;
xl = xlim; xlim([0, xl(2)*1.18]);

% Shade spatial-coordinate rows to indicate they are not uncertain inputs
for j = 1:n_feat
    if ~is_physical(sort_idx(n_feat-j+1))
        patch([0 xl(2)*1.18 xl(2)*1.18 0], [j-0.5 j-0.5 j+0.5 j+0.5], ...
              [0.93 0.93 0.93], 'FaceAlpha',0.4,'EdgeColor','none');
    end
end

% Right panel: normalised bar chart (% of total)
ax2 = subplot(1,2,2);
b2  = barh(n_feat:-1:1, pct_sorted, 'FaceColor','flat');
b2.CData = flipud(sorted_colors);
hold on;
for j = 1:n_feat
    text(pct_sorted(j) + 0.5, n_feat-j+1, ...
         sprintf('%.1f%%', pct_sorted(j)), ...
         'VerticalAlignment','middle','FontSize',9,'FontWeight','bold', ...
         'Color', sorted_colors(j,:));
end
set(ax2,'YTick',1:n_feat,'YTickLabel',fliplr(names_sorted),'FontSize',9);
xlabel('% of total OOB importance','FontSize',9);
title({'\bfNormalised importance', '(sums to 100%)'}, ...
    'FontSize',9,'Interpreter','tex');
grid on; axis tight; xlim([0 115]);
xline(100/n_feat, 'k:', 'LineWidth',0.9,'HandleVisibility','off');
text(100/n_feat + 0.3, 0.55, 'Equal weight', 'FontSize',7,'Color',[0.5 0.5 0.5]);

% Shade spatial rows
for j = 1:n_feat
    if ~is_physical(sort_idx(n_feat-j+1))
        patch([0 115 115 0], [j-0.5 j-0.5 j+0.5 j+0.5], ...
              [0.93 0.93 0.93], 'FaceAlpha',0.4,'EdgeColor','none');
    end
end

% Shared legend for colour coding
annotation(fig,'textbox',[0.01 0.01 0.98 0.06], ...
    'String', ['Colour: Blue/Orange/Green = physical uncertain inputs  |  ' ...
               'Grey shading = spatial coordinates (not uncertain inputs)'], ...
    'FontSize',7,'EdgeColor','none','HorizontalAlignment','center', ...
    'Color',[0.3 0.3 0.3]);

sgtitle({ ...
    '\bfRF OOB Permutation Importance — which features does the RF rely on?', ...
    sprintf('N=%d LHS samples | %d trees | Test accuracy=%.1f%% | AUC=%.4f', ...
            100, RF.NumTrees, 100*res.metrics.accuracy, res.metrics.AUC)}, ...
    'FontSize',10,'Interpreter','tex');

saveFig(fig, fullfile('figures','RF_permutation_importance'));
close(fig);
fprintf('\nSaved RF_permutation_importance\n');

%% ── TABLE FILES ───────────────────────────────────────────────────────────
% CSV
T = table(feat_names(:), oob_imp, oob_pct, ...
    'VariableNames', {'Parameter','OOB_Importance','OOB_Pct'});
csv_path = fullfile('results','feature_importance_table.csv');
writetable(T, csv_path);
fprintf('Saved %s\n', csv_path);

% Formatted text
txt_path = fullfile('results','feature_importance_table.txt');
fid = fopen(txt_path,'w');
SEP  = [repmat('=',1,70) '\n'];
DASH = [repmat('-',1,56) '\n'];

fprintf(fid, SEP);
fprintf(fid, '  RF Permutation Feature Importance\n');
fprintf(fid, '  N=100 LHS samples | %d trees | Accuracy=%.1f%% | AUC=%.4f\n', ...
        RF.NumTrees, 100*res.metrics.accuracy, res.metrics.AUC);
fprintf(fid, SEP);
fprintf(fid, '\n  METHOD: OOB Permutation Importance\n');
fprintf(fid, DASH);
fprintf(fid, '  Each feature is randomly shuffled (one at a time) across\n');
fprintf(fid, '  out-of-bag samples.  The drop in accuracy = importance.\n');
fprintf(fid, '  Larger drop = RF relied on this feature more.\n\n');
fprintf(fid, '  Rank  %-28s  %11s  %7s  Bar\n','Feature','OOB Imp.','%% total');
fprintf(fid, DASH);
for j = 1:n_feat
    bw = max(1, round(30 * pct_sorted(j)/100));
    tag = '';
    if ~is_physical(sort_idx(j)), tag = '  [spatial coord]'; end
    fprintf(fid,'  [%d]   %-28s  %11.5f  %6.1f%%  %s%s\n', ...
            j, names_sorted{j}, oob_sorted(j), pct_sorted(j), ...
            repmat('|',1,bw), tag);
end
fprintf(fid,'\n  Note: Range and Depth are spatial coordinates, not uncertain\n');
fprintf(fid,'  inputs — their high importance reflects that acoustic shadow\n');
fprintf(fid,'  zones have fixed spatial structure, not parameter sensitivity.\n');
fprintf(fid,'\n  Among uncertain physical inputs only (Freq + zS + Temp = %.1f%%):\n', ...
        sum(oob_pct(1:3)));
phys_pct = 100 * oob_imp(1:3) / sum(oob_imp(1:3));
[~, ps] = sort(oob_imp(1:3),'descend');
for j=1:3
    k = ps(j);
    fprintf(fid,'    [%d] %-26s  %.1f%% of physical-param importance\n', ...
            j, feat_names{k}, phys_pct(k));
end
fprintf(fid,'\n');
fprintf(fid, SEP);
fclose(fid);
fprintf('Saved %s\n', txt_path);

fprintf('\n=== Done ===\n');

%% ── Local functions ───────────────────────────────────────────────────────
function saveFig(fig, base_path)
    print(fig, base_path, '-dpng', '-r200');
    try
        exportgraphics(fig,[base_path '.pdf'],'ContentType','image','Resolution',200);
    catch
    end
end
