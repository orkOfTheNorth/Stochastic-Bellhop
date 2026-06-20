%% regenerate_pce_summary.m
% Regenerates Methods/PCE/pce_summary_all.png from saved pce_results.mat files.
% Applies y-axis clipping to [0, 30%] so the high-order numerical blow-up
% (K=18-20, values reaching 10^16%) does not destroy the plot scale.

ROOT = fileparts(fileparts(mfilename('fullpath')));
try; cd(ROOT); catch; end
addpath(genpath(fullfile(ROOT, 'core')));

PARAMS = {'freq', 'zS', 'svp'};
YLIM   = 30;   % clip y-axis at 30% relative L2 — the 10% convergence line is clearly visible

% Collect all pce_results.mat files
results_files = dir(fullfile('Methods','PCE','*','*','results','pce_results.mat'));
n_combos = numel(results_files);
if n_combos == 0
    error('No pce_results.mat files found. Run run_pce.m first.');
end

sum_L2    = {cell(n_combos,1), cell(n_combos,1), cell(n_combos,1)};
sum_label = cell(n_combos, 1);
MAX_ORDER_all = zeros(n_combos,1);

for ci = 1:n_combos
    r = load(results_files(ci).folder + "/" + results_files(ci).name);
    sum_label{ci} = sprintf('%s / %s', r.sc_name, r.dist_name);
    MAX_ORDER_all(ci) = r.MAX_ORDER;
    for pi = 1:3
        sum_L2{pi}{ci} = r.L2_mat(:, pi);
    end
end

MAX_ORDER = max(MAX_ORDER_all);
cmap = lines(n_combos);

fig_sum = figure('Position', [50 50 1400 900], 'Visible', 'off');
for pi = 1:3
    ax = subplot(3, 1, pi);
    for ci = 1:n_combos
        l2 = sum_L2{pi}{ci};
        if isempty(l2) || all(isnan(l2)), continue; end
        mo = numel(l2);
        plot(ax, 1:mo, l2*100, '-o', 'Color', cmap(ci,:), ...
             'LineWidth', 1.2, 'MarkerSize', 4, 'DisplayName', sum_label{ci});
        hold(ax,'on');
    end
    yline(ax, 10, 'k--', 'LineWidth', 1.5, 'DisplayName', '10% threshold');
    ylim(ax, [0, YLIM]);
    set(ax, 'XTick', max(1, 1:2:MAX_ORDER));
    ylabel(ax, 'Relative L2 (%)');
    title(ax, sprintf('param = %s  [y-axis clipped at %d%% — blow-up at K>17 omitted]', ...
                      PARAMS{pi}, YLIM), 'FontSize', 9, 'Interpreter', 'none');
    if pi == 1
        legend(ax, 'Location', 'northeast', 'FontSize', 6, 'NumColumns', 2);
    end
    if pi == 3, xlabel(ax, 'PCE Order K'); end
    grid(ax, 'on');
end
sgtitle('PCE Convergence — All Scenarios × Distributions (y-axis clipped at 30%)', ...
        'FontSize', 11);

out_path = fullfile('Methods', 'PCE', 'pce_summary_all');
saveas(fig_sum, [out_path '.png']);
close(fig_sum);
fprintf('Saved: %s.png\n', out_path);
