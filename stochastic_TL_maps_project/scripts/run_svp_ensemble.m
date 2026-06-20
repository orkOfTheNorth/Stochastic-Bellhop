%% run_svp_ensemble.m — Generate SVP ensemble figures for all completed MC combos
%
% For each completed scenario × distribution, plots all N SVP profiles
% overlaid on one graph (lightning-strike style) with the nominal in bold.
% Saves to Methods/MC/<scen>/<dist>/figures/svp_ensemble.png
% Replaces the old per-sample individual SVP PNGs.

clear; close all; clc; warning('off');
ROOT = fileparts(fileparts(mfilename('fullpath')));
try, cd(ROOT); catch; end
set(0, 'DefaultFigureVisible', 'off');
addpath(genpath(fullfile(ROOT, 'core')));
addpath(fullfile(ROOT, 'binaries'));

cfg   = loadConfig();
N     = cfg.MC.N;

fprintf('=== SVP Ensemble Figures ===\n');

for si = 1:numel(cfg.scenarios)
    sc = cfg.scenarios(si);
    for di = 1:numel(cfg.distributions)
        dist = cfg.distributions(di);

        fig_dir = fullfile('Methods', 'MC', sc.name, dist.name, 'figures');
        out_path = fullfile(fig_dir, 'svp_ensemble');

        if isfile([out_path '.png'])
            fprintf('  [SKIP] %s / %s — already exists\n', sc.name, dist.name);
            continue;
        end

        % Check that this combo has MC results (i.e. was run)
        res_dir = fullfile('Methods', 'MC', sc.name, dist.name, 'results');
        if ~isfile(fullfile(res_dir, 'MC_freq.mat'))
            fprintf('  [SKIP] %s / %s — MC not done yet\n', sc.name, dist.name);
            continue;
        end

        fprintf('  Generating: %s / %s\n', sc.name, dist.name);
        if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

        B = computeVarianceBounds(cfg, dist);
        S = lhsSample(N, cfg, dist, B);

        svp_nom = makeSVPNoise(0, sc.maxDepth_m, cfg);
        z_vec   = svp_nom(:,1);

        fig = figure('Position', [50 50 420 560]);
        ax  = axes;
        hold(ax, 'on');
        for i = 1:N
            svp_i = makeSVPNoise(S.svp(i), sc.maxDepth_m, cfg);
            plot(ax, svp_i(:,2), z_vec, '-', ...
                 'Color', [0.4 0.6 0.9 0.12], 'LineWidth', 0.8);
        end
        plot(ax, svp_nom(:,2), z_vec, 'k-', 'LineWidth', 2.5, 'DisplayName', 'Nominal');
        set(ax, 'YDir', 'reverse');
        xlabel(ax, 'Sound Speed (m/s)');
        ylabel(ax, 'Depth (m)');
        title(ax, sprintf('SVP Ensemble | %s | %s | N=%d', sc.name, dist.name, N), ...
              'Interpreter', 'none');
        grid(ax, 'on');
        saveFigPNG(fig, out_path);
        close(fig);
        fprintf('    Saved: %s\n', [out_path '.png']);
    end
end

fprintf('=== Done ===\n');
