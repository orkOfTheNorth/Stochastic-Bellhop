%% fill_findings.m — Extract computed values for findings.tex dagger placeholders.
%
% Reads cached .mat result files and prints LaTeX-ready numbers for every
% dagger marker (†) left in findings.tex.  Run once after a full pipeline
% pass; copy-paste the printed snippets into the corresponding tables.
%
% Usage (from stochastic_TL_maps_project/):
%   matlab -batch "fill_findings"
%   — or interactively: fill_findings

try
    cd(fileparts(mfilename('fullpath')));
catch
end

addpath(genpath('Shared_Utils'));

cfg    = loadConfig();
FOM    = cfg.nominal.FOM_dB;
PARAMS = {'freq', 'zS', 'svp'};

fprintf('\n==================================================================\n');
fprintf('  fill_findings.m — LaTeX values for findings.tex dagger markers\n');
fprintf('==================================================================\n\n');

% ─────────────────────────────────────────────────────────────────────────────
%  TABLE 1: tab:ks_ln3  — Median KS statistic per scenario × distribution
%  Recomputes from MC TL cubes using ln3moments (no new LN3 run needed).
% ─────────────────────────────────────────────────────────────────────────────
TARGET_DISTS = {'Normal_5pct', 'Normal_10pct'};
scen_names   = {cfg.scenarios.name};

fprintf('--- TABLE: tab:ks_ln3  (Median KS for LN3 fit) ---\n');
fprintf('%-18s  %16s  %16s\n', 'Scenario', 'Normal 5% (KS)', 'Normal 10% (KS)');
fprintf('%s\n', repmat('-', 1, 56));

for si = 1:numel(cfg.scenarios)
    sc = cfg.scenarios(si);
    row_vals = nan(1, numel(TARGET_DISTS));

    for di = 1:numel(TARGET_DISTS)
        % TL samples always come from the IS base distribution cache (Normal_10pct).
        % IS reweighting only affects statistics, not the underlying sample set.
        tl_file = fullfile('Cache', sc.name, cfg.IS_base_dist, sprintf('TL_%s.mat', PARAMS{2}));
        if ~isfile(tl_file)
            tl_file = fullfile('Cache', sc.name, cfg.IS_base_dist, 'TL_zS.mat');
        end
        if ~isfile(tl_file)
            continue;
        end
        try
            d   = load(tl_file, 'TL_save');
            TL  = double(d.TL_save);
            [Nz, Nr, N] = size(TL);
            TL_pix = reshape(permute(TL, [3 1 2]), N, Nz*Nr);
            % KS statistic: compare empirical CDF to fitted LN3 CDF at each pixel
            ks_vals = zeros(1, Nz*Nr);
            for px = 1:Nz*Nr
                x = TL_pix(:, px);
                x = x(isfinite(x));
                if numel(x) < 4, continue; end
                gamma = 0.95 * min(x);
                y = log(x - gamma);
                y = y(isfinite(y));
                if numel(y) < 4, continue; end
                mu_y  = mean(y);
                sig_y = std(y);
                if sig_y < 1e-12, continue; end
                [~, ~, ks_stat] = kstest((y - mu_y) / sig_y);
                ks_vals(px) = ks_stat;
            end
            row_vals(di) = median(ks_vals(ks_vals > 0));
        catch me
            fprintf('  [WARN] %s / %s : %s\n', sc.name, TARGET_DISTS{di}, me.message);
        end
    end

    v1 = row_vals(1); v2 = row_vals(2);
    s1 = ternary(isnan(v1), 'N/A', sprintf('%.3f', v1));
    s2 = ternary(isnan(v2), 'N/A', sprintf('%.3f', v2));
    fprintf('%-18s  %16s  %16s\n', sc.name, s1, s2);
end

fprintf('\nLaTeX snippet (copy into tab:ks_ln3):\n');
fprintf('%%  (re-run fill_findings.m and paste values below)\n\n');

% ─────────────────────────────────────────────────────────────────────────────
%  TABLE 2: tab:analytical_l1  — L1 errors from analytical waveguide validation
% ─────────────────────────────────────────────────────────────────────────────
fprintf('\n--- TABLE: tab:analytical_l1  (Analytical Delta validation L1) ---\n');
anal_dir = fullfile('analytical', 'results');
dist_labels = {'Normal_1pct', 'Normal_5pct', 'Normal_10pct'};
sigma_vals  = [0.5, 2.5, 5.0];

fprintf('%-14s  %8s  %22s  %22s\n', 'Distribution', 'σ_zS (m)', 'L1 (Analytical J)', 'L1 (FD J)');
fprintf('%s\n', repmat('-', 1, 72));

for k = 1:3
    % run_analytical_delta_validation.m saves to validation_<dist>.mat
    afile = fullfile(anal_dir, sprintf('validation_%s.mat', dist_labels{k}));
    if ~isfile(afile)
        fprintf('%-14s  %8.1f  %22s  %22s\n', dist_labels{k}, sigma_vals(k), 'file missing', 'file missing');
        continue;
    end
    try
        a = load(afile);
        % Variable names from run_analytical_delta_validation.m
        null_mask = a.TL_nom(:) <= 120;
        mean_VM   = max(mean(a.Var_MC(null_mask)), 1e-10);
        l1_anal = mean(abs(a.Var_delta_anal(null_mask) - a.Var_MC(null_mask))) / mean_VM;
        l1_fd   = mean(abs(a.Var_delta_fd(null_mask)   - a.Var_MC(null_mask))) / mean_VM;
        fprintf('%-14s  %8.1f  %22.2f%%  %22.2f%%\n', dist_labels{k}, sigma_vals(k), ...
                100*l1_anal, 100*l1_fd);
    catch me
        fprintf('%-14s  %8.1f  [error: %s]\n', dist_labels{k}, sigma_vals(k), me.message);
    end
end

fprintf('\n==================================================================\n');
fprintf('  Done.  Paste the values above into the dagger cells in findings.tex.\n');
fprintf('==================================================================\n\n');


function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
