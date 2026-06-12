%% run_MC_convergence.m — MC estimator convergence diagnostics
%
% Uses the IS framework to sub-sample the 300 Normal_10pct base runs and
% show how EX and Var converge as N increases — for ALL THREE distributions
% (Normal_1pct, Normal_5pct, Normal_10pct) without extra Bellhop runs.
%
% For each scenario × param × distribution, evaluates statistics at
%   N_vec = [10 20 30 50 75 100 150 200 300]
% sub-samples drawn WITHOUT replacement from the 300 base samples.
%
% Outputs per scenario:
%   Methods/MC_convergence/<scen>/figures/
%     mc_conv_scatter_<param>_<dist>.png  — pixel scatter vs N (all 9 N values)
%     mc_conv_mean_<param>.png            — mean EX and Var vs N, all 3 dists
%     mc_conv_deriv_<param>.png           — d(mean)/dN vs N, all 3 dists

function run_MC_convergence(sc_target)
if nargin < 1, sc_target = ''; end
close all;
warning('off', 'MATLAB:unknownObjectIEEE');
warning('off', 'MATLAB:singularMatrix');
warning('off', 'MATLAB:rankDeficientMatrix');
try cd(fileparts(mfilename('fullpath'))); catch; end

addpath(genpath('Shared_Utils'));

cfg = loadConfig();

N_VEC  = [10 20 30 50 75 100 150 200 300];
PARAMS = {'freq', 'zS', 'svp'};
PARAM_LBLS = {'Frequency', 'Source depth z_S', 'SVP shift'};
base_name  = cfg.IS_base_dist;   % 'Normal_10pct'

base_idx = find(strcmp({cfg.distributions.name}, base_name), 1);
B_base   = computeVarianceBounds(cfg, cfg.distributions(base_idx));
sig_base = [B_base.sig_freq, B_base.sig_zS, B_base.sig_svp];

rng(cfg.MC.rng_seed);   % reproducible sub-sampling

for si = 1:numel(cfg.scenarios)
    sc = cfg.scenarios(si);
    if ~isempty(sc_target) && ~strcmp(sc.name, sc_target), continue; end

    % FOM not used here (convergence is stats-only, no detection threshold needed)
    fig_dir = fullfile('Methods','MC_convergence', sc.name, 'figures');
    if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

    % Check for completion marker
    if isfile(fullfile(fig_dir, sprintf('mc_conv_mean_%s.png', PARAMS{end})))
        fprintf('  SKIP MC convergence %s (already done)\n', sc.name);
        continue;
    end

    fprintf('\n=== MC Convergence | %s | FOM=%ddB ===\n', sc.name, FOM);

    for pi = 1:numel(PARAMS)
        param = PARAMS{pi};
        clear TL_all x_samples;

        %% Load 300-sample base TL cube
        tl_file = fullfile('Cache', sc.name, base_name, sprintf('TL_%s.mat', param));
        if ~isfile(tl_file)
            fprintf('  SKIP %s — TL cache missing\n', param);
            continue;
        end
        tmp      = load(tl_file, 'TL_save');
        TL_all   = double(tmp.TL_save);   % [Nz × Nr × 300]
        [Nz, Nr, N_full] = size(TL_all);
        clear tmp;

        % Load x_samples (perturbation values drawn from base distribution)
        x_file = fullfile('Cache', sc.name, base_name, sprintf('x_%s.mat', param));
        if isfile(x_file)
            xs      = load(x_file, 'x_samples');
            x_samples = xs.x_samples(:);   % [N_full × 1]
        else
            % Regenerate from LHS with fixed seed (reproducible)
            S_base = lhsSample(N_full, cfg, cfg.distributions(base_idx), B_base);
            switch pi
                case 1, x_samples = S_base.freq - cfg.nominal.freq_Hz;
                case 2, x_samples = S_base.zS   - cfg.nominal.zS_m;
                case 3, x_samples = S_base.svp;
            end
        end

        fprintf('  [%s] loaded %d×%d×%d, computing convergence...\n', ...
                param, Nz, Nr, N_full);

        % Pre-compute IS log-weights for all distributions (reuse per sub-sample)
        n_dists   = numel(cfg.distributions);
        log_w_all = zeros(N_full, n_dists);
        for di = 1:n_dists
            B_di = computeVarianceBounds(cfg, cfg.distributions(di));
            switch pi
                case 1, sig_di = B_di.sig_freq;
                case 2, sig_di = B_di.sig_zS;
                case 3, sig_di = B_di.sig_svp;
            end
            lw = (x_samples.^2) / 2 * (1/sig_base(pi)^2 - 1/sig_di^2);
            log_w_all(:, di) = lw - max(lw);   % stabilise
        end

        % Storage: [numel(N_VEC) × n_dists] for mean_EX and mean_Var across pixels
        mean_EX_conv  = zeros(numel(N_VEC), n_dists);
        mean_Var_conv = zeros(numel(N_VEC), n_dists);
        % [Nz × Nr × numel(N_VEC)] for scatter (one distribution: Normal_10pct)
        Var_scatter   = zeros(Nz, Nr, numel(N_VEC));

        for ni = 1:numel(N_VEC)
            Ni    = N_VEC(ni);
            idx_s = sort(randperm(N_full, min(Ni, N_full)));   % sub-sample indices

            TL_sub    = TL_all(:,:,idx_s);     % [Nz × Nr × Ni]

            for di = 1:n_dists
                lw_sub = log_w_all(idx_s, di);
                w      = exp(lw_sub - max(lw_sub));
                w_norm = w / sum(w);
                w3     = reshape(w_norm, 1, 1, Ni);

                EX_ni   = sum(TL_sub .* w3, 3);
                dTL     = TL_sub - EX_ni;
                Var_ni  = sum(w3 .* dTL.^2, 3) / (1 - sum(w_norm.^2));

                mean_EX_conv(ni, di)  = mean(EX_ni(:));
                mean_Var_conv(ni, di) = mean(Var_ni(:));

                % Store scatter data for the base distribution only
                if di == base_idx
                    Var_scatter(:,:,ni) = Var_ni;
                end
            end
        end

        % Load MC reference (full 300 sample) for scatter reference line
        mc_ref_file = fullfile('Methods','MC', sc.name, base_name, 'results', ...
                               sprintf('MC_%s.mat', param));
        if isfile(mc_ref_file)
            mc_ref   = load(mc_ref_file, 'MC_Var');
            Var_true = mc_ref.MC_Var;
            have_ref = true;
        else
            have_ref = false;
        end

        dist_lbls = {cfg.distributions.name};
        colors    = lines(n_dists);

        %% Fig 1: Scatter of pixel Var vs N (for base distribution, 9 subplots)
        fig_sc = figure('Position', [20 20 numel(N_VEC)*280 400]);
        for ni = 1:numel(N_VEC)
            ax = subplot(1, numel(N_VEC), ni);
            Vest = reshape(Var_scatter(:,:,ni), [], 1);
            if have_ref
                Vref = Var_true(:);
                scatter(ax, Vref, Vest, 2, [0.2 0.4 0.8], 'filled', ...
                        'MarkerFaceAlpha', 0.2);
                hold(ax,'on');
                vmax = max(max(Vref), max(Vest));
                plot(ax, [0 vmax], [0 vmax], 'r-', 'LineWidth',1.2);
                xlabel(ax, 'Var_{N=300}','FontSize',7);
                ylabel(ax, 'Var_{sub}','FontSize',7);
            else
                histogram(ax, Vest, 20, 'Normalization','pdf');
                xlabel(ax, 'Var','FontSize',7);
            end
            title(ax, sprintf('N=%d', N_VEC(ni)), 'FontSize',8);
            axis(ax,'tight');
        end
        sgtitle(sprintf('MC Var scatter | %s | %s | %s (base)', ...
                        sc.name, param, base_name), 'Interpreter','none','FontSize',10);
        saveFigPNG(fig_sc, fullfile(fig_dir, sprintf('mc_conv_scatter_%s_%s', param, base_name)));
        drawnow; close(fig_sc);

        %% Fig 2: Mean EX and Var vs N for all 3 distributions
        fig_mv = figure('Position', [50 50 900 380]);
        ax_ex  = subplot(1,2,1);
        ax_var = subplot(1,2,2);
        for di = 1:n_dists
            plot(ax_ex,  N_VEC, mean_EX_conv(:,di),  '-o', 'Color', colors(di,:), ...
                 'LineWidth',1.5,'MarkerSize',5,'DisplayName', strrep(dist_lbls{di},'_',' '));
            hold(ax_ex,'on');
            plot(ax_var, N_VEC, mean_Var_conv(:,di), '-o', 'Color', colors(di,:), ...
                 'LineWidth',1.5,'MarkerSize',5,'DisplayName', strrep(dist_lbls{di},'_',' '));
            hold(ax_var,'on');
        end
        xlabel(ax_ex,  'N samples'); ylabel(ax_ex,  'Mean EX (dB)');
        xlabel(ax_var, 'N samples'); ylabel(ax_var, 'Mean Var (dB²)');
        title(ax_ex,  sprintf('E[TL] vs N | %s | %s', sc.name, PARAM_LBLS{pi}), ...
              'FontSize',9,'Interpreter','none');
        title(ax_var, sprintf('Var[TL] vs N | %s | %s', sc.name, PARAM_LBLS{pi}), ...
              'FontSize',9,'Interpreter','none');
        legend(ax_ex,  'Location','best','FontSize',7);
        legend(ax_var, 'Location','best','FontSize',7);
        grid(ax_ex,'on'); grid(ax_var,'on');
        sgtitle(sprintf('MC Convergence — mean statistics | %s | %s', sc.name, param), ...
                'Interpreter','none','FontSize',11);
        saveFigPNG(fig_mv, fullfile(fig_dir, sprintf('mc_conv_mean_%s', param)));
        drawnow; close(fig_mv);

        %% Fig 3: Derivative d(mean)/dN — shows when statistics plateau
        fig_dv = figure('Position', [50 50 900 380]);
        ax_dex = subplot(1,2,1);
        ax_dvr = subplot(1,2,2);
        for di = 1:n_dists
            dEX  = diff(mean_EX_conv(:,di))  ./ diff(N_VEC(:));
            dVar = diff(mean_Var_conv(:,di))  ./ diff(N_VEC(:));
            N_mid = (N_VEC(1:end-1) + N_VEC(2:end)) / 2;
            plot(ax_dex, N_mid, dEX,  '-o', 'Color', colors(di,:), ...
                 'LineWidth',1.5,'MarkerSize',5,'DisplayName', strrep(dist_lbls{di},'_',' '));
            hold(ax_dex,'on');
            plot(ax_dvr, N_mid, dVar, '-o', 'Color', colors(di,:), ...
                 'LineWidth',1.5,'MarkerSize',5,'DisplayName', strrep(dist_lbls{di},'_',' '));
            hold(ax_dvr,'on');
        end
        yline(ax_dex, 0, 'k--', 'LineWidth',1);
        yline(ax_dvr, 0, 'k--', 'LineWidth',1);
        xlabel(ax_dex,  'N (midpoint)'); ylabel(ax_dex,  'd(EX)/dN (dB/sample)');
        xlabel(ax_dvr,  'N (midpoint)'); ylabel(ax_dvr,  'd(Var)/dN (dB²/sample)');
        title(ax_dex, sprintf('dE[TL]/dN | %s | %s', sc.name, PARAM_LBLS{pi}), ...
              'FontSize',9,'Interpreter','none');
        title(ax_dvr, sprintf('dVar/dN | %s | %s', sc.name, PARAM_LBLS{pi}), ...
              'FontSize',9,'Interpreter','none');
        legend(ax_dex, 'Location','best','FontSize',7);
        legend(ax_dvr, 'Location','best','FontSize',7);
        grid(ax_dex,'on'); grid(ax_dvr,'on');
        sgtitle(sprintf('MC Convergence — derivative | %s | %s', sc.name, param), ...
                'Interpreter','none','FontSize',11);
        saveFigPNG(fig_dv, fullfile(fig_dir, sprintf('mc_conv_deriv_%s', param)));
        drawnow; close(fig_dv);

        fprintf('  Figures saved for %s\n', param);
    end
end

fprintf('\n=== run_MC_convergence.m complete. ===\n');
end
