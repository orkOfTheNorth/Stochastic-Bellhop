%% run_MC_convergence.m — MC estimator convergence diagnostics
%
% For each scenario × param × distribution, sub-samples the N=150 direct
% Bellhop runs and shows how E[TL] and Var[TL] converge as N increases.
%
%   N_vec = [10 20 30 50 75 100 150]
% sub-samples drawn WITHOUT replacement from the 150 base samples.
%
% Outputs per scenario:
%   Methods/MC_convergence/<scen>/figures/
%     mc_conv_scatter_<param>_<dist>.png  — pixel scatter vs N (all N values)
%     mc_conv_mean_<param>.png            — mean EX and Var vs N, all 3 dists
%     mc_conv_deriv_<param>.png           — d(mean)/dN vs N, all 3 dists

function run_MC_convergence(sc_target)
if nargin < 1, sc_target = ''; end
close all;
warning('off', 'MATLAB:unknownObjectIEEE');
warning('off', 'MATLAB:singularMatrix');
warning('off', 'MATLAB:rankDeficientMatrix');
ROOT = fileparts(fileparts(mfilename('fullpath')));
try; cd(ROOT); catch; end

addpath(genpath(fullfile(ROOT, 'core')));
addpath(fullfile(ROOT, 'binaries'));

cfg = loadConfig();

N_VEC  = [10 20 30 50 75 100 150];
PARAMS = {'freq', 'zS', 'svp'};
PARAM_LBLS = {'Frequency', 'Source depth z_S', 'SVP shift'};

rng(cfg.MC.rng_seed);   % reproducible sub-sampling

for si = 1:numel(cfg.scenarios)
    sc = cfg.scenarios(si);
    if ~isempty(sc_target) && ~strcmp(sc.name, sc_target), continue; end

    fig_dir = fullfile('Methods','MC_convergence', sc.name, 'figures');
    if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

    if isfile(fullfile(fig_dir, sprintf('mc_conv_mean_%s.png', PARAMS{end})))
        fprintf('  SKIP MC convergence %s (already done)\n', sc.name);
        continue;
    end

    fprintf('\n=== MC Convergence | %s ===\n', sc.name);

    n_dists  = numel(cfg.distributions);
    colors   = lines(n_dists);

    for pi = 1:numel(PARAMS)
        param = PARAMS{pi};

        % Storage across all distributions
        mean_EX_conv  = NaN(numel(N_VEC), n_dists);
        mean_Var_conv = NaN(numel(N_VEC), n_dists);

        Nz_ref = 0;  Nr_ref = 0;

        for di = 1:n_dists
            dist_name = cfg.distributions(di).name;
            tl_file = fullfile('Cache', sc.name, dist_name, sprintf('TL_%s.mat', param));
            if ~isfile(tl_file)
                fprintf('  SKIP %s/%s — TL cache missing\n', dist_name, param);
                continue;
            end
            tmp      = load(tl_file, 'TL_save');
            TL_all   = double(tmp.TL_save);   % [Nz × Nr × N]
            [Nz, Nr, N_full] = size(TL_all);
            if Nz_ref == 0, Nz_ref = Nz; Nr_ref = Nr; end
            clear tmp;

            fprintf('  [%s/%s] %d×%d×%d samples\n', dist_name, param, Nz, Nr, N_full);

            for ni = 1:numel(N_VEC)
                Ni    = min(N_VEC(ni), N_full);
                idx_s = sort(randperm(N_full, Ni));
                TL_sub = TL_all(:,:,idx_s);

                EX_ni  = mean(TL_sub, 3);
                Var_ni = var(TL_sub, 0, 3);

                mean_EX_conv(ni, di)  = mean(EX_ni(:));
                mean_Var_conv(ni, di) = mean(Var_ni(:));
            end

            %% Fig: Var scatter vs N (reference = full N from MC results)
            mc_ref_file = fullfile('Methods','MC', sc.name, dist_name, 'results', ...
                                   sprintf('MC_%s.mat', param));
            have_ref = isfile(mc_ref_file);
            if have_ref
                mc_ref   = load(mc_ref_file, 'MC_Var');
                Var_true = mc_ref.MC_Var;
            end

            fig_sc = figure('Position', [20 20 numel(N_VEC)*280 400]);
            for ni = 1:numel(N_VEC)
                Ni    = min(N_VEC(ni), N_full);
                idx_s = sort(randperm(N_full, Ni));
                TL_sub = TL_all(:,:,idx_s);
                Var_est = var(TL_sub, 0, 3);

                ax = subplot(1, numel(N_VEC), ni);
                Vest = reshape(Var_est, [], 1);
                if have_ref
                    Vref = Var_true(:);
                    scatter(ax, Vref, Vest, 2, [0.2 0.4 0.8], 'filled', ...
                            'MarkerFaceAlpha', 0.2);
                    hold(ax,'on');
                    vmax = max(max(Vref), max(Vest));
                    plot(ax, [0 vmax], [0 vmax], 'r-', 'LineWidth',1.2);
                    xlabel(ax, 'Var_{full}','FontSize',7);
                    ylabel(ax, 'Var_{sub}','FontSize',7);
                else
                    histogram(ax, Vest, 20, 'Normalization','pdf');
                    xlabel(ax, 'Var','FontSize',7);
                end
                title(ax, sprintf('N=%d', N_VEC(ni)), 'FontSize',8);
                axis(ax,'tight');
            end
            sgtitle(sprintf('MC Var scatter | %s | %s | %s', sc.name, param, dist_name), ...
                    'Interpreter','none','FontSize',10);
            saveFigPNG(fig_sc, fullfile(fig_dir, sprintf('mc_conv_scatter_%s_%s', param, dist_name)));
            drawnow; close(fig_sc);

            clear TL_all;
        end

        dist_lbls = {cfg.distributions.name};

        %% Fig: Mean EX and Var vs N for all distributions
        fig_mv = figure('Position', [50 50 900 380]);
        ax_ex  = subplot(1,2,1);
        ax_var = subplot(1,2,2);
        for di = 1:n_dists
            if all(isnan(mean_EX_conv(:,di))), continue; end
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

        %% Fig: Derivative d(mean)/dN
        fig_dv = figure('Position', [50 50 900 380]);
        ax_dex = subplot(1,2,1);
        ax_dvr = subplot(1,2,2);
        for di = 1:n_dists
            if all(isnan(mean_EX_conv(:,di))), continue; end
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
