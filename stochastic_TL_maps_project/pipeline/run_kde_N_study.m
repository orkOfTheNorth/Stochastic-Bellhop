function run_kde_N_study(sc_target, dist_target)
%% run_kde_N_study  KDE vs empirical P_shadow accuracy as a function of N.
%
% Loads the 1000-sample LHS cache and sub-samples at N_VEC.
% At each N compares two P_shadow estimators:
%   Empirical: fraction of samples with TL < FOM
%   KDE:       Gaussian-kernel CDF with Silverman bandwidth
%
% Both are evaluated against the N=1000 empirical reference.
% The KDE estimator has lower variance (acts as kernel interpolation),
% motivating its use as a cost-effective alternative to direct MC.
%
% Output: Methods/MC_convergence/<sc>/figures/kde_vs_N.png

if nargin < 1, sc_target  = 'deep_water';   end
if nargin < 2, dist_target = 'Normal_5pct'; end

close all;
ROOT = fileparts(fileparts(mfilename('fullpath')));
try; cd(ROOT); catch; end
addpath(genpath(fullfile(ROOT,'core')));
set(0,'DefaultFigureVisible','off');

cfg    = loadConfig();
N_VEC  = [10 20 30 50 75 100 150 200 300 500 750 1000];
PARAMS = {'freq','zS','svp'};
PARAM_LBL = {'תדר','z_S','SVP'};
COLORS = {[0.2 0.5 0.9], [0.9 0.3 0.2], [0.2 0.7 0.3]};  % blue, red, green

for si = 1:numel(cfg.scenarios)
    sc = cfg.scenarios(si);
    if ~strcmp(sc.name, sc_target), continue; end
    FOM = getFOM(cfg, sc.name);

    for di = 1:numel(cfg.distributions)
        dist = cfg.distributions(di);
        if ~strcmp(dist.name, dist_target), continue; end

        fig_dir = fullfile('Methods','MC_convergence', sc.name, 'figures');
        if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
        out_png = fullfile(fig_dir, 'kde_vs_N.png');
        if isfile(out_png)
            fprintf('[SKIP] %s\n', out_png); return;
        end

        fig = figure('Position',[50 50 1100 370]);

        for pi = 1:numel(PARAMS)
            param    = PARAMS{pi};
            lhs_file = fullfile('Cache', sc.name, dist.name, ...
                                sprintf('TL_%s.mat', param));
            if ~isfile(lhs_file)
                fprintf('[MISS] %s\n', lhs_file); continue;
            end

            fprintf('Loading %s: %s\n', param, lhs_file);
            d    = load(lhs_file,'TL_save');
            TL   = double(d.TL_save);    % [Nz x Nr x 1000]
            N_MAX = size(TL,3);
            [Nz, Nr, ~] = size(TL);

            % ── reference (N=1000 empirical) ─────────────────────────────
            ref = mean(TL < FOM, 3);     % [Nz x Nr]

            nv       = numel(N_VEC);
            rmse_emp = zeros(1,nv);
            rmse_kde = zeros(1,nv);

            rng(42);
            for ni = 1:nv
                N = N_VEC(ni);
                if N >= N_MAX
                    rmse_emp(ni) = 0; rmse_kde(ni) = 0; continue;
                end

                idx    = randperm(N_MAX, N);
                TL_sub = TL(:,:,idx);           % [Nz x Nr x N]

                % Empirical
                pd_emp     = mean(TL_sub < FOM, 3);
                rmse_emp(ni) = sqrt(mean((pd_emp(:)-ref(:)).^2));

                % KDE with Silverman's bandwidth h = 1.06*sigma*N^{-1/5}
                sig_px = std(TL_sub, 0, 3);     % [Nz x Nr]
                h      = 1.06 .* sig_px .* N^(-0.2);
                h      = max(h, 0.1);            % floor to avoid /0
                pd_kde = zeros(Nz, Nr);
                for k = 1:N
                    pd_kde = pd_kde + normcdf((FOM - TL_sub(:,:,k)) ./ h);
                end
                pd_kde = pd_kde / N;
                rmse_kde(ni) = sqrt(mean((pd_kde(:)-ref(:)).^2));
            end

            % ── subplot ──────────────────────────────────────────────────
            ax = subplot(1,3,pi);
            loglog(ax, N_VEC, max(rmse_emp,1e-8), '--o', ...
                   'Color', COLORS{pi}, 'LineWidth',1.8, 'MarkerSize',5, ...
                   'DisplayName','אמפירי (ספירה)'); hold(ax,'on');
            loglog(ax, N_VEC, max(rmse_kde,1e-8), '-s', ...
                   'Color', COLORS{pi}*0.6, 'LineWidth',2.0, 'MarkerSize',5, ...
                   'DisplayName','KDE (גאוסי)');
            grid(ax,'on'); legend(ax,'Location','northeast','FontSize',8);
            xlabel(ax,'N (גודל מדגם)');
            ylabel(ax,'RMSE של P_{shadow}');
            title(ax, sprintf('פרמטר: %s', PARAM_LBL{pi}), 'Interpreter','tex');
            xlim(ax,[N_VEC(1) N_MAX]);
        end

        sgtitle(sprintf('KDE vs. ספירה אמפירית | %s | %s\n(קו מלא = KDE, מקווקוו = אמפירי)', ...
            sc.name, dist.name), 'Interpreter','none', 'FontSize',11);

        print(fig, out_png, '-dpng', '-r150');
        close(fig);
        fprintf('Saved: %s\n', out_png);
    end
end
fprintf('run_kde_N_study done.\n');
end
