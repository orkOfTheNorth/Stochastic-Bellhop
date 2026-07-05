function run_lhs_iid_convergence(sc_target, dist_target)
%% run_lhs_iid_convergence  LHS vs IID convergence — sub-sample from cache.
%
% Loads pre-existing 1000-sample caches (no new Bellhop runs):
%   LHS: Cache/<sc>/<dist>/TL_<param>.mat
%   IID: Cache/<sc>/<dist>_iid/TL_<param>.mat
%
% Sub-samples at N_VEC without replacement.
% Reference = N=1000 statistics from each cube.
% Computes spatial-mean RMSE vs reference for E[TL], Var[TL], P_shadow.
%
% Output: Methods/MC_convergence/<sc>/figures/lhs_vs_iid_<param>.png
%         (3-panel log-log: LHS blue solid, IID red dashed; NO reference line)

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
PARAM_LBL = {'תדר (freq)','עומק מקור (z_S)','SVP'};

for si = 1:numel(cfg.scenarios)
    sc = cfg.scenarios(si);
    if ~strcmp(sc.name, sc_target), continue; end
    FOM = getFOM(cfg, sc.name);

    for di = 1:numel(cfg.distributions)
        dist = cfg.distributions(di);
        if ~strcmp(dist.name, dist_target), continue; end

        fig_dir = fullfile('Methods','MC_convergence', sc.name, 'figures');
        if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

        for pi = 1:numel(PARAMS)
            param   = PARAMS{pi};
            out_png = fullfile(fig_dir, sprintf('lhs_vs_iid_%s.png', param));
            if isfile(out_png)
                fprintf('[SKIP] %s\n', out_png); continue;
            end

            % ── load caches ──────────────────────────────────────────────
            lhs_file = fullfile('Cache', sc.name, dist.name, ...
                                sprintf('TL_%s.mat', param));
            iid_file = fullfile('Cache', sc.name, [dist.name '_iid'], ...
                                sprintf('TL_%s.mat', param));

            if ~isfile(lhs_file)
                fprintf('[MISS] LHS cache: %s\n', lhs_file); continue;
            end
            if ~isfile(iid_file)
                fprintf('[MISS] IID cache: %s\n', iid_file); continue;
            end

            fprintf('Loading LHS: %s\n', lhs_file);
            d = load(lhs_file,'TL_save');
            TL_lhs = double(d.TL_save);   % [Nz x Nr x 1000]

            fprintf('Loading IID: %s\n', iid_file);
            d = load(iid_file,'TL_save');
            TL_iid = double(d.TL_save);

            N_MAX = size(TL_lhs,3);
            assert(size(TL_iid,3) == N_MAX, 'Cache size mismatch');

            % ── N=1000 references ────────────────────────────────────────
            ref_EX_lhs  = mean(TL_lhs, 3);
            ref_Var_lhs = var(TL_lhs,  0, 3);
            ref_Pd_lhs  = mean(TL_lhs < FOM, 3);

            ref_EX_iid  = mean(TL_iid, 3);
            ref_Var_iid = var(TL_iid,  0, 3);
            ref_Pd_iid  = mean(TL_iid < FOM, 3);

            nv = numel(N_VEC);
            rmse_EX_lhs  = zeros(1,nv);  rmse_EX_iid  = zeros(1,nv);
            rmse_Var_lhs = zeros(1,nv);  rmse_Var_iid = zeros(1,nv);
            rmse_Pd_lhs  = zeros(1,nv);  rmse_Pd_iid  = zeros(1,nv);

            rng(42);
            for ni = 1:nv
                N = N_VEC(ni);
                if N >= N_MAX
                    % full set → RMSE = 0 by definition
                    rmse_EX_lhs(ni)  = 0; rmse_EX_iid(ni)  = 0;
                    rmse_Var_lhs(ni) = 0; rmse_Var_iid(ni) = 0;
                    rmse_Pd_lhs(ni)  = 0; rmse_Pd_iid(ni)  = 0;
                    continue;
                end

                idx_l = randperm(N_MAX, N);
                sl    = TL_lhs(:,:,idx_l);
                rmse_EX_lhs(ni)  = rmse2(mean(sl,3),     ref_EX_lhs);
                rmse_Var_lhs(ni) = rmse2(var(sl,0,3),    ref_Var_lhs);
                rmse_Pd_lhs(ni)  = rmse2(mean(sl<FOM,3), ref_Pd_lhs);

                idx_i = randperm(N_MAX, N);
                si2   = TL_iid(:,:,idx_i);
                rmse_EX_iid(ni)  = rmse2(mean(si2,3),     ref_EX_iid);
                rmse_Var_iid(ni) = rmse2(var(si2,0,3),    ref_Var_iid);
                rmse_Pd_iid(ni)  = rmse2(mean(si2<FOM,3), ref_Pd_iid);
            end

            % ── figure ───────────────────────────────────────────────────
            fig = figure('Position',[50 50 1200 380]);
            TITLES = {'RMSE של E[TL] (dB)','RMSE של Var[TL] (dB^2)', ...
                      'RMSE של P_{shadow}'};
            lhs_data = {rmse_EX_lhs, rmse_Var_lhs, rmse_Pd_lhs};
            iid_data = {rmse_EX_iid, rmse_Var_iid, rmse_Pd_iid};

            for p = 1:3
                ax = subplot(1,3,p);
                loglog(ax, N_VEC, max(lhs_data{p},1e-8), 'b-o', ...
                       'LineWidth',1.8, 'MarkerSize',5, ...
                       'DisplayName','LHS (סטרטיפיקציה)'); hold(ax,'on');
                loglog(ax, N_VEC, max(iid_data{p},1e-8), 'r--s', ...
                       'LineWidth',1.8, 'MarkerSize',5, ...
                       'DisplayName','IID (עצמאי)');
                grid(ax,'on'); legend(ax,'Location','northeast','FontSize',8);
                xlabel(ax,'N (גודל מדגם)'); ylabel(ax, TITLES{p});
                title(ax, TITLES{p}, 'Interpreter','tex');
                xlim(ax,[N_VEC(1) N_MAX]);
            end
            sgtitle(sprintf('LHS vs IID | %s | %s | %s', ...
                sc.name, dist.name, PARAM_LBL{pi}), ...
                'Interpreter','none', 'FontSize',11);

            print(fig, out_png, '-dpng', '-r150');
            close(fig);
            fprintf('Saved: %s\n', out_png);
        end
    end
end
fprintf('run_lhs_iid_convergence done.\n');
end

%% helper
function v = rmse2(A, B)
    v = sqrt(mean((A(:)-B(:)).^2));
end
