%% run_LN3.m — LN3 fit diagnostic + pixel-wise CDF RMSE maps
%
% For each scenario × distribution × parameter-subset:
%
%  1. 4×4 histogram grid (existing): 16 spatial sample locations,
%     each showing MC histogram + LN3 PDF + FOM line with KS statistic.
%
%  2. CDF RMSE map (new): pixel-by-pixel L2 distance between fitted
%     LN3 CDF and empirical MC CDF across 50 evaluation points.
%     Saved to results/LN3_RMSE_<param>.mat.
%
%  3. RMSE summary figure (new): per-scenario figure with 3-param ×
%     3-distribution grid of RMSE maps (LN3_RMSE_all.png).
%
% Detection convention: P(TL < FOM) = P(detect) throughout.
%
% Outputs per scenario × distribution:
%   Methods/LN3/<scen>/<dist>/figures/LN3_<param>.png
%   Methods/LN3/<scen>/<dist>/results/LN3_RMSE_<param>.mat
% Per scenario (across distributions):
%   Methods/LN3/<scen>/figures/LN3_RMSE_all.png

function run_LN3(sc_target, dist_target)
%% run_LN3 — LN3 fit diagnostic + RMSE maps.  No args = all combos.
if nargin < 2, sc_target = ''; dist_target = ''; end
close all;
warning('off', 'MATLAB:unknownObjectIEEE');
warning('off', 'MATLAB:singularMatrix');
warning('off', 'MATLAB:rankDeficientMatrix');
try
    cd(fileparts(mfilename('fullpath')));
catch
end

addpath(genpath('Shared_Utils'));

cfg = loadConfig();

PARAMS  = {'freq', 'zS', 'svp'};
N_BINS  = 20;
GRID_R  = 4;
GRID_Z  = 4;
N_EVAL  = 50;   % CDF evaluation points for RMSE

for si = 1:numel(cfg.scenarios)
    sc = cfg.scenarios(si);

    for di = 1:numel(cfg.distributions)
        dist = cfg.distributions(di);
        if ~isempty(sc_target) && ~(strcmp(sc.name,sc_target) && strcmp(dist.name,dist_target))
            continue;
        end

        FOM = getFOM(cfg, sc.name);
        fprintf('\n=== LN3 | %s | %s | FOM=%ddB ===\n', sc.name, dist.name, FOM);

        mc_res   = fullfile('Methods','MC',  sc.name, dist.name, 'results');
        fig_dir  = fullfile('Methods','LN3', sc.name, dist.name, 'figures');
        res_dir  = fullfile('Methods','LN3', sc.name, dist.name, 'results');
        for d = {fig_dir, res_dir}
            if ~exist(d{1},'dir'), mkdir(d{1}); end
        end

        for ki = 1:numel(PARAMS)
            clear TL_all r_km z_m;   % ensure each param loads its own data
            sn      = PARAMS{ki};
            mc_file = fullfile(mc_res, sprintf('MC_%s.mat', sn));

            if ~isfile(mc_file)
                fprintf('  SKIP %s — MC file not found\n', sn);
                continue;
            end

            %% ── Histogram figure (skip-if-done) ────────────────────────────
            out_fig = fullfile(fig_dir, sprintf('LN3_%s.png', sn));
            if skipIfDone(out_fig, sprintf('%s / %s / %s', sc.name, dist.name, sn))
                % histogram already done; TL_all loaded below only if RMSE needed
            else
                try
                    D = load(mc_file, 'TL_all','r_km','z_m','MC_P_detect');
                catch ME
                    fprintf('  SKIP %s — could not load file (%s)\n', sn, ME.message);
                    continue;
                end
                TL_all = D.TL_all;
                r_km   = D.r_km;
                z_m    = D.z_m;
                [Nz, Nr, N] = size(TL_all);

                ri_vec = round(linspace(Nr*0.10, Nr*0.90, GRID_R));
                zi_vec = round(linspace(Nz*0.10, Nz*0.85, GRID_Z));
                ri_vec = max(1, min(Nr, ri_vec));
                zi_vec = max(1, min(Nz, zi_vec));

                fig = figure('Position', [50 50 1600 1000]);
                t = tiledlayout(GRID_Z, GRID_R, 'TileSpacing','compact','Padding','compact');
                title(t, sprintf('LN3 Fit Quality | %s | %s | %s | N=%d | FOM=%ddB', ...
                                 sc.name, dist.name, sn, N, FOM), ...
                      'FontSize', 11, 'Interpreter','none');

                for zi_i = 1:GRID_Z
                    for ri_i = 1:GRID_R
                        zi = zi_vec(zi_i);
                        ri = ri_vec(ri_i);

                        samps  = squeeze(TL_all(zi, ri, :));
                        p_emp  = D.MC_P_detect(zi, ri);   % P(detect) from MC
                        [p_ln3, gam, mu_ln, sig_ln] = ln3fit(samps, FOM);
                        ks     = ksStatLN3(samps, gam, mu_ln, sig_ln);

                        ax = nexttile;
                        histogram(ax, samps, N_BINS, 'Normalization','pdf', ...
                                  'FaceColor',[0.7 0.85 1.0], 'EdgeColor','none');
                        hold(ax,'on');

                        x_pdf = linspace(min(samps)*0.98, max(samps)*1.02, 300);
                        sh    = x_pdf - gam;
                        sh(sh<=0) = NaN;
                        y_pdf = lognpdf(sh, mu_ln, sig_ln);
                        plot(ax, x_pdf, y_pdf, 'b-', 'LineWidth', 1.5);
                        xline(ax, FOM, 'k--', 'LineWidth', 1.2);
                        hold(ax,'off');

                        xlabel(ax,'TL (dB)','FontSize',7);
                        ylabel(ax,'pdf','FontSize',7);

                        if ks < 0.10,     tc = [0.0 0.5 0.0];
                        elseif ks < 0.15, tc = [0.8 0.5 0.0];
                        else,             tc = [0.8 0.0 0.0];
                        end
                        title(ax, sprintf('r=%.1fkm z=%.0fm\nP_{det}(LN3)=%.2f P_{det}(MC)=%.2f KS=%.3f', ...
                                          r_km(ri), z_m(zi), p_ln3, p_emp, ks), ...
                              'Color', tc, 'FontSize', 7, 'Interpreter','none');
                    end
                end

                saveFigPNG(fig, fullfile(fig_dir, sprintf('LN3_%s', sn)));
                drawnow; close(fig);
                fprintf('  Saved: LN3_%s.png\n', sn);
            end

            %% ── CDF RMSE map ────────────────────────────────────────────────
            rmse_file = fullfile(res_dir, sprintf('LN3_RMSE_%s.mat', sn));
            if isfile(rmse_file)
                fprintf('  SKIP RMSE %s (cached)\n', sn);
                continue;
            end

            % Load TL_all if not already in memory from histogram step
            if ~exist('TL_all','var')
                try
                    D = load(mc_file, 'TL_all','r_km','z_m');
                catch ME
                    fprintf('  SKIP RMSE %s — could not load (%s)\n', sn, ME.message);
                    continue;
                end
                TL_all = D.TL_all;
                r_km   = D.r_km;
                z_m    = D.z_m;
            end

            [Nz, Nr, N] = size(TL_all);
            Npix = Nz * Nr;
            fprintf('  Computing LN3 RMSE map [%s] (%d×%d pixels)...\n', sn, Nz, Nr);

            % Reshape to [N × Npix]
            TL_mat = reshape(permute(TL_all, [3 1 2]), N, Npix);

            % Vectorized LN3 fit for all pixels
            gam_v = 0.95 * min(TL_mat, [], 1);          % [1 × Npix]
            sh_v  = max(TL_mat - gam_v, 1e-6);           % [N × Npix]
            lsh_v = log(sh_v);
            mu_v  = mean(lsh_v, 1);                      % [1 × Npix]
            sig_v = max(std(lsh_v, 0, 1), 1e-6);         % [1 × Npix]

            % CDF evaluation grid: 2nd–98th percentile of all TL values
            tl_flat  = TL_mat(:);
            eval_pts = linspace(prctile(tl_flat, 2), prctile(tl_flat, 98), N_EVAL);
            clear tl_flat;

            % Empirical and LN3 CDF at each eval point — [N_EVAL × Npix]
            F_emp = zeros(N_EVAL, Npix, 'single');
            F_ln3 = zeros(N_EVAL, Npix, 'single');

            for ti = 1:N_EVAL
                t = eval_pts(ti);
                F_emp(ti,:) = single(mean(TL_mat <= t, 1));
                sh_eval = t - gam_v;                     % [1 × Npix]
                fcdf = logncdf(sh_eval, mu_v, sig_v);
                fcdf(sh_eval <= 0) = 0;
                F_ln3(ti,:) = single(fcdf);
            end

            RMSE_map = reshape(sqrt(mean((F_ln3 - F_emp).^2, 1)), Nz, Nr);
            clear TL_mat F_emp F_ln3 sh_v lsh_v;

            save(rmse_file, 'RMSE_map','r_km','z_m','FOM','N_EVAL','eval_pts');
            fprintf('  RMSE map saved → %s\n', rmse_file);

            clear TL_all;  % free memory between params
        end
    end

    %% ── Per-scenario RMSE summary figure (3 params × 3 distributions) ──────
    sc_fig_dir = fullfile('Methods','LN3', sc.name, 'figures');
    if ~exist(sc_fig_dir,'dir'), mkdir(sc_fig_dir); end
    rmse_all_fig = fullfile(sc_fig_dir, 'LN3_RMSE_all.png');

    if ~isempty(sc_target) || isfile(rmse_all_fig)
        % regenerate when running per-scenario or first time
    end

    % Collect RMSE maps across distributions and params
    n_dists  = numel(cfg.distributions);
    n_params = numel(PARAMS);
    rmse_maps = cell(n_dists, n_params);
    have_any  = false;

    for di = 1:n_dists
        for ki = 1:n_params
            rfile = fullfile('Methods','LN3', sc.name, ...
                             cfg.distributions(di).name, 'results', ...
                             sprintf('LN3_RMSE_%s.mat', PARAMS{ki}));
            if isfile(rfile)
                tmp = load(rfile, 'RMSE_map','r_km','z_m');
                rmse_maps{di,ki} = tmp.RMSE_map;
                r_km_ref = tmp.r_km;
                z_m_ref  = tmp.z_m;
                have_any = true;
            end
        end
    end

    if ~have_any
        fprintf('  [RMSE all] No RMSE results yet for %s — skipping summary figure.\n', sc.name);
        continue;
    end

    bathy_m = bathymetryMaker(sc.bathy_type, sc.maxR_m);

    % Shared colorscale across all panels
    all_rmse = cellfun(@(m) max(m(:)), rmse_maps(~cellfun(@isempty, rmse_maps)));
    rmse_max = max([all_rmse; 0.01]);

    param_lbls = {'Frequency', 'Source depth z_S', 'SVP shift'};
    dist_names = {cfg.distributions.name};

    fig_rmse = figure('Position',[50 50 n_params*380 n_dists*280]);
    for di = 1:n_dists
        for ki = 1:n_params
            ax = subplot(n_dists, n_params, (di-1)*n_params + ki);
            if isempty(rmse_maps{di,ki})
                axis(ax,'off');
                text(ax, 0.5, 0.5, 'N/A', 'HorizontalAlignment','center');
                continue;
            end
            pcolor(ax, r_km_ref, z_m_ref, rmse_maps{di,ki});
            shading(ax,'interp'); set(ax,'YDir','reverse');
            colormap(ax, hot(256));
            clim(ax, [0, rmse_max]);
            if ki == n_params
                cb = colorbar(ax,'eastoutside');
                cb.Label.String = 'CDF RMSE';
                cb.FontSize = 7;
            end
            overlayBathymetry(ax, bathy_m, sc.maxDepth_m);
            if di == 1
                title(ax, param_lbls{ki}, 'FontSize',9,'Interpreter','none');
            end
            if ki == 1
                ylabel(ax, strrep(dist_names{di},'_',' '), 'FontSize',8,'Interpreter','none');
            else
                set(ax,'YTickLabel',[]);
            end
            if di == n_dists
                xlabel(ax,'Range (km)','FontSize',8);
            else
                set(ax,'XTickLabel',[]);
            end
        end
    end
    sgtitle(sprintf('LN3 CDF RMSE | %s | rows=distribution, cols=param', sc.name), ...
            'Interpreter','none', 'FontSize',11);
    saveFigPNG(fig_rmse, strrep(rmse_all_fig, '.png', ''));
    drawnow; close(fig_rmse);
    fprintf('  LN3_RMSE_all.png saved → %s\n', sc_fig_dir);
end

fprintf('\n=== run_LN3.m complete. ===\n');
end   % function run_LN3

%% ── KS statistic helper ──────────────────────────────────────────────────────
function ks = ksStatLN3(samps, gam, mu_ln, sig_ln)
samps  = sort(double(samps(:)));
n      = numel(samps);
F_emp  = (1:n)' / n;
sh     = samps - gam;
sh(sh <= 0) = 1e-9;
F_ln3  = logncdf(sh, mu_ln, sig_ln);
ks     = max(abs(F_emp - F_ln3));
end
