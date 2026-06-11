%% run_Comparison.m — Compare MC vs Delta and MC vs LN3
%
% Must run AFTER run_Delta.m, run_MC.m, and run_LN3.m.
%
% For each scenario × distribution, produces 10 PNGs:
%
%   EX_diff_MC_Delta.png     — MC_EX − Delta E[TL]  (redblue colormap)
%   EX_diff_MC_LN3.png       — MC_EX − LN3_EX       (redblue)
%   Var_diff_MC_Delta.png    — MC_Var − Delta Var    (redblue)
%   Var_diff_MC_LN3.png      — MC_Var − LN3_Var      (redblue)
%   shadow_95_MC_vs_Delta.png — side-by-side MC empirical | Delta Cheb at 95%
%   shadow_95_MC_vs_LN3.png  — side-by-side MC empirical | LN3 at 95%
%   shadow_90_MC_vs_Delta.png — same at 90%
%   shadow_90_MC_vs_LN3.png
%   shadow_80_MC_vs_Delta.png — same at 80%
%   shadow_80_MC_vs_LN3.png
%
% Note: No combined all-methods figure is generated.

function run_Comparison(sc_target, dist_target)
%% run_Comparison — MC vs Delta comparison.  No args = all combos.
if nargin < 2, sc_target = ''; dist_target = ''; end
close all; clc; warning('off');
try cd(fileparts(mfilename('fullpath'))); catch; end

addpath(genpath('Shared_Utils'));
addpath(genpath('Bellhop'));

cfg = loadConfig();

FOM = cfg.nominal.FOM_dB;
THRESHOLDS = [0.50, 0.60, 0.70, 0.80, 0.90, 0.95];



for si = 1:numel(cfg.scenarios)
    sc = cfg.scenarios(si);

    for di = 1:numel(cfg.distributions)
        dist = cfg.distributions(di);
        if ~isempty(sc_target) && ~(strcmp(sc.name,sc_target) && strcmp(dist.name,dist_target))
            continue;
        end

        mc_dir    = fullfile('Methods','MC',          sc.name, dist.name, 'results');
        delta_dir = fullfile('Methods','Delta',       sc.name, dist.name, 'results');
        fig_dir   = fullfile('Methods','Comparison',  sc.name, dist.name, 'figures');
        res_dir   = fullfile('Methods','Comparison',  sc.name, dist.name, 'results');

        % Skip if primary output already exists
        if isfile(fullfile(fig_dir, 'EX_diff_MC_Delta.png')) && ...
           isfile(fullfile(fig_dir, 'EX_diff_MC_LN3.png'))
            fprintf('  SKIP Comparison (already done): %s / %s\n', sc.name, dist.name);
            continue;
        end

        fprintf('\n=== Comparison | %s | %s ===\n', sc.name, dist.name);

        for d = {fig_dir, res_dir}
            if ~exist(d{1},'dir'), mkdir(d{1}); end
        end

        %% Determine which subsets to compare:
        %   Prefer zS_freq_svp (most informative).
        %   Fall back to each available single-variable subset when all-3 is missing.
        PREFERRED = 'zS_freq_svp';
        FALLBACK  = {'zS', 'freq', 'svp'};

        pref_mc    = fullfile(mc_dir,    sprintf('MC_%s.mat',    PREFERRED));
        pref_delta = fullfile(delta_dir, sprintf('delta_%s.mat', PREFERRED));

        if isfile(pref_mc) && isfile(pref_delta)
            compare_subsets = {PREFERRED};
        else
            compare_subsets = {};
            for ki = 1:numel(FALLBACK)
                fb = FALLBACK{ki};
                if isfile(fullfile(mc_dir,    sprintf('MC_%s.mat',    fb))) && ...
                   isfile(fullfile(delta_dir, sprintf('delta_%s.mat', fb)))
                    compare_subsets{end+1} = fb; %#ok<AGROW>
                end
            end
        end

        if isempty(compare_subsets)
            fprintf('  SKIP — required files missing for %s / %s\n', sc.name, dist.name);
            continue;
        end

        n_figs_total = 0;

        for csi = 1:numel(compare_subsets)
            sn         = compare_subsets{csi};
            mc_file    = fullfile(mc_dir,    sprintf('MC_%s.mat',    sn));
            delta_file = fullfile(delta_dir, sprintf('delta_%s.mat', sn));

            M = load(mc_file,    'MC_EX','MC_Var','MC_PrFOM','Cheb_lb','LN3_prob', ...
                                 'r_km','z_m','FOM','N');   % TL_all loaded lazily below
            D = load(delta_file, 'TL_expected','Var_TL','Cheb_lb_s','r_km','z_m','DeltaLN3_prob');

            r_km = M.r_km;  z_m = M.z_m;
            N    = M.N;

            %% LN3 derived moments — load cache if available, else compute once and save
            ln3_cache = fullfile(res_dir, sprintf('LN3_moments_%s.mat', sn));
            if isfile(ln3_cache)
                tmp = load(ln3_cache, 'LN3_EX', 'LN3_Var');
                LN3_EX  = tmp.LN3_EX;
                LN3_Var = tmp.LN3_Var;
                clear tmp;
                fprintf('  Loaded LN3 cache [%s]\n', sn);
            else
                fprintf('  Computing LN3 moments map [%s]...\n', sn);
                tmp_tl = load(mc_file, 'TL_all');
                [Nz, Nr, ~] = size(tmp_tl.TL_all);
                TL_pix = reshape(permute(tmp_tl.TL_all, [3 1 2]), N, Nz*Nr);
                clear tmp_tl;
                [ex_v, var_v] = ln3moments(TL_pix, FOM);
                clear TL_pix;
                LN3_EX  = reshape(ex_v,  Nz, Nr);
                LN3_Var = reshape(var_v, Nz, Nr);
                save(ln3_cache, 'LN3_EX', 'LN3_Var', '-v7.3');
            end

            tag = sprintf('%s | %s | %s', sc.name, dist.name, sn);

            % Suffix for filenames — blank for preferred subset (backwards compat)
            if strcmp(sn, PREFERRED)
                sfx = '';
            else
                sfx = sprintf('_%s', sn);
            end

            %% ── DIFFERENCE MAPS ─────────────────────────────────────────────
            pairings = {
                ['EX_diff_MC_Delta'  sfx],  M.MC_EX,  D.TL_expected, 'E[TL]: MC − Delta';
                ['EX_diff_MC_LN3'   sfx],   M.MC_EX,  LN3_EX,        'E[TL]: MC − LN3';
                ['Var_diff_MC_Delta' sfx],  M.MC_Var, D.Var_TL,      'Var[TL]: MC − Delta';
                ['Var_diff_MC_LN3'  sfx],   M.MC_Var, LN3_Var,       'Var[TL]: MC − LN3';
            };

            for pi = 1:size(pairings,1)
                fname  = pairings{pi,1};
                map_a  = pairings{pi,2};
                map_b  = pairings{pi,3};
                ttl    = pairings{pi,4};

                diff_map = map_a - map_b;
                mx = max(abs(diff_map(:)));
                if mx == 0, mx = 1; end

                fig = figure('Position',[50 50 800 520]);
                pcolor(r_km, z_m, diff_map);
                shading interp; set(gca,'YDir','reverse');
                colormap(redblue(256)); colorbar; clim([-mx mx]);
                xlabel('Range (km)'); ylabel('Depth (m)');
                title(sprintf('%s | %s', ttl, tag), 'Interpreter','none');
                saveFigPNG(fig, fullfile(fig_dir, fname));
                drawnow; close(fig);
            end

            %% ── SHADOW ZONE COMPARISONS (at 80%, 90%, 95%) ──────────────────
            % Include Delta-LN3 if available in the delta results file
            if isfield(D, 'DeltaLN3_prob')
                shadow_pairings = {
                    ['MC_vs_Delta'    sfx], M.MC_PrFOM, D.Cheb_lb_s,    'MC Empirical', 'Delta Chebyshev';
                    ['MC_vs_LN3'     sfx], M.MC_PrFOM, M.LN3_prob,     'MC Empirical', 'LN3 (MC-fit)';
                    ['MC_vs_DeltaLN3' sfx], M.MC_PrFOM, D.DeltaLN3_prob,'MC Empirical', 'Delta-LN3';
                };
            else
                shadow_pairings = {
                    ['MC_vs_Delta' sfx], M.MC_PrFOM, D.Cheb_lb_s, 'MC Empirical', 'Delta Chebyshev';
                    ['MC_vs_LN3'  sfx], M.MC_PrFOM, M.LN3_prob,  'MC Empirical', 'LN3 (MC-fit)';
                };
            end

            for thr_idx = 1:numel(THRESHOLDS)
                thr     = THRESHOLDS(thr_idx);
                thr_pct = round(thr * 100);

                for spi = 1:size(shadow_pairings,1)
                    pair_name = shadow_pairings{spi,1};
                    prob_A    = shadow_pairings{spi,2};
                    prob_B    = shadow_pairings{spi,3};
                    lbl_A     = shadow_pairings{spi,4};
                    lbl_B     = shadow_pairings{spi,5};

                    fig = figure('Position',[50 50 1300 520]);

                    ax1 = subplot(1,2,1);
                    shadowCategoryMap(ax1, r_km, z_m, prob_A, ...
                        sprintf('%s — %d%%', lbl_A, thr_pct), THRESHOLDS);

                    ax2 = subplot(1,2,2);
                    shadowCategoryMap(ax2, r_km, z_m, prob_B, ...
                        sprintf('%s — %d%%', lbl_B, thr_pct), THRESHOLDS);

                    sgtitle(sprintf('Shadow Zones %d%% | %s | %s', thr_pct, tag, pair_name), ...
                            'Interpreter','none');

                    fname = sprintf('shadow_%d_%s', thr_pct, pair_name);
                    saveFigPNG(fig, fullfile(fig_dir, fname));
                    drawnow; close(fig);
                end
            end

            n_figs_total = n_figs_total + numel(THRESHOLDS) * size(shadow_pairings,1) + 4;
        end

        fprintf('  %d figures saved to %s\n', n_figs_total, fig_dir);

        %% ── 3×2 Variance grid: MC (row 1) vs Delta (row 2), cols = zS / freq / svp ──
        % Only generated when all 3 single-variable subsets exist for both methods.
        SINGLE_SUBS  = {'zS',   'freq',  'svp'};
        SINGLE_LBLS  = {'z_S only', 'Freq only', 'SVP only'};
        mc_var_maps    = cell(1,3);
        delta_var_maps = cell(1,3);
        have_all = true;
        for ki = 1:3
            f_mc    = fullfile(mc_dir,    sprintf('MC_%s.mat',    SINGLE_SUBS{ki}));
            f_delta = fullfile(delta_dir, sprintf('delta_%s.mat', SINGLE_SUBS{ki}));
            if ~isfile(f_mc) || ~isfile(f_delta)
                have_all = false; break;
            end
            try
                Vm = load(f_mc,    'MC_Var','r_km','z_m');
                Vd = load(f_delta, 'Var_TL');
            catch
                have_all = false; break;
            end
            mc_var_maps{ki}    = Vm.MC_Var;
            delta_var_maps{ki} = Vd.Var_TL;
            r_km_s = Vm.r_km;  z_m_s = Vm.z_m;
        end

        if have_all
            % shared colorscale across all 6 panels
            all_v = [cellfun(@(m) max(m(:)), mc_var_maps), ...
                     cellfun(@(m) max(m(:)), delta_var_maps)];
            v_max_s = max(all_v);
            if v_max_s == 0, v_max_s = 1; end

            bathy_s = bathymetryMaker(sc.bathy_type, sc.maxR_m);

            fig3x2 = figure('Position',[50 50 1200 640]);
            row_lbls = {'MC  Var[TL]', 'Delta  Var[TL]'};
            for row = 1:2
                maps = {mc_var_maps, delta_var_maps};
                for col = 1:3
                    ax = subplot(2, 3, (row-1)*3 + col);
                    pcolor(ax, r_km_s, z_m_s, maps{row}{col});
                    shading(ax,'interp');
                    set(ax,'YDir','reverse');
                    colormap(ax, hot(256));
                    clim(ax, [0, v_max_s]);
                    axis(ax,'tight');
                    if col == 3
                        cb = colorbar(ax,'eastoutside');
                        cb.Label.String = 'dB²';
                        cb.FontSize = 7;
                    end
                    overlayBathymetry(ax, bathy_s, sc.maxDepth_m);
                    if row == 1
                        title(ax, SINGLE_LBLS{col}, 'FontSize',9,'Interpreter','none');
                    end
                    if col == 1
                        ylabel(ax, row_lbls{row}, 'FontSize',9,'Interpreter','none');
                    else
                        set(ax,'YTickLabel',[]);
                    end
                    if row == 2
                        xlabel(ax,'Range (km)','FontSize',8);
                    else
                        set(ax,'XTickLabel',[]);
                    end
                end
            end
            sgtitle(sprintf('Var[TL] — MC vs Delta | %s | %s | shared colorscale', ...
                            sc.name, dist.name), 'FontSize',11, 'Interpreter','none');
            saveFigPNG(fig3x2, fullfile(fig_dir, 'var_MC_vs_Delta_3x2'));
            drawnow; close(fig3x2);
            fprintf('  var_MC_vs_Delta_3x2.png saved\n');
        end
    end
end


fprintf('\n=== run_Comparison.m complete. ===\n');
end   % function run_Comparison
