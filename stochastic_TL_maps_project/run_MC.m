function run_MC(sc_target, dist_target)
%% run_MC — Monte Carlo UQ.  Call with no args to run all combos,
%  or run_MC('baseline','Normal_10pct') to process one combo only.
if nargin < 2, sc_target = ''; dist_target = ''; end

close all; clc;
warning('off', 'MATLAB:unknownObjectIEEE');
warning('off', 'MATLAB:singularMatrix');
warning('off', 'MATLAB:rankDeficientMatrix');
try
    cd(fileparts(mfilename('fullpath')));
catch
end

set(0, 'DefaultFigureVisible', 'off');
try
    if isempty(gcp('nocreate'))
        parpool('local', 10);
    end
catch
    % Parallel Computing Toolbox unavailable — parfor runs sequentially
end

addpath(genpath('Shared_Utils'));
addpath(genpath('Bellhop'));

cfg = loadConfig();
[subset_names, subset_labels, active] = subsetDefs();

FOM   = cfg.nominal.FOM_dB;
freq0 = cfg.nominal.freq_Hz;
zS0   = cfg.nominal.zS_m;
geo   = cfg.nominal.geo(:)';
N     = cfg.MC.N;

THRESHOLDS = cfg.thresholds(:)';

%% ── LOOP: scenarios × distributions ───────────────────────────────────────
for si = 1:numel(cfg.scenarios)
    sc = cfg.scenarios(si);

    for di = 1:numel(cfg.distributions)
        dist = cfg.distributions(di);
        if ~isempty(sc_target) && ~(strcmp(sc.name,sc_target) && strcmp(dist.name,dist_target))
            continue;
        end
        B    = computeVarianceBounds(cfg, dist);
        S    = lhsSample(N, cfg, dist, B);    % struct: S.freq, S.zS, S.svp

        fprintf('\n════════════════════════════════════════\n');
        fprintf('  MC | %s | %s | N=%d\n', sc.name, dist.name, N);
        fprintf('════════════════════════════════════════\n');

        res_dir   = fullfile('Methods','MC', sc.name, dist.name, 'results');
        fig_dir   = fullfile('Methods','MC', sc.name, dist.name, 'figures');
        tl_cache  = fullfile('Cache', sc.name, dist.name);
        raw_cache = fullfile('Cache', sc.name, 'bellhop_raw');
        for d = {res_dir, fig_dir, tl_cache, raw_cache}
            if ~exist(d{1},'dir'), mkdir(d{1}); end
        end

        max_depth = sc.maxDepth_m;
        bathy_m   = bathymetryMaker(sc.bathy_type, sc.maxR_m);

        %% ── SVP ensemble figure — all N samples overlaid (one PNG per dist) ──
        svp_fig_path = fullfile(fig_dir, 'svp_ensemble');
        if ~isfile([svp_fig_path '.png'])
            svp_nom = makeSVPNoise(0, max_depth, cfg);
            z_vec   = svp_nom(:,1);
            fig_svp = figure('Position', [50 50 420 560]);
            ax_svp  = axes;
            hold(ax_svp, 'on');
            for i = 1:N
                svp_i = makeSVPNoise(S.svp(i), max_depth, cfg);
                plot(ax_svp, svp_i(:,2), z_vec, '-', ...
                     'Color', [0.4 0.6 0.9 0.12], 'LineWidth', 0.8);
            end
            plot(ax_svp, svp_nom(:,2), z_vec, 'k-', 'LineWidth', 2.5, ...
                 'DisplayName', 'Nominal');
            set(ax_svp, 'YDir', 'reverse');
            xlabel(ax_svp, 'Sound Speed (m/s)');
            ylabel(ax_svp, 'Depth (m)');
            title(ax_svp, sprintf('SVP Ensemble | %s | %s | N=%d', sc.name, dist.name, N), ...
                  'Interpreter', 'none');
            grid(ax_svp, 'on');
            saveFigPNG(fig_svp, svp_fig_path);
            drawnow; close(fig_svp);
        end

        %% ── Dummy run for grid size ────────────────────────────────────────
        sim_nom = {freq0, sc.maxR_m, zS0, 0, 0, 'summer', sc.bathy_type, geo, FOM};
        [TL_dummy, r_grid, z_grid] = bellhopCached(sim_nom, raw_cache);
        r_km = r_grid/1000;  z_m = z_grid;
        [Nz, Nr] = size(TL_dummy);
        clear TL_dummy;

        all_stats = struct();

        %% ── MAIN SUBSET LOOP ───────────────────────────────────────────────
        for s = 1:7
            fl = active{s};
            sn = subset_names{s};
            if sum(fl) > 1   % skip mixed-perturbation subsets
                fprintf('  [SKIP mixed] %s\n', subset_labels{s});
                continue;
            end
            fprintf('\n=== Subset %d/7: [%s] ===\n', s, subset_labels{s});

            tl_file = fullfile(tl_cache, sprintf('TL_%s.mat', sn));

            %% ── SKIP IF RESULT ALREADY EXISTS ──────────────────────────────
            mc_file = fullfile(res_dir, sprintf('MC_%s.mat', sn));
            if isfile(mc_file)
                fprintf('  [SKIP] %s — result already exists\n', sn);
                tmp_mc = load(mc_file, 'MC_EX','MC_Var','MC_PrFOM','Cheb_lb','LN3_prob');
                all_stats.(sn) = struct('MC_EX',tmp_mc.MC_EX, 'MC_Var',tmp_mc.MC_Var, ...
                                        'MC_PrFOM',tmp_mc.MC_PrFOM, 'Cheb_lb',tmp_mc.Cheb_lb, ...
                                        'LN3_prob',tmp_mc.LN3_prob);
                continue;
            end

            %% ── TL CACHE CHECK ─────────────────────────────────────────────
            TL_all    = zeros(Nz, Nr, N);
            use_cache = false;
            runs_todo = 1:N;

            if isfile(tl_file)
                tc = load(tl_file, 'N_saved','dist_name','level_pct');
                bounds_ok = isfield(tc,'dist_name') && ...
                            strcmp(tc.dist_name, dist.name) && ...
                            isfield(tc,'level_pct') && tc.level_pct == dist.level_pct;
                if bounds_ok && tc.N_saved >= N
                    tmp    = load(tl_file, 'TL_save');
                    TL_all = double(tmp.TL_save(:,:,1:N));
                    use_cache = true;
                    fprintf('  [CACHE HIT] %d runs loaded\n', N);
                elseif bounds_ok && tc.N_saved > 0
                    rd = load(tl_file, 'TL_save', 'runs_done');
                    if isfield(rd,'runs_done') && numel(rd.runs_done) >= N
                        runs_todo = find(~rd.runs_done(1:N));
                        TL_all    = double(rd.TL_save);
                        fprintf('  [RESUMING] %d/%d done, %d remaining\n', ...
                                tc.N_saved, N, numel(runs_todo));
                    end
                else
                    fprintf('  [CACHE MISS] bounds changed — fresh run\n');
                end
            end

            %% ── RUN BELLHOP (parallel, 10 workers) ─────────────────────────
            if ~use_cache
                freq_s  = S.freq;  zS_s = S.zS;  svp_s = S.svp;
                maxR_m  = sc.maxR_m;  btype = sc.bathy_type;
                do_freq = logical(fl(1));
                do_zS   = logical(fl(2));
                do_svp  = logical(fl(3));
                TL_col  = cell(N, 1);
                t_start = tic;

                parfor i = runs_todo
                    fi = freq0; if do_freq, fi = freq_s(i); end
                    zi = zS0;   if do_zS,   zi = zS_s(i);   end
                    sv = 0;     if do_svp,  sv = svp_s(i);  end
                    sim = {fi, maxR_m, zi, 0, 0, 'summer', btype, geo, FOM};
                    tmp  = [tempname '_blhp'];  mkdir(tmp);  prev = cd(tmp);
                    if do_svp
                        svp_i = makeSVPNoise(sv, max_depth, cfg);
                        sim{6} = 'custom';
                        TL_col{i} = bellhopCached(sim, raw_cache, 'CustomSVP', svp_i);
                    else
                        TL_col{i} = bellhopCached(sim, raw_cache);
                    end
                    cd(prev);  rmdir(tmp, 's');
                end

                for i = runs_todo, TL_all(:,:,i) = TL_col{i}; end
                fprintf('  [ALL DONE] %d runs in %.0fs\n', N, toc(t_start));
                TL_save   = single(TL_all);
                N_saved   = int32(N);
                dist_name = dist.name;
                level_pct = dist.level_pct;
                save(tl_file, 'TL_save', 'N_saved', 'dist_name', 'level_pct', '-v7.3');
                fprintf('  [SAVED] %s\n', tl_file);
            end

            %% ── STATISTICS ─────────────────────────────────────────────────
            MC_EX    = mean(TL_all, 3);
            MC_Var   = var(TL_all, 0, 3);
            MC_PrFOM = mean(TL_all > FOM, 3);
            Cheb_lb  = chebyshevBound(MC_EX, MC_Var, FOM);

            %% LN3 probability map — vectorized over all pixels
            fprintf('  Computing LN3 map...\n');
            t_ln3 = tic;
            TL_pix = reshape(permute(TL_all, [3 1 2]), N, Nz*Nr);
            [~, ~, prob_vec] = ln3moments(TL_pix, FOM);
            LN3_prob = reshape(prob_vec, Nz, Nr);
            clear TL_pix prob_vec;
            fprintf('  LN3 done in %.1fs\n', toc(t_ln3));

            %% Save results
            fname = fullfile(res_dir, sprintf('MC_%s.mat', sn));
            save(fname, 'MC_EX','MC_Var','MC_PrFOM','Cheb_lb','LN3_prob', ...
                 'TL_all','r_grid','z_grid','r_km','z_m','FOM','N', ...
                 'B', '-v7.3');
            fprintf('  Saved: %s\n', fname);

            all_stats.(sn) = struct('MC_EX',MC_EX,'MC_Var',MC_Var, ...
                                    'MC_PrFOM',MC_PrFOM,'Cheb_lb',Cheb_lb, ...
                                    'LN3_prob',LN3_prob);
        end

        %% ── FIGURES ───────────────────────────────────────────────────────

        % ---- TL maps — computed subsets only ----
        for s = 1:7
            sn  = subset_names{s};
            lbl = subset_labels{s};
            if ~isfield(all_stats, sn), continue; end
            st  = all_stats.(sn);

            fig = figure('Position',[50 50 780 500]);
            plotTLmap(gca, r_km, z_m, st.MC_EX, FOM);
            overlayBathymetry(gca, bathy_m, max_depth);
            title(sprintf('MC E[TL] | %s | %s | %s | N=%d', ...
                          sc.name, dist.name, lbl, N), 'Interpreter','none');

            % Datatip showing EX and Var
            dcm = datacursormode(fig);
            set(dcm, 'UpdateFcn', @(~,evt) tlDatatip(evt, r_km, z_m, ...
                                                       st.MC_EX, st.MC_Var));
            saveFigPNG(fig, fullfile(fig_dir, sprintf('TL_%s', sn)));
            drawnow; close(fig);
        end

        % ---- Combined E[TL] — computed subsets only ----
        avail_sn = fieldnames(all_stats);
        n_avail  = numel(avail_sn);
        fw = max(800, 220 * n_avail);

        fig_ex = figure('Position',[50 50 fw 420]);
        all_ex = cellfun(@(k) all_stats.(k).MC_EX, avail_sn, 'UniformOutput', false);
        ex_min = min(cellfun(@(m) min(m(:)), all_ex));
        ex_max = max(cellfun(@(m) max(m(:)), all_ex));
        ex_clim = [max(50, ex_min), min(150, ex_max)];

        for s = 1:n_avail
            sn  = avail_sn{s};
            idx = find(strcmp(subset_names, sn), 1);
            ax  = subplot(1, n_avail, s);
            pcolor(r_km, z_m, all_stats.(sn).MC_EX);
            shading interp; set(ax,'YDir','reverse');
            colormap(ax, jet); colorbar(ax); clim(ex_clim);
            hold(ax,'on');
            TL_sm = movmean(movmean(all_stats.(sn).MC_EX, 20, 2), 20, 1);
            contour(r_km, z_m, TL_sm, [FOM FOM], 'w-', 'LineWidth', 1.0);
            hold(ax,'off');
            overlayBathymetry(ax, bathy_m, max_depth);
            xlabel(ax,'Range (km)'); if s==1, ylabel(ax,'Depth (m)'); end
            title(ax, subset_labels{idx}, 'FontSize',8, 'Interpreter','none');
        end
        sgtitle(sprintf('MC E[TL] — %s | %s | N=%d | FOM=%ddB', ...
                        sc.name, dist.name, N, FOM));
        saveFigPNG(fig_ex, fullfile(fig_dir,'combined_EX'));
        close(fig_ex);

        % ---- Combined Var[TL] — computed subsets only ----
        fig_var = figure('Position',[50 50 fw 420]);
        all_var = cellfun(@(k) all_stats.(k).MC_Var, avail_sn, 'UniformOutput', false);
        v_max = max(cellfun(@(m) max(m(:)), all_var));
        v_max = max(v_max, 1e-6);

        for s = 1:n_avail
            sn  = avail_sn{s};
            idx = find(strcmp(subset_names, sn), 1);
            ax  = subplot(1, n_avail, s);
            pcolor(r_km, z_m, all_stats.(sn).MC_Var);
            shading interp; set(ax,'YDir','reverse');
            colormap(ax, hot); colorbar(ax); clim([0, v_max]);
            overlayBathymetry(ax, bathy_m, max_depth);
            xlabel(ax,'Range (km)'); if s==1, ylabel(ax,'Depth (m)'); end
            title(ax, subset_labels{idx}, 'FontSize',8, 'Interpreter','none');
        end
        sgtitle(sprintf('MC Var[TL] — %s | %s | N=%d', sc.name, dist.name, N));
        saveFigPNG(fig_var, fullfile(fig_dir,'combined_Var'));
        close(fig_var);

        % ---- Chebyshev + Empirical bound maps — computed subsets only ----
        for s = 1:7
            sn  = subset_names{s};
            lbl = subset_labels{s};
            if ~isfield(all_stats, sn), continue; end
            st  = all_stats.(sn);

            figC = figure('Position',[50 50 700 480]);
            shadowCategoryMap(gca, r_km, z_m, st.Cheb_lb, ...
                sprintf('Cheb | %s | %s | %s', sc.name, dist.name, lbl), ...
                THRESHOLDS);
            overlayBathymetry(gca, bathy_m, max_depth);
            saveFigPNG(figC, fullfile(fig_dir, sprintf('cheb_%s', sn)));
            drawnow; close(figC);

            figE = figure('Position',[50 50 700 480]);
            shadowCategoryMap(gca, r_km, z_m, st.MC_PrFOM, ...
                sprintf('Empirical | %s | %s | %s', sc.name, dist.name, lbl), ...
                THRESHOLDS);
            overlayBathymetry(gca, bathy_m, max_depth);
            saveFigPNG(figE, fullfile(fig_dir, sprintf('empirical_%s', sn)));
            drawnow; close(figE);
        end

        fprintf('  All figures saved: %s\n', fig_dir);
    end
end


fprintf('\n=== run_MC.m complete. ===\n');
end   % function run_MC

%% ── Datatip helper ────────────────────────────────────────────────────────
function txt = tlDatatip(evt, r_km, z_m, EX_map, Var_map)
pos = evt.Position;
[~, ri] = min(abs(r_km - pos(1)));
[~, zi] = min(abs(z_m  - pos(2)));
ex  = EX_map(zi, ri);
va  = Var_map(zi, ri);
txt = {sprintf('Range: %.2f km', r_km(ri)), ...
       sprintf('Depth: %.1f m',  z_m(zi)), ...
       '──────────────', ...
       sprintf('E[TL] = %.2f dB', ex), ...
       sprintf('Var[TL] = %.4f',  va), ...
       sprintf('σ[TL]  = %.4f dB', sqrt(va))};
end
