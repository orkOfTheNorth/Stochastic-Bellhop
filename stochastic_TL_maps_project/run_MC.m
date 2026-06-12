function run_MC(sc_target, dist_target)
%% run_MC — Monte Carlo UQ with Importance Sampling recycling.
%
%  Strategy:
%   - Run 300 Bellhop samples from Normal_10pct (broadest Normal support).
%   - Derive Normal_5pct and Normal_1pct statistics for FREE via IS reweighting.
%   - Detection probability = P(TL < FOM)  [inverted from previous convention].
%
%  Usage:
%   run_MC()                        — all scenarios
%   run_MC('deep_water','')         — one scenario, all distributions
%   run_MC('deep_water','Normal_5pct') — one combo (IS-derived; fast)
%
if nargin < 2, sc_target = ''; dist_target = ''; end

close all; clc;
warning('off', 'MATLAB:unknownObjectIEEE');
warning('off', 'MATLAB:singularMatrix');
warning('off', 'MATLAB:rankDeficientMatrix');
try, cd(fileparts(mfilename('fullpath'))); catch; end

set(0, 'DefaultFigureVisible', 'off');
try
    if isempty(gcp('nocreate')), parpool('local', 10); end
catch; end

addpath(genpath('Shared_Utils'));
addpath(genpath('Bellhop'));

cfg            = loadConfig();
[snames, slbls, active] = subsetDefs();
IS_BASE        = cfg.IS_base_dist;      % 'Normal_10pct'
N              = cfg.N_MC;              % 300
THRESHOLDS     = cfg.thresholds(:)';
freq0          = cfg.nominal.freq_Hz;
zS0            = cfg.nominal.zS_m;
geo            = cfg.nominal.geo(:)';

%% ── OUTER LOOP: scenarios ────────────────────────────────────────────────────
for si = 1:numel(cfg.scenarios)
    sc  = cfg.scenarios(si);
    FOM = getFOM(cfg, sc.name);

    if ~isempty(sc_target) && ~strcmp(sc.name, sc_target), continue; end

    fprintf('\n╔══════════════════════════════════════════╗\n');
    fprintf('  Scenario: %s  |  FOM = %d dB\n', sc.name, FOM);
    fprintf('╚══════════════════════════════════════════╝\n');

    raw_cache = fullfile('Cache', sc.name, 'bellhop_raw');
    if ~exist(raw_cache,'dir'), mkdir(raw_cache); end

    %% ── Get grid size from one nominal run ───────────────────────────────────
    sim_nom = {freq0, sc.maxR_m, zS0, 0, 0, 'summer', sc.bathy_type, geo, FOM};
    [TL_dummy, r_grid, z_grid] = bellhopCached(sim_nom, raw_cache);
    r_km = r_grid / 1000;  z_m = z_grid;
    [Nz, Nr] = size(TL_dummy);  clear TL_dummy;
    bathy_m   = bathymetryMaker(sc.bathy_type, sc.maxR_m);
    max_depth = sc.maxDepth_m;

    %% ── STEP 1: Run Bellhop for IS base distribution (Normal_10pct, 300 runs) ─
    base_dist = cfg.distributions(strcmp({cfg.distributions.name}, IS_BASE));
    if isempty(base_dist)
        error('IS base distribution "%s" not found in config.', IS_BASE);
    end
    B_base = computeVarianceBounds(cfg, base_dist);
    S_base = lhsSample(N, cfg, base_dist, B_base);   % S_base.freq/zS/svp [N×1]

    % Perturbation values (zero-mean): used for IS weight computation
    x_base.freq = S_base.freq - freq0;
    x_base.zS   = S_base.zS   - zS0;
    x_base.svp  = S_base.svp;

    base_cache_dir = fullfile('Cache', sc.name, IS_BASE);
    if ~exist(base_cache_dir,'dir'), mkdir(base_cache_dir); end

    fprintf('\n--- Running IS base: %s ---\n', IS_BASE);
    [TL_base, x_base] = runOrLoadTLcache(sc, base_cache_dir, raw_cache, ...
        snames, active, S_base, x_base, freq0, zS0, geo, FOM, N, Nz, Nr, ...
        max_depth, cfg, base_dist);

    %% ── STEP 2: Compute & save stats for ALL distributions ───────────────────
    for di = 1:numel(cfg.distributions)
        dist = cfg.distributions(di);
        if ~isempty(sc_target) && ~isempty(dist_target)
            if ~strcmp(sc.name,sc_target) || ~strcmp(dist.name,dist_target)
                continue;
            end
        elseif ~isempty(sc_target) && ~strcmp(sc.name,sc_target)
            continue;
        end

        res_dir = fullfile('Methods','MC', sc.name, dist.name, 'results');
        fig_dir = fullfile('Methods','MC', sc.name, dist.name, 'figures');
        for d = {res_dir, fig_dir}
            if ~exist(d{1},'dir'), mkdir(d{1}); end
        end

        B_dist    = computeVarianceBounds(cfg, dist);
        is_base   = strcmp(dist.name, IS_BASE);

        fprintf('\n  --- Distribution: %s (IS=%d) ---\n', dist.name, ~is_base);

        all_stats = struct();

        for s = 1:numel(snames)
            fl = active{s};  sn = snames{s};
            if sum(fl) > 1, continue; end   % skip mixed subsets

            mc_file = fullfile(res_dir, sprintf('MC_%s.mat', sn));
            if isfile(mc_file)
                fprintf('    [SKIP] %s already done\n', sn);
                tmp = load(mc_file, 'MC_EX','MC_Var','MC_P_detect','Cheb_detect','LN3_prob');
                all_stats.(sn) = tmp;
                continue;
            end

            %% Determine which parameter is perturbed
            if fl(1),     param = 'freq'; sig_large = B_base.sig_freq; sig_tgt = B_dist.sig_freq; x_samp = x_base.freq;
            elseif fl(2), param = 'zS';   sig_large = B_base.sig_zS;   sig_tgt = B_dist.sig_zS;   x_samp = x_base.zS;
            else,         param = 'svp';  sig_large = B_base.sig_svp;  sig_tgt = B_dist.sig_svp;  x_samp = x_base.svp;
            end

            if ~isfield(TL_base, sn) || isempty(TL_base.(sn))
                fprintf('    [WARN] No TL cache for %s/%s\n', IS_BASE, sn);
                continue;
            end
            TL_all = double(TL_base.(sn));   % [Nz × Nr × N]

            %% Compute stats (uniform weights for base, IS weights for others)
            if is_base
                MC_EX       = mean(TL_all, 3);
                MC_Var      = var(TL_all, 0, 3);
                MC_P_detect = mean(TL_all < FOM, 3);
                ESS         = N;
                w_norm      = ones(N,1) / N;
            else
                [st, ESS] = importanceSample(TL_all, x_samp, sig_large, sig_tgt, FOM);
                MC_EX       = st.EX;
                MC_Var      = st.Var;
                MC_P_detect = st.P_detect;
                w_norm      = st.w_norm;
                fprintf('    %s: ESS=%.0f/%.0f (%.0f%%)\n', ...
                        sn, ESS, N, 100*ESS/N);
            end

            %% Chebyshev detection upper bound: P_detect ≤ 1 - Cheb_lb_shadow
            Cheb_shadow = chebyshevBound(MC_EX, MC_Var, FOM);
            Cheb_detect = 1 - Cheb_shadow;   % upper bound on P_detect

            %% LN3 probability (using IS-weighted samples via weighted KDE trick)
            fprintf('    Computing LN3 map (%s)...\n', sn);
            t_ln3 = tic;
            TL_pix   = reshape(permute(TL_all, [3 1 2]), N, Nz*Nr);
            [~, ~, prob_vec] = ln3moments(TL_pix, FOM, w_norm);
            LN3_prob = reshape(prob_vec, Nz, Nr);
            clear TL_pix prob_vec;
            fprintf('    LN3 done in %.1fs\n', toc(t_ln3));

            %% Save result
            save(mc_file, 'MC_EX','MC_Var','MC_P_detect','Cheb_detect', ...
                 'LN3_prob','r_grid','z_grid','r_km','z_m','FOM','N','ESS', ...
                 'B_dist','w_norm','-v7.3');
            fprintf('    Saved: %s\n', mc_file);

            all_stats.(sn) = struct('MC_EX',MC_EX,'MC_Var',MC_Var, ...
                'MC_P_detect',MC_P_detect,'Cheb_detect',Cheb_detect, ...
                'LN3_prob',LN3_prob);
        end

        %% ── Figures ──────────────────────────────────────────────────────────
        plotMCFigures(all_stats, snames, slbls, active, fig_dir, ...
                      r_km, z_m, bathy_m, max_depth, sc, dist, N, FOM, THRESHOLDS);

        %% ── Param-distribution figure ────────────────────────────────────────
        plotParamDistributions(S_base, x_base, B_base, B_dist, dist, IS_BASE, ...
                               is_base, N, fig_dir);
    end
end

fprintf('\n=== run_MC complete ===\n');
end   % run_MC


%% ═══════════════════════════════════════════════════════════════════════════
function [TL_base, x_base] = runOrLoadTLcache(sc, cache_dir, raw_cache, ...
    snames, active, S_base, x_base, freq0, zS0, geo, FOM, N, Nz, Nr, ...
    max_depth, cfg, dist)
%RUNORLOADTLCACHE  Run or load the N=300 Bellhop cache for the IS base dist.
    TL_base = struct();
    for s = 1:numel(snames)
        fl = active{s};  sn = snames{s};
        if sum(fl) > 1, continue; end

        tl_file = fullfile(cache_dir, sprintf('TL_%s.mat', sn));
        x_file  = fullfile(cache_dir, sprintf('x_%s.mat',  sn));

        if isfile(tl_file) && isfile(x_file)
            tc = load(tl_file, 'N_saved');
            if tc.N_saved >= N
                tmp = load(tl_file, 'TL_save');
                TL_base.(sn) = double(tmp.TL_save(:,:,1:N));
                xc = load(x_file, 'x_samples');
                if fl(1),     x_base.freq = xc.x_samples;
                elseif fl(2), x_base.zS   = xc.x_samples;
                else,         x_base.svp  = xc.x_samples;
                end
                fprintf('  [CACHE HIT] %s/%s\n', dist.name, sn);
                continue;
            end
        end

        %% Run Bellhop
        do_freq = logical(fl(1));  do_zS = logical(fl(2));  do_svp = logical(fl(3));
        freq_s  = S_base.freq;  zS_s = S_base.zS;  svp_s = S_base.svp;
        btype   = sc.bathy_type;  maxR_m = sc.maxR_m;

        TL_col  = cell(N, 1);
        t_start = tic;
        fprintf('  Running N=%d for %s/%s...\n', N, dist.name, sn);

        parfor i = 1:N
            fi = freq0; if do_freq, fi = freq_s(i); end
            zi = zS0;   if do_zS,   zi = zS_s(i);   end
            sv = 0;     if do_svp,  sv = svp_s(i);  end
            sim = {fi, maxR_m, zi, 0, 0, 'summer', btype, geo, FOM};
            tmp_dir = [tempname '_blhp'];  mkdir(tmp_dir);  prev = cd(tmp_dir);
            if do_svp
                svp_i = makeSVPNoise(sv, max_depth, cfg);
                sim{6} = 'custom';
                TL_col{i} = bellhopCached(sim, raw_cache, 'CustomSVP', svp_i);
            else
                TL_col{i} = bellhopCached(sim, raw_cache);
            end
            cd(prev);  rmdir(tmp_dir, 's');
        end

        TL_all = zeros(Nz, Nr, N);
        for i = 1:N, TL_all(:,:,i) = TL_col{i}; end
        fprintf('  Done in %.0fs\n', toc(t_start));

        %% Save TL cache
        TL_save = single(TL_all);
        N_saved = int32(N);
        dist_name = dist.name;
        save(tl_file, 'TL_save','N_saved','dist_name','-v7.3');

        %% Save perturbation samples (x_samples)
        if fl(1),     x_samples = S_base.freq - freq0;
        elseif fl(2), x_samples = S_base.zS   - zS0;
        else,         x_samples = S_base.svp;
        end
        save(x_file, 'x_samples', '-v7');
        fprintf('  Saved: %s\n', tl_file);

        TL_base.(sn) = TL_all;
        if fl(1),     x_base.freq = x_samples;
        elseif fl(2), x_base.zS   = x_samples;
        else,         x_base.svp  = x_samples;
        end
    end
end


%% ═══════════════════════════════════════════════════════════════════════════
function plotMCFigures(all_stats, snames, slbls, active, fig_dir, ...
    r_km, z_m, bathy_m, max_depth, sc, dist, N, FOM, THRESHOLDS)
%PLOTMCFIGURES  Generate all MC output figures.
    avail = fieldnames(all_stats);
    if isempty(avail), return; end

    for s = 1:numel(snames)
        sn = snames{s};
        if ~isfield(all_stats,sn), continue; end
        st = all_stats.(sn);

        %% TL map
        fig = figure('Position',[50 50 780 500]);
        plotTLmap(gca, r_km, z_m, st.MC_EX, FOM);
        overlayBathymetry(gca, bathy_m, max_depth);
        title(sprintf('MC E[TL] | %s | %s | %s | N=%d', ...
              sc.name, dist.name, sn, N), 'Interpreter','none');
        saveFigPNG(fig, fullfile(fig_dir, sprintf('TL_%s', sn)));
        close(fig);

        %% Detection probability map (empirical and Chebyshev upper bound)
        figD = figure('Position',[50 50 700 480]);
        detectionCategoryMap(gca, r_km, z_m, st.MC_P_detect, ...
            sprintf('P(detect) MC | %s | %s | %s', sc.name, dist.name, sn), ...
            THRESHOLDS);
        overlayBathymetry(gca, bathy_m, max_depth);
        saveFigPNG(figD, fullfile(fig_dir, sprintf('detect_empirical_%s', sn)));
        close(figD);

        figC = figure('Position',[50 50 700 480]);
        detectionCategoryMap(gca, r_km, z_m, st.Cheb_detect, ...
            sprintf('P(detect) Cheb UB | %s | %s | %s', sc.name, dist.name, sn), ...
            THRESHOLDS);
        overlayBathymetry(gca, bathy_m, max_depth);
        saveFigPNG(figC, fullfile(fig_dir, sprintf('detect_cheb_%s', sn)));
        close(figC);
    end

    n_avail = numel(avail);
    fw = max(800, 220 * n_avail);

    %% Combined E[TL]
    fig_ex = figure('Position',[50 50 fw 420]);
    for k = 1:n_avail
        sn  = avail{k};
        idx = find(strcmp(snames, sn), 1);
        ax  = subplot(1, n_avail, k);
        pcolor(r_km, z_m, all_stats.(sn).MC_EX);
        shading interp; set(ax,'YDir','reverse');
        colormap(ax, jet); colorbar(ax); clim([50 150]);
        overlayBathymetry(ax, bathy_m, max_depth);
        xlabel(ax,'Range (km)'); if k==1, ylabel(ax,'Depth (m)'); end
        title(ax, slbls{idx}, 'FontSize',8,'Interpreter','none');
    end
    sgtitle(sprintf('MC E[TL] — %s | %s | N=%d | FOM=%ddB', sc.name, dist.name, N, FOM));
    saveFigPNG(fig_ex, fullfile(fig_dir,'combined_EX'));  close(fig_ex);

    %% Combined Var[TL]
    fig_var = figure('Position',[50 50 fw 420]);
    v_max = max(cellfun(@(k) max(all_stats.(k).MC_Var(:)), avail));
    v_max = max(v_max, 1e-6);
    for k = 1:n_avail
        sn  = avail{k};
        idx = find(strcmp(snames, sn), 1);
        ax  = subplot(1, n_avail, k);
        pcolor(r_km, z_m, all_stats.(sn).MC_Var);
        shading interp; set(ax,'YDir','reverse');
        colormap(ax, hot); colorbar(ax); clim([0 v_max]);
        overlayBathymetry(ax, bathy_m, max_depth);
        xlabel(ax,'Range (km)'); if k==1, ylabel(ax,'Depth (m)'); end
        title(ax, slbls{idx}, 'FontSize',8,'Interpreter','none');
    end
    sgtitle(sprintf('MC Var[TL] — %s | %s | N=%d', sc.name, dist.name, N));
    saveFigPNG(fig_var, fullfile(fig_dir,'combined_Var'));  close(fig_var);
end


%% ═══════════════════════════════════════════════════════════════════════════
function plotParamDistributions(S_base, x_base, B_base, B_dist, dist, IS_BASE, ...
                                is_base, N, fig_dir)
%PLOTPARAMDISTRIBUTIONS  Histogram of the 300 drawn samples vs theoretical PDFs.
    fig = figure('Position',[50 50 1200 380]);
    params     = {'freq', 'zS',  'svp'};
    param_lbls = {'Frequency perturbation (Hz)', 'Source depth perturbation (m)', 'SVP perturbation (°C)'};
    sigs_base  = [B_base.sig_freq, B_base.sig_zS, B_base.sig_svp];
    sigs_tgt   = [B_dist.sig_freq, B_dist.sig_zS, B_dist.sig_svp];
    x_samps    = {x_base.freq, x_base.zS, x_base.svp};

    for p = 1:3
        ax  = subplot(1, 3, p);
        xs  = x_samps{p};
        sb  = sigs_base(p);
        st  = sigs_tgt(p);
        xi  = linspace(-4*sb, 4*sb, 300);

        histogram(ax, xs, 30, 'Normalization','pdf', ...
                  'FaceColor',[0.6 0.8 1.0], 'FaceAlpha',0.7, 'EdgeColor','none');
        hold(ax,'on');
        plot(ax, xi, normpdf(xi, 0, sb), 'b-', 'LineWidth', 2, ...
             'DisplayName', sprintf('P: %s (σ=%.2g)', IS_BASE, sb));
        if ~is_base
            plot(ax, xi, normpdf(xi, 0, st), 'r--', 'LineWidth', 2, ...
                 'DisplayName', sprintf('Q: %s (σ=%.2g)', dist.name, st));
        end
        hold(ax,'off');
        xlabel(ax, param_lbls{p});
        ylabel(ax, 'Density');
        legend(ax, 'Location','northeast','FontSize',7);
        title(ax, sprintf('%s | N=%d', params{p}, N), 'Interpreter','none');
        grid(ax,'on');
    end
    sgtitle(sprintf('Sampled parameter distributions — %s (base: %s)', dist.name, IS_BASE));
    saveFigPNG(fig, fullfile(fig_dir, 'MC_param_distributions'));
    close(fig);
end
