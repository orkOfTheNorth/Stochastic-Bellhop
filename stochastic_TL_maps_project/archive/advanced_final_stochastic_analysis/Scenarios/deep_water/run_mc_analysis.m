function run_mc_analysis(config)
% run_mc_analysis  –  Core Monte Carlo UQ analysis function.
%   Accepts a config struct and performs the full analysis for that scenario.

    %% ── UNPACK CONFIG AND SET UP PATHS ───────────────────────────────────────
    SCN_LABEL       = config.label;
    bathy_type      = config.bathy_type;
    maxR            = config.maxR;
    output_base_dir = config.output_base_dir;
    N               = config.N;
    base_dir        = config.base_dir;

    results_dir = fullfile(output_base_dir, 'results');
    figures_dir = fullfile(output_base_dir, 'figures');
    tl_cache_dir = fullfile(output_base_dir, 'TL_cache');
    if ~exist(results_dir,'dir'),  mkdir(results_dir);  end
    if ~exist(figures_dir,'dir'),  mkdir(figures_dir);  end
    if ~exist(tl_cache_dir,'dir'), mkdir(tl_cache_dir); end
    cache_dir = fullfile(base_dir, 'cache');

    %% ── PARAMETERS ───────────────────────────────────────────────────────────
    freq0 = 10000;   zS0 = 5;   FOM = 100;
    geo   = [0.989*1500, 1.63, 0.07];
    force_rerun = isfield(config, 'force_rerun') && config.force_rerun;

    freq_bnd = [freq0*0.99, freq0*1.01];
    zS_bnd   = [zS0*0.99,  zS0*1.01];
    svp_bnd  = [-0.25, 0.25];   % °C

    %% ── SUBSET DEFINITIONS ───────────────────────────────────────────────────
    subset_names  = {'zS','freq','svp','zS_freq','zS_svp','freq_svp','zS_freq_svp'};
    subset_labels = {'z_S only','Freq only','SVP only', ...
                     'z_S+Freq','z_S+SVP','Freq+SVP','z_S+Freq+SVP'};
    active = {[0 1 0],[1 0 0],[0 0 1],[1 1 0],[0 1 1],[1 0 1],[1 1 1]};

    %% ── LHS SAMPLES ──────────────────────────────────────────────────────────
    rng(42);
    lhs = lhsdesign(N, 3);
    samp_freq = freq_bnd(1) + lhs(:,1) * diff(freq_bnd);
    samp_zS   = zS_bnd(1)   + lhs(:,2) * diff(zS_bnd);
    samp_svp  = svp_bnd(1)  + lhs(:,3) * diff(svp_bnd);

    %% ── DUMMY RUN FOR GRID ───────────────────────────────────────────────────
    set(0,'DefaultFigureVisible','off');
    sim_nom = {freq0, maxR, zS0, 0, 0, "summer", bathy_type, geo, FOM};
    [TL_dummy, r_grid, z_grid] = bellhopCached(sim_nom, cache_dir);
    close all;
    r_km = r_grid / 1000;   z_m = z_grid;
    [Nz, Nr] = size(TL_dummy);
    max_depth = max(z_grid);

    %% ── SVP PROFILE INSPECTION ─────────────────────────────────────────────
    fig_svp = figure('Name','MC SVP Profiles','Visible','off','Position',[50 50 520 640]);
    svp_nom = makeSVPNoise(0, max_depth);
    hold on;
    for ii = 1:N
        svp_i = makeSVPNoise(samp_svp(ii), max_depth);
        plot(svp_i(:,2), svp_i(:,1), 'Color', [0.65 0.65 0.65], 'LineWidth', 0.6);
    end
    plot(svp_nom(:,2), svp_nom(:,1), 'k-', 'LineWidth', 2.5);
    hold off;
    set(gca, 'YDir', 'reverse');
    xlabel('Sound Speed (m/s)'); ylabel('Depth (m)');
    title(sprintf('SVP Profiles  N=%d  svp bnd=[%.2f, %.2f] deg C', N, svp_bnd(1), svp_bnd(2)));
    grid on;
    saveFig(fig_svp, fullfile(figures_dir,'svp_profiles'));
    close(fig_svp);

    %% ── MAIN LOOP ────────────────────────────────────────────────────────────
    all_stats = struct();
    for s = 1:7
        fl = active{s};
        sn = subset_names{s};
        fprintf('\n=== Subset %d/7: [%s] ===\n', s, subset_labels{s});

        tl_file   = fullfile(tl_cache_dir, sprintf('TL_%s.mat', sn));
        use_cache = false;
        runs_todo = 1:N;
        TL_all    = zeros(Nz, Nr, N);

        if isfile(tl_file) && ~force_rerun
            tc = load(tl_file, 'N_saved','freq_bnd','zS_bnd','svp_bnd');
            bounds_ok = isfield(tc,'svp_bnd') && isequal(tc.freq_bnd,freq_bnd) && ...
                        isequal(tc.zS_bnd,zS_bnd) && isequal(tc.svp_bnd,svp_bnd);
            if bounds_ok
                if tc.N_saved >= N
                    tmp    = load(tl_file,'TL_save');
                    TL_all = double(tmp.TL_save(:,:,1:N));
                    use_cache = true;
                    fprintf('  [CACHE HIT] %d/%d runs loaded\n', N, tc.N_saved);
                elseif tc.N_saved > 0
                    rd  = load(tl_file,'TL_save','runs_done');
                    if isfield(rd,'runs_done') && numel(rd.runs_done) >= N
                        runs_todo = find(~rd.runs_done(1:N));
                        TL_all    = double(rd.TL_save);
                        fprintf('  [RESUMING] %d/%d done, %d remaining\n', ...
                                tc.N_saved, N, numel(runs_todo));
                    end
                end
            else
                fprintf('  [CACHE MISS] bounds changed — fresh run\n');
            end
        end

        if ~use_cache
            m_tl = matfile(tl_file, 'Writable', true);
            if numel(runs_todo) == N
                m_tl.TL_save   = zeros(Nz, Nr, N, 'single');
                m_tl.runs_done = false(1, N);
                m_tl.freq_bnd  = freq_bnd;
                m_tl.zS_bnd    = zS_bnd;
                m_tl.svp_bnd   = svp_bnd;
                m_tl.N_saved   = int32(0);
            end
            rd_vec = m_tl.runs_done;
            n_done  = sum(rd_vec);
            t_start = tic;
            for iter = 1:numel(runs_todo)
                i = runs_todo(iter);
                fi = freq0;  if fl(1), fi = samp_freq(i); end
                zi = zS0;    if fl(2), zi = samp_zS(i);   end
                si = 0;      if fl(3), si = samp_svp(i);  end
                sim = {fi, maxR, zi, 0, 0, "summer", bathy_type, geo, FOM};
                if fl(3)
                    sim{6} = "custom";
                    [TL,~,~] = bellhopCached(sim, cache_dir, 'CustomSVP', makeSVPNoise(si, max_depth));
                else
                    [TL,~,~] = bellhopCached(sim, cache_dir);
                end
                TL_all(:,:,i) = TL;
                m_tl.TL_save(:,:,i) = single(TL);
                rd_vec(i)           = true;
                m_tl.runs_done      = rd_vec;
                n_done              = n_done + 1;
                m_tl.N_saved        = int32(n_done);
                close all;
                if mod(n_done,10)==0 || n_done==N
                    progressBar(n_done, N, sn, t_start);
                end
            end
            fprintf('  [ALL SAVED] %s  (N=%d, %.0fMB)\n', tl_file, N, Nz*Nr*N*4/1e6);
        end

        % ── Statistics ───────────────────────────────────────────────────
        MC_EX    = mean(TL_all, 3);
        MC_Var   = var(TL_all, 0, 3);
        MC_PrFOM = mean(TL_all > FOM, 3);
        d  = MC_EX - FOM;
        lb = zeros(Nz, Nr);
        mask = d > 0;
        lb(mask) = max(0, 1 - MC_Var(mask) ./ (MC_Var(mask) + d(mask).^2));
        Cheb_lb = lb;

        fprintf('  Computing LN3 shadow map...\n');
        t_ln3 = tic;
        num_pixels = Nz * Nr;
        TL_pixels = reshape(permute(TL_all, [3 1 2]), N, num_pixels);
        LN3_prob_vec = zeros(1, num_pixels);
        parfor p = 1:num_pixels
            LN3_prob_vec(p) = ln3prob(TL_pixels(:, p), FOM);
        end
        LN3_prob = reshape(LN3_prob_vec, Nz, Nr);
        fprintf('  LN3 done in %.1fs\n', toc(t_ln3));
        
        fname = fullfile(results_dir, sprintf('MC_%s.mat', sn));
        save(fname, 'MC_EX','MC_Var','MC_PrFOM','Cheb_lb','LN3_prob', ...
                    'TL_all','r_grid','z_grid','r_km','z_m','FOM','N', ...
                    'freq_bnd','zS_bnd','svp_bnd','freq0','zS0', '-v7.3');
        fprintf('  Saved %s\n', fname);

        all_stats.(sn) = struct('MC_EX',MC_EX,'MC_Var',MC_Var,'MC_PrFOM',MC_PrFOM, ...
                                'Cheb_lb',Cheb_lb,'LN3_prob',LN3_prob);
    end

    %% ── COMPOSITE SUMMARY FIGURES ────────────────────────────────────────────
    nom_file = fullfile(output_base_dir, '..', 'Delta_Method', 'results', 'delta_jacobians.mat');
    if ~isfile(nom_file) % Fallback for scenario runs
        nom_file = fullfile(base_dir, 'Delta_Method', 'results', 'delta_jacobians.mat');
    end

    if isfile(nom_file)
        tmp = load(nom_file,'TL_nom'); raw_ref = tmp.TL_nom;
    else
        raw_ref = all_stats.(subset_names{7}).MC_EX;
    end
    TL_ref = movmean(movmean(raw_ref, 20, 2), 20, 1);

    create_summary_figure('MC E[TL]', 'MC_all_subsets_mean_TL', all_stats, 'MC_EX', 'mean', r_km, z_m, N, TL_ref, FOM, SCN_LABEL, figures_dir);
    create_summary_figure('MC Var[TL]', 'MC_all_subsets_variance', all_stats, 'MC_Var', 'variance', r_km, z_m, N, [], [], SCN_LABEL, figures_dir);
    create_summary_figure('MC Chebyshev P(shadow)', 'MC_all_subsets_cheb', all_stats, 'Cheb_lb', 'shadow', r_km, z_m, N, [], [], SCN_LABEL, figures_dir);
    create_summary_figure('MC Empirical P(shadow)', 'MC_all_subsets_empirical', all_stats, 'MC_PrFOM', 'shadow', r_km, z_m, N, [], [], SCN_LABEL, figures_dir);
    create_summary_figure('MC LN3 P(shadow)', 'MC_all_subsets_ln3', all_stats, 'LN3_prob', 'shadow', r_km, z_m, N, [], [], SCN_LABEL, figures_dir);

    set(0,'DefaultFigureVisible','on');
    fprintf('\n=== Monte Carlo complete: %s ===\n', SCN_LABEL);
end

%% ── LOCAL HELPER FUNCTIONS ───────────────────────────────────────────────────
function create_summary_figure(fig_title, file_name, all_stats, field_name, plot_type, r_km, z_m, N, TL_ref, FOM, SCN_LABEL, figures_dir)
    subset_names  = {'zS','freq','svp','zS_freq','zS_svp','freq_svp','zS_freq_svp'};
    subset_labels = {'z_S only','Freq only','SVP only', ...
                     'z_S+Freq','z_S+SVP','Freq+SVP','z_S+Freq+SVP'};
    N_s = numel(subset_names);
    fw = max(1400, 200 * N_s);
    fig = figure('Name', fig_title, 'Visible', 'off', 'Position', [50 50 fw 380]);

    for s = 1:N_s
        sn = subset_names{s};
        if ~isfield(all_stats, sn), continue; end
        ax = subplot(1, N_s, s);
        data = all_stats.(sn).(field_name);

        switch plot_type
            case 'mean'
                pcolor(r_km, z_m, data); shading interp;
                colormap(ax, jet); colorbar; clim([50 150]);
                hold on; contour(r_km, z_m, TL_ref, [FOM FOM], 'w-', 'LineWidth', 0.8); hold off;
            case 'variance'
                pcolor(r_km, z_m, data); shading interp;
                colormap(ax, hot); colorbar;
            case 'shadow'
                shadowThresholdMaps(ax, r_km, z_m, data, subset_labels{s}, [0.70 0.80 0.90 0.95]);
        end
        
        if ~strcmp(plot_type, 'shadow')
            set(ax, 'YDir', 'reverse');
            xlabel('Range (km)');
            if s == 1, ylabel('Depth (m)'); end
            title(subset_labels{s}, 'FontSize', 8);
        end
    end
    
    if isempty(SCN_LABEL)
        sgtitle(sprintf('%s  N=%d', fig_title, N));
    else
        sgtitle(sprintf('%s  N=%d  |  %s', fig_title, N, SCN_LABEL));
    end
    saveFig(fig, fullfile(figures_dir, file_name));
    close(fig);
end