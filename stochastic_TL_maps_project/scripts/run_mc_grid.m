%% run_MC_grid.m — Pre-computed parameter grid for all scenarios
%
% Runs Bellhop once on a dense 1-D grid per parameter per scenario, then
% samples TL for any distribution via interpolation — no re-runs needed.
%
% Grid layout (per scenario, per parameter):
%   Cache/<scen>/grid/TL_freq_grid.mat   — [Nz × Nr × N_GRID]
%   Cache/<scen>/grid/TL_zS_grid.mat
%   Cache/<scen>/grid/TL_svp_grid.mat
%   Cache/<scen>/grid/grid_axes.mat      — freq_grid, zS_grid, svp_grid vectors
%
% MC results from the grid (per scenario × distribution):
%   Methods/MC_grid/<scen>/<dist>/results/MC_<param>.mat
%   Methods/MC_grid/<scen>/<dist>/figures/  (same figures as run_MC)
%
% Usage: run this script ONCE per scenario to build the grids.
% After that, any distribution / N combination is ~instant.

clear; close all; clc; warning('off');
ROOT = fileparts(fileparts(mfilename('fullpath')));
try, cd(ROOT); catch; end

set(0, 'DefaultFigureVisible', 'off');
try
    if isempty(gcp('nocreate'))
        parpool('local', 10);
    end
catch
end

addpath(genpath(fullfile(ROOT, 'core')));
addpath(fullfile(ROOT, 'binaries'));

cfg   = loadConfig();
FOM   = cfg.nominal.FOM_dB;
freq0 = cfg.nominal.freq_Hz;
zS0   = cfg.nominal.zS_m;
geo   = cfg.nominal.geo(:)';
N_MC  = cfg.MC.N;

% Grid resolution — wide enough to cover Normal_10pct ±3σ
N_GRID = 200;

% Widest distribution bounds (Normal_10pct = ±3σ where σ = 5% of nominal)
WIDE_DIST = cfg.distributions(end);   % Normal_10pct
B_wide    = computeVarianceBounds(cfg, WIDE_DIST);

% Extend grid 20% beyond widest distribution to avoid extrapolation
MARGIN = 1.2;
freq_grid = linspace(freq0 + MARGIN*( B_wide.freq_bnd(1)-freq0 ), ...
                     freq0 + MARGIN*( B_wide.freq_bnd(2)-freq0 ), N_GRID)';
zS_grid   = linspace(zS0   + MARGIN*( B_wide.zS_bnd(1)  -zS0   ), ...
                     zS0   + MARGIN*( B_wide.zS_bnd(2)  -zS0   ), N_GRID)';
svp_grid  = linspace(MARGIN* B_wide.svp_bnd(1), ...
                     MARGIN* B_wide.svp_bnd(2),  N_GRID)';

fprintf('=== run_MC_grid.m ===\n');
fprintf('Grid: %d pts  freq=[%.0f,%.0f]Hz  zS=[%.2f,%.2f]m  svp=[%.3f,%.3f]C\n', ...
        N_GRID, freq_grid(1), freq_grid(end), ...
        zS_grid(1), zS_grid(end), svp_grid(1), svp_grid(end));

%% ── Build grids per scenario ─────────────────────────────────────────────
for si = 1:numel(cfg.scenarios)
    sc        = cfg.scenarios(si);
    grid_dir  = fullfile('Cache', sc.name, 'grid');
    raw_cache = fullfile('Cache', sc.name, 'bellhop_raw');
    for d = {grid_dir, raw_cache}
        if ~exist(d{1},'dir'), mkdir(d{1}); end
    end

    axes_file = fullfile(grid_dir, 'grid_axes.mat');
    maxR_m    = sc.maxR_m;
    btype     = sc.bathy_type;
    max_depth = sc.maxDepth_m;

    % Save axes regardless (cheap)
    save(axes_file, 'freq_grid', 'zS_grid', 'svp_grid', 'N_GRID');

    fprintf('\n── Scenario: %s ──\n', sc.name);

    %% freq grid
    freq_file = fullfile(grid_dir, 'TL_freq_grid.mat');
    if ~isfile(freq_file)
        fprintf('  [FREQ] Running %d Bellhop calls...\n', N_GRID);
        TL_col = cell(N_GRID, 1);
        t0 = tic;
        parfor i = 1:N_GRID
            sim = {freq_grid(i), maxR_m, zS0, 0, 0, 'summer', btype, geo, FOM};
            tmp = [tempname '_blhp']; mkdir(tmp); prev = cd(tmp);
            TL_col{i} = bellhopCached(sim, raw_cache);
            cd(prev); rmdir(tmp, 's');
        end
        TL_freq_grid = cat(3, TL_col{:});   % [Nz × Nr × N_GRID]
        TL_freq_grid = single(TL_freq_grid);
        save(freq_file, 'TL_freq_grid', 'freq_grid', '-v7.3');
        fprintf('  [FREQ] Done in %.0fs\n', toc(t0));
    else
        fprintf('  [FREQ] Grid exists — skipping\n');
    end

    %% zS grid
    zS_file = fullfile(grid_dir, 'TL_zS_grid.mat');
    if ~isfile(zS_file)
        fprintf('  [ZS] Running %d Bellhop calls...\n', N_GRID);
        TL_col = cell(N_GRID, 1);
        t0 = tic;
        parfor i = 1:N_GRID
            sim = {freq0, maxR_m, zS_grid(i), 0, 0, 'summer', btype, geo, FOM};
            tmp = [tempname '_blhp']; mkdir(tmp); prev = cd(tmp);
            TL_col{i} = bellhopCached(sim, raw_cache);
            cd(prev); rmdir(tmp, 's');
        end
        TL_zS_grid = cat(3, TL_col{:});
        TL_zS_grid = single(TL_zS_grid);
        save(zS_file, 'TL_zS_grid', 'zS_grid', '-v7.3');
        fprintf('  [ZS] Done in %.0fs\n', toc(t0));
    else
        fprintf('  [ZS] Grid exists — skipping\n');
    end

    %% svp grid
    svp_file = fullfile(grid_dir, 'TL_svp_grid.mat');
    if ~isfile(svp_file)
        fprintf('  [SVP] Running %d Bellhop calls...\n', N_GRID);
        TL_col = cell(N_GRID, 1);
        t0 = tic;
        parfor i = 1:N_GRID
            svp_i = makeSVPNoise(svp_grid(i), max_depth, cfg);
            sim   = {freq0, maxR_m, zS0, 0, 0, 'custom', btype, geo, FOM};
            tmp   = [tempname '_blhp']; mkdir(tmp); prev = cd(tmp);
            TL_col{i} = bellhopCached(sim, raw_cache, 'CustomSVP', svp_i);
            cd(prev); rmdir(tmp, 's');
        end
        TL_svp_grid = cat(3, TL_col{:});
        TL_svp_grid = single(TL_svp_grid);
        save(svp_file, 'TL_svp_grid', 'svp_grid', '-v7.3');
        fprintf('  [SVP] Done in %.0fs\n', toc(t0));
    else
        fprintf('  [SVP] Grid exists — skipping\n');
    end
end

fprintf('\n=== Grid build complete. Running MC sampling from grids... ===\n');

%% ── Sample MC statistics from grids (all dists, all scenarios) ───────────
for si = 1:numel(cfg.scenarios)
    sc       = cfg.scenarios(si);
    grid_dir = fullfile('Cache', sc.name, 'grid');

    % Load grids once per scenario
    Gf = load(fullfile(grid_dir, 'TL_freq_grid.mat'), 'TL_freq_grid', 'freq_grid');
    Gz = load(fullfile(grid_dir, 'TL_zS_grid.mat'),   'TL_zS_grid',   'zS_grid');
    Gs = load(fullfile(grid_dir, 'TL_svp_grid.mat'),  'TL_svp_grid',  'svp_grid');

    [Nz, Nr, ~] = size(Gf.TL_freq_grid);
    TLf = double(reshape(Gf.TL_freq_grid, Nz*Nr, N_GRID)');  % [N_GRID × Npix]
    TLz = double(reshape(Gz.TL_zS_grid,   Nz*Nr, N_GRID)');
    TLs = double(reshape(Gs.TL_svp_grid,  Nz*Nr, N_GRID)');

    % Load r/z axes from any existing MC result for this scenario
    ax_file = fullfile('Methods','MC', sc.name, cfg.distributions(1).name, ...
                       'results', 'MC_freq.mat');
    if isfile(ax_file)
        ax = load(ax_file, 'r_km', 'z_m');
        r_km = ax.r_km;  z_m = ax.z_m;
    else
        r_km = [];  z_m = [];
    end

    for di = 1:numel(cfg.distributions)
        dist    = cfg.distributions(di);
        B       = computeVarianceBounds(cfg, dist);
        res_dir = fullfile('Methods','MC_grid', sc.name, dist.name, 'results');
        fig_dir = fullfile('Methods','MC_grid', sc.name, dist.name, 'figures');
        for d = {res_dir, fig_dir}
            if ~exist(d{1},'dir'), mkdir(d{1}); end
        end

        fprintf('  Sampling: %s / %s\n', sc.name, dist.name);

        params = {'freq', 'zS', 'svp'};
        grids  = {Gf.freq_grid, Gz.zS_grid, Gs.svp_grid};
        TLmats = {TLf, TLz, TLs};

        for pi = 1:3
            param   = params{pi};
            gx      = grids{pi};
            TLmat   = TLmats{pi};   % [N_GRID × Npix]
            mc_file = fullfile(res_dir, sprintf('MC_%s.mat', param));
            if isfile(mc_file), continue; end

            % Draw N_MC samples from the distribution
            rng(cfg.MC.rng_seed + pi);
            switch param
                case 'freq'
                    if strcmp(dist.type,'normal')
                        smp = freq0 + B.sig_freq * randn(N_MC, 1);
                    else
                        smp = B.freq_bnd(1) + diff(B.freq_bnd) * rand(N_MC, 1);
                    end
                case 'zS'
                    if strcmp(dist.type,'normal')
                        smp = zS0 + B.sig_zS * randn(N_MC, 1);
                    else
                        smp = B.zS_bnd(1) + diff(B.zS_bnd) * rand(N_MC, 1);
                    end
                case 'svp'
                    if strcmp(dist.type,'normal')
                        smp = B.sig_svp * randn(N_MC, 1);
                    else
                        smp = B.svp_bnd(1) + diff(B.svp_bnd) * rand(N_MC, 1);
                    end
            end
            % Clamp to grid range
            smp = max(min(smp, gx(end)), gx(1));

            % Interpolate TL at each sample point — [N_MC × Npix]
            TL_samp = interp1(gx, TLmat, smp, 'pchip');   % [N_MC × Npix]
            TL_all  = reshape(TL_samp', Nz, Nr, N_MC);     % [Nz × Nr × N_MC]

            MC_EX  = mean(TL_all, 3);
            MC_Var = var(TL_all,  0, 3);

            % Empirical P(shadow)
            MC_PrFOM = mean(TL_all >= FOM, 3);

            % Chebyshev lower bound on P(shadow)
            Cheb_lb = zeros(Nz, Nr);
            shadow  = MC_EX > FOM;
            d2      = (MC_EX(shadow) - FOM).^2;
            v       = MC_Var(shadow);
            Cheb_lb(shadow) = max(0, 1 - v ./ (v + d2));

            % LN3 probability from Delta moments
            [LN3_prob, ~, ~, ~] = deltaLN3prob(MC_EX, MC_Var, FOM);

            save(mc_file, 'MC_EX','MC_Var','MC_PrFOM','Cheb_lb','LN3_prob', ...
                 'r_km','z_m', '-v7.3');
        end
    end
end

fprintf('\n=== run_MC_grid.m complete ===\n');
