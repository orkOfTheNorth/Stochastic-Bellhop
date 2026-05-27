%% run_MC.m  –  Monte Carlo UQ  (2026-05-26)
%
% Runs N Monte Carlo samples for each of the 7 parameter subsets.
% Default N=100 (override: N=5; run_MC).
%
% TWO-LEVEL CACHING
%   Level 1 – bellhopCached   : individual Bellhop runs saved in cache/
%   Level 2 – TL_cache        : full TL_all cube per subset (single precision)
%     On re-run with same N and bounds → loads TL_cache, skips all Bellhop.
%     On re-run with smaller N         → uses first N slices from TL_cache.
%     On re-run with larger N          → runs fresh, overwrites TL_cache.
%     Set force_rerun=true to ignore TL_cache entirely.
%
% Saved per subset:  results/MC_<subset>.mat
%   MC_EX / MC_Var / MC_PrFOM / Cheb_lb / r_grid / z_grid / r_km / z_m / FOM / N
% Full subset also saves TL_all (for show_MC.m).
% TL_cache/ stores TL_save (single) + metadata per subset.

clearvars -except N; close all; clc; warning('off');
try, cd(fileparts(mfilename('fullpath'))); catch; end
addpath(genpath('../../uniform source depth analasis 13.5/code/Functions'));
addpath('..');   % for bellhopCached.m
if ~exist('results','dir'),  mkdir('results');  end
if ~exist('figures','dir'),  mkdir('figures');  end
cache_dir    = fullfile('..','cache');   % shared with Delta_Method
tl_cache_dir = 'TL_cache';
if ~exist(tl_cache_dir,'dir'), mkdir(tl_cache_dir); end

%% ── PARAMETERS ──────────────────────────────────────────────────────────────
freq0 = 10000;   maxR = 50000;   zS0 = 5;   FOM = 100;
geo   = [0.989*1500, 1.63, 0.07];
if ~exist('N','var'), N = 100; end   % override: N=5; run_MC
force_rerun      = false;   % set true to ignore TL_cache for ALL subsets
force_rerun_temp = true;    % re-run only temp-containing subsets (fixes SVP depth bug)

% Uniform bounds
freq_bnd = [freq0*0.99, freq0*1.01];
zS_bnd   = [zS0*0.99,  zS0*1.01];
temp_bnd = [-1, 1];

%% ── SUBSET DEFINITIONS ───────────────────────────────────────────────────────
% active flags = [freq_active, zS_active, temp_active]
subset_names  = {'zS','freq','temp','zS_freq','zS_temp','freq_temp','zS_freq_temp'};
subset_labels = {'z_S only','Freq only','Temp only', ...
                 'z_S + Freq','z_S + Temp','Freq + Temp','z_S + Freq + Temp'};
active = {[0 1 0],[1 0 0],[0 0 1],[1 1 0],[0 1 1],[1 0 1],[1 1 1]};

%% ── PRE-GENERATE SAMPLES (Latin Hypercube over all 3 dims) ──────────────────
rng(42);  % reproducibility
lhs = lhsdesign(N, 3);   % N×3  in [0,1]
% Map to physical values
samp_freq = freq_bnd(1) + lhs(:,1) * diff(freq_bnd);   % Hz
samp_zS   = zS_bnd(1)   + lhs(:,2) * diff(zS_bnd);    % m
samp_temp = temp_bnd(1)  + lhs(:,3) * diff(temp_bnd); % °C shift

%% ── DUMMY RUN TO GET GRID DIMENSIONS ────────────────────────────────────────
set(0,'DefaultFigureVisible','off');
sim_nom = {freq0, maxR, zS0, 0, 0, "summer", "const_35", geo, FOM};
[TL_dummy, r_grid, z_grid] = bellhopCached(sim_nom, cache_dir);
close all;
r_km = r_grid / 1000;
z_m  = z_grid;
[Nz, Nr] = size(TL_dummy);

%% ── MAIN LOOP OVER SUBSETS ───────────────────────────────────────────────────
for s = 1:7
    fl  = active{s};
    sn  = subset_names{s};
    fprintf('\n=== Subset %d/7: [%s] ===\n', s, subset_labels{s});

    % ── TL cache check ───────────────────────────────────────────────────────
    tl_file   = fullfile(tl_cache_dir, sprintf('TL_%s.mat', sn));
    use_cache = false;

    rerun_this = force_rerun || (force_rerun_temp && fl(3));
    if isfile(tl_file) && ~rerun_this
        tc = load(tl_file, 'N_saved','freq_bnd','zS_bnd','temp_bnd');
        bounds_ok = isequal(tc.freq_bnd, freq_bnd) && ...
                    isequal(tc.zS_bnd,   zS_bnd)   && ...
                    isequal(tc.temp_bnd, temp_bnd);
        if bounds_ok && tc.N_saved >= N
            fprintf('  [TL CACHE HIT]  using %d/%d saved samples from %s\n', ...
                    N, tc.N_saved, tl_file);
            tmp    = load(tl_file, 'TL_save');
            TL_all = double(tmp.TL_save(:,:,1:N));
            use_cache = true;
        elseif bounds_ok
            fprintf('  [TL CACHE] cached N=%d < requested N=%d — running fresh.\n', ...
                    tc.N_saved, N);
        else
            fprintf('  [TL CACHE] parameter bounds changed — running fresh.\n');
        end
    end

    if ~use_cache
        TL_all = zeros(Nz, Nr, N);

        for i = 1:N
            if mod(i,20)==0 || i==1
                fprintf('  Run %d/%d ...\n', i, N);
            end

            % Frequency (active or nominal)
            fi = freq0;
            if fl(1), fi = samp_freq(i); end

            % Source depth (active or nominal)
            zi = zS0;
            if fl(2), zi = samp_zS(i); end

            % Temperature shift (active or nominal=0)
            ti = 0;
            if fl(3), ti = samp_temp(i); end

            % Build sim_pars
            sim = {fi, maxR, zi, 0, 0, "summer", "const_35", geo, FOM};
            if fl(3)  % custom SVP
                sim{6} = "custom";
                [TL,~,~] = bellhopCached(sim, cache_dir, 'CustomSVP', makeSVP(ti, max(z_grid)));
            else
                [TL,~,~] = bellhopCached(sim, cache_dir);
            end
            TL_all(:,:,i) = TL;
            close all;
        end

        % Save TL cache (single precision to halve disk usage)
        TL_save = single(TL_all);
        N_saved = N;
        save(tl_file, 'TL_save','N_saved','freq_bnd','zS_bnd','temp_bnd', '-v7.3');
        fprintf('  [TL CACHE SAVED]  %s  (N=%d, %.0f MB)\n', ...
                tl_file, N, numel(TL_save)*4/1e6);
    end

    % Statistics
    MC_EX    = mean(TL_all, 3);
    MC_Var   = var(TL_all, 0, 3);
    MC_PrFOM = mean(TL_all > FOM, 3);   % empirical P(TL > FOM)

    % Cantelli lower bound (same formula as Delta method)
    d  = MC_EX - FOM;
    lb = zeros(Nz, Nr);
    mask = d > 0;
    lb(mask) = max(0,  1 - MC_Var(mask) ./ (MC_Var(mask) + d(mask).^2));
    Cheb_lb = lb;

    % Save (full subset: also save TL_all for interactive viewer)
    fname = fullfile('results', sprintf('MC_%s.mat', sn));
    if s == 7
        save(fname, 'MC_EX','MC_Var','MC_PrFOM','Cheb_lb', ...
                    'TL_all','r_grid','z_grid','r_km','z_m','FOM','N', ...
                    'freq_bnd','zS_bnd','temp_bnd','freq0','zS0');
    else
        save(fname, 'MC_EX','MC_Var','MC_PrFOM','Cheb_lb', ...
                    'r_grid','z_grid','r_km','z_m','FOM','N');
    end
    fprintf('  Saved %s\n', fname);

    % Quick summary figure (saved, not interactive)
    % Load nominal TL for FOM contour (from Delta results if available)
    nom_file = fullfile('..','Delta_Method','results','delta_jacobians.mat');
    if isfile(nom_file)
        tmp = load(nom_file,'TL_nom'); raw_ref = tmp.TL_nom;
    else
        raw_ref = MC_EX;
    end
    % Smooth for contour drawing only (no toolbox needed)
    TL_nom_ref = movmean(movmean(raw_ref, 20, 2), 20, 1);

    fig = figure('Name',sprintf('MC %s',subset_labels{s}),'Position',[50 50 1100 460]);

    subplot(1,2,1);
    pcolor(r_km, z_m, MC_EX); shading interp; set(gca,'YDir','reverse');
    colormap(gca,jet); colorbar; clim([50 150]);
    hold on;
    contour(r_km, z_m, TL_nom_ref,                      [FOM FOM],   'w-', 'LineWidth',0.8);
    contour(r_km, z_m, movmean(movmean(Cheb_lb,20,2),20,1), [0.95 0.95], 'w--','LineWidth',0.8);
    hold off;
    xlabel('Range (km)'); ylabel('Depth (m)');
    title(sprintf('MC E[TL] | %s | N=%d | f=%dHz zS=%.1fm',subset_labels{s},N,freq0,zS0));

    subplot(1,2,2);
    pcolor(r_km, z_m, MC_Var); shading interp; set(gca,'YDir','reverse');
    colormap(gca,hot); colorbar;
    hold on;
    contour(r_km, z_m, TL_nom_ref,                      [FOM FOM],   'w-', 'LineWidth',0.8);
    contour(r_km, z_m, movmean(movmean(Cheb_lb,20,2),20,1), [0.95 0.95], 'w--','LineWidth',0.8);
    hold off;
    xlabel('Range (km)'); ylabel('Depth (m)');
    title(sprintf('MC Var[TL] | %s | solid=FOM  dashed=Cheb95%%',subset_labels{s}));

    saveFig(fig, fullfile('figures', sprintf('MC_%s_summary', sn)));
    close(fig);
end

set(0,'DefaultFigureVisible','on');
fprintf('\n=== Monte Carlo complete. Open show_MC.m for interactive viewer. ===\n');

%% ── LOCAL FUNCTIONS ──────────────────────────────────────────────────────────
function saveFig(fig, base_path)
    print(fig, base_path, '-dpng', '-r300');
    try
        exportgraphics(fig, [base_path '.pdf'], 'ContentType','image', 'Resolution',300);
    catch
    end
end

function svp = makeSVP(temp_shift, max_depth)
    % max_depth: actual water column depth (m) — halfspace must start here,
    % not at 5000 m, so Bellhop places the rigid bottom at the correct depth.
    if nargin < 2, max_depth = 35; end
    depths  = linspace(0, max_depth, 200);
    base_td = [0 30 180 400 5000; 25 25 17 13.6 13.6]';
    t  = interp1(base_td(:,1), base_td(:,2) + temp_shift, depths, 'linear', 'extrap');
    S  = 37;
    sv = 1499.2 + 4.6*t - 0.055*t.^2 + 0.00029*t.^3 + ...
         (1.34 - 0.01*t).*(S-35) + 0.016*depths;
    svp = [depths.' sv.'];
end
