%% regenerate_pce_dist.m
% Regenerates pce_dist_*.png figures for all 12 scenario/distribution combos.
%
% WHY THE OLD VERSION SHOWED A SPIKE:
%   run_pce.m saves C_all using ridge lambda = 1e-10 * trace(Phi'*Phi).
%   trace(Phi'*Phi) is dominated by the He_20 column (norm^2 ~ N * 20! ~ 2.4e21),
%   making lambda ~ 2.4e11.  This completely overwhelms all PCE coefficients
%   (even c_1, whose column norm^2 ~ N * 1! = 1000), forcing C_all -> 0.
%   As a result, sampling from C_all gives y_pce ~ constant = spike.
%
%   The convergence figures and VarK_best look correct because they use
%   QR-based coefficients (no over-regularization), computed separately
%   inside run_pce.m.  Only the distribution sampling was broken.
%
% FIX (this script):
%   Re-fit the K*-order PCE using QR least squares WITH NO RIDGE on the
%   K*+1 columns only.  The well-conditioned K*-column Vandermonde does NOT
%   need regularization; the ridge in run_pce.m was only needed to stabilise
%   the MAX_ORDER=20 full system.  This gives correct coefficients c_0..c_K*
%   from which we can sample valid PCE distributions.
%
% FIGURE:
%   3 rows (one per parameter: freq, zS, svp)
%   3 cols (range slices: 25%, 50%, 75% of max range)
%   At mid-depth (50% of water column)
%   Blue = MC OAT histogram  |  Orange = PCE (K*) sampled distribution
%   OAT = One-At-a-Time: TL_pi.mat was generated with ONLY parameter pi varying.
%   Comparing PCE_pi vs MC_OAT_pi is the correct 1-D comparison.
%
% RUN:
%   cd stochastic_TL_maps_project
%   addpath(genpath('core')); addpath('scripts'); addpath('pipeline');
%   regenerate_pce_dist

%% Configuration

ROOT = fileparts(fileparts(mfilename('fullpath')));
try; cd(ROOT); catch; end
addpath(genpath(fullfile(ROOT,'core')));
addpath(fullfile(ROOT,'binaries'));

cfg    = loadConfig();
PARAMS = {'freq','zS','svp'};

M         = 5000;           % PCE samples drawn for surrogate distribution
N_BINS    = 40;             % histogram bins per panel
DIST_TYPE = 'normal';

freq0 = cfg.nominal.freq_Hz;
zS0   = cfg.nominal.zS_m;

rng(42);

%% Main loop

results_files = dir(fullfile('Methods','PCE','*','*','results','pce_results.mat'));
fprintf('Found %d pce_results.mat files.\n', numel(results_files));
n_done = 0;

for ci = 1:numel(results_files)

    %% Load PCE results (Kstar only — C_all is NOT used)
    mat_path = fullfile(results_files(ci).folder, results_files(ci).name);
    r = load(mat_path, 'sc_name','dist_name','Kstar');

    sc_name   = r.sc_name;
    dist_name = r.dist_name;
    Kstar     = r.Kstar;      % [1x3] LOO-selected PCE order per parameter

    % N from cache (how many TL samples to use)
    if isfield(r,'N_MC'),  N = r.N_MC;
    elseif isfield(r,'N'), N = r.N;
    else,                  N = cfg.N_MC;  end

    fig_dir = fullfile('Methods','PCE', sc_name, dist_name, 'figures');
    if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
    out_path = fullfile(fig_dir, sprintf('pce_dist_%s_%s', sc_name, dist_name));

    fprintf('[%d/%d] %s / %s  K*=[%d %d %d]  N=%d\n', ci, numel(results_files), ...
            sc_name, dist_name, Kstar(1), Kstar(2), Kstar(3), N);

    %% Get grid dimensions from first available TL cache
    Nz = 0;  Nr = 0;  r_km = [];  z_m = [];
    for pi = 1:3
        tf = fullfile('Cache', sc_name, dist_name, sprintf('TL_%s.mat', PARAMS{pi}));
        if ~isfile(tf), continue; end
        tmp = load(tf, 'TL_save');
        [Nz, Nr, ~] = size(tmp.TL_save);
        if isfield(tmp,'r_vec'), r_km = tmp.r_vec / 1e3;
        else,                    r_km = linspace(0, 50, Nr);  end
        if isfield(tmp,'z_vec'), z_m = tmp.z_vec;
        else,                    z_m = linspace(0, 2500, Nz); end
        clear tmp;
        break;
    end
    if Nz == 0
        fprintf('  SKIP -- no TL cache\n');
        continue;
    end
    Npix = Nz * Nr;

    %% Regenerate xi_train using same lhsSample call as run_pce.m
    % lhsSample(N, cfg, dist, B) uses rng(cfg.MC.rng_seed) internally →
    % deterministic output for given (N, dist).  Same call as in run_pce.m.
    dist_idx = find(strcmp({cfg.distributions.name}, dist_name), 1);
    if isempty(dist_idx)
        fprintf('  SKIP -- distribution not found in cfg\n');
        continue;
    end
    dist_cfg = cfg.distributions(dist_idx);
    B_target = computeVarianceBounds(cfg, dist_cfg);

    S_train  = lhsSample(N, cfg, dist_cfg, B_target);   % same call as run_pce.m
    xi_all   = [(S_train.freq - freq0) / B_target.sig_freq, ...
                (S_train.zS   - zS0)   / B_target.sig_zS, ...
                 S_train.svp           / B_target.sig_svp];   % [N x 3]

    %% Panel positions
    z_mid  = max(1, min(Nz, round(0.50 * Nz)));   % mid depth for all rows
    r_idxs = max(1, min(Nr, round([0.25 0.50 0.75] * Nr)));

    %% Draw figure: 3 rows (params) x 3 cols (range slices)
    colors = lines(3);
    fig = figure('Position', [50 50 1200 900], 'Visible', 'off');

    legend_drawn = false;
    any_panel = false;

    for pi = 1:3

        %% Load TL cache for parameter pi (OAT: only this param varies)
        tf = fullfile('Cache', sc_name, dist_name, sprintf('TL_%s.mat', PARAMS{pi}));
        if ~isfile(tf)
            fprintf('  WARN: no TL cache for %s\n', PARAMS{pi});
            continue;
        end
        tmp    = load(tf, 'TL_save');
        TL_mat = reshape(permute(double(tmp.TL_save(:,:,1:N)), [3 1 2]), N, Npix);
        clear tmp;

        %% Re-fit PCE at K* using QR (no ridge -- K*-column system is well-conditioned)
        Kp       = max(1, Kstar(pi));
        xi_pi    = xi_all(:, pi);                          % [N x 1]
        Phi_tr   = pce_basis(xi_pi, Kp, DIST_TYPE);       % [N x (Kp+1)]
        c_best   = Phi_tr \ TL_mat;                        % [Kp+1 x Npix], QR backslash

        %% Sample from PCE at M new independent xi
        xi_new   = randn(M, 1);
        Phi_new  = pce_basis(xi_new,  Kp, DIST_TYPE);     % [M x (Kp+1)]
        y_pce    = Phi_new * c_best;                       % [M x Npix]

        %% Plot 3 panels (one per range slice) in row pi
        for ri_i = 1:3
            pix = (z_mid - 1) * Nr + r_idxs(ri_i);

            mc_samps  = TL_mat(:, pix);   % OAT MC for this parameter
            pce_samps = y_pce(:, pix);

            all_v = [mc_samps; pce_samps];
            lo    = prctile(all_v, 0.5);
            hi    = prctile(all_v, 99.5);
            if hi <= lo, hi = lo + 1; end
            edges = linspace(lo, hi, N_BINS);

            ax = subplot(3, 3, (pi-1)*3 + ri_i);

            histogram(ax, mc_samps, edges, 'Normalization','pdf', ...
                      'FaceColor', colors(pi,:), 'FaceAlpha', 0.70, ...
                      'DisplayName', sprintf('MC OAT (N=%d)', N));
            hold(ax,'on');
            histogram(ax, pce_samps, edges, 'Normalization','pdf', ...
                      'FaceColor', [0.92 0.50 0.10], 'FaceAlpha', 0.55, ...
                      'DisplayName', sprintf('PCE K*=%d', Kp));

            r_str = '';  z_str = '';
            if ~isempty(r_km), r_str = sprintf('r=%.0fkm', r_km(r_idxs(ri_i))); end
            if ~isempty(z_m),  z_str = sprintf('z=%.0fm',  z_m(z_mid)); end
            title(ax, sprintf('%s  |  %s  %s', PARAMS{pi}, r_str, z_str), ...
                  'FontSize', 8, 'Interpreter','none');
            xlabel(ax, 'TL (dB)', 'FontSize', 7);
            ylabel(ax, 'PDF',     'FontSize', 7);
            grid(ax, 'on');

            if ~legend_drawn
                legend(ax, 'Location','best', 'FontSize', 7);
                legend_drawn = true;
            end
            any_panel = true;
        end

        fprintf('  [%s] K*=%d  c_1_range=[%.3f, %.3f]\n', ...
                PARAMS{pi}, Kp, min(c_best(2,:)), max(c_best(2,:)));
    end

    if ~any_panel
        close(fig);
        fprintf('  SKIP -- no panels drawn\n');
        continue;
    end

    sgtitle(sprintf('PCE (K*) vs MC OAT  |  %s  |  %s\nRows: freq / zS / svp   Cols: near / mid / far range   z=mid-depth', ...
                    sc_name, dist_name), 'FontSize', 10, 'Interpreter','none');

    saveFigPNG(fig, out_path);
    close(fig);
    n_done = n_done + 1;
    fprintf('  Saved: %s.png\n', out_path);

end

fprintf('\nDone. Regenerated %d pce_dist figures.\n', n_done);
