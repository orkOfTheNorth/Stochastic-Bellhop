%% plot_taylor_order_study.m — PCE-based variance convergence study
%
% Uses Polynomial Chaos Expansion (PCE) to assess how many expansion terms
% are needed to accurately capture Var[TL] for each parameter.
%
% METHOD: Gauss quadrature projection
% ─────────────────────────────────────────────────────────────────────────
% TL(θ) ≈ Σ_{n=0}^{K} c_n Φ_n(ξ)
%
%   ξ = (θ−θ₀)/σ         (standardised variable)
%   Φ_n = Hermite He_n    (Normal distributions)
%   Φ_n = Legendre P_n    (Uniform distributions)
%
% PCE coefficients via Gauss quadrature (K+1 points → exact for deg 2K+1):
%   c_n = Σ_i w_i TL(θ_i) Φ_n(ξ_i) / γ_n
%
% Order-K variance (orthogonality eliminates cross-terms):
%   Var_K = Σ_{n=1}^{K} c_n² γ_n
%
% For each order K=1..MAX_ORDER, K+1 Bellhop runs are needed.
% To avoid re-running: evaluate at MAX_ORDER+1 points once, reuse for all K.
%
% Quadrature nodes live WITHIN the distribution support (no ±10σ blowup).
% Numerically stable — no catastrophic cancellation.
%
% Outputs: Methods/Delta/taylor_order_analysis/taylor_order_study_<scen>_<dist>.png

clear; close all; clc; warning('off');
ROOT = fileparts(fileparts(mfilename('fullpath')));
try, cd(ROOT); catch; end

addpath(genpath(fullfile(ROOT, 'core')));
addpath(fullfile(ROOT, 'binaries'));

set(0, 'DefaultFigureVisible', 'off');
try
    if isempty(gcp('nocreate'))
        parpool('local', 10);
    end
catch
end

cfg   = loadConfig();
FOM   = cfg.nominal.FOM_dB;
freq0 = cfg.nominal.freq_Hz;
zS0   = cfg.nominal.zS_m;
geo   = cfg.nominal.geo(:)';

STUDY_CASES = {
    'downslope',  'Uniform_5pct';
    'baseline',   'Normal_10pct';
    'deep_water', 'Normal_10pct';
};

MAX_ORDER  = 10;   % PCE order to assess
N_EVAL     = 50;   % fixed quadrature points for coefficient computation (>> MAX_ORDER)
PARAM_LIST = {'freq', 'zS', 'svp'};

fig_dir = fullfile('Methods','Delta','taylor_order_analysis');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

fprintf('=== PCE variance convergence study (orders 1–%d) ===\n', MAX_ORDER);
fprintf('    %d Bellhop runs per parameter per case.\n', MAX_ORDER+1);

for ci = 1:size(STUDY_CASES,1)
    sc_name   = STUDY_CASES{ci,1};
    dist_name = STUDY_CASES{ci,2};

    fprintf('\n────────────────────────────────────────\n');
    fprintf('  Case %d/3: %s / %s\n', ci, sc_name, dist_name);
    fprintf('────────────────────────────────────────\n');

    sc        = cfg.scenarios(strcmp({cfg.scenarios.name}, sc_name));
    dist_idx  = find(strcmp({cfg.distributions.name}, dist_name), 1);
    dist_cfg  = cfg.distributions(dist_idx);
    B         = computeVarianceBounds(cfg, dist_cfg);

    % N_EVAL Gauss quadrature nodes — enough to accurately compute all c_n up to MAX_ORDER
    % (need N_EVAL >= 2*MAX_ORDER+1 for exact integration of degree-2*MAX_ORDER polynomials)
    [xi_nodes, w_quad, gamma_n] = gaussQuadNodes(N_EVAL, dist_cfg.type);

    % Map standardised nodes to physical parameter values
    sig_vals = [B.sig_freq, B.sig_zS, B.sig_svp];
    a_vals   = [B.freq_bnd(2)-freq0, B.zS_bnd(2)-zS0, B.svp_bnd(2)];  % half-width for uniform

    mc_res    = fullfile('Methods','MC', sc_name, dist_name, 'results');
    raw_cache = fullfile(pwd, 'Cache', sc_name, 'bellhop_raw');
    if ~exist(raw_cache,'dir'), mkdir(raw_cache); end

    L1_cells  = cell(1,3);
    Vk_cells  = cell(1,3);
    Vmc_cells = cell(1,3);
    rk_cell   = cell(1,3);
    zm_cell   = cell(1,3);

    for pi = 1:3
        param   = PARAM_LIST{pi};
        mc_file = fullfile(mc_res, sprintf('MC_%s.mat', param));
        if ~isfile(mc_file)
            fprintf('  SKIP param=%s — MC data missing\n', param);
            L1_cells{pi}  = nan(MAX_ORDER,1);
            Vk_cells{pi}  = cell(MAX_ORDER,1);
            Vmc_cells{pi} = [];
            continue;
        end

        Mc = load(mc_file, 'MC_Var', 'r_km', 'z_m');
        [Nz, Nr] = size(Mc.MC_Var);
        rk_cell{pi} = Mc.r_km;
        zm_cell{pi} = Mc.z_m;

        % Physical quadrature abscissae
        switch param
            case 'freq'
                if strcmp(dist_cfg.type,'normal')
                    theta_q = freq0 + sig_vals(1) * xi_nodes;
                else
                    theta_q = freq0 + a_vals(1) * xi_nodes;
                end
            case 'zS'
                if strcmp(dist_cfg.type,'normal')
                    theta_q = zS0 + sig_vals(2) * xi_nodes;
                else
                    theta_q = zS0 + a_vals(2) * xi_nodes;
                end
            case 'svp'
                if strcmp(dist_cfg.type,'normal')
                    theta_q = sig_vals(3) * xi_nodes;
                else
                    theta_q = a_vals(3) * xi_nodes;
                end
        end

        Nq        = N_EVAL;   % fixed large set — same for all orders
        maxR_m    = sc.maxR_m;
        btype     = sc.bathy_type;
        max_depth = sc.maxDepth_m;

        pce_cache_dir  = fullfile(pwd, 'Cache', sc_name, 'pce');
        if ~exist(pce_cache_dir,'dir'), mkdir(pce_cache_dir); end
        pce_cache_file = fullfile(pce_cache_dir, ...
            sprintf('TL_pce_%s_%s_N%d.mat', param, dist_name, N_EVAL));

        if isfile(pce_cache_file)
            fprintf('  [%s/%s] param=%s: loading PCE cache...\n', sc_name, dist_name, param);
            tmp = load(pce_cache_file, 'TL_quad_mat', 'xi_nodes_saved', 'theta_q_saved');
            TL_quad_mat = tmp.TL_quad_mat;   % [Nq x Nz*Nr] already stacked
        else
            fprintf('  [%s/%s] param=%s: %d quadrature Bellhop calls...\n', ...
                    sc_name, dist_name, param, Nq);
            TL_quad_mat = [];   % will be filled below
        end

        if isempty(TL_quad_mat)
            TL_quad = cell(Nq, 1);
            t0 = tic;
            parfor qi = 1:Nq
            switch param
                case 'freq'
                    fi = theta_q(qi); zi = zS0; sv = 0;
                    sim = {fi, maxR_m, zi, 0, 0, 'summer', btype, geo, FOM};
                    tmp = [tempname '_blhp']; mkdir(tmp); prev = cd(tmp);
                    TL_quad{qi} = bellhopCached(sim, raw_cache);
                    cd(prev); rmdir(tmp,'s');
                case 'zS'
                    zi = max(theta_q(qi), 0.1); fi = freq0; sv = 0;
                    sim = {fi, maxR_m, zi, 0, 0, 'summer', btype, geo, FOM};
                    tmp = [tempname '_blhp']; mkdir(tmp); prev = cd(tmp);
                    TL_quad{qi} = bellhopCached(sim, raw_cache);
                    cd(prev); rmdir(tmp,'s');
                case 'svp'
                    sv = theta_q(qi); fi = freq0; zi = zS0;
                    svp_i = makeSVPNoise(sv, max_depth, cfg);
                    sim = {fi, maxR_m, zi, 0, 0, 'custom', btype, geo, FOM};
                    tmp = [tempname '_blhp']; mkdir(tmp); prev = cd(tmp);
                    TL_quad{qi} = bellhopCached(sim, raw_cache, 'CustomSVP', svp_i);
                    cd(prev); rmdir(tmp,'s');
            end
        end
            fprintf('  [%s] Bellhop done in %.0fs\n', param, toc(t0));

            % Stack into [Nq x Npix] and save cache
            TL_q        = cat(3, TL_quad{:});
            TL_quad_mat = reshape(permute(TL_q,[3 1 2]), Nq, Nz*Nr);
            xi_nodes_saved = xi_nodes;
            theta_q_saved  = theta_q;
            save(pce_cache_file, 'TL_quad_mat','xi_nodes_saved','theta_q_saved','-v7.3');
            fprintf('  PCE cache saved: %s\n', pce_cache_file);
        end

        % TL at quadrature points: [Nq x Npix]
        TL_mat = TL_quad_mat;

        % Evaluate PCE basis Phi_0..Phi_MAX_ORDER at all N_EVAL nodes [N_EVAL x MAX_ORDER+1]
        Phi = pce_basis(xi_nodes, MAX_ORDER, dist_cfg.type);

        % PCE coefficients via Gauss quadrature inner product (exact up to aliasing):
        %   c_n = (1/gamma_n) * Σ_i w_i * TL(xi_i) * Phi_n(xi_i)
        % With N_EVAL=50 >> MAX_ORDER=10, coefficients are accurately resolved.
        % C [MAX_ORDER+1 x Npix]
        WPhi = diag(w_quad) * Phi;        % [N_EVAL x MAX_ORDER+1], weighted basis
        C    = WPhi' * TL_mat;                        % [MAX_ORDER+1 x Npix]
        C    = C ./ gamma_n(1:MAX_ORDER+1);           % divide each row by gamma_n

        % ── Compute Var_K for each order K ───────────────────────────────
        Var_MC_pix = Mc.MC_Var(:);
        mean_VM    = max(mean(Var_MC_pix), 1e-10);

        l1_vec     = zeros(MAX_ORDER, 1);
        vk_cell_p  = cell(MAX_ORDER, 1);

        for K = 1:MAX_ORDER
            % Var_K = Σ_{n=1}^{K} c_n² γ_n  (orthogonality — no cross terms)
            Var_pix = zeros(1, Nz*Nr);
            for n = 1:K
                Var_pix = Var_pix + C(n+1,:).^2 * gamma_n(n+1);
            end
            Var_pix = max(Var_pix, 0);

            vk_cell_p{K}  = reshape(Var_pix, Nz, Nr);
            l1_vec(K)     = mean(abs(Var_pix(:) - Var_MC_pix)) / mean_VM;
        end

        L1_cells{pi}  = l1_vec;
        Vk_cells{pi}  = vk_cell_p;
        Vmc_cells{pi} = Mc.MC_Var;

        fprintf('  [%s/%s] param=%s done.\n', sc_name, dist_name, param);
    end

    % ── Print L1 table ────────────────────────────────────────────────────
    L1_mat = cell2mat(L1_cells);
    fprintf('\n  Relative L1 error per PCE order:\n');
    fprintf('  Order |  freq   |   zS    |   svp\n');
    for K = 1:MAX_ORDER
        fprintf('    %2d  |  %.4f  |  %.4f  |  %.4f\n', K, L1_mat(K,1), L1_mat(K,2), L1_mat(K,3));
    end

    good_ord = zeros(1,3);
    for pi = 1:3
        idx = find(L1_mat(:,pi) < 0.10, 1);
        if isempty(idx), good_ord(pi) = NaN; else, good_ord(pi) = idx; end
    end
    fprintf('  First PCE order with L1<10%%: freq=%s  zS=%s  svp=%s\n', ...
            fmt_ord(good_ord(1)), fmt_ord(good_ord(2)), fmt_ord(good_ord(3)));

    % ── Figure ───────────────────────────────────────────────────────────
    r_km = rk_cell{find(~cellfun(@isempty,rk_cell),1)};
    z_m  = zm_cell{find(~cellfun(@isempty,zm_cell),1)};
    clrs = {[0 0.4 0.8],[0.8 0.2 0],[0.1 0.6 0.1]};

    % Show maps for the WORST (most nonlinear) parameter — most informative
    [~, best_pi] = max(mean(L1_mat, 1, 'omitnan'));
    [~, best_k]  = min(L1_mat(:,best_pi));

    fig = figure('Position',[50 50 1600 500]);

    ax1 = subplot(1,5,1); hold(ax1,'on');
    for pi = 1:3
        if ~any(isnan(L1_cells{pi}))
            plot(ax1, 1:MAX_ORDER, L1_mat(:,pi), '-o', 'Color',clrs{pi}, ...
                 'LineWidth',1.8, 'DisplayName',PARAM_LIST{pi});
        end
    end
    yline(ax1,0.10,'k--','LineWidth',1.2,'DisplayName','10%');
    xlabel(ax1,'PCE Order K'); ylabel(ax1,'Relative L1');
    title(ax1, sprintf('PCE L1 vs order\n%s | %s',sc_name,dist_name), ...
          'Interpreter','none','FontSize',8);
    legend(ax1,'Location','best','FontSize',7);
    ylim(ax1,[0, min(max(L1_mat(:))*1.1,3)]);
    xticks(ax1,1:MAX_ORDER); grid(ax1,'on'); hold(ax1,'off');

    show_orders = unique([1 5 min(best_k,MAX_ORDER)]);
    for col = 1:3
        ax = subplot(1,5,col+1);
        k_show = show_orders(min(col,numel(show_orders)));
        if ~isempty(Vk_cells{best_pi}{k_show}) && ~isempty(Vmc_cells{best_pi})
            diff_m = Vk_cells{best_pi}{k_show} - Vmc_cells{best_pi};
            mx = max(abs(diff_m(:)));
            if mx==0, mx=1; end
            pcolor(ax, r_km, z_m, diff_m); shading(ax,'interp');
            set(ax,'YDir','reverse'); colormap(ax,redblue(256));
            clim(ax,[-mx mx]); colorbar(ax);
        end
        xlabel(ax,'Range (km)');
        title(ax, sprintf('Var_{K=%d} − Var_{MC}\n%s', k_show, PARAM_LIST{best_pi}), ...
              'Interpreter','none','FontSize',8);
    end

    ax5 = subplot(1,5,5);
    if ~isempty(Vk_cells{best_pi}{best_k}) && ~isempty(Vmc_cells{best_pi})
        Vk  = Vk_cells{best_pi}{best_k}(:);
        Vmc = Vmc_cells{best_pi}(:);
        ok  = Vmc > 0;
        scatter(ax5, Vmc(ok), Vk(ok), 3, 'filled', ...
                'MarkerFaceAlpha',0.25,'MarkerFaceColor',clrs{best_pi});
        hold(ax5,'on');
        lv = max([Vmc(ok); Vk(ok)]);
        plot(ax5,[0 lv],[0 lv],'r-','LineWidth',1.5); hold(ax5,'off');
        R2 = corr(Vmc(ok),Vk(ok))^2;
        title(ax5,sprintf('Scatter | %s\nR²=%.3f',PARAM_LIST{best_pi},R2), ...
              'FontSize',8,'Interpreter','none');
    end
    xlabel(ax5,'Var_{MC}'); ylabel(ax5,sprintf('Var_{K=%d}',best_k));
    grid(ax5,'on'); axis(ax5,'equal');

    sgtitle(sprintf('PCE Variance Convergence | %s | %s', sc_name, dist_name), ...
            'FontSize',11,'Interpreter','none');

    out_path = fullfile(fig_dir, sprintf('taylor_order_study_%s_%s', sc_name, dist_name));
    saveFigPNG(fig, out_path);
    drawnow; close(fig);
    fprintf('  Saved: %s.png\n', out_path);
end

fprintf('\n=== plot_taylor_order_study.m complete ===\n');


%% ── Local helpers ────────────────────────────────────────────────────────────
function Phi = pce_basis(xi, K, dist_type)
% Evaluate PCE basis polynomials Phi_0..Phi_K at nodes xi [N x 1].
% Returns Phi [N x K+1].
N   = numel(xi);
Phi = zeros(N, K+1);
switch lower(dist_type)
    case 'normal'
        % Probabilist Hermite: He_0=1, He_1=x, He_{n+1}=x*He_n - n*He_{n-1}
        Phi(:,1) = ones(N,1);
        if K >= 1, Phi(:,2) = xi; end
        for n = 1:K-1
            Phi(:,n+2) = xi .* Phi(:,n+1) - n * Phi(:,n);
        end
    case 'uniform'
        % Legendre: P_0=1, P_1=x, P_{n+1}=((2n+1)*x*P_n - n*P_{n-1})/(n+1)
        Phi(:,1) = ones(N,1);
        if K >= 1, Phi(:,2) = xi; end
        for n = 1:K-1
            Phi(:,n+2) = ((2*n+1)*xi.*Phi(:,n+1) - n*Phi(:,n)) / (n+1);
        end
end
end

function s = fmt_ord(n)
if isnan(n), s='>10'; else, s=num2str(n); end
end
