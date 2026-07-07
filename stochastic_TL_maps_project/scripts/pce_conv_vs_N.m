function pce_conv_vs_N(sc_target, dist_target)
%% pce_conv_vs_N  PCE(K*) fit-quality vs training sample size N.
%
% For the LOO-CV-selected order K* (loaded from
% Methods/PCE/<sc>/<dist>/results/pce_results.mat, field Kstar [1x3]),
% the PCE coefficients are refitted on N training samples sub-sampled
% from the existing N=1000 LHS cache Cache/<sc>/<dist>/TL_<param>.mat
% (same cube run_pce.m fits on; cfg.N_MC = 1000).  No new Bellhop runs.
%
% Fit: same Hermite basis (core/math/pce_basis.m) + QR least squares:
%   Phi_n = pce_basis(xi(idx), K*);  [Q,R] = qr(Phi_n,0);  C = R\(Q'*TL);
% Fitted on RAW (uncentred) TL, exactly like scripts/regenerate_pce_dist.m,
% so the reconstructed distribution has the correct mean (avoids the
% centred-mean bug documented in README Known Issues).
%
% PCE-predicted CDF: Monte-Carlo evaluation of the fitted polynomial at
% M=2000 fresh xi ~ N(0,1) draws (K* coefficients only — never K_max,
% to avoid Hermite blow-up at |xi|>2).
%
% Metric per pixel, on a grid x spanning the pixel's full-sample TL range:
%   KS = max(abs(Fpce - Femp));  W1 = trapz(x, abs(Fpce - Femp))
% where Femp = empirical CDF of the FULL N=1000 LHS MC samples (fixed ref).
% Reported value per N = spatial median over pixels, POOLED across the 3
% params (freq, zS, svp) — one combined KS line and one W1 line.
% Pixels are decimated with stride 10 in z and r (~4400 pixels) for runtime.
%
% N values below the per-param regression threshold N >= 2*(K*+1) are
% skipped for that param (regression underdetermined margin); the vertical
% dashed line marks N = 2*(max(Kstar)+1), the strictest threshold among
% the 3 params.
%
% Output: Methods/PCE/<sc>/<dist>/figures/pce_conv_vs_N.png

if nargin < 1 || isempty(sc_target), sc_target = 'deep_water'; end
if nargin < 2 || isempty(dist_target)
    dist_list = {'Normal_5pct','Normal_10pct'};
else
    dist_list = {dist_target};
end

close all;
ROOT = fileparts(fileparts(mfilename('fullpath')));
try; cd(ROOT); catch; end
addpath(genpath(fullfile(ROOT,'core')));
set(0,'DefaultFigureVisible','off');

cfg     = loadConfig();
N_MAX   = cfg.N_MC;                            % 1000
% 50..1000 per the report, extended down to 20/30 so the region just above
% the regression threshold N = 2*(K*+1) (typically 12-16 here) is visible;
% N values below the per-param threshold are dropped automatically.
N_VEC   = [20 30 50 75 100 150 200 300 500 750 1000];
PARAMS  = {'freq','zS','svp'};
N_GRID  = cfg.ln3.N_EVAL;                      % 50 CDF evaluation points
STRIDE  = 10;                                  % spatial pixel decimation
M_EVAL  = 2000;                                % PCE reconstruction draws
freq0   = cfg.nominal.freq_Hz;
zS0     = cfg.nominal.zS_m;

for dd = 1:numel(dist_list)
    dist_name = dist_list{dd};

    fig_dir = fullfile('Methods','PCE', sc_target, dist_name, 'figures');
    if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
    out_png = fullfile(fig_dir, 'pce_conv_vs_N.png');
    if isfile(out_png)
        fprintf('[SKIP] %s\n', out_png); continue;
    end

    res_file = fullfile('Methods','PCE', sc_target, dist_name, ...
                        'results','pce_results.mat');
    if ~isfile(res_file)
        fprintf('[MISS] %s\n', res_file); continue;
    end
    r = load(res_file, 'Kstar');
    Kstar = r.Kstar;                          % [1 x 3] LOO-selected order
    fprintf('%s / %s : K* = [%d %d %d]\n', sc_target, dist_name, ...
            Kstar(1), Kstar(2), Kstar(3));

    % ── regenerate the same standardised xi used by run_pce.m ────────────────
    di = find(strcmp({cfg.distributions.name}, dist_name), 1);
    dist_cfg = cfg.distributions(di);
    B = computeVarianceBounds(cfg, dist_cfg);
    S = lhsSample(N_MAX, cfg, dist_cfg, B);   % deterministic (rng inside)
    xi_all = [(S.freq - freq0) / B.sig_freq, ...
              (S.zS   - zS0)   / B.sig_zS, ...
               S.svp           / B.sig_svp];   % [N_MAX x 3]

    nv = numel(N_VEC);
    KS_all = cell(1,3);  W1_all = cell(1,3);

    for pi = 1:numel(PARAMS)
        param = PARAMS{pi};
        Kp    = max(1, Kstar(pi));
        N_min = 2 * (Kp + 1);                 % regression threshold

        tl_file = fullfile('Cache', sc_target, dist_name, ...
                           sprintf('TL_%s.mat', param));
        if ~isfile(tl_file), fprintf('[MISS] %s\n', tl_file); continue; end

        fprintf('Loading: %s\n', tl_file);
        d = load(tl_file,'TL_save');
        [Nz, Nr, Nc] = size(d.TL_save);
        assert(Nc >= N_MAX, 'Cache smaller than N_MC');
        iz = 1:STRIDE:Nz;  ir = 1:STRIDE:Nr;
        TL = double(reshape(d.TL_save(iz,ir,1:N_MAX), [], N_MAX));  % [P x 1000]
        clear d;
        ok = all(isfinite(TL),2);
        TL = TL(ok,:);
        P  = size(TL,1);
        fprintf('  [%s] K*=%d  N_min=%d  %d pixels (stride %d)\n', ...
                param, Kp, N_min, P, STRIDE);

        % fixed reference: full-N=1000 empirical CDF on per-pixel grids
        [x, Femp, dx] = ecdfGrid(TL, N_GRID);

        % fixed reconstruction draws (same for every N → smooth curves)
        rng(42);
        xi_new  = randn(M_EVAL, 1);
        Phi_new = pce_basis(xi_new, Kp, 'normal');    % [M x Kp+1]

        xi = xi_all(:, pi);                           % [N_MAX x 1]

        KS_all{pi} = NaN(nv, P);  W1_all{pi} = NaN(nv, P);

        rng(42);
        for ni = 1:nv
            N = N_VEC(ni);
            if N < N_min, continue; end               % underdetermined margin
            if N >= N_MAX, idx = 1:N_MAX;
            else,          idx = sort(randperm(N_MAX, N));
            end

            % QR least-squares refit at K* on N training samples
            Phi_n  = pce_basis(xi(idx), Kp, 'normal');   % [N x Kp+1]
            [Q, R_qr] = qr(Phi_n, 0);
            C = R_qr \ (Q' * TL(:, idx)');               % [Kp+1 x P]

            % PCE-predicted CDF via MC evaluation of the polynomial
            y_pce = Phi_new * C;                          % [M x P]
            Fpce = zeros(N_GRID, P);
            for j = 1:N_GRID
                Fpce(j,:) = mean(y_pce <= x(j,:), 1);
            end

            D = abs(Fpce - Femp);
            KS_all{pi}(ni,:) = max(D, [], 1);
            W1_all{pi}(ni,:) = dx .* (sum(D,1) - 0.5*(D(1,:) + D(end,:)));
        end
        fprintf('  [%s] done.\n', param);
        clear TL x Femp;
    end

    % ── pool pixels across params, spatial median per N ─────────────────────
    ks_med = median(cat(2, KS_all{:}), 2, 'omitnan');
    w1_med = median(cat(2, W1_all{:}), 2, 'omitnan');
    N_thresh = 2 * (max(Kstar) + 1);   % strictest per-param threshold

    % ── figure ────────────────────────────────────────────────────────────
    fig = figure('Position',[50 50 900 380]);
    TITLES = {'median KS = max|F_{PCE} - F_{emp}|', ...
              'median W_1 = \int|F_{PCE} - F_{emp}|dx  (dB)'};
    data = {ks_med, w1_med};

    for p = 1:2
        ax = subplot(1,2,p);
        loglog(ax, N_VEC, max(data{p},1e-8), 'b-o', ...
               'LineWidth',1.8, 'MarkerSize',5, ...
               'DisplayName','PCE(K*) vs MC'); hold(ax,'on');
        xline(ax, N_thresh, 'k--', sprintf('N=2(K*+1)=%d', N_thresh), ...
              'LineWidth',1.2, 'LabelVerticalAlignment','bottom', ...
              'DisplayName','regression threshold');
        grid(ax,'on'); legend(ax,'Location','northeast','FontSize',8);
        xlabel(ax,'N (training samples)');
        ylabel(ax, TITLES{p});
        title(ax, TITLES{p}, 'Interpreter','tex');
        xlim(ax,[min(N_VEC(1), 0.8*N_thresh) N_VEC(end)]);
    end
    sgtitle(sprintf('PCE(K*) fit convergence vs N | %s | %s | ref: empirical CDF (N=1000)', ...
        sc_target, dist_name), 'Interpreter','none', 'FontSize',11);

    print(fig, out_png, '-dpng', '-r150');
    close(fig);
    fprintf('Saved: %s\n', out_png);
end
fprintf('pce_conv_vs_N done.\n');
end

%% ── helper ──────────────────────────────────────────────────────────────────

function [x, Femp, dx] = ecdfGrid(TL, N_GRID)
% Per-pixel uniform grid over the full-sample TL range + empirical CDF on it.
%   TL [P x N]  →  x [N_GRID x P], Femp [N_GRID x P], dx [1 x P]
    lo = min(TL, [], 2)';
    hi = max(TL, [], 2)';
    hi = max(hi, lo + 1e-6);
    dx = (hi - lo) / (N_GRID - 1);
    x  = lo + (0:N_GRID-1)' .* dx;
    Femp = zeros(N_GRID, size(TL,1));
    for j = 1:N_GRID
        Femp(j,:) = mean(TL <= x(j,:)', 2)';
    end
end
