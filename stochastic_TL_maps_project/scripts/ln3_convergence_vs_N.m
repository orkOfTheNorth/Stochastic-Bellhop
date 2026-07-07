function ln3_convergence_vs_N(sc_target, dist_target)
%% ln3_convergence_vs_N  LN3 fit-quality convergence vs sample size N.
%
% For each N in N_VEC, sub-samples N TL realisations from the existing
% N=1000 caches (no new Bellhop runs):
%   MOM fit (ln3FitMOM equations, vectorised) on LHS sub-samples
%     from Cache/<sc>/<dist>/TL_<param>.mat
%   MLE fit (ln3Fit / lognfit equations, vectorised)  on IID sub-samples
%     from Cache/<sc>/<dist>_iid/TL_<param>.mat
%
% Metric per pixel (evaluated on a common grid x spanning the pixel's
% full-sample TL range, N_EVAL = cfg.ln3.N_EVAL points):
%   KS = max(abs(Fln3 - Femp));  W1 = trapz(x, abs(Fln3 - Femp))
% where Femp is the FIXED empirical CDF of the full N=1000 MC samples
% (LHS cube for the MOM curve, IID cube for the MLE curve).
%
% Reported value per N = SPATIAL MEDIAN over pixels, POOLED across the
% 3 physical parameters (freq, zS, svp).  All 3 params are used; the
% median is taken over the concatenated pixel populations of all params
% (this is the "one combined line per fit-method" the report caption asks
% for).  Pixels are spatially subsampled with a stride of 10 in both z and
% r (~4400 of ~440k pixels) to keep runtime at a few minutes; the spatial
% median is insensitive to this decimation.
%
% Output: Methods/LN3/<sc>/<dist>/figures/ln3_conv_vs_N.png
%         (2-panel loglog: KS left, W1 right;
%          MOM/LHS blue solid 'o', MLE/IID red dashed 's')
%
% Usage:
%   ln3_convergence_vs_N                              % deep_water, 5% + 10%
%   ln3_convergence_vs_N('deep_water','Normal_5pct')  % single combo

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
N_VEC   = [10 20 30 50 75 100 150 200 300 500 750 1000];
PARAMS  = {'freq','zS','svp'};
N_GRID  = cfg.ln3.N_EVAL;          % 50 CDF evaluation points per pixel
STRIDE  = 10;                      % spatial pixel decimation (z and r)

for dd = 1:numel(dist_list)
    dist_name = dist_list{dd};

    fig_dir = fullfile('Methods','LN3', sc_target, dist_name, 'figures');
    if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
    out_png = fullfile(fig_dir, 'ln3_conv_vs_N.png');
    if isfile(out_png)
        fprintf('[SKIP] %s\n', out_png); continue;
    end

    nv = numel(N_VEC);
    % per-param cell arrays of [nv x P] metric matrices (pooled at the end)
    KS_mom = cell(1,3);  W1_mom = cell(1,3);
    KS_mle = cell(1,3);  W1_mle = cell(1,3);

    for pi = 1:numel(PARAMS)
        param = PARAMS{pi};

        lhs_file = fullfile('Cache', sc_target, dist_name, ...
                            sprintf('TL_%s.mat', param));
        iid_file = fullfile('Cache', sc_target, [dist_name '_iid'], ...
                            sprintf('TL_%s.mat', param));
        if ~isfile(lhs_file), fprintf('[MISS] %s\n', lhs_file); continue; end
        if ~isfile(iid_file), fprintf('[MISS] %s\n', iid_file); continue; end

        fprintf('Loading LHS: %s\n', lhs_file);
        d = load(lhs_file,'TL_save');
        [Nz, Nr, N_MAX] = size(d.TL_save);
        iz = 1:STRIDE:Nz;  ir = 1:STRIDE:Nr;
        TL_lhs = double(reshape(d.TL_save(iz,ir,:), [], N_MAX));  % [P x 1000]
        clear d;

        fprintf('Loading IID: %s\n', iid_file);
        d = load(iid_file,'TL_save');
        TL_iid = double(reshape(d.TL_save(iz,ir,:), [], N_MAX));
        clear d;
        assert(size(TL_iid,2) == N_MAX, 'Cache size mismatch');

        % keep only pixels with finite samples in both cubes
        ok = all(isfinite(TL_lhs),2) & all(isfinite(TL_iid),2);
        TL_lhs = TL_lhs(ok,:);  TL_iid = TL_iid(ok,:);
        P = size(TL_lhs,1);
        fprintf('  [%s] %d pixels (stride %d), N_MAX=%d\n', param, P, STRIDE, N_MAX);

        % ── fixed reference: full-N empirical CDF on per-pixel grids ────────
        [x_l, Femp_l, dx_l] = ecdfGrid(TL_lhs, N_GRID);   % LHS reference
        [x_i, Femp_i, dx_i] = ecdfGrid(TL_iid, N_GRID);   % IID reference

        KS_mom{pi} = NaN(nv, P);  W1_mom{pi} = NaN(nv, P);
        KS_mle{pi} = NaN(nv, P);  W1_mle{pi} = NaN(nv, P);

        rng(42);
        for ni = 1:nv
            N = N_VEC(ni);
            if N >= N_MAX
                idx_l = 1:N_MAX;  idx_i = 1:N_MAX;
            else
                idx_l = sort(randperm(N_MAX, N));
                idx_i = sort(randperm(N_MAX, N));
            end

            % MOM fit on LHS sub-sample (vectorised ln3FitMOM equations)
            [gam, mu, sig] = fitMOMvec(TL_lhs(:, idx_l));
            F = ln3cdfGrid(x_l, gam, mu, sig);
            D = abs(F - Femp_l);
            KS_mom{pi}(ni,:) = max(D, [], 1);
            W1_mom{pi}(ni,:) = trapzUniform(D, dx_l);

            % MLE fit on IID sub-sample (vectorised ln3Fit / lognfit equations)
            [gam, mu, sig] = fitMLEvec(TL_iid(:, idx_i));
            F = ln3cdfGrid(x_i, gam, mu, sig);
            D = abs(F - Femp_i);
            KS_mle{pi}(ni,:) = max(D, [], 1);
            W1_mle{pi}(ni,:) = trapzUniform(D, dx_i);
        end
        fprintf('  [%s] done.\n', param);
        clear TL_lhs TL_iid x_l x_i Femp_l Femp_i;
    end

    % ── pool pixels across the 3 params, spatial median per N ────────────────
    ks_mom = median(cat(2, KS_mom{:}), 2, 'omitnan');
    w1_mom = median(cat(2, W1_mom{:}), 2, 'omitnan');
    ks_mle = median(cat(2, KS_mle{:}), 2, 'omitnan');
    w1_mle = median(cat(2, W1_mle{:}), 2, 'omitnan');

    % ── figure (visual style follows run_lhs_iid_convergence.m) ─────────────
    fig = figure('Position',[50 50 900 380]);
    TITLES = {'median KS = max|F_{LN3} - F_{emp}|', ...
              'median W_1 = \int|F_{LN3} - F_{emp}|dx  (dB)'};
    mom_data = {ks_mom, w1_mom};
    mle_data = {ks_mle, w1_mle};

    for p = 1:2
        ax = subplot(1,2,p);
        loglog(ax, N_VEC, max(mom_data{p},1e-8), 'b-o', ...
               'LineWidth',1.8, 'MarkerSize',5, ...
               'DisplayName','MOM (LHS)'); hold(ax,'on');
        loglog(ax, N_VEC, max(mle_data{p},1e-8), 'r--s', ...
               'LineWidth',1.8, 'MarkerSize',5, ...
               'DisplayName','MLE (IID)');
        grid(ax,'on'); legend(ax,'Location','northeast','FontSize',8);
        xlabel(ax,'N (sample size)');
        ylabel(ax, TITLES{p});
        title(ax, TITLES{p}, 'Interpreter','tex');
        xlim(ax,[N_VEC(1) N_VEC(end)]);
    end
    sgtitle(sprintf('LN3 fit convergence vs N | %s | %s | ref: empirical CDF (N=1000)', ...
        sc_target, dist_name), 'Interpreter','none', 'FontSize',11);

    print(fig, out_png, '-dpng', '-r150');
    close(fig);
    fprintf('Saved: %s\n', out_png);
end
fprintf('ln3_convergence_vs_N done.\n');
end

%% ── helpers ──────────────────────────────────────────────────────────────────

function [x, Femp, dx] = ecdfGrid(TL, N_GRID)
% Per-pixel uniform grid over the full-sample TL range + empirical CDF on it.
%   TL   [P x N]  samples per pixel
%   x    [N_GRID x P], Femp [N_GRID x P], dx [1 x P]
    lo = min(TL, [], 2)';                 % [1 x P]
    hi = max(TL, [], 2)';
    hi = max(hi, lo + 1e-6);
    dx = (hi - lo) / (N_GRID - 1);
    x  = lo + (0:N_GRID-1)' .* dx;        % [N_GRID x P]
    Femp = zeros(N_GRID, size(TL,1));
    for j = 1:N_GRID
        Femp(j,:) = mean(TL <= x(j,:)', 2)';
    end
end

function [gam, mu, sig] = fitMOMvec(TL)
% Vectorised ln3FitMOM: gamma = 0.95*min, then MOM for LN2 on shifted samples.
%   TL [P x n]  →  gam, mu, sig  [1 x P]
    gam = 0.95 * min(TL, [], 2)';         % [1 x P]
    sh  = TL - gam';
    sh(sh <= 0) = 1e-6;
    M1 = mean(sh, 2)';
    M2 = mean(sh.^2, 2)';
    sig2 = log(max(M2 ./ max(M1.^2, eps), 1 + 1e-12));
    sig  = sqrt(sig2);
    mu   = log(max(M1, eps)) - 0.5 * sig2;
    % degenerate pixels (M2<=M1^2): fall back to log-moment (MLE) estimates,
    % mirroring the lognfit fallback in ln3FitMOM.m
    bad = (M2 <= M1.^2) | (M1 <= 0);
    if any(bad)
        L = log(sh(bad,:));
        mu(bad)  = mean(L, 2)';
        sig(bad) = sqrt(mean((L - mean(L,2)).^2, 2))';
    end
end

function [gam, mu, sig] = fitMLEvec(TL)
% Vectorised ln3Fit: gamma = 0.95*min, then lognormal MLE on shifted samples
% (lognfit MLE = sample mean / sqrt(biased variance) of log-samples).
    gam = 0.95 * min(TL, [], 2)';
    sh  = TL - gam';
    sh(sh <= 0) = 1e-6;
    L   = log(sh);
    mu  = mean(L, 2)';
    sig = sqrt(mean((L - mean(L,2)).^2, 2))';
    sig = max(sig, 1e-8);
end

function F = ln3cdfGrid(x, gam, mu, sig)
% LN3 CDF on per-pixel grid: F = logncdf(x - gam, mu, sig), 0 for x <= gam.
%   x [G x P]; gam, mu, sig [1 x P]  →  F [G x P]
    y = x - gam;
    F = zeros(size(x));
    m = y > 0;
    Z = (log(max(y,realmin)) - mu) ./ (sqrt(2) * sig);
    Fall = 0.5 * erfc(-Z);
    F(m) = Fall(m);
end

function w = trapzUniform(D, dx)
% Trapezoidal integral along dim 1 with per-column uniform spacing dx [1 x P].
    w = dx .* (sum(D,1) - 0.5*(D(1,:) + D(end,:)));
end
