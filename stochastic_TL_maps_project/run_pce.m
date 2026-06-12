%% run_pce.m — PCE variance estimation from cached MC TL cubes
%
% For each completed scenario × distribution, fits Polynomial Chaos Expansion
% coefficients using the Normal_10pct base-distribution TL samples (no new Bellhop runs).
%
% ARCHITECTURE (IS-aware)
% ───────────────────────
%   PCE is fitted ONCE per scenario/param on Normal_10pct samples.
%   Variance under other distributions uses analytic rescaling:
%
%     Var_K(σ_target) = Σ_{n=1}^K  c_n² · n! · (σ_target/σ_base)^{2n}
%
%   This is exact for polynomial TL functions and avoids running Bellhop
%   for Normal_1pct / Normal_5pct.
%
% COMBINED LOO-CV CRITERION
% ──────────────────────────
%   LOO_combined(K) = α · LOO_Var_norm(K) + (1−α) · LOO_EX_norm(K)
%
%   where α = cfg.LOO_alpha (= 0.5), and:
%     LOO_EX(K)  = analytic leave-one-out MSE (hat matrix, standard)
%     LOO_Var(K) = LOO variance prediction error vs MC Var
%
%   K* = argmin_K LOO_combined(K)    [single order for both EX and Var]
%
% Outputs per scenario × distribution:
%   Methods/PCE/<scen>/<dist>/results/pce_results.mat
%   Methods/PCE/<scen>/<dist>/figures/pce_convergence_<scen>_<dist>.png
%   Methods/PCE/<scen>/<dist>/figures/pce_order_relevance.png
%   Methods/PCE/<scen>/<dist>/figures/pce_detect_scatter.png
% Summary:
%   Methods/PCE/pce_summary_all.png

function run_pce(sc_target, dist_target)
%% run_pce — PCE from cached MC data.  No args = all combos.
if nargin < 2, sc_target = ''; dist_target = ''; end
close all;
warning('off', 'MATLAB:unknownObjectIEEE');
warning('off', 'MATLAB:singularMatrix');
warning('off', 'MATLAB:rankDeficientMatrix');
try cd(fileparts(mfilename('fullpath'))); catch; end

set(0, 'DefaultFigureVisible', 'off');
addpath(genpath('Shared_Utils'));

cfg   = loadConfig();
freq0 = cfg.nominal.freq_Hz;
zS0   = cfg.nominal.zS_m;
N     = cfg.N_MC;

MAX_ORDER = cfg.max_pce_order;   % 20 (from config)
alpha_loo = cfg.LOO_alpha;       % 0.5 — weight for LOO_Var vs LOO_EX
PARAMS    = {'freq', 'zS', 'svp'};

% Base distribution (Normal_10pct): PCE is always fitted on this
base_name = cfg.IS_base_dist;    % 'Normal_10pct'
base_idx  = find(strcmp({cfg.distributions.name}, base_name), 1);
if isempty(base_idx)
    error('IS base distribution "%s" not found in config.', base_name);
end
B_base = computeVarianceBounds(cfg, cfg.distributions(base_idx));
sig_base = [B_base.sig_freq, B_base.sig_zS, B_base.sig_svp];   % [1 × 3]

%% ── Identify completed combos ────────────────────────────────────────────────
fprintf('=== PCE (order 1..%d, combined LOO) ===\n', MAX_ORDER);

completed = {};
for si = 1:numel(cfg.scenarios)
    sc = cfg.scenarios(si);
    for di = 1:numel(cfg.distributions)
        dist = cfg.distributions(di);
        res_dir = fullfile('Methods','MC', sc.name, dist.name, 'results');
        all_done = true;
        for pi = 1:3
            if ~isfile(fullfile(res_dir, sprintf('MC_%s.mat', PARAMS{pi})))
                all_done = false; break;
            end
        end
        if all_done, completed{end+1} = {sc, dist}; end %#ok<AGROW>
    end
end
fprintf('Found %d completed scenario/dist combos.\n\n', numel(completed));

%% ── Storage for summary figure ───────────────────────────────────────────────
sum_L2   = {cell(numel(completed),1), cell(numel(completed),1), cell(numel(completed),1)};
sum_label = cell(numel(completed), 1);

%% ── Main loop ────────────────────────────────────────────────────────────────
for ci = 1:numel(completed)
    sc   = completed{ci}{1};
    dist = completed{ci}{2};
    if ~isempty(sc_target) && ~(strcmp(sc.name,sc_target) && strcmp(dist.name,dist_target))
        continue;
    end
    sc_name   = sc.name;
    dist_name = dist.name;
    dist_type = lower(dist.type);

    fprintf('────────────────────────────────────────\n');
    fprintf('  PCE | %s | %s\n', sc_name, dist_name);
    fprintf('────────────────────────────────────────\n');

    mc_res_dir = fullfile('Methods','MC',  sc_name, dist_name, 'results');
    pce_dir    = fullfile('Methods','PCE', sc_name, dist_name);
    if skipIfDone(fullfile(pce_dir,'results','pce_results.mat'), ...
                  sprintf('%s / %s', sc_name, dist_name))
        continue;
    end
    fig_dir = fullfile(pce_dir, 'figures');
    out_dir = fullfile(pce_dir, 'results');
    for d = {fig_dir, out_dir}
        if ~exist(d{1},'dir'), mkdir(d{1}); end
    end

    FOM = getFOM(cfg, sc_name);

    %% Variance bounds for the TARGET distribution
    B_target = computeVarianceBounds(cfg, dist);
    sig_tgt  = [B_target.sig_freq, B_target.sig_zS, B_target.sig_svp];

    %% Regenerate LHS samples for the BASE distribution (fixed seed → reproducible)
    S_base = lhsSample(N, cfg, cfg.distributions(base_idx), B_base);

    % Standardised xi ∈ N(0,1) for base distribution
    xi_all_base = [(S_base.freq - freq0) / B_base.sig_freq, ...
                   (S_base.zS   - zS0)   / B_base.sig_zS, ...
                    S_base.svp           / B_base.sig_svp];   % [N × 3]

    % Gamma norms for analytic Var: <He_n, He_n>_{N(0,1)} = n!
    gamma_vec = factorial(0:MAX_ORDER)';   % [MAX_ORDER+1 × 1]

    %% Per-parameter PCE fit ──────────────────────────────────────────────────
    L2_mat       = zeros(MAX_ORDER, 3);
    LOO_EX_mat   = zeros(MAX_ORDER, 3);
    LOO_Var_mat  = zeros(MAX_ORDER, 3);
    LOO_comb_mat = zeros(MAX_ORDER, 3);
    kappa_vec    = zeros(1, 3);
    Kstar        = zeros(1, 3);
    C_all        = cell(1, 3);
    VarK_best    = cell(1, 3);
    Pdetect_best = cell(1, 3);

    Nz = []; Nr = [];

    for pi = 1:3
        param = PARAMS{pi};

        %% Load TL cube from Normal_10pct cache
        tl_file = fullfile('Cache', sc_name, base_name, sprintf('TL_%s.mat', param));
        if ~isfile(tl_file)
            % fallback: try old per-distribution path for backwards compat
            tl_file = fullfile('Cache', sc_name, dist_name, sprintf('TL_%s.mat', param));
        end
        if ~isfile(tl_file)
            fprintf('  [SKIP] TL cache missing for %s: %s\n', param, tl_file);
            L2_mat(:,pi) = NaN;
            continue;
        end
        tmp    = load(tl_file, 'TL_save');
        TL_all = double(tmp.TL_save(:,:,1:N));
        if isempty(Nz), [Nz, Nr, ~] = size(TL_all); end
        Npix = Nz * Nr;

        %% Load MC reference: Var and P_detect for TARGET distribution
        mc_file = fullfile(mc_res_dir, sprintf('MC_%s.mat', param));
        mc_data = load(mc_file, 'MC_Var');
        Var_MC  = mc_data.MC_Var(:);   % [Npix × 1]
        std_VM  = max(std(Var_MC), 1e-10);
        clear mc_data;

        %% Standardised xi for base distribution
        xi = xi_all_base(:, pi);   % [N × 1]  ~  N(0,1)

        %% Build basis and QR decomposition
        Phi    = pce_basis(xi, MAX_ORDER, dist_type);           % [N × MAX_ORDER+1]
        TL_mat = reshape(permute(TL_all, [3 1 2]), N, Npix);   % [N × Npix]
        kappa_vec(pi) = cond(Phi);

        [Q, R_qr] = qr(Phi, 0);   % thin QR: Q [N×MAX_ORDER+1], R [MAX_ORDER+1×MAX_ORDER+1]
        QTL       = Q' * TL_mat;  % [MAX_ORDER+1 × Npix]

        % Ridge-regularised coefficients (full order, for dist figure)
        lambda    = 1e-10 * trace(Phi'*Phi);
        C_all{pi} = (Phi'*Phi + lambda*eye(MAX_ORDER+1)) \ (Phi'*TL_mat);

        % Distribution rescaling ratio:  σ_target / σ_base
        ratio = sig_tgt(pi) / sig_base(pi);

        %% Convergence metrics (nested QR projections, rank-1 updates) ─────────
        % Analytic variance: Var_K = Σ_{n=1}^K c_n² · n! · ratio^{2n}
        %
        % PCE mean (zeroth coeff): EX_K = c₀ = QTL(1,:)/sqrt(N) = mean(TL_mat)
        % → constant for all K; EX L2 metric is distribution-specific
        %
        h_cumsum = cumsum(Q.^2, 2);      % [N × MAX_ORDER+1]
        y_hat_K  = zeros(N, Npix);       % accumulated PCE predictions (base dist)

        for K = 1:MAX_ORDER
            % Rank-1 update: add K-th order contribution
            y_hat_K    = y_hat_K + Q(:,K+1) * QTL(K+1,:);
            c_K_accum(K+1,:) = QTL(K+1,:);   % = Q_{K+1}' TL (unnormalized)

            % Recover PCE coefficient c_n (in original He_n basis) via R:
            % c = R^{-1} * QTL(1:K+1,:)
            c_K = R_qr(1:K+1,1:K+1) \ QTL(1:K+1,:);   % [K+1 × Npix]

            % Analytic Var under TARGET distribution (rescaling formula)
            ratio_n2 = (ratio .^ (2*(1:K)'))';            % [K × 1] scaled factors
            Var_K_pix = sum((c_K(2:end,:).^2) .* (gamma_vec(2:K+1) .* ratio_n2), 1); % [1 × Npix]
            Var_K_pix = max(Var_K_pix, 0);

            % L2 Var error vs MC
            L2_mat(K, pi) = sqrt(mean((Var_K_pix' - Var_MC).^2)) / std_VM;

            % LOO-CV for EX prediction (analytic hat matrix)
            h_K    = h_cumsum(:, K+1);
            r_K    = TL_mat - y_hat_K;
            loo_K  = mean((r_K ./ max(1 - h_K, 1e-10)).^2, 1);
            LOO_EX_mat(K, pi) = mean(loo_K);

            % LOO Var: variance of LOO predictions, analytically scaled to target dist
            y_hat_loo   = y_hat_K + (h_K ./ max(1 - h_K, 1e-10)) .* r_K;  % [N × Npix]
            % Scale LOO variance by ratio² (first-order approximation)
            Var_loo_K   = max(var(y_hat_loo, 0, 1) * ratio^2, 0);          % [1 × Npix]
            LOO_Var_mat(K, pi) = mean((Var_loo_K' - Var_MC).^2) / max(var(Var_MC), eps);
        end

        % Normalise and combine LOO criteria
        LOO_EX_n  = LOO_EX_mat(:,pi) / max(LOO_EX_mat(1,pi), eps);
        LOO_Var_n = LOO_Var_mat(:,pi) / max(LOO_Var_mat(1,pi), eps);
        LOO_comb_mat(:,pi) = alpha_loo * LOO_Var_n + (1-alpha_loo) * LOO_EX_n;

        [~, Kstar(pi)] = min(LOO_comb_mat(:,pi));

        %% Best-order Var map and detection map ────────────────────────────────
        Ks = Kstar(pi);
        c_best  = R_qr(1:Ks+1,1:Ks+1) \ QTL(1:Ks+1,:);
        ratio_n2_best = (ratio .^ (2*(1:Ks)'))';
        Var_best_pix = sum((c_best(2:end,:).^2) .* (gamma_vec(2:Ks+1) .* ratio_n2_best), 1);
        Var_best_pix = max(Var_best_pix, 0);
        VarK_best{pi}    = reshape(Var_best_pix, Nz, Nr);

        % P(detect) = P(TL < FOM) via Gaussian approximation
        EX_best = c_best(1,:);   % [1 × Npix]
        Pdetect_best{pi} = reshape( ...
            normcdf(FOM, EX_best, sqrt(Var_best_pix + eps)), Nz, Nr);

        fprintf('  [%s] K*(LOO-comb)=%d  cond=%.1e  L2_Var@K=1: %.3f  ratio=%.2f\n', ...
                param, Kstar(pi), kappa_vec(pi), L2_mat(1,pi), ratio);
    end

    %% Save results ────────────────────────────────────────────────────────────
    mc_ref = load(fullfile(mc_res_dir, 'MC_freq.mat'), 'r_km','z_m');
    r_km  = mc_ref.r_km;
    z_m   = mc_ref.z_m;
    save(fullfile(out_dir, 'pce_results.mat'), ...
         'C_all','L2_mat','LOO_EX_mat','LOO_Var_mat','LOO_comb_mat', ...
         'Kstar','kappa_vec','VarK_best','Pdetect_best', ...
         'r_km','z_m','MAX_ORDER','dist_name','sc_name','FOM', '-v7.3');

    %% ── Figure 1: convergence panels (LOO-CV + L2 + scatter) ────────────────
    colors = lines(3);
    fig    = figure('Position', [50 50 1600 640]);

    % ── (1,1): Combined LOO-CV vs order ──────────────────────────────────────
    ax_loo = subplot(2,4,1);
    for pi = 1:3
        plot(ax_loo, 1:MAX_ORDER, LOO_comb_mat(:,pi), '-o', ...
             'Color', colors(pi,:), 'LineWidth', 1.6, 'MarkerSize', 4, ...
             'DisplayName', PARAMS{pi});
        hold(ax_loo, 'on');
        xline(ax_loo, Kstar(pi), '--', 'Color', colors(pi,:)*0.7, 'LineWidth', 1, ...
              'HandleVisibility','off');
    end
    set(ax_loo, 'XTick', max(1, 1:2:MAX_ORDER));
    xlabel(ax_loo, 'PCE Order K');
    ylabel(ax_loo, 'LOO_{combined} (norm.)');
    title(ax_loo, sprintf('LOO Combined  [K*=vert]\n%s | %s', sc_name, dist_name), ...
          'FontSize', 8, 'Interpreter', 'none');
    legend(ax_loo, 'Location','northeast','FontSize',7);
    grid(ax_loo, 'on');

    % ── (1,2): L2 Var error vs order ─────────────────────────────────────────
    ax_l2 = subplot(2,4,2);
    for pi = 1:3
        plot(ax_l2, 1:MAX_ORDER, L2_mat(:,pi)*100, '-o', ...
             'Color', colors(pi,:), 'LineWidth', 1.5, 'MarkerSize', 4, ...
             'DisplayName', PARAMS{pi});
        hold(ax_l2, 'on');
    end
    yline(ax_l2, 10, 'k--', '10%', 'LineWidth', 1.2, 'LabelHorizontalAlignment','left');
    set(ax_l2, 'XTick', max(1, 1:2:MAX_ORDER));
    xlabel(ax_l2, 'PCE Order K');
    ylabel(ax_l2, 'Relative L2 Var error (%)');
    title(ax_l2, sprintf('L2 Variance\n%s | %s', sc_name, dist_name), ...
          'FontSize', 8, 'Interpreter', 'none');
    legend(ax_l2, 'Location','northeast','FontSize',7);
    grid(ax_l2, 'on');

    % ── (1,3): LOO_EX vs order ────────────────────────────────────────────────
    ax_lex = subplot(2,4,3);
    for pi = 1:3
        loo_ex_n = LOO_EX_mat(:,pi) / max(LOO_EX_mat(1,pi),eps);
        plot(ax_lex, 1:MAX_ORDER, loo_ex_n, '-o', ...
             'Color', colors(pi,:), 'LineWidth', 1.5, 'MarkerSize', 4, ...
             'DisplayName', PARAMS{pi});
        hold(ax_lex, 'on');
    end
    set(ax_lex, 'XTick', max(1, 1:2:MAX_ORDER));
    xlabel(ax_lex, 'PCE Order K');
    ylabel(ax_lex, 'LOO_{EX} (norm.)');
    title(ax_lex, sprintf('LOO EX (prediction)\n%s | %s', sc_name, dist_name), ...
          'FontSize', 8, 'Interpreter', 'none');
    legend(ax_lex, 'Location','northeast','FontSize',7);
    grid(ax_lex, 'on');

    % ── (1,4): LOO_Var vs order ───────────────────────────────────────────────
    ax_lv = subplot(2,4,4);
    for pi = 1:3
        loo_var_n = LOO_Var_mat(:,pi) / max(LOO_Var_mat(1,pi),eps);
        plot(ax_lv, 1:MAX_ORDER, loo_var_n, '-o', ...
             'Color', colors(pi,:), 'LineWidth', 1.5, 'MarkerSize', 4, ...
             'DisplayName', PARAMS{pi});
        hold(ax_lv, 'on');
    end
    set(ax_lv, 'XTick', max(1, 1:2:MAX_ORDER));
    xlabel(ax_lv, 'PCE Order K');
    ylabel(ax_lv, 'LOO_{Var} (norm.)');
    title(ax_lv, sprintf('LOO Var (variance)\n%s | %s', sc_name, dist_name), ...
          'FontSize', 8, 'Interpreter', 'none');
    legend(ax_lv, 'Location','northeast','FontSize',7);
    grid(ax_lv, 'on');

    % ── (2,1-3): scatter Var_{K*} vs Var_MC per param ────────────────────────
    for pi = 1:3
        ax_r = subplot(2, 4, pi+4);
        if isempty(VarK_best{pi})
            title(ax_r, sprintf('%s (no data)', PARAMS{pi}), 'FontSize',8,'Interpreter','none');
            continue;
        end
        mc_pi = load(fullfile(mc_res_dir, sprintf('MC_%s.mat', PARAMS{pi})), 'MC_Var');
        vmc   = mc_pi.MC_Var(:);
        vpce  = VarK_best{pi}(:);

        ok       = vmc > prctile(vmc, 5);
        xlim_hi  = prctile(vmc(ok),  99);
        ylim_hi  = prctile(vpce(ok), 99);
        lim      = max(xlim_hi, ylim_hi);

        scatter(ax_r, vmc(ok), vpce(ok), 3, 'filled', ...
                'MarkerFaceAlpha', 0.25, 'MarkerFaceColor', colors(pi,:));
        hold(ax_r,'on');
        plot(ax_r, [0 lim], [0 lim], 'k-', 'LineWidth', 1.5);
        xlim(ax_r, [0 xlim_hi]); ylim(ax_r, [0 ylim_hi]);

        ok2 = ok & isfinite(vmc) & isfinite(vpce);
        R2  = 0;
        if sum(ok2) > 2, R2 = corr(vmc(ok2), vpce(ok2))^2; end

        xlabel(ax_r, 'Var_{MC}'); ylabel(ax_r, sprintf('Var_{K=%d}', Kstar(pi)));
        title(ax_r, sprintf('Var Scatter | %s\nR^2=%.3f K*=%d', PARAMS{pi}, R2, Kstar(pi)), ...
              'FontSize', 8, 'Interpreter','none');
        grid(ax_r,'on'); axis(ax_r,'equal');
    end

    % ── (2,4): scatter P_detect PCE vs MC ────────────────────────────────────
    ax_pd = subplot(2, 4, 8);
    hold(ax_pd,'on');
    for pi = 1:3
        if isempty(Pdetect_best{pi}), continue; end
        mc_pd = load(fullfile(mc_res_dir, sprintf('MC_%s.mat', PARAMS{pi})), 'MC_P_detect');
        scatter(ax_pd, mc_pd.MC_P_detect(:), Pdetect_best{pi}(:), 3, 'filled', ...
                'MarkerFaceAlpha', 0.2, 'MarkerFaceColor', colors(pi,:), ...
                'DisplayName', PARAMS{pi});
    end
    plot(ax_pd, [0 1], [0 1], 'k-', 'LineWidth', 1.5, 'HandleVisibility','off');
    xlim(ax_pd, [0 1]); ylim(ax_pd, [0 1]);
    xlabel(ax_pd, 'P_{detect} MC'); ylabel(ax_pd, 'P_{detect} PCE');
    title(ax_pd, sprintf('P(detect) scatter\nall params @ K*'), ...
          'FontSize', 8, 'Interpreter','none');
    legend(ax_pd, 'Location','northwest','FontSize',7);
    grid(ax_pd,'on'); axis(ax_pd,'equal');

    sgtitle(sprintf('PCE Variance Convergence | %s | %s', sc_name, dist_name), ...
            'FontSize', 10, 'Interpreter','none');
    saveFigPNG(fig, fullfile(fig_dir, sprintf('pce_convergence_%s_%s', sc_name, dist_name)));
    close(fig);
    fprintf('  Convergence figure saved.\n');

    %% ── Figure 2: Order-relevance (cumulative Var improvement) ──────────────
    fig_rel = figure('Position', [50 50 1000 380]);
    for pi = 1:3
        ax_rel = subplot(1, 3, pi);
        l2_v   = L2_mat(:, pi);
        if all(isnan(l2_v)), axis(ax_rel,'off'); continue; end
        l2_start = l2_v(1);
        l2_end   = min(l2_v);
        total_imp = max(l2_start - l2_end, eps);
        relevance = max(l2_start - l2_v, 0) / total_imp * 100;

        plot(ax_rel, 1:MAX_ORDER, relevance, 'b-o', 'LineWidth',1.5,'MarkerSize',5);
        hold(ax_rel,'on');
        yline(ax_rel, 90, 'r--', '90%', 'LineWidth',1);
        xline(ax_rel, Kstar(pi), 'g--', sprintf('K*=%d', Kstar(pi)), 'LineWidth',1.2);
        hold(ax_rel,'off');
        set(ax_rel, 'XTick', max(1, 1:2:MAX_ORDER));
        xlabel(ax_rel, 'PCE Order K');
        ylabel(ax_rel, '% Var improvement (cumul.)');
        title(ax_rel, sprintf('%s\nKstar=%d', PARAMS{pi}, Kstar(pi)), ...
              'FontSize',9,'Interpreter','none');
        grid(ax_rel,'on');
    end
    sgtitle(sprintf('PCE Order Relevance — cumulative Var improvement | %s | %s', ...
                    sc_name, dist_name), 'Interpreter','none','FontSize',10);
    saveFigPNG(fig_rel, fullfile(fig_dir, 'pce_order_relevance'));
    close(fig_rel);

    %% ── Figure 3: Distribution comparison (Normal_10pct only) ───────────────
    if strcmp(dist_name, base_name)
        dist_fig_path = fullfile(fig_dir, sprintf('pce_dist_%s_%s.png', sc_name, dist_name));
        if ~isfile(dist_fig_path)
            M = 2000;  rng(42);
            xi_new = randn(M, 1);

            r_idx = max(1, min(Nr, round([0.25 0.50 0.75] * Nr)));
            z_idx = max(1, round(Nz / 2));

            fig_d = figure('Position', [50 50 1100 800]);
            for pi = 1:3
                if isempty(C_all{pi}), continue; end
                tl_f  = fullfile('Cache', sc_name, base_name, sprintf('TL_%s.mat', PARAMS{pi}));
                if ~isfile(tl_f), tl_f = fullfile('Cache', sc_name, dist_name, sprintf('TL_%s.mat', PARAMS{pi})); end
                if ~isfile(tl_f), continue; end
                tmp_pi  = load(tl_f, 'TL_save');
                TL_pi   = double(tmp_pi.TL_save(:,:,1:N));
                TL_pi_m = reshape(permute(TL_pi, [3 1 2]), N, Nz*Nr);
                clear TL_pi tmp_pi;

                Phi_new = pce_basis(xi_new, MAX_ORDER, dist_type);
                y_pce   = Phi_new * C_all{pi};

                for ci_px = 1:3
                    ax_d = subplot(3, 3, (pi-1)*3 + ci_px);
                    pix  = (z_idx - 1)*Nr + r_idx(ci_px);

                    tl_mc  = TL_pi_m(:, pix);
                    tl_pce = y_pce(:, pix);
                    all_v  = [tl_mc; tl_pce];
                    edges  = linspace(prctile(all_v,1), prctile(all_v,99), 22);

                    histogram(ax_d, tl_mc,  edges, 'Normalization','pdf', ...
                              'FaceColor', colors(pi,:), 'FaceAlpha', 0.65, ...
                              'DisplayName', sprintf('MC (N=%d)', N));
                    hold(ax_d,'on');
                    histogram(ax_d, tl_pce, edges, 'Normalization','pdf', ...
                              'FaceColor', [0.9 0.5 0.1], 'FaceAlpha', 0.45, ...
                              'DisplayName', sprintf('PCE (M=%d)', M));
                    xlabel(ax_d,'TL (dB)','FontSize',7);
                    if ci_px == 1, ylabel(ax_d,'PDF','FontSize',7); end
                    title(ax_d, sprintf('%s  r=%.0fkm', PARAMS{pi}, r_km(r_idx(ci_px))), ...
                          'FontSize',7,'Interpreter','none');
                    if pi == 1 && ci_px == 1
                        legend(ax_d,'Location','best','FontSize',6);
                    end
                    grid(ax_d,'on');
                end
            end
            sgtitle(sprintf('PCE vs MC Distribution | %s | %s\nrows=params, cols=range slices z=%dm', ...
                            sc_name, dist_name, z_m(z_idx)), 'FontSize',9);
            saveFigPNG(fig_d, strrep(dist_fig_path,'.png',''));
            close(fig_d);
            fprintf('  Dist figure saved.\n');
        end
    end

    %% Store for summary ───────────────────────────────────────────────────────
    for pi = 1:3
        sum_L2{pi}{ci} = L2_mat(:, pi);
    end
    sum_label{ci} = sprintf('%s / %s', sc_name, dist_name);
end

%% ── Summary figure: all combos × all params ─────────────────────────────────
n_combos = numel(completed);
cmap     = lines(n_combos);

fig_sum = figure('Position', [50 50 1400 900]);
for pi = 1:3
    ax = subplot(3, 1, pi);
    for ci = 1:n_combos
        if isempty(sum_L2{pi}{ci}), continue; end
        plot(ax, 1:MAX_ORDER, sum_L2{pi}{ci}*100, '-o', 'Color', cmap(ci,:), ...
             'LineWidth', 1.2, 'MarkerSize', 4, 'DisplayName', sum_label{ci});
        hold(ax,'on');
    end
    yline(ax, 10, 'k--', 'LineWidth', 1.2);
    set(ax, 'XTick', max(1, 1:2:MAX_ORDER));
    ylabel(ax, 'Relative L2 (%)');
    title(ax, sprintf('param = %s', PARAMS{pi}), 'FontSize',9,'Interpreter','none');
    if pi == 1
        legend(ax,'Location','northeast','FontSize',6,'NumColumns',2);
    end
    if pi == 3, xlabel(ax,'PCE Order K'); end
    grid(ax,'on');
end
sgtitle('PCE Convergence — All Scenarios × Distributions', 'FontSize',11);

pce_root = fullfile('Methods','PCE');
if ~exist(pce_root,'dir'), mkdir(pce_root); end
saveFigPNG(fig_sum, fullfile(pce_root, 'pce_summary_all'));
close(fig_sum);
fprintf('\nSaved summary: Methods/PCE/pce_summary_all.png\n');

%% ── Console K* table ─────────────────────────────────────────────────────────
fprintf('\n  K* (argmin combined LOO):\n');
fprintf('  %-35s  freq   zS    svp\n', 'Scenario / Distribution');
fprintf('  %s\n', repmat('-',1,65));
for ci = 1:numel(completed)
    sc   = completed{ci}{1};
    dist = completed{ci}{2};
    pf   = fullfile('Methods','PCE', sc.name, dist.name, 'results', 'pce_results.mat');
    if ~isfile(pf), continue; end
    r = load(pf, 'Kstar');
    fprintf('  %-35s  %-5d  %-5d  %-5d\n', ...
            sprintf('%s / %s', sc.name, dist.name), r.Kstar(1), r.Kstar(2), r.Kstar(3));
end
fprintf('\n=== run_pce.m complete ===\n');
end   % function run_pce
