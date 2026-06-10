%% run_pce.m — PCE variance estimation from cached MC TL cubes
%
% For each completed scenario × distribution, fits Polynomial Chaos Expansion
% coefficients using the existing MC TL samples (no new Bellhop runs).
% Produces convergence figures and saves PCE coefficients as reusable outputs.
%
% Method:
%   TL(xi) ≈ Σ_{n=0}^K c_n Φ_n(xi)   (1-D PCE per parameter)
%   Coefficients via least-squares regression on N LHS samples.
%   Var_K = Σ_{n=1}^K c_n^2 γ_n   (Sobol variance decomposition)
%
% Outputs per scenario × distribution:
%   Methods/PCE/<scen>/<dist>/results/pce_results.mat
%   Methods/PCE/<scen>/<dist>/figures/pce_convergence_<scen>_<dist>.png
% Summary:
%   Methods/PCE/pce_summary_all.png

clear; close all; clc; warning('off');
try, cd(fileparts(mfilename('fullpath'))); catch; end

set(0, 'DefaultFigureVisible', 'off');
addpath(genpath('Shared_Utils'));

cfg   = loadConfig();
freq0 = cfg.nominal.freq_Hz;
zS0   = cfg.nominal.zS_m;
N     = cfg.MC.N;

MAX_ORDER = 10;
PARAMS    = {'freq', 'zS', 'svp'};
PARAM_COL = [1, 2, 3];   % lhsSample column index for each param

%% ── Identify completed combos ────────────────────────────────────────────────
fprintf('=== PCE from cached MC data (order 1..%d) ===\n', MAX_ORDER);

completed = {};
for si = 1:numel(cfg.scenarios)
    sc = cfg.scenarios(si);
    for di = 1:numel(cfg.distributions)
        dist = cfg.distributions(di);
        % Check that all 3 single-param MC result files exist
        res_dir = fullfile('Methods','MC', sc.name, dist.name, 'results');
        all_done = true;
        for pi = 1:3
            if ~isfile(fullfile(res_dir, sprintf('MC_%s.mat', PARAMS{pi})))
                all_done = false; break;
            end
        end
        if all_done
            completed{end+1} = {sc, dist}; %#ok<AGROW>
        end
    end
end
fprintf('Found %d completed scenario/dist combos.\n\n', numel(completed));

%% ── Storage for summary figure ───────────────────────────────────────────────
% sum_L1{pi}{combo_idx} = L1 vector [MAX_ORDER x 1]
sum_L1    = {cell(numel(completed),1), cell(numel(completed),1), cell(numel(completed),1)};
sum_label = cell(numel(completed), 1);
sum_scen  = cell(numel(completed), 1);

%% ── Main loop ────────────────────────────────────────────────────────────────
for ci = 1:numel(completed)
    sc   = completed{ci}{1};
    dist = completed{ci}{2};
    sc_name   = sc.name;
    dist_name = dist.name;
    dist_type = lower(dist.type);   % 'normal' or 'uniform'

    fprintf('────────────────────────────────────────\n');
    fprintf('  PCE | %s | %s\n', sc_name, dist_name);
    fprintf('────────────────────────────────────────\n');

    res_dir = fullfile('Methods','MC',  sc_name, dist_name, 'results');
    pce_dir = fullfile('Methods','PCE', sc_name, dist_name);
    fig_dir = fullfile(pce_dir, 'figures');
    out_dir = fullfile(pce_dir, 'results');
    for d = {fig_dir, out_dir}
        if ~exist(d{1},'dir'), mkdir(d{1}); end
    end

    %% Compute variance bounds (needed to regenerate sample values)
    B = computeVarianceBounds(cfg, dist);

    %% Regenerate LHS samples (same seed → identical to original MC run)
    S = lhsSample(N, cfg, dist, B);   % S.freq [N×1], S.zS [N×1], S.svp [N×1]

    % gamma not needed — variance computed directly from fitted values (see below)

    %% Standardised sample values per parameter
    switch dist_type
        case 'normal'
            xi_all = [(S.freq - freq0) / B.sig_freq, ...
                      (S.zS   - zS0)   / B.sig_zS, ...
                       S.svp           / B.sig_svp];   % [N × 3]
        case 'uniform'
            xi_all = [(S.freq - freq0) / (diff(B.freq_bnd)/2), ...
                      (S.zS   - zS0)   / (diff(B.zS_bnd)/2), ...
                       S.svp           / (diff(B.svp_bnd)/2)];   % [N × 3]
    end

    %% Per-parameter PCE fit ──────────────────────────────────────────────────
    L1_mat    = zeros(MAX_ORDER, 3);   % variance relative error vs MC
    LOO_mat   = zeros(MAX_ORDER, 3);   % analytic LOO-CV MSE (main K* criterion)
    kappa_vec = zeros(1, 3);           % condition number of Phi
    Kstar     = zeros(1, 3);           % K* from LOO argmin
    Kstar_L1  = zeros(1, 3);           % K* from old L1<10% rule (kept for comparison)
    C_all     = cell(1, 3);
    VarK_best = cell(1, 3);

    Nz = []; Nr = [];   % filled on first param

    for pi = 1:3
        param = PARAMS{pi};

        %% Load TL cube
        tl_file = fullfile('Cache', sc_name, dist_name, sprintf('TL_%s.mat', param));
        if ~isfile(tl_file)
            fprintf('  [SKIP] TL cache missing: %s\n', tl_file);
            L1_mat(:,pi) = NaN;
            continue;
        end
        tmp    = load(tl_file, 'TL_save');
        TL_all = double(tmp.TL_save(:,:,1:N));   % [Nz × Nr × N]
        if isempty(Nz), [Nz, Nr, ~] = size(TL_all); end

        %% Load MC reference variance
        mc_file = fullfile(res_dir, sprintf('MC_%s.mat', param));
        mc_res  = load(mc_file, 'MC_Var');
        Var_MC  = mc_res.MC_Var(:);   % [Nz*Nr × 1]
        mean_VM = max(mean(Var_MC), 1e-10);

        %% Standardised xi for this param
        xi = xi_all(:, pi);   % [N × 1]

        %% Build basis, QR decomposition, fit PCE coefficients
        Phi    = pce_basis(xi, MAX_ORDER, dist_type);          % [N × MAX_ORDER+1]
        TL_mat = reshape(permute(TL_all, [3 1 2]), N, Nz*Nr);  % [N × Nz*Nr]

        % Condition number check — cond(Phi) > 1e8 means basis is numerically unreliable
        kappa_vec(pi) = cond(Phi);

        % QR for nested LOO/variance — Q(:,1:K+1)*QTL(1:K+1,:) is the K-th order projection.
        [Q, R_qr] = qr(Phi, 0);   % thin QR: Q [N×MAX_ORDER+1], R [MAX_ORDER+1×MAX_ORDER+1]
        QTL = Q' * TL_mat;         % [MAX_ORDER+1 × Nz*Nr]

        % Ridge-regularized coefficients for saving (stabilises high-order terms).
        % lambda is ~1e-10 relative to operator scale — negligible bias, reduces variance.
        lambda   = 1e-10 * trace(Phi'*Phi);
        C_all{pi} = (Phi'*Phi + lambda*eye(MAX_ORDER+1)) \ (Phi'*TL_mat);

        %% Convergence metrics via nested QR projections
        for K = 1:MAX_ORDER
            y_hat_K   = Q(:,1:K+1) * QTL(1:K+1,:);    % K-th order fit [N × Nz*Nr]
            Var_K_pix = var(y_hat_K, 0, 1);

            % L1 variance metric (kept for reference)
            L1_mat(K, pi) = mean(abs(Var_K_pix - Var_MC')) / mean_VM;

            % Analytic LOO-CV: e_LOO_i = r_i / (1 - h_ii), h_ii = ||Q(i,1:K+1)||^2
            % No refitting required — exact for OLS, excellent approximation overall.
            h_K   = sum(Q(:,1:K+1).^2, 2);              % hat-matrix diagonal [N×1]
            r_K   = TL_mat - y_hat_K;                    % residuals [N × Nz*Nr]
            loo_K = mean((r_K ./ (1 - h_K)).^2, 1);     % per-pixel LOO MSE [1×Nz*Nr]
            LOO_mat(K, pi) = mean(loo_K);                % scalar summary
        end

        %% K* — argmin LOO (primary) and old L1<10% (secondary, for comparison)
        [~, Kstar(pi)]    = min(LOO_mat(:,pi));
        idx_l1            = find(L1_mat(:,pi) < 0.10, 1);
        Kstar_L1(pi)      = idx_l1;  % may be empty → 0

        %% Store best-order Var map (using LOO-selected K*)
        y_hat_best = Q(:,1:Kstar(pi)+1) * QTL(1:Kstar(pi)+1,:);
        VarK_best{pi} = reshape(var(y_hat_best, 0, 1), Nz, Nr);

        fprintf('  [%s] K*(LOO)=%d  K*(L1)=%s  cond=%.1e  LOO@K=1: %.3f  L1@K=1: %.3f\n', ...
                param, Kstar(pi), fmt_kstar(Kstar_L1(pi)), kappa_vec(pi), ...
                LOO_mat(1,pi), L1_mat(1,pi));
    end

    %% Save results ────────────────────────────────────────────────────────────
    mc_ref = load(fullfile(res_dir, 'MC_freq.mat'), 'r_grid','z_grid','r_km','z_m');
    r_km  = mc_ref.r_km;
    z_m   = mc_ref.z_m;
    save(fullfile(out_dir, 'pce_results.mat'), ...
         'C_all','L1_mat','LOO_mat','Kstar','Kstar_L1','kappa_vec', ...
         'VarK_best','r_km','z_m','MAX_ORDER','dist_name','sc_name', '-v7.3');

    %% Figure: L1 convergence (left) + 3 scatter plots (right) ───────────────
    fig = figure('Position', [50 50 1400 400]);

    % ── Left panel: LOO-CV (primary) + L1 (dashed reference) ────────────────
    ax_l   = subplot(1,4,1);
    colors = lines(3);
    for pi = 1:3
        % Normalise LOO to K=1 so curves start at 1 regardless of units
        loo_norm = LOO_mat(:,pi) / max(LOO_mat(1,pi), 1e-30);
        plot(ax_l, 1:MAX_ORDER, loo_norm, '-o', ...
             'Color', colors(pi,:), 'LineWidth', 1.8, 'MarkerSize', 5, ...
             'DisplayName', PARAMS{pi});
        hold(ax_l, 'on');
        % Mark LOO-selected K* with vertical tick
        xline(ax_l, Kstar(pi), '--', 'Color', colors(pi,:) * 0.7, 'LineWidth', 1, ...
              'HandleVisibility', 'off');
    end
    % L1 curves as thin dashed reference (secondary y not needed — both 0-1 after norm)
    for pi = 1:3
        l1_norm = L1_mat(:,pi) / max(L1_mat(1,pi), 1e-30);
        plot(ax_l, 1:MAX_ORDER, l1_norm, '--', ...
             'Color', colors(pi,:), 'LineWidth', 0.8, ...
             'HandleVisibility', 'off');
    end
    yline(ax_l, 0, 'k-', 'LineWidth', 0.5);
    set(ax_l, 'XTick', 1:MAX_ORDER);
    xlabel(ax_l, 'PCE Order K');
    ylabel(ax_l, 'Normalised error (rel. K=1)');
    title(ax_l, sprintf('LOO-CV (solid) + L1 (dash)\n%s | %s', sc_name, dist_name), ...
          'FontSize', 8, 'Interpreter', 'none');
    legend(ax_l, 'Location', 'northeast', 'FontSize', 7);
    grid(ax_l, 'on');

    % ── Right panels: scatter Var_{K*} vs Var_MC for each param ─────────────
    for pi = 1:3
        ax_r = subplot(1, 4, pi+1);
        if isempty(VarK_best{pi})
            title(ax_r, sprintf('%s\n(no data)', PARAMS{pi}), 'FontSize', 8, 'Interpreter','none');
            continue;
        end
        mc_pi = load(fullfile(res_dir, sprintf('MC_%s.mat', PARAMS{pi})), 'MC_Var');
        vmc   = mc_pi.MC_Var(:);
        vpce  = VarK_best{pi}(:);

        ok   = vmc > prctile(vmc, 5);
        xlim_hi = prctile(vmc(ok),  99);
        ylim_hi = prctile(vpce(ok), 99);
        lim  = max(xlim_hi, ylim_hi);

        scatter(ax_r, vmc(ok), vpce(ok), 3, 'filled', ...
                'MarkerFaceAlpha', 0.25, 'MarkerFaceColor', colors(pi,:));
        hold(ax_r, 'on');
        plot(ax_r, [0 lim], [0 lim], 'k-', 'LineWidth', 1.5);
        xlim(ax_r, [0 xlim_hi]); ylim(ax_r, [0 ylim_hi]);

        ok2 = ok & isfinite(vmc) & isfinite(vpce);
        R2  = corr(vmc(ok2), vpce(ok2))^2;

        xlabel(ax_r, 'Var_{MC}');
        ylabel(ax_r, sprintf('Var_{K=%s}', fmt_kstar(Kstar(pi))));
        title(ax_r, sprintf('Scatter | %s\nR^2=%.3f', PARAMS{pi}, R2), ...
              'FontSize', 8, 'Interpreter', 'none');
        grid(ax_r, 'on');
        axis(ax_r, 'equal');
    end

    sgtitle(sprintf('PCE Variance Convergence | %s | %s', sc_name, dist_name), ...
            'FontSize', 10, 'Interpreter', 'none');

    out_fig = fullfile(fig_dir, sprintf('pce_convergence_%s_%s.png', sc_name, dist_name));
    saveas(fig, out_fig);
    close(fig);
    fprintf('  Saved: %s\n', out_fig);

    %% Store for summary ───────────────────────────────────────────────────────
    for pi = 1:3
        sum_L1{pi}{ci} = L1_mat(:, pi);
    end
    sum_label{ci} = sprintf('%s / %s', sc_name, dist_name);
    sum_scen{ci}  = sc_name;
end

%% ── Summary figure: all combos × all params ─────────────────────────────────
scen_names = unique(sum_scen(~cellfun(@isempty, sum_scen)));
n_combos   = numel(completed);
cmap       = lines(n_combos);

fig_sum = figure('Position', [50 50 1400 900]);
for pi = 1:3
    ax = subplot(3, 1, pi);
    for ci = 1:n_combos
        if isempty(sum_L1{pi}{ci}), continue; end
        l1v = sum_L1{pi}{ci} * 100;
        plot(ax, 1:MAX_ORDER, l1v, '-o', 'Color', cmap(ci,:), ...
             'LineWidth', 1.2, 'MarkerSize', 4, 'DisplayName', sum_label{ci});
        hold(ax, 'on');
    end
    yline(ax, 10, 'k--', 'LineWidth', 1.2);
    set(ax, 'XTick', 1:MAX_ORDER);
    ylabel(ax, 'Relative L1 (%)');
    title(ax, sprintf('param = %s', PARAMS{pi}), 'FontSize', 9, 'Interpreter','none');
    if pi == 1
        legend(ax, 'Location', 'northeast', 'FontSize', 6, 'NumColumns', 2);
    end
    if pi == 3, xlabel(ax, 'PCE Order K'); end
    grid(ax, 'on');
end
sgtitle('PCE Convergence — All Scenarios × Distributions', 'FontSize', 11);

pce_root = fullfile('Methods', 'PCE');
if ~exist(pce_root,'dir'), mkdir(pce_root); end
saveas(fig_sum, fullfile(pce_root, 'pce_summary_all.png'));
close(fig_sum);
fprintf('\nSaved summary: Methods/PCE/pce_summary_all.png\n');

%% ── Console K* table ─────────────────────────────────────────────────────────
fprintf('\n  K* (argmin LOO-CV | old L1<10%% in brackets):\n');
fprintf('  %-35s  freq   zS    svp\n', 'Scenario / Distribution');
fprintf('  %s\n', repmat('-', 1, 65));
for ci = 1:numel(completed)
    sc   = completed{ci}{1};
    dist = completed{ci}{2};
    pce_file = fullfile('Methods','PCE', sc.name, dist.name, 'results', 'pce_results.mat');
    if ~isfile(pce_file), continue; end
    r = load(pce_file, 'Kstar', 'Kstar_L1');
    fprintf('  %-35s  %-9s  %-9s  %-9s\n', ...
            sprintf('%s / %s', sc.name, dist.name), ...
            sprintf('%d[%s]', r.Kstar(1), fmt_kstar(r.Kstar_L1(1))), ...
            sprintf('%d[%s]', r.Kstar(2), fmt_kstar(r.Kstar_L1(2))), ...
            sprintf('%d[%s]', r.Kstar(3), fmt_kstar(r.Kstar_L1(3))));
end
fprintf('\n=== run_pce.m complete ===\n');

%% ── Helper ───────────────────────────────────────────────────────────────────
function s = fmt_kstar(k)
if isnan(k) || k > 10, s = '>10'; else, s = num2str(k); end
end
