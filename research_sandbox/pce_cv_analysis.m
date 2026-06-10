%% pce_cv_analysis.m  —  K-fold cross-validation for PCE order selection
%
% Investigates PCE overfitting: with N=50 samples and order up to K=10,
% the design matrix has up to 11 columns. This script compares training L1
% vs. K-fold CV L1 to detect where overfitting begins.
%
% Loads: Cache/deep_water/Normal_10pct/TL_zS.mat  (and TL_freq.mat)
%        Methods/MC/deep_water/Normal_10pct/results/MC_zS.mat
%
% Outputs: figures saved in research_sandbox/figures/
%          results saved in research_sandbox/cv_results.mat

clear; close all; clc; warning('off','all');

%% Paths -----------------------------------------------------------------
PROJECT = 'C:\Users\orind\Stochastic-Bellhop\stochastic_TL_maps_project';
SANDBOX = 'C:\Users\orind\Stochastic-Bellhop\research_sandbox';
addpath(genpath(fullfile(PROJECT, 'Shared_Utils')));

fig_dir = fullfile(SANDBOX, 'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

set(0, 'DefaultFigureVisible', 'off');

%% Config ----------------------------------------------------------------
cfg    = loadConfig();
freq0  = cfg.nominal.freq_Hz;
zS0    = cfg.nominal.zS_m;
N      = cfg.MC.N;          % 50
K_FOLD = 5;
MAX_K  = 10;
PARAMS = {'freq','zS','svp'};

dist_name = 'Normal_10pct';
sc_name   = 'deep_water';

dist      = cfg.distributions(strcmp({cfg.distributions.name}, dist_name));
B         = computeVarianceBounds(cfg, dist);
S         = lhsSample(N, cfg, dist, B);

dist_type = lower(dist.type);   % 'normal'
xi_all = [(S.freq - freq0) / B.sig_freq, ...
          (S.zS   - zS0)   / B.sig_zS, ...
           S.svp           / B.sig_svp];   % [50 x 3]

fprintf('=== PCE K-fold CV Analysis ===\n');
fprintf('Scenario: %s  |  Distribution: %s\n', sc_name, dist_name);
fprintf('N=%d samples,  K_FOLD=%d,  MAX_K=%d\n\n', N, K_FOLD, MAX_K);

%% Storage ---------------------------------------------------------------
results = struct();

for pi = 1:3
    param = PARAMS{pi};

    %% Load TL cube
    tl_file = fullfile(PROJECT, 'Cache', sc_name, dist_name, sprintf('TL_%s.mat', param));
    if ~isfile(tl_file)
        fprintf('[SKIP] missing %s\n', tl_file);
        continue;
    end
    tmp    = load(tl_file, 'TL_save');
    TL_all = double(tmp.TL_save(:,:,1:N));
    [Nz, Nr, ~] = size(TL_all);
    TL_mat = reshape(permute(TL_all, [3 1 2]), N, Nz*Nr);  % [N x pixels]

    %% Load MC reference variance
    mc_file = fullfile(PROJECT, 'Methods', 'MC', sc_name, dist_name, 'results', ...
                       sprintf('MC_%s.mat', param));
    mc_res  = load(mc_file, 'MC_Var');
    Var_MC  = mc_res.MC_Var(:)';   % [1 x pixels]
    mean_VM = max(mean(Var_MC), 1e-10);

    xi = xi_all(:, pi);

    %% Full training L1 (same as run_pce.m) --------------------------------
    Phi = pce_basis(xi, MAX_K, dist_type);  % [N x MAX_K+1]
    [Q, ~] = qr(Phi, 0);
    QTL = Q' * TL_mat;

    L1_train = zeros(MAX_K, 1);
    for K = 1:MAX_K
        yhat_K = Q(:,1:K+1) * QTL(1:K+1,:);
        Var_K  = var(yhat_K, 0, 1);
        L1_train(K) = mean(abs(Var_K - Var_MC)) / mean_VM;
    end

    %% K-fold cross-validation ---------------------------------------------
    % Partition N=50 samples into K_FOLD folds
    idx      = randperm(N, N);   % shuffle with fixed seed implicitly
    fold_ids = reshape(idx(1:floor(N/K_FOLD)*K_FOLD), K_FOLD, []);
    % Each row of fold_ids is one fold's sample indices

    Var_oof_sum = zeros(K_FOLD, Nz*Nr);   % out-of-fold variance accumulator
    % We store per-fold, per-K predictions
    oof_pred_all = zeros(K_FOLD, MAX_K, Nz*Nr);  % [folds x orders x pixels]

    for f = 1:K_FOLD
        test_idx  = fold_ids(f,:);
        train_idx = reshape(fold_ids(setdiff(1:K_FOLD, f),:), 1, []);

        xi_train  = xi(train_idx);
        xi_test   = xi(test_idx);
        TL_train  = TL_mat(train_idx, :);
        TL_test   = TL_mat(test_idx,  :);

        Phi_train = pce_basis(xi_train, MAX_K, dist_type);
        Phi_test  = pce_basis(xi_test,  MAX_K, dist_type);

        [Q_tr, R_tr] = qr(Phi_train, 0);

        for K = 1:MAX_K
            % Fit on training set using order K
            c_K = R_tr(1:K+1, 1:K+1) \ (Q_tr(:,1:K+1)' * TL_train);
            % Predict on test set
            yhat_test = Phi_test(:, 1:K+1) * c_K;   % [n_test x pixels]
            oof_pred_all(f, K, :) = mean(yhat_test, 1);
        end
        Var_oof_sum(f,:) = var(TL_test, 0, 1);
    end

    %% CV L1: use out-of-fold predicted mean to reconstruct variance --------
    % Strategy: aggregate all OOF predictions; compute variance across folds
    % For each K: CV variance = variance of {fold mean predictions over test sets}
    %   A simpler proxy: CV-L1(K) = mean over folds of
    %   mean|Var(yhat_test,K) - Var_MC| / mean(Var_MC)
    L1_cv = zeros(MAX_K, 1);
    for K = 1:MAX_K
        fold_L1s = zeros(K_FOLD, 1);
        for f = 1:K_FOLD
            test_idx = fold_ids(f,:);
            xi_train = xi(setdiff(1:N, fold_ids(f,:)));
            xi_test  = xi(test_idx);
            TL_train = TL_mat(setdiff(1:N, fold_ids(f,:)), :);
            TL_test  = TL_mat(test_idx, :);

            Phi_tr   = pce_basis(xi_train, MAX_K, dist_type);
            Phi_te   = pce_basis(xi_test,  MAX_K, dist_type);
            c_K = (Phi_tr(:,1:K+1)' * Phi_tr(:,1:K+1)) \ (Phi_tr(:,1:K+1)' * TL_train);
            yhat_te  = Phi_te(:,1:K+1) * c_K;
            Var_te   = var(yhat_te, 0, 1);
            fold_L1s(f) = mean(abs(Var_te - Var_MC)) / mean_VM;
        end
        L1_cv(K) = mean(fold_L1s);
    end

    %% Find K* from train and from CV ---------------------------------------
    kstar_train = find(L1_train < 0.10, 1);
    kstar_cv    = find(L1_cv    < 0.10, 1);
    if isempty(kstar_train), kstar_train = NaN; end
    if isempty(kstar_cv),    kstar_cv    = NaN; end

    %% Print summary --------------------------------------------------------
    fprintf('Param: %s\n', param);
    fprintf('  K:        '); fprintf('%4d ', 1:MAX_K); fprintf('\n');
    fprintf('  Train L1: '); fprintf('%4.2f ', L1_train*100); fprintf('%%\n');
    fprintf('  CV L1:    '); fprintf('%4.2f ', L1_cv*100);    fprintf('%%\n');
    fprintf('  K*(train)=%s  K*(CV)=%s\n\n', ...
            fmt_kstar(kstar_train), fmt_kstar(kstar_cv));

    %% Detect overfitting ---------------------------------------------------
    overfit_K = find(L1_cv > L1_train * 1.5 & L1_cv > 0.05, 1);
    if ~isempty(overfit_K)
        fprintf('  --> Overfitting detected at K=%d (CV/Train ratio=%.2fx)\n\n', ...
                overfit_K, L1_cv(overfit_K)/L1_train(overfit_K));
    else
        fprintf('  --> No clear overfitting (CV and train track closely)\n\n');
    end

    results.(param).L1_train    = L1_train;
    results.(param).L1_cv       = L1_cv;
    results.(param).kstar_train = kstar_train;
    results.(param).kstar_cv    = kstar_cv;
    results.(param).Var_MC      = Var_MC;
    results.(param).mean_VM     = mean_VM;

    %% Plot -----------------------------------------------------------------
    fig = figure('Position', [50 50 700 450]);
    colors = [0.2 0.4 0.8; 0.8 0.2 0.2];
    plot(1:MAX_K, L1_train*100, '-o', 'Color', colors(1,:), 'LineWidth', 2, ...
         'MarkerSize', 6, 'DisplayName', 'Training L1');
    hold on;
    plot(1:MAX_K, L1_cv*100,    '-s', 'Color', colors(2,:), 'LineWidth', 2, ...
         'MarkerSize', 6, 'DisplayName', sprintf('%d-fold CV L1', K_FOLD));
    yline(10, 'k--', '10\% threshold', 'LineWidth', 1.5, ...
          'LabelHorizontalAlignment', 'left');

    % Mark K*
    if ~isnan(kstar_train)
        plot(kstar_train, L1_train(kstar_train)*100, 'v', ...
             'Color', colors(1,:), 'MarkerSize', 10, 'MarkerFaceColor', colors(1,:), ...
             'HandleVisibility','off');
        text(kstar_train+0.15, L1_train(kstar_train)*100+1.5, ...
             sprintf('K^*_{train}=%d', kstar_train), 'Color', colors(1,:), 'FontSize', 9);
    end
    if ~isnan(kstar_cv)
        plot(kstar_cv, L1_cv(kstar_cv)*100, '^', ...
             'Color', colors(2,:), 'MarkerSize', 10, 'MarkerFaceColor', colors(2,:), ...
             'HandleVisibility','off');
        text(kstar_cv+0.15, L1_cv(kstar_cv)*100+1.5, ...
             sprintf('K^*_{CV}=%d', kstar_cv), 'Color', colors(2,:), 'FontSize', 9);
    end

    set(gca, 'XTick', 1:MAX_K);
    xlabel('PCE Order K');
    ylabel('Relative L_1 error (%)');
    title(sprintf('PCE Overfitting Diagnostic: %s | %s | param=%s', ...
                  sc_name, dist_name, param), 'Interpreter', 'none');
    legend('Location', 'northeast', 'FontSize', 9);
    grid on;
    ylim([0, max([L1_train; L1_cv])*100 * 1.2 + 5]);

    saveas(fig, fullfile(fig_dir, sprintf('cv_L1_%s_%s_%s.png', sc_name, dist_name, param)));
    close(fig);
    fprintf('  Saved figure: cv_L1_%s_%s_%s.png\n', sc_name, dist_name, param);
end

%% Save all results
save(fullfile(SANDBOX, 'cv_results.mat'), 'results', 'sc_name', 'dist_name', ...
     'K_FOLD', 'MAX_K', '-v7.3');
fprintf('\nSaved: research_sandbox/cv_results.mat\n');
fprintf('=== pce_cv_analysis.m complete ===\n');

function s = fmt_kstar(k)
if isnan(k) || (isnumeric(k) && k > 10)
    s = '>10';
else
    s = num2str(k);
end
end
