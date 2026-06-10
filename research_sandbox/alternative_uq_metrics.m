%% alternative_uq_metrics.m  —  Task C: Beyond-variance UQ metrics
%
% Computes from existing MC TL cubes:
%   1. P95 envelope  (95th percentile TL at each range-depth pixel)
%   2. P(TL > threshold) detection failure probability map
%   3. 90% CI width map (P95 - P5)
%   4. Pixel-level TL histograms + normality check (Lilliefors test)
%
% Uses: Cache/baseline/Normal_10pct/TL_freq.mat  (most samples, interesting scenario)
%       Cache/deep_water/Normal_10pct/TL_freq.mat
%       Cache/deep_water/Normal_10pct/TL_zS.mat
%       Cache/deep_water/Normal_10pct/TL_svp.mat
%
% Outputs: figures in research_sandbox/figures/
%          results in research_sandbox/uq_metrics_results.mat

clear; close all; clc; warning('off','all');

PROJECT = 'C:\Users\orind\Stochastic-Bellhop\stochastic_TL_maps_project';
SANDBOX = 'C:\Users\orind\Stochastic-Bellhop\research_sandbox';
addpath(genpath(fullfile(PROJECT, 'Shared_Utils')));

fig_dir = fullfile(SANDBOX, 'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
set(0, 'DefaultFigureVisible', 'off');

cfg   = loadConfig();
N     = cfg.MC.N;     % 50
TH    = 100;          % detection threshold in dB (TL > TH => target not detectable)

COMBOS = {
    'baseline',   'Normal_10pct', 'freq';
    'deep_water', 'Normal_10pct', 'freq';
    'deep_water', 'Normal_10pct', 'zS';
    'deep_water', 'Normal_10pct', 'svp';
};

all_results = struct();

for ci = 1:size(COMBOS,1)
    sc_name   = COMBOS{ci,1};
    dist_name = COMBOS{ci,2};
    param     = COMBOS{ci,3};

    tl_file = fullfile(PROJECT, 'Cache', sc_name, dist_name, sprintf('TL_%s.mat', param));
    if ~isfile(tl_file)
        fprintf('[SKIP] %s\n', tl_file);
        continue;
    end

    fprintf('\n=== %s | %s | %s ===\n', sc_name, dist_name, param);

    tmp    = load(tl_file, 'TL_save');
    TL_all = double(tmp.TL_save(:,:,1:N));   % [Nz x Nr x N]
    [Nz, Nr, ~] = size(TL_all);

    % Load r_km, z_m from MC results
    mc_file = fullfile(PROJECT, 'Methods', 'MC', sc_name, dist_name, 'results', ...
                       sprintf('MC_%s.mat', param));
    mc_res  = load(mc_file, 'r_km', 'z_m', 'MC_Mean', 'MC_Var');
    r_km    = mc_res.r_km;
    z_m     = mc_res.z_m;
    MC_Mean = mc_res.MC_Mean;
    MC_Var  = mc_res.MC_Var;

    %% 1. P5, P50 (median), P95 maps
    P05 = prctile(TL_all, 5,  3);   % [Nz x Nr]
    P50 = prctile(TL_all, 50, 3);
    P95 = prctile(TL_all, 95, 3);
    CI90_width = P95 - P05;          % 90% CI width (dB)

    %% 2. Detection failure probability: P(TL > TH)
    P_fail = mean(TL_all > TH, 3);   % [Nz x Nr]

    %% 3. Normality check at sampled pixels
    % Sample 200 random pixels and run Lilliefors test
    rng(42);
    n_test = min(200, Nz*Nr);
    pix_idx = randperm(Nz*Nr, n_test);
    [z_flat, r_flat] = ind2sub([Nz Nr], pix_idx);  % pixel coordinates
    TL_flat = reshape(TL_all, Nz*Nr, N)';  % [N x Nz*Nr]

    lillie_h = zeros(n_test, 1);
    lillie_p = zeros(n_test, 1);
    skew_pix = zeros(n_test, 1);
    kurt_pix = zeros(n_test, 1);

    for k = 1:n_test
        samples = TL_flat(:, pix_idx(k));
        if var(samples) < 1e-12, continue; end
        try
            [lillie_h(k), lillie_p(k)] = lillietest(samples, 'Alpha', 0.05);
        catch
            lillie_h(k) = NaN; lillie_p(k) = NaN;
        end
        skew_pix(k) = skewness(samples);
        kurt_pix(k) = kurtosis(samples) - 3;  % excess kurtosis
    end

    pct_normal = 100 * mean(lillie_h(~isnan(lillie_h)) == 0);
    fprintf('  Normality (Lilliefors 5%%): %.1f%% of pixels accept H0 (normal)\n', pct_normal);
    fprintf('  Mean |skewness| = %.3f,  Mean excess kurtosis = %.3f\n', ...
            mean(abs(skew_pix)), mean(kurt_pix));
    fprintf('  CI90 width: mean=%.2f dB, max=%.2f dB\n', ...
            mean(CI90_width(:)), max(CI90_width(:)));
    fprintf('  P(TL>%ddB): mean=%.3f, max=%.3f\n', TH, mean(P_fail(:)), max(P_fail(:)));

    key = sprintf('%s_%s_%s', sc_name, dist_name, param);
    all_results.(key).P05       = P05;
    all_results.(key).P50       = P50;
    all_results.(key).P95       = P95;
    all_results.(key).CI90      = CI90_width;
    all_results.(key).P_fail    = P_fail;
    all_results.(key).pct_norm  = pct_normal;
    all_results.(key).skew_pix  = skew_pix;
    all_results.(key).kurt_pix  = kurt_pix;
    all_results.(key).r_km      = r_km;
    all_results.(key).z_m       = z_m;
    all_results.(key).MC_Mean   = MC_Mean;
    all_results.(key).MC_Var    = MC_Var;

    %% Plot: 4-panel map figure
    fig = figure('Position', [50 50 1300 700]);

    subplot(2,3,1);
    imagesc(r_km, z_m, MC_Mean); colorbar; axis xy;
    xlabel('Range (km)'); ylabel('Depth (m)');
    title(sprintf('E[TL] (dB)\n%s|%s|%s', sc_name, dist_name, param), ...
          'Interpreter','none','FontSize',8);
    colormap(gca, parula);

    subplot(2,3,2);
    imagesc(r_km, z_m, P95); colorbar; axis xy;
    xlabel('Range (km)'); ylabel('Depth (m)');
    title('P95 TL (dB)', 'FontSize', 8);
    colormap(gca, parula);

    subplot(2,3,3);
    imagesc(r_km, z_m, CI90_width); colorbar; axis xy;
    xlabel('Range (km)'); ylabel('Depth (m)');
    title('90% CI width (dB)', 'FontSize', 8);
    colormap(gca, hot);

    subplot(2,3,4);
    imagesc(r_km, z_m, sqrt(MC_Var)); colorbar; axis xy;
    xlabel('Range (km)'); ylabel('Depth (m)');
    title('\sigma_{TL} (dB)', 'FontSize', 8);
    colormap(gca, hot);

    subplot(2,3,5);
    imagesc(r_km, z_m, P_fail); colorbar; axis xy; caxis([0 1]);
    xlabel('Range (km)'); ylabel('Depth (m)');
    title(sprintf('P(TL > %d dB)', TH), 'FontSize', 8);
    colormap(gca, hot);

    % Panel 6: histogram of 5 representative pixels
    subplot(2,3,6);
    depths_idx = round(linspace(1, Nz, 5));
    range_mid  = round(Nr/3);
    colors6 = lines(5);
    for kk = 1:5
        samples = squeeze(TL_all(depths_idx(kk), range_mid, :));
        [h_cnt, h_edge] = histcounts(samples, 12, 'Normalization','pdf');
        h_centers = (h_edge(1:end-1) + h_edge(2:end))/2;
        plot(h_centers, h_cnt, '-', 'Color', colors6(kk,:), 'LineWidth', 1.5, ...
             'DisplayName', sprintf('z=%.0fm', z_m(depths_idx(kk))));
        hold on;
        % Overlay normal fit
        mu_k = mean(samples); sg_k = std(samples);
        if sg_k > 0.01
            xg = linspace(min(samples)-1, max(samples)+1, 100);
            plot(xg, normpdf(xg, mu_k, sg_k), '--', 'Color', colors6(kk,:), ...
                 'LineWidth', 0.8, 'HandleVisibility','off');
        end
    end
    xlabel('TL (dB)'); ylabel('PDF');
    title(sprintf('Pixel histograms + Normal fits\n(range=%.1fkm)', r_km(range_mid)), ...
          'FontSize', 8);
    legend('Location','best','FontSize',6);
    grid on;

    sgtitle(sprintf('Alternative UQ Metrics: %s | %s | param=%s', sc_name, dist_name, param), ...
            'Interpreter','none', 'FontSize', 10);

    saveas(fig, fullfile(fig_dir, sprintf('uq_metrics_%s_%s_%s.png', sc_name, dist_name, param)));
    close(fig);
    fprintf('  Saved: uq_metrics_%s_%s_%s.png\n', sc_name, dist_name, param);
end

%% Save
save(fullfile(SANDBOX, 'uq_metrics_results.mat'), 'all_results', '-v7.3');
fprintf('\nSaved: research_sandbox/uq_metrics_results.mat\n');
fprintf('=== alternative_uq_metrics.m complete ===\n');
