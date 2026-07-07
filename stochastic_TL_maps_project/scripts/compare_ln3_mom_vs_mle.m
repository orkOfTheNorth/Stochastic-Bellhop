%% compare_ln3_mom_vs_mle.m
% Generates LN3 MOM vs MLE comparison figures for each scenario.
% Shows:
%   - Panel grid: 4 spatial locations, each showing MC histogram + MLE fit + MOM fit
%   - Summary scatter: P(detect) from MOM vs MLE vs MC
%   - KS statistic comparison: MOM vs MLE
%
% MLE = lognfit (gamma = 0.95*min, then MLE for mu/sigma) [current pipeline]
% MOM = method of moments (gamma = 0.95*min, then MOM for mu/sigma)
%
% Saves to Methods/LN3/<scen>/figures/LN3_MOM_vs_MLE_<param>.png

ROOT = fileparts(fileparts(mfilename('fullpath')));
try; cd(ROOT); catch; end
addpath(genpath(fullfile(ROOT,'core')));
addpath(fullfile(ROOT,'binaries'));

cfg = loadConfig();
PARAMS = {'freq','zS','svp'};
FOM_DEFAULT = 70;

for si = 1:numel(cfg.scenarios)
    sc = cfg.scenarios(si);
    sc_name = sc.name;

    for di = 1:numel(cfg.distributions)
        dist = cfg.distributions(di);
        dist_name = dist.name;
        FOM = getFOM(cfg, sc_name);

        fig_dir = fullfile('Methods','LN3', sc_name, dist_name, 'figures');
        if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

        for ki = 1:numel(PARAMS)
            sn = PARAMS{ki};
            out_fig = fullfile(fig_dir, sprintf('LN3_MOM_vs_MLE_%s.png', sn));
            if isfile(out_fig)
                fprintf('  SKIP %s/%s/%s (exists)\n', sc_name, dist_name, sn);
                continue;
            end

            tl_cache = fullfile('Cache', sc_name, dist_name, sprintf('TL_%s.mat', sn));
            mc_file  = fullfile('Methods','MC', sc_name, dist_name, 'results', ...
                                sprintf('MC_%s.mat', sn));
            if ~isfile(tl_cache) || ~isfile(mc_file)
                fprintf('  SKIP %s/%s/%s — missing data\n', sc_name, dist_name, sn);
                continue;
            end

            fprintf('[%s/%s/%s] Computing MOM vs MLE...\n', sc_name, dist_name, sn);

            D_tl = load(tl_cache, 'TL_save');
            D_mc = load(mc_file, 'r_km','z_m','MC_P_detect');
            TL_all = double(D_tl.TL_save);
            [Nz, Nr, N] = size(TL_all);
            r_km = D_mc.r_km;
            z_m  = D_mc.z_m;
            MC_P = D_mc.MC_P_detect;
            clear D_tl D_mc;

            % Pixel grid: 4 sample locations
            ri_vec = round(linspace(Nr*0.15, Nr*0.85, 4));
            zi_vec = round(linspace(Nz*0.25, Nz*0.75, 4));
            ri_vec = max(1, min(Nr, ri_vec));
            zi_vec = max(1, min(Nz, zi_vec));

            % Figure: histogram grid (4×4 = 4 range × 4 depth)
            fig = figure('Position',[50 50 1400 900]);
            t   = tiledlayout(4, 4, 'TileSpacing','compact','Padding','compact');
            title(t, sprintf('LN3 MLE vs MOM | %s | %s | %s | FOM=%ddB', ...
                             sc_name, dist_name, sn, FOM), ...
                  'FontSize',11,'Interpreter','none');

            ks_mle_all = NaN(Nz,Nr);
            ks_mom_all = NaN(Nz,Nr);
            p_mle_all  = NaN(Nz,Nr);
            p_mom_all  = NaN(Nz,Nr);

            for zi_i = 1:4
                for ri_i = 1:4
                    zi = zi_vec(zi_i);
                    ri = ri_vec(ri_i);
                    samps = squeeze(TL_all(zi, ri, :));

                    % MLE fit (current pipeline)
                    [p_mle, gam_mle, mu_mle, sig_mle] = ln3Fit(samps, FOM);
                    ks_mle = ksStatLN3(samps, gam_mle, mu_mle, sig_mle);

                    % MOM fit (new)
                    [p_mom, gam_mom, mu_mom, sig_mom] = ln3FitMOM(samps, FOM);
                    ks_mom = ksStatLN3(samps, gam_mom, mu_mom, sig_mom);

                    ax = nexttile;
                    histogram(ax, samps, 25, 'Normalization','pdf', ...
                              'FaceColor',[0.75 0.85 0.95],'EdgeColor','none', ...
                              'DisplayName','MC');
                    hold(ax,'on');

                    x_pdf = linspace(min(samps)*0.98, max(samps)*1.02, 300);

                    % MLE curve
                    sh_mle = x_pdf - gam_mle;
                    sh_mle(sh_mle<=0) = NaN;
                    y_mle  = lognpdf(sh_mle, mu_mle, sig_mle);
                    plot(ax, x_pdf, y_mle, 'b-', 'LineWidth',1.8, 'DisplayName', ...
                         sprintf('MLE (KS=%.3f)', ks_mle));

                    % MOM curve
                    sh_mom = x_pdf - gam_mom;
                    sh_mom(sh_mom<=0) = NaN;
                    y_mom  = lognpdf(sh_mom, mu_mom, sig_mom);
                    plot(ax, x_pdf, y_mom, 'r--', 'LineWidth',1.8, 'DisplayName', ...
                         sprintf('MOM (KS=%.3f)', ks_mom));

                    xline(ax, FOM, 'k:', 'LineWidth',1.2, 'DisplayName','FOM');
                    hold(ax,'off');

                    xlabel(ax,'TL (dB)','FontSize',7);
                    ylabel(ax,'pdf','FontSize',7);
                    title(ax, sprintf('r=%.1fkm z=%.0fm\nP_{MC}=%.2f P_{MLE}=%.2f P_{MOM}=%.2f', ...
                                      r_km(ri), z_m(zi), MC_P(zi,ri), p_mle, p_mom), ...
                          'FontSize',6,'Interpreter','none');
                    if zi_i==1 && ri_i==1, legend(ax,'Location','northwest','FontSize',6); end
                end
            end

            saveFigPNG(fig, fullfile(fig_dir, sprintf('LN3_MOM_vs_MLE_%s', sn)));
            close(fig);
            fprintf('  Saved: LN3_MOM_vs_MLE_%s.png\n', sn);

            % Also compute pixel-wise scatter: MLE vs MOM P(detect)
            fprintf('  Computing pixel-wise MLE/MOM maps...\n');
            p_mle_map = zeros(Nz,Nr);
            p_mom_map = zeros(Nz,Nr);
            ks_mle_map = zeros(Nz,Nr);
            ks_mom_map = zeros(Nz,Nr);
            nll_mle_map = zeros(Nz,Nr);
            nll_mom_map = zeros(Nz,Nr);

            for ri = 1:Nr
                for zi = 1:Nz
                    samps = squeeze(TL_all(zi,ri,:));
                    [pm, gm, mu_m, sg_m] = ln3Fit(samps, FOM);
                    p_mle_map(zi,ri)  = pm;
                    ks_mle_map(zi,ri) = ksStatLN3(samps, gm, mu_m, sg_m);
                    nll_mle_map(zi,ri) = ln3NLL(samps, gm, mu_m, sg_m);
                    [po, go, mu_o, sg_o] = ln3FitMOM(samps, FOM);
                    p_mom_map(zi,ri)  = po;
                    ks_mom_map(zi,ri) = ksStatLN3(samps, go, mu_o, sg_o);
                    nll_mom_map(zi,ri) = ln3NLL(samps, go, mu_o, sg_o);
                end
            end

            % Scatter plot: MLE vs MOM P(detect)
            fig2 = figure('Position',[50 50 1600 450]);
            subplot(1,4,1);
            scatter(MC_P(:), p_mle_map(:), 2, [0.2 0.4 0.8],'filled','MarkerFaceAlpha',0.2);
            hold on; plot([0 1],[0 1],'r-','LineWidth',1.2); hold off;
            xlabel('MC P(detect)'); ylabel('MLE P(detect)');
            title(sprintf('MLE vs MC\nRMSE=%.4f', rms(p_mle_map(:)-MC_P(:))));
            grid on;

            subplot(1,4,2);
            scatter(MC_P(:), p_mom_map(:), 2, [0.8 0.2 0.2],'filled','MarkerFaceAlpha',0.2);
            hold on; plot([0 1],[0 1],'r-','LineWidth',1.2); hold off;
            xlabel('MC P(detect)'); ylabel('MOM P(detect)');
            title(sprintf('MOM vs MC\nRMSE=%.4f', rms(p_mom_map(:)-MC_P(:))));
            grid on;

            subplot(1,4,3);
            scatter(ks_mle_map(:), ks_mom_map(:), 2, [0.2 0.7 0.2],'filled','MarkerFaceAlpha',0.2);
            hold on;
            mx = max([ks_mle_map(:); ks_mom_map(:)]);
            plot([0 mx],[0 mx],'k--','LineWidth',1);
            xline(0.10,'r:','LineWidth',1); yline(0.10,'r:','LineWidth',1);
            hold off;
            xlabel('KS (MLE)'); ylabel('KS (MOM)');
            title(sprintf('KS: MLE vs MOM\nMLE<0.10: %.1f%% | MOM<0.10: %.1f%%', ...
                          100*mean(ks_mle_map(:)<0.10), 100*mean(ks_mom_map(:)<0.10)));
            grid on;

            % NLL comparison — the LN3-native metric: MLE is by
            % construction the fit minimizing in-sample NLL, so this is
            % the theoretically "fair" comparison (replaces the previous
            % unverified "MLE 10-30% lower RMSE" claim with real numbers).
            subplot(1,4,4);
            scatter(nll_mle_map(:), nll_mom_map(:), 2, [0.6 0.3 0.7],'filled','MarkerFaceAlpha',0.2);
            hold on;
            mx = max([nll_mle_map(:); nll_mom_map(:)]);
            mn = min([nll_mle_map(:); nll_mom_map(:)]);
            plot([mn mx],[mn mx],'k--','LineWidth',1);
            hold off;
            xlabel('mean NLL (MLE)'); ylabel('mean NLL (MOM)');
            pct_mle_better = 100*mean(nll_mle_map(:) < nll_mom_map(:));
            title(sprintf('NLL: MLE vs MOM\nMLE lower-NLL in %.1f%% of pixels', pct_mle_better));
            grid on;

            sgtitle(sprintf('MLE vs MOM | %s | %s | %s', sc_name, dist_name, sn), ...
                    'Interpreter','none','FontSize',11);

            saveFigPNG(fig2, fullfile(fig_dir, sprintf('LN3_MOM_vs_MLE_scatter_%s', sn)));
            close(fig2);
            fprintf('  Saved: LN3_MOM_vs_MLE_scatter_%s.png\n', sn);

            clear TL_all;
        end
    end
end
fprintf('\nDone. MOM vs MLE comparison complete.\n');
