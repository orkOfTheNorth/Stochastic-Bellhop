function regen_pce_relevance()
%% regen_pce_relevance  Re-draw pce_order_relevance.png from saved results.
% Loads pce_results.mat for every scenario/dist and regenerates the
% order-relevance figure showing both K_Var (peak of curve) and K* (LOO).

close all;
ROOT = fileparts(fileparts(mfilename('fullpath')));
try; cd(ROOT); catch; end
addpath(genpath(fullfile(ROOT,'core')));
set(0,'DefaultFigureVisible','off');

cfg    = loadConfig();
PARAMS = {'freq','zS','svp'};

for si = 1:numel(cfg.scenarios)
    sc = cfg.scenarios(si);
    for di = 1:numel(cfg.distributions)
        dist = cfg.distributions(di);
        res_file = fullfile('Methods','PCE', sc.name, dist.name, ...
                            'results','pce_results.mat');
        if ~isfile(res_file)
            fprintf('[MISS] %s\n', res_file); continue;
        end

        d = load(res_file, 'L2_mat','LOO_comb_mat','Kstar','MAX_ORDER', ...
                           'sc_name','dist_name');
        L2_mat       = d.L2_mat;
        LOO_comb_mat = d.LOO_comb_mat;
        Kstar        = d.Kstar;
        MAX_ORDER    = d.MAX_ORDER;
        sc_name      = d.sc_name;
        dist_name    = d.dist_name;

        fig_dir = fullfile('Methods','PCE', sc.name, dist.name, 'figures');
        if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
        out_png = fullfile(fig_dir, 'pce_order_relevance');

        fig_rel = figure('Position', [50 50 1000 380]);
        for pi = 1:3
            ax_rel = subplot(1, 3, pi);
            l2_v   = L2_mat(:, pi);
            if all(isnan(l2_v)) || all(l2_v == 0)
                axis(ax_rel,'off'); continue;
            end

            [~, K_var_star] = min(l2_v);   % K at max relevance
            l2_start  = l2_v(1);
            l2_end    = min(l2_v);
            total_imp = max(l2_start - l2_end, eps);
            relevance = max(l2_start - l2_v, 0) / total_imp * 100;

            plot(ax_rel, 1:MAX_ORDER, relevance, 'b-o', ...
                 'LineWidth',1.5,'MarkerSize',5);
            hold(ax_rel,'on');
            yline(ax_rel, 90, 'r--', '90%', 'LineWidth',1);
            xline(ax_rel, K_var_star, 'g--', ...
                  sprintf('K_{Var}=%d', K_var_star), 'LineWidth',1.2);
            xline(ax_rel, Kstar(pi), 'm--', ...
                  sprintf('K*=%d (LOO)', Kstar(pi)), 'LineWidth',1.2);
            hold(ax_rel,'off');
            set(ax_rel, 'XTick', max(1, 1:2:MAX_ORDER));
            xlabel(ax_rel, 'PCE Order K');
            ylabel(ax_rel, '% Var improvement (cumul.)');
            title(ax_rel, sprintf('%s\nK_{Var}=%d  K*=%d', ...
                  PARAMS{pi}, K_var_star, Kstar(pi)), ...
                  'FontSize',9,'Interpreter','none');
            grid(ax_rel,'on');
        end
        sgtitle(sprintf('PCE Order Relevance | %s | %s', sc_name, dist_name), ...
                'Interpreter','none','FontSize',10);
        saveFigPNG(fig_rel, out_png);
        close(fig_rel);
        fprintf('Saved: %s\n', out_png);
    end
end
fprintf('regen_pce_relevance done.\n');
end
