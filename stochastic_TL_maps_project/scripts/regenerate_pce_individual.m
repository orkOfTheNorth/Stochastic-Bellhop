%% regenerate_pce_individual.m
% Regenerates all individual PCE convergence figures (pce_convergence_*.png)
% from saved pce_results.mat files, with y-axis clipping on LOO and L2 panels
% to prevent K>17 numerical blow-up from obscuring the scientifically meaningful
% convergence region (K=1-17).
%
% This script does NOT re-run any Bellhop or PCE computations; it only
% regenerates the figures from already-saved results.

ROOT = fileparts(fileparts(mfilename('fullpath')));
try; cd(ROOT); catch; end

YLIM_L2  = 30;    % clip L2 Var error panel at 30%
YLIM_LOO =  3;    % clip normalized LOO panels at 3x initial value

PARAMS = {'freq', 'zS', 'svp'};
colors = lines(3);

results_files = dir(fullfile('Methods','PCE','*','*','results','pce_results.mat'));
fprintf('Found %d pce_results.mat files.\n', numel(results_files));

for ci = 1:numel(results_files)
    mat_path = fullfile(results_files(ci).folder, results_files(ci).name);
    r = load(mat_path);

    sc_name   = r.sc_name;
    dist_name = r.dist_name;
    MAX_ORDER = r.MAX_ORDER;

    fig_dir = fullfile('Methods','PCE', sc_name, dist_name, 'figures');
    if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

    mc_res_dir = fullfile('Methods','MC', sc_name, dist_name, 'results');
    fprintf('[%d/%d] %s / %s\n', ci, numel(results_files), sc_name, dist_name);

    fig = figure('Position', [50 50 1600 640], 'Visible', 'off');

    %% ── (1,1): Combined LOO-CV vs order ────────────────────────────────────
    ax_loo = subplot(2,4,1);
    for pi = 1:3
        v = r.LOO_comb_mat(:,pi);
        if all(isnan(v)), continue; end
        plot(ax_loo, 1:MAX_ORDER, v, '-o', ...
             'Color', colors(pi,:), 'LineWidth', 1.6, 'MarkerSize', 4, ...
             'DisplayName', PARAMS{pi});
        hold(ax_loo, 'on');
        xline(ax_loo, r.Kstar(pi), '--', 'Color', colors(pi,:)*0.7, 'LineWidth', 1, ...
              'HandleVisibility','off');
    end
    ylim(ax_loo, [0, YLIM_LOO]);
    set(ax_loo, 'XTick', max(1, 1:2:MAX_ORDER));
    xlabel(ax_loo, 'PCE Order K');
    ylabel(ax_loo, 'LOO_{combined} (norm.)');
    title(ax_loo, sprintf('LOO Combined  [K*=vert]\n%s | %s', sc_name, dist_name), ...
          'FontSize', 8, 'Interpreter', 'none');
    legend(ax_loo, 'Location','northeast','FontSize',7);
    grid(ax_loo, 'on');

    %% ── (1,2): L2 Var error vs order ───────────────────────────────────────
    ax_l2 = subplot(2,4,2);
    for pi = 1:3
        v = r.L2_mat(:,pi);
        if all(isnan(v)), continue; end
        plot(ax_l2, 1:MAX_ORDER, v*100, '-o', ...
             'Color', colors(pi,:), 'LineWidth', 1.5, 'MarkerSize', 4, ...
             'DisplayName', PARAMS{pi});
        hold(ax_l2, 'on');
    end
    yline(ax_l2, 10, 'k--', '10%', 'LineWidth', 1.2, 'LabelHorizontalAlignment','left');
    ylim(ax_l2, [0, YLIM_L2]);
    set(ax_l2, 'XTick', max(1, 1:2:MAX_ORDER));
    xlabel(ax_l2, 'PCE Order K');
    ylabel(ax_l2, sprintf('Relative L2 Var error %%  [clipped at %d%%]', YLIM_L2));
    title(ax_l2, sprintf('L2 Variance [y-clipped]\n%s | %s', sc_name, dist_name), ...
          'FontSize', 8, 'Interpreter', 'none');
    legend(ax_l2, 'Location','northeast','FontSize',7);
    grid(ax_l2, 'on');

    %% ── (1,3): LOO_EX vs order ──────────────────────────────────────────────
    ax_lex = subplot(2,4,3);
    for pi = 1:3
        v = r.LOO_EX_mat(:,pi);
        if all(isnan(v)), continue; end
        loo_ex_n = v / max(v(1), eps);
        plot(ax_lex, 1:MAX_ORDER, loo_ex_n, '-o', ...
             'Color', colors(pi,:), 'LineWidth', 1.5, 'MarkerSize', 4, ...
             'DisplayName', PARAMS{pi});
        hold(ax_lex, 'on');
    end
    ylim(ax_lex, [0, YLIM_LOO]);
    set(ax_lex, 'XTick', max(1, 1:2:MAX_ORDER));
    xlabel(ax_lex, 'PCE Order K');
    ylabel(ax_lex, 'LOO_{EX} (norm.)');
    title(ax_lex, sprintf('LOO EX (prediction)\n%s | %s', sc_name, dist_name), ...
          'FontSize', 8, 'Interpreter', 'none');
    legend(ax_lex, 'Location','northeast','FontSize',7);
    grid(ax_lex, 'on');

    %% ── (1,4): LOO_Var vs order ─────────────────────────────────────────────
    ax_lv = subplot(2,4,4);
    for pi = 1:3
        v = r.LOO_Var_mat(:,pi);
        if all(isnan(v)), continue; end
        loo_var_n = v / max(v(1), eps);
        plot(ax_lv, 1:MAX_ORDER, loo_var_n, '-o', ...
             'Color', colors(pi,:), 'LineWidth', 1.5, 'MarkerSize', 4, ...
             'DisplayName', PARAMS{pi});
        hold(ax_lv, 'on');
    end
    ylim(ax_lv, [0, YLIM_LOO]);
    set(ax_lv, 'XTick', max(1, 1:2:MAX_ORDER));
    xlabel(ax_lv, 'PCE Order K');
    ylabel(ax_lv, 'LOO_{Var} (norm.)');
    title(ax_lv, sprintf('LOO Var (variance)\n%s | %s', sc_name, dist_name), ...
          'FontSize', 8, 'Interpreter', 'none');
    legend(ax_lv, 'Location','northeast','FontSize',7);
    grid(ax_lv, 'on');

    %% ── (2,1-3): scatter Var_{K*} vs Var_MC per param ──────────────────────
    for pi = 1:3
        ax_r = subplot(2, 4, pi+4);
        if isempty(r.VarK_best{pi})
            title(ax_r, sprintf('%s (no data)', PARAMS{pi}), 'FontSize',8,'Interpreter','none');
            continue;
        end
        mc_file = fullfile(mc_res_dir, sprintf('MC_%s.mat', PARAMS{pi}));
        if ~isfile(mc_file)
            title(ax_r, sprintf('%s (MC missing)', PARAMS{pi}), 'FontSize',8,'Interpreter','none');
            continue;
        end
        mc_pi = load(mc_file, 'MC_Var');
        vmc  = mc_pi.MC_Var(:);
        vpce = r.VarK_best{pi}(:);

        ok      = vmc > prctile(vmc, 5);
        xlim_hi = prctile(vmc(ok),  99);
        ylim_hi = prctile(vpce(ok), 99);
        lim     = max(xlim_hi, ylim_hi);

        scatter(ax_r, vmc(ok), vpce(ok), 3, 'filled', ...
                'MarkerFaceAlpha', 0.25, 'MarkerFaceColor', colors(pi,:));
        hold(ax_r,'on');
        plot(ax_r, [0 lim], [0 lim], 'k-', 'LineWidth', 1.5);
        xlim(ax_r, [0 xlim_hi]); ylim(ax_r, [0 ylim_hi]);

        ok2 = ok & isfinite(vmc) & isfinite(vpce);
        R2  = 0;
        if sum(ok2) > 2, R2 = corr(vmc(ok2), vpce(ok2))^2; end

        xlabel(ax_r, 'Var_{MC}');
        ylabel(ax_r, sprintf('Var_{K=%d}', r.Kstar(pi)));
        title(ax_r, sprintf('Var Scatter | %s\nR^2=%.3f K*=%d', PARAMS{pi}, R2, r.Kstar(pi)), ...
              'FontSize', 8, 'Interpreter','none');
        grid(ax_r,'on'); axis(ax_r,'equal');
    end

    %% ── (2,4): scatter P_detect PCE vs MC ──────────────────────────────────
    ax_pd = subplot(2, 4, 8);
    hold(ax_pd,'on');
    for pi = 1:3
        if isempty(r.Pdetect_best{pi}), continue; end
        mc_pd_file = fullfile(mc_res_dir, sprintf('MC_%s.mat', PARAMS{pi}));
        if ~isfile(mc_pd_file), continue; end
        mc_pd = load(mc_pd_file, 'MC_P_detect');
        scatter(ax_pd, mc_pd.MC_P_detect(:), r.Pdetect_best{pi}(:), 3, 'filled', ...
                'MarkerFaceAlpha', 0.2, 'MarkerFaceColor', colors(pi,:), ...
                'DisplayName', PARAMS{pi});
    end
    plot(ax_pd, [0 1], [0 1], 'k-', 'LineWidth', 1.5, 'HandleVisibility','off');
    xlim(ax_pd, [0 1]); ylim(ax_pd, [0 1]);
    xlabel(ax_pd, 'P_{detect} MC');
    ylabel(ax_pd, 'P_{detect} PCE');
    title(ax_pd, sprintf('P(detect) scatter\nall params @ K*'), ...
          'FontSize', 8, 'Interpreter','none');
    legend(ax_pd, 'Location','northwest','FontSize',7);
    grid(ax_pd,'on'); axis(ax_pd,'equal');

    sgtitle(sprintf('PCE Variance Convergence | %s | %s  [LOO/L2 y-axes clipped]', ...
                    sc_name, dist_name), 'FontSize', 10, 'Interpreter','none');

    out_name = fullfile(fig_dir, sprintf('pce_convergence_%s_%s', sc_name, dist_name));
    saveas(fig, [out_name '.png']);
    close(fig);
    fprintf('  Saved: %s.png\n', out_name);
end

fprintf('\nDone. Regenerated %d individual PCE convergence figures.\n', numel(results_files));
