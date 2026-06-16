%% plot_delta_order_comparison.m
%
% Visualises the importance of the 2nd-order Delta correction for
% VARIANCE across every subset combination × distribution — ALL scenarios.
%
% Two output figures per scenario:
%
%  Fig 1  —  ABSOLUTE correction   ΔVar = Var_2nd − Var_1st   [dB²]
%             Layout: n_dist rows × 7 columns
%             Colorscale shared within each row (distribution) so
%             the relative magnitude across subsets is clear.
%
%  Fig 2  —  RELATIVE correction   100 × ΔVar / Var_1st       [%]
%             Same layout.  Each subplot has its own colorscale so
%             spatial patterns show up even when magnitudes differ.
%             Saturates at 100 % — cells where 1st-order ≈ 0 are masked.
%
% Saved under Methods/Delta/<scenario>/delta_order_diff_absolute.png
%                                      delta_order_diff_relative.png

clear; close all;
try, cd(fileparts(mfilename('fullpath'))); catch; end

addpath(genpath('Shared_Utils'));

cfg = loadConfig();
[subset_names, subset_labels, ~] = subsetDefs();

N_DIST = numel(cfg.distributions);
N_SUB  = 7;
PW = 210;  PH = 170;   % pixels per panel

for si = 1:numel(cfg.scenarios)
    sc      = cfg.scenarios(si);
    SC_NAME = sc.name;

    fprintf('\n── %s ──\n', SC_NAME);

    % ── Pre-load all data ─────────────────────────────────────────────────
    data     = cell(N_DIST, N_SUB);
    r_km_ref = [];  z_m_ref = [];

    for di = 1:N_DIST
        dist    = cfg.distributions(di);
        res_dir = fullfile('Methods','Delta', SC_NAME, dist.name, 'results');

        for s = 1:N_SUB
            sfile = fullfile(res_dir, sprintf('delta_%s.mat', subset_names{s}));
            if ~isfile(sfile), continue; end
            try
                D = load(sfile, 'Var_TL','Vs_1st','r_km','z_m');
            catch
                continue;
            end

            if isempty(r_km_ref)
                r_km_ref = D.r_km;
                z_m_ref  = D.z_m;
            end

            diff_abs = D.Var_TL - D.Vs_1st;   % always ≥ 0
            rel = NaN(size(diff_abs));
            ok  = D.Vs_1st > 1e-8;
            rel(ok) = 100 * diff_abs(ok) ./ D.Vs_1st(ok);
            rel = min(rel, 100);

            data{di,s} = struct('diff_abs', diff_abs, 'diff_rel', rel, ...
                                'dist_name', dist.name);
        end
    end

    if isempty(r_km_ref)
        fprintf('  No data found — skip\n');
        continue;
    end

    out_dir = fullfile('Methods','Delta', SC_NAME);
    if ~exist(out_dir,'dir'), mkdir(out_dir); end

    FIG_W = PW * N_SUB;
    FIG_H = PH * N_DIST;

    % ── Figure 1: Absolute correction ─────────────────────────────────────
    fig1 = figure('Name', sprintf('Abs — %s', SC_NAME), ...
                  'Position',[30 30 FIG_W FIG_H]);

    for di = 1:N_DIST
        row_max = max(cellfun(@(d) max(d.diff_abs(:)), data(di, ~cellfun(@isempty, data(di,:)))));
        if isempty(row_max), row_max = 1e-10; end
        row_max = max(row_max, 1e-10);

        for s = 1:N_SUB
            ax = subplot(N_DIST, N_SUB, (di-1)*N_SUB + s);
            if isempty(data{di,s}), axis(ax,'off'); continue; end

            imagesc(ax, r_km_ref, z_m_ref, data{di,s}.diff_abs);
            colormap(ax, hot(256));
            clim(ax, [0, row_max]);
            set(ax,'YDir','reverse'); axis(ax,'tight');

            if s == 1
                ylabel(ax, strrep(data{di,s}.dist_name,'_',' '),'FontSize',7);
            else
                set(ax,'YTickLabel',[]);
            end
            if di == N_DIST, xlabel(ax,'Range (km)','FontSize',7);
            else,            set(ax,'XTickLabel',[]); end
            title(ax, subset_labels{s},'FontSize',7,'Interpreter','none');
            if s == N_SUB
                cb = colorbar(ax,'eastoutside');
                cb.Label.String = 'dB²'; cb.FontSize = 6;
            end
        end
    end

    sgtitle(fig1, sprintf('\\DeltaVar = Var_{2nd} − Var_{1st}  [dB²]  |  %s', SC_NAME), 'FontSize',11);
    out1 = fullfile(out_dir, 'delta_order_diff_absolute');
    saveFigPNG(fig1, out1);
    fprintf('  Saved: %s.png\n', out1);
    close(fig1);

    % ── Figure 2: Relative correction ─────────────────────────────────────
    fig2 = figure('Name', sprintf('Rel%% — %s', SC_NAME), ...
                  'Position',[60 60 FIG_W FIG_H]);

    for di = 1:N_DIST
        for s = 1:N_SUB
            ax = subplot(N_DIST, N_SUB, (di-1)*N_SUB + s);
            if isempty(data{di,s}), axis(ax,'off'); continue; end

            rel = data{di,s}.diff_rel;
            pmax = max(rel(~isnan(rel)));
            if isempty(pmax) || pmax == 0, pmax = 1; end

            imagesc(ax, r_km_ref, z_m_ref, rel);
            colormap(ax, parula(256));
            clim(ax, [0, pmax]);
            set(ax,'YDir','reverse'); axis(ax,'tight');

            if s == 1
                ylabel(ax, strrep(data{di,s}.dist_name,'_',' '),'FontSize',7);
            else
                set(ax,'YTickLabel',[]);
            end
            if di == N_DIST, xlabel(ax,'Range (km)','FontSize',7);
            else,            set(ax,'XTickLabel',[]); end
            title(ax, subset_labels{s},'FontSize',7,'Interpreter','none');
            if s == N_SUB
                cb = colorbar(ax,'eastoutside');
                cb.Label.String = '%'; cb.FontSize = 6;
            end
        end
    end

    sgtitle(fig2, sprintf('100 × ΔVar / Var_{1st}  [%%]  |  %s', SC_NAME), 'FontSize',11);
    out2 = fullfile(out_dir, 'delta_order_diff_relative');
    saveFigPNG(fig2, out2);
    fprintf('  Saved: %s.png\n', out2);
    close(fig2);

    % ── Console summary ───────────────────────────────────────────────────
    fprintf('\n  %-20s', 'Distribution');
    for s = 1:N_SUB, fprintf('  %9s', subset_labels{s}); end
    fprintf('\n');
    for di = 1:N_DIST
        dist = cfg.distributions(di);
        fprintf('  %-20s', dist.name);
        for s = 1:N_SUB
            if isempty(data{di,s})
                fprintf('  %9s', '—');
            else
                rel = data{di,s}.diff_rel;
                med = median(rel(~isnan(rel) & rel > 0));
                if isnan(med), med = 0; end
                fprintf('  %8.1f%%', med);
            end
        end
        fprintf('\n');
    end

    r_km_ref = [];  % reset for next scenario
end

fprintf('\n=== plot_delta_order_comparison.m complete ===\n');
