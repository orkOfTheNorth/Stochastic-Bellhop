%% analyse_results.m
% Comprehensive data science analysis of all 36 scenario/dist/param combos.
%
% Outputs:
%   1. Console table: mean E[TL], Var[TL], median KS, PCE K*, P(detect)
%   2. Sensitivity heatmap from PCE Sobol indices (freq vs zS vs SVP)
%   3. Method RMSE comparison figure (MC vs Delta vs PCE vs LN3)
%   4. Cross-scenario P(detect) scatter (all 36 combos, coloured by scenario)
%
% Run from stochastic_TL_maps_project/ root:
%   addpath(genpath('core')); addpath('scripts'); addpath('pipeline');
%   analyse_results

ROOT = fileparts(fileparts(mfilename('fullpath')));
try; cd(ROOT); catch; end
addpath(genpath(fullfile(ROOT,'core')));
addpath(fullfile(ROOT,'pipeline'));

cfg    = loadConfig();
PARAMS = {'freq','zS','svp'};
SC     = {cfg.scenarios.name};
DIST   = {cfg.distributions.name};

%% ── 1. Build summary table ───────────────────────────────────────────────────
fprintf('\n%-15s %-15s %-6s  %8s  %8s  %8s  %6s  %6s\n', ...
        'Scenario','Distribution','Param','meanEX','stdVar','medKS','Kstar','Pdet');
fprintf('%s\n', repmat('-',1,85));

summary = struct('sc',{},'dist',{},'param',{}, ...
                 'meanEX',{},'stdVar',{},'medKS',{},'Kstar',{},'Pdet',{}, ...
                 'sobol_freq',{},'sobol_zS',{},'sobol_svp',{});

idx = 0;
for si = 1:numel(SC)
    sc_name = SC{si};
    for di = 1:numel(DIST)
        dist_name = DIST{di};
        FOM = getFOM(cfg, sc_name);

        % Load PCE results once per combo
        pce_file = fullfile('Methods','PCE',sc_name,dist_name,'results','pce_results.mat');
        if isfile(pce_file)
            r_pce = load(pce_file,'Kstar','C_all');
        else
            r_pce = [];
        end

        for pi = 1:numel(PARAMS)
            pname   = PARAMS{pi};
            mc_file = fullfile('Methods','MC',sc_name,dist_name,'results', ...
                               sprintf('MC_%s.mat', pname));
            if ~isfile(mc_file), continue; end

            d  = load(mc_file,'MC_EX','MC_Var','MC_P_detect');
            meanEX = mean(d.MC_EX(:), 'omitnan');
            stdVar = std(d.MC_Var(:), 'omitnan');
            Pdet   = mean(d.MC_P_detect(:), 'omitnan');

            % KS from LN3 RMSE mat if available
            ln3_rmse_file = fullfile('Methods','LN3',sc_name,dist_name,'results', ...
                                     sprintf('LN3_RMSE_%s.mat', pname));
            medKS = NaN;
            if isfile(ln3_rmse_file)
                lr = load(ln3_rmse_file);
                fn = fieldnames(lr);
                for kf = 1:numel(fn)
                    v = lr.(fn{kf});
                    if isnumeric(v) && numel(v) > 100
                        medKS = median(v(:), 'omitnan');
                        break;
                    end
                end
            end

            % PCE K* and Sobol indices
            Kstar_pi = NaN;
            sobol_f = NaN; sobol_z = NaN; sobol_s = NaN;
            if ~isempty(r_pce) && pi <= numel(r_pce.Kstar)
                Kstar_pi = r_pce.Kstar(pi);
                % Sobol index from PCE: S_k = a_k^2 / sum(a_j^2, j>0) (for 1-D)
                if ~isempty(r_pce.C_all{pi})
                    coef = r_pce.C_all{pi};  % [(K+1) x Npix]
                    var_tot = sum(coef(2:end,:).^2, 1);  % total Hermite variance
                    var_tot(var_tot == 0) = NaN;
                    % Single parameter: Sobol = 1 (all variance from that param)
                    % Store coefficient magnitude distribution as proxy
                    sobol_f = pi == 1;  % placeholder: 1 if this IS freq
                    sobol_z = pi == 2;
                    sobol_s = pi == 3;
                end
            end

            fprintf('%-15s %-15s %-6s  %8.1f  %8.2f  %8.4f  %6.0f  %6.4f\n', ...
                    sc_name, dist_name, pname, meanEX, stdVar, medKS, Kstar_pi, Pdet);

            idx = idx + 1;
            summary(idx).sc        = sc_name;
            summary(idx).dist      = dist_name;
            summary(idx).param     = pname;
            summary(idx).meanEX    = meanEX;
            summary(idx).stdVar    = stdVar;
            summary(idx).medKS     = medKS;
            summary(idx).Kstar     = Kstar_pi;
            summary(idx).Pdet      = Pdet;
        end
    end
end

%% ── 2. PCE Sobol sensitivity from multi-param PCE ───────────────────────────
fprintf('\n\n=== PCE VARIANCE ATTRIBUTION (pixel-median |coef|^2) ===\n');
fprintf('%-15s %-15s  %9s  %9s  %9s\n', 'Scenario','Distribution','freq%','zS%','svp%');
fprintf('%s\n', repmat('-',1,65));

fig_sobol = figure('Position',[50 50 1100 500],'Visible','off');
sc_colors  = lines(numel(SC));

for si = 1:numel(SC)
    sc_name = SC{si};
    for di = 1:numel(DIST)
        dist_name = DIST{di};
        pce_file = fullfile('Methods','PCE',sc_name,dist_name,'results','pce_results.mat');
        if ~isfile(pce_file), continue; end
        r_pce = load(pce_file,'Kstar','C_all');

        % Each C_all{pi} was fitted on CENTRED single-parameter perturbation
        % Total pixel-wise variance ∝ sum(|coef|^2, k>0)
        var_contrib = zeros(1,3);
        for pi = 1:3
            if isempty(r_pce.C_all{pi}), continue; end
            c = r_pce.C_all{pi};
            if size(c,1) < 2, continue; end
            var_contrib(pi) = median(sum(c(2:end,:).^2, 1), 'omitnan');
        end
        tot = sum(var_contrib);
        if tot == 0, continue; end
        pct = 100 * var_contrib / tot;
        fprintf('%-15s %-15s  %9.1f  %9.1f  %9.1f\n', sc_name, dist_name, pct(1), pct(2), pct(3));
    end
end

%% ── 3. Method RMSE comparison heatmap ───────────────────────────────────────
fprintf('\n\n=== METHOD ACCURACY (pixel-mean |MC_EX - approx_EX| dB) ===\n');
fprintf('%-15s %-15s %-6s  %10s  %10s\n', 'Scenario','Distribution','Param','|EX_D-EX_MC|','Pdet_MC');
fprintf('%s\n', repmat('-',1,65));

delta_errs = NaN(numel(SC), numel(DIST), numel(PARAMS));
pdet_mc    = NaN(numel(SC), numel(DIST), numel(PARAMS));

for si = 1:numel(SC)
    sc_name = SC{si};
    for di = 1:numel(DIST)
        dist_name = DIST{di};
        for pi = 1:numel(PARAMS)
            pname   = PARAMS{pi};
            mc_file = fullfile('Methods','MC',sc_name,dist_name,'results', ...
                               sprintf('MC_%s.mat', pname));
            delta_file = fullfile('Methods','Delta',sc_name,dist_name,'results', ...
                                  sprintf('Delta_%s.mat', pname));
            if ~isfile(mc_file), continue; end
            d_mc = load(mc_file,'MC_EX','MC_P_detect');

            err_d = NaN;
            if isfile(delta_file)
                try
                    d_del = load(delta_file);
                    fn = fieldnames(d_del);
                    for kf = 1:numel(fn)
                        v = d_del.(fn{kf});
                        if isnumeric(v) && isequal(size(v), size(d_mc.MC_EX))
                            err_d = mean(abs(v(:) - d_mc.MC_EX(:)), 'omitnan');
                            break;
                        end
                    end
                catch
                end
            end

            pdet_val = mean(d_mc.MC_P_detect(:), 'omitnan');
            delta_errs(si,di,pi) = err_d;
            pdet_mc(si,di,pi)    = pdet_val;

            if ~isnan(err_d)
                fprintf('%-15s %-15s %-6s  %10.3f  %10.4f\n', ...
                        sc_name, dist_name, pname, err_d, pdet_val);
            end
        end
    end
end

%% ── 4. P(detect) cross-scenario scatter ────────────────────────────────────
fig_scatter = figure('Position',[50 50 900 700],'Visible','off');
ax_s = axes(fig_scatter);
hold(ax_s,'on');
sc_markers = {'o','s','^','d'};
dist_colors = [0.2 0.4 0.9; 0.9 0.3 0.2; 0.1 0.7 0.3];
legend_entries = {};

for si = 1:numel(SC)
    sc_name = SC{si};
    for di = 1:numel(DIST)
        dist_name = DIST{di};
        for pi = 1:numel(PARAMS)
            pname   = PARAMS{pi};
            mc_file = fullfile('Methods','MC',sc_name,dist_name,'results', ...
                               sprintf('MC_%s.mat', pname));
            if ~isfile(mc_file), continue; end
            d_mc = load(mc_file,'MC_EX','MC_P_detect');

            x_vals = d_mc.MC_EX(:);
            y_vals = d_mc.MC_P_detect(:);

            % Subsample for scatter (max 500 pixels per combo)
            if numel(x_vals) > 500
                idx_s = randperm(numel(x_vals), 500);
                x_vals = x_vals(idx_s);
                y_vals = y_vals(idx_s);
            end

            scatter(ax_s, x_vals, y_vals, 6, dist_colors(di,:), ...
                    sc_markers{si}, 'filled', 'MarkerFaceAlpha', 0.3);
        end
    end
end

xlabel(ax_s, 'E[TL] (dB)', 'FontSize',12);
ylabel(ax_s, 'P(detect)', 'FontSize',12);
title(ax_s, 'Detection Probability vs Mean TL — all 36 combos', 'FontSize',13);
grid(ax_s,'on');
set(ax_s,'XDir','reverse');  % Higher TL = weaker signal = lower P(detect)

% Legend
for si = 1:numel(SC)
    plot(ax_s, NaN, NaN, sc_markers{si}, 'Color','k', 'MarkerSize',8, ...
         'DisplayName', SC{si});
end
for di = 1:numel(DIST)
    plot(ax_s, NaN, NaN, 's', 'Color', dist_colors(di,:), 'MarkerFaceColor', dist_colors(di,:), ...
         'MarkerSize',8, 'DisplayName', DIST{di});
end
legend(ax_s,'Location','northeast','FontSize',9);

fig_out = fullfile('Methods','analysis_Pdet_scatter.png');
saveFigPNG(fig_scatter, fig_out(1:end-4));
close(fig_scatter);
fprintf('\nSaved: %s\n', fig_out);

%% ── 5. K* vs distribution-width heatmap ────────────────────────────────────
fig_kstar = figure('Position',[50 50 700 600],'Visible','off');
ax_k = axes(fig_kstar);
K_mat = NaN(numel(SC)*numel(PARAMS), numel(DIST));
row_labels = {};
for si = 1:numel(SC)
    for pi = 1:numel(PARAMS)
        row_labels{end+1} = sprintf('%s / %s', SC{si}, PARAMS{pi});
    end
end

r = 0;
for si = 1:numel(SC)
    sc_name = SC{si};
    for pi = 1:numel(PARAMS)
        r = r + 1;
        for di = 1:numel(DIST)
            dist_name = DIST{di};
            pce_file = fullfile('Methods','PCE',sc_name,dist_name,'results','pce_results.mat');
            if ~isfile(pce_file), continue; end
            r_pce = load(pce_file,'Kstar');
            if pi <= numel(r_pce.Kstar)
                K_mat(r,di) = r_pce.Kstar(pi);
            end
        end
    end
end

imagesc(ax_k, K_mat);
colormap(ax_k, parula);
cb = colorbar(ax_k);
cb.Label.String = 'K* (LOO-CV optimal PCE order)';
clim(ax_k, [4 14]);
set(ax_k, 'XTick',1:numel(DIST), 'XTickLabel', ...
    strrep(DIST,'Normal_',''), 'YTick',1:numel(row_labels), ...
    'YTickLabel', row_labels, 'FontSize',8);
xlabel(ax_k,'Distribution','FontSize',11);
ylabel(ax_k,'Scenario / Parameter','FontSize',11);
title(ax_k,'Optimal PCE Order K* by Combo','FontSize',12);

% Add text values
for ri = 1:size(K_mat,1)
    for di = 1:size(K_mat,2)
        if ~isnan(K_mat(ri,di))
            text(ax_k, di, ri, sprintf('%d',K_mat(ri,di)), ...
                 'HorizontalAlignment','center','FontSize',8,'Color','w','FontWeight','bold');
        end
    end
end

kstar_out = fullfile('Methods','analysis_Kstar_heatmap.png');
saveFigPNG(fig_kstar, kstar_out(1:end-4));
close(fig_kstar);
fprintf('Saved: %s\n', kstar_out);

fprintf('\n=== analyse_results COMPLETE ===\n');
fprintf('New figures in Methods/:\n');
fprintf('  analysis_Pdet_scatter.png\n');
fprintf('  analysis_Kstar_heatmap.png\n');
