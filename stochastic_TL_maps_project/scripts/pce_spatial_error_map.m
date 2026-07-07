function pce_spatial_error_map(dist_name, param)
%% pce_spatial_error_map  PCE(K*) spatial relative-error map, all 4 scenarios.
%
% Shows WHERE PCE(K*) reconstruction agrees/disagrees with MC — the
% question "does PCE actually perform well near the convergence zones
% (CZ), where the physics matters most, or just on average?" — answered
% spatially instead of with a single scalar.
%
% Uses the already-fitted PCE reconstruction (Methods/PCE/<sc>/<dist>/
% results/pce_results.mat, field VarK_best{param_idx} — the Var[TL] field
% reconstructed at the LOO-CV-selected order K*, no refit needed) against
% the true MC variance field (Methods/MC/<sc>/<dist>/results/MC_<param>.mat,
% field MC_Var).
%
% Metric: relative error |Var_PCE - Var_MC| / Var_MC per pixel, %.
% Seafloor overlaid (critical for upslope/downslope — see run_comparison.m).
%
% Output: Methods/PCE/pce_spatial_error_<dist>_<param>.png (2x2 grid,
%         one panel per scenario: deep_water, shallow_water, upslope, downslope)
%
% Usage: pce_spatial_error_map                      % Normal_5pct, svp (default)
%        pce_spatial_error_map('Normal_10pct','zS')

if nargin < 1 || isempty(dist_name), dist_name = 'Normal_5pct'; end
if nargin < 2 || isempty(param),     param     = 'svp';         end

ROOT = fileparts(fileparts(mfilename('fullpath')));
try; cd(ROOT); catch; end
addpath(genpath(fullfile(ROOT,'core')));
set(0,'DefaultFigureVisible','off');

cfg = loadConfig();
PARAMS = {'freq','zS','svp'};
pi_idx = find(strcmp(PARAMS, param), 1);
assert(~isempty(pi_idx), 'Unknown param: %s', param);

out_dir = fullfile('Methods','PCE');
if ~exist(out_dir,'dir'), mkdir(out_dir); end
out_png = fullfile(out_dir, sprintf('pce_spatial_error_%s_%s.png', dist_name, param));

sc_names = {cfg.scenarios.name};
fig = figure('Position',[50 50 1100 850]);

all_relerr = {};
for si = 1:numel(cfg.scenarios)
    sc = cfg.scenarios(si);
    pce_file = fullfile('Methods','PCE', sc.name, dist_name, 'results','pce_results.mat');
    mc_file  = fullfile('Methods','MC',  sc.name, dist_name, 'results', sprintf('MC_%s.mat', param));
    if ~isfile(pce_file) || ~isfile(mc_file)
        fprintf('[MISS] %s / %s\n', pce_file, mc_file);
        continue;
    end
    P = load(pce_file, 'VarK_best', 'r_km', 'z_m', 'Kstar');
    M = load(mc_file, 'MC_Var');

    Var_pce = P.VarK_best{pi_idx};
    Var_mc  = M.MC_Var;
    relerr  = 100 * abs(Var_pce - Var_mc) ./ max(Var_mc, eps);   % %
    all_relerr{end+1} = relerr; %#ok<AGROW>

    ax = subplot(2,2,si);
    pcolor(ax, P.r_km, P.z_m, relerr);
    shading(ax,'interp'); set(ax,'YDir','reverse');
    colormap(ax, hot(256));
    clim(ax, [0, 30]);   % 0-30% relative error range; >30% saturates red
    cb = colorbar(ax); cb.Label.String = 'rel. error (%)'; cb.FontSize = 7;

    bathy = bathymetryMaker(sc.bathy_type, sc.maxR_m);
    overlayBathymetry(ax, bathy, sc.maxDepth_m);

    title(ax, sprintf('%s (K^*=%d)', strrep(sc.name,'_',' '), P.Kstar(pi_idx)), ...
          'FontSize',10,'Interpreter','tex');
    xlabel(ax,'Range (km)'); ylabel(ax,'Depth (m)');
    fprintf('%s: median rel. error = %.2f%%, mean = %.2f%%\n', ...
            sc.name, median(relerr(:)), mean(relerr(:)));
end

sgtitle(sprintf('PCE(K^*) relative error in Var[TL] vs MC | %s | param=%s | seafloor overlaid', ...
        dist_name, param), 'Interpreter','tex', 'FontSize',12);

print(fig, out_png, '-dpng', '-r150');
close(fig);
fprintf('Saved: %s\n', out_png);
end
