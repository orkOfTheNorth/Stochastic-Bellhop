%% run_LN3.m — LN3 fit diagnostic for 3 unmixed subsets
%
% Loads MC TL_all samples for the 3 single-parameter subsets (freq, zS, svp)
% and produces one fit-quality figure per subset:
%
%   4×4 = 16 spatial sample locations (evenly spaced in range and depth)
%   Each cell: 20-bin histogram + LN3 PDF + FOM line
%   Title per cell: r, z, P_LN3, P_emp, KS statistic
%
% Outputs per scenario × distribution (3 PNGs):
%   Methods/LN3/<scen>/<dist>/figures/LN3_<subset>.png

function run_LN3(sc_target, dist_target)
%% run_LN3 — LN3 fit diagnostic.  No args = all combos.
if nargin < 2, sc_target = ''; dist_target = ''; end
close all; clc;
warning('off', 'MATLAB:unknownObjectIEEE');
warning('off', 'MATLAB:singularMatrix');
warning('off', 'MATLAB:rankDeficientMatrix');
try
    cd(fileparts(mfilename('fullpath')));
catch
end

addpath(genpath('Shared_Utils'));

cfg = loadConfig();

UNMIXED = {'freq', 'zS', 'svp'};   % the 3 single-parameter subsets only
N_BINS  = 20;
GRID_R  = 4;   % sampling grid in range
GRID_Z  = 4;   % sampling grid in depth  →  4x4 = 16 points



for si = 1:numel(cfg.scenarios)
    sc = cfg.scenarios(si);

    for di = 1:numel(cfg.distributions)
        dist = cfg.distributions(di);
        if ~isempty(sc_target) && ~(strcmp(sc.name,sc_target) && strcmp(dist.name,dist_target))
            continue;
        end

        fprintf('\n=== LN3 | %s | %s ===\n', sc.name, dist.name);

        mc_res  = fullfile('Methods','MC',  sc.name, dist.name, 'results');
        fig_dir = fullfile('Methods','LN3', sc.name, dist.name, 'figures');
        if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

        for ki = 1:numel(UNMIXED)
            sn = UNMIXED{ki};
            mc_file = fullfile(mc_res, sprintf('MC_%s.mat', sn));

            if ~isfile(mc_file)
                fprintf('  SKIP %s — MC file not found\n', sn);
                continue;
            end
            out_fig = fullfile(fig_dir, sprintf('LN3_%s.png', sn));
            if skipIfDone(out_fig, sprintf('%s / %s / %s', sc.name, dist.name, sn))
                continue;
            end

            try
                D = load(mc_file, 'TL_all','r_km','z_m','FOM','MC_PrFOM');
            catch ME
                fprintf('  SKIP %s — could not load file (%s)\n', sn, ME.message);
                continue;
            end
            TL_all = D.TL_all;
            r_km   = D.r_km;
            z_m    = D.z_m;
            FOM    = D.FOM;
            [Nz, Nr, N] = size(TL_all);

            %% Choose 16 sample locations (4x4 grid, avoid edges)
            ri_vec = round(linspace(Nr*0.10, Nr*0.90, GRID_R));
            zi_vec = round(linspace(Nz*0.10, Nz*0.85, GRID_Z));
            ri_vec = max(1, min(Nr, ri_vec));
            zi_vec = max(1, min(Nz, zi_vec));

            n_pts = GRID_R * GRID_Z;   % 16

            fig = figure('Position', [50 50 1600 1000]);
            t = tiledlayout(GRID_Z, GRID_R, 'TileSpacing','compact','Padding','compact');
            title(t, sprintf('LN3 Fit Quality | %s | %s | %s | N=%d | FOM=%ddB', ...
                             sc.name, dist.name, sn, N, FOM), ...
                  'FontSize', 11, 'Interpreter','none');

            pt_idx = 0;
            for zi_i = 1:GRID_Z
                for ri_i = 1:GRID_R
                    pt_idx = pt_idx + 1;
                    zi = zi_vec(zi_i);
                    ri = ri_vec(ri_i);

                    samps  = squeeze(TL_all(zi, ri, :));
                    p_emp  = D.MC_PrFOM(zi, ri);
                    [p_ln3, gam, mu_ln, sig_ln] = ln3fit(samps, FOM);
                    ks     = ksStatLN3(samps, gam, mu_ln, sig_ln);

                    ax = nexttile;
                    histogram(ax, samps, N_BINS, 'Normalization','pdf', ...
                              'FaceColor',[0.7 0.85 1.0], 'EdgeColor','none');
                    hold(ax,'on');

                    % LN3 PDF overlay
                    x_pdf = linspace(min(samps)*0.98, max(samps)*1.02, 300);
                    sh    = x_pdf - gam;
                    sh(sh<=0) = NaN;
                    y_pdf = lognpdf(sh, mu_ln, sig_ln);
                    plot(ax, x_pdf, y_pdf, 'b-', 'LineWidth', 1.5);

                    % FOM line
                    xline(ax, FOM, 'k--', 'LineWidth', 1.2);

                    hold(ax,'off');
                    xlabel(ax,'TL (dB)','FontSize',7);
                    ylabel(ax,'pdf','FontSize',7);

                    % Color-code title by KS quality
                    if ks < 0.10,     tc = [0.0 0.5 0.0];   % green — good
                    elseif ks < 0.15, tc = [0.8 0.5 0.0];   % orange — ok
                    else,             tc = [0.8 0.0 0.0];    % red — poor
                    end
                    title(ax, sprintf('r=%.1fkm z=%.0fm\nP_{LN3}=%.2f P_{emp}=%.2f KS=%.3f', ...
                                      r_km(ri), z_m(zi), p_ln3, p_emp, ks), ...
                          'Color', tc, 'FontSize', 7, 'Interpreter','none');
                end
            end

            out_path = fullfile(fig_dir, sprintf('LN3_%s', sn));
            saveFigPNG(fig, out_path);
            drawnow; close(fig);
            fprintf('  Saved: LN3_%s.png  (%d points)\n', sn, n_pts);
        end
    end
end


fprintf('\n=== run_LN3.m complete. ===\n');
end   % function run_LN3

%% ── KS statistic helper ───────────────────────────────────────────────────
function ks = ksStatLN3(samps, gam, mu_ln, sig_ln)
% Kolmogorov-Smirnov distance between empirical CDF and fitted LN3 CDF.
samps  = sort(double(samps(:)));
n      = numel(samps);
F_emp  = (1:n)' / n;
sh     = samps - gam;
sh(sh <= 0) = 1e-9;
F_ln3  = logncdf(sh, mu_ln, sig_ln);
ks     = max(abs(F_emp - F_ln3));
end
