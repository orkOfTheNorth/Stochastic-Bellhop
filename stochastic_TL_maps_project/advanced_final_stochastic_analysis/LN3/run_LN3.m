%% run_LN3.m  –  3-Parameter Lognormal Fit Quality  (per subset)
%
% For each of the 7 subsets, loads the MC TL cube and shows the LN3 fit
% at a 4×4 grid of representative (range, depth) points:
%
%   Columns = 4 ranges:  ~10%, 30%, 60%, 90% of max range
%   Rows    = 4 depths:  ~10%, 30%, 60%, 85% of max depth
%
% Each panel: histogram of N TL samples + fitted LN3 PDF + FOM line
%
% Output:  figures/LN3_grid_fits_<subset>.png   (7 files)
%          figures/LN3_probability_maps.png      (composite P(shadow) all subsets)
%
% Requires: MC results in ../Monte_Carlo/results/MC_<subset>.mat
%           (produced by ../Monte_Carlo/run_MC.m)

clear; close all; clc; warning('off');
try, cd(fileparts(mfilename('fullpath'))); catch; end

base_dir = fullfile('..');
addpath(genpath(fullfile(base_dir, 'Functions')));

if ~exist('figures','dir'), mkdir('figures'); end

subset_names  = {'zS','freq','svp','zS_freq','zS_svp','freq_svp','zS_freq_svp'};
subset_labels = {'z_S only','Freq only','SVP only', ...
                 'z_S+Freq','z_S+SVP','Freq+SVP','z_S+Freq+SVP'};

% Fractional positions along range and depth axes
r_fracs = [0.10, 0.30, 0.60, 0.90];
z_fracs = [0.10, 0.30, 0.60, 0.85];

Nr_pts = numel(r_fracs);
Nz_pts = numel(z_fracs);

all_ln3 = struct();

for s = 1:7
    sn    = subset_names{s};
    fname = fullfile('..','Monte_Carlo','results', sprintf('MC_%s.mat', sn));

    if ~isfile(fname)
        fprintf('[SKIP] %s — results not found (%s)\n', sn, fname);
        continue;
    end
    fprintf('Loading %s ...\n', fname);
    d = load(fname, 'TL_all','r_grid','z_grid','r_km','z_m','FOM','N');

    TL_all = double(d.TL_all);
    r_km   = d.r_km;
    z_m    = d.z_m;
    FOM    = d.FOM;
    N      = d.N;
    Nr     = numel(r_km);
    Nz     = numel(z_m);

    % ── Select sample points ─────────────────────────────────────────────
    ri_pts = round(r_fracs * (Nr-1)) + 1;
    zi_pts = round(z_fracs * (Nz-1)) + 1;
    ri_pts = min(max(ri_pts, 1), Nr);
    zi_pts = min(max(zi_pts, 1), Nz);

    % ── Build grid-fit figure ─────────────────────────────────────────────
    fig = figure('Name', sprintf('LN3 Grid Fits – %s', subset_labels{s}), ...
                 'Visible', 'off', 'Position', [50 50 1200 900]);

    for zi_idx = 1:Nz_pts
        for ri_idx = 1:Nr_pts
            zi = zi_pts(zi_idx);
            ri = ri_pts(ri_idx);

            samps = squeeze(TL_all(zi, ri, :));
            samps = samps(isfinite(samps));
            if numel(samps) < 5
                continue;
            end

            % LN3 fit
            gam  = min(samps) * 0.95;
            sh   = samps - gam;
            sh(sh <= 0) = 1e-6;
            p    = lognfit(sh);
            mu_ln = p(1);  sig_ln = p(2);

            % P(TL > FOM)
            if FOM > gam
                prob_ln3 = 1 - logncdf(FOM - gam, mu_ln, sig_ln);
            else
                prob_ln3 = 1.0;
            end
            prob_emp = mean(samps > FOM);

            % Plot
            ax = subplot(Nz_pts, Nr_pts, (zi_idx-1)*Nr_pts + ri_idx);
            histogram(ax, samps, min(N, 15), 'Normalization', 'pdf', ...
                      'FaceColor', [0.55 0.75 0.95], 'EdgeColor', 'w');
            hold(ax, 'on');

            % LN3 PDF
            x_lo = max(gam + 0.01, min(samps) - 2);
            x_hi = max(samps) + 2;
            x_pdf = linspace(x_lo, x_hi, 300);
            y_pdf = lognpdf(x_pdf - gam, mu_ln, sig_ln);
            y_pdf(x_pdf <= gam) = 0;
            plot(ax, x_pdf, y_pdf, 'r-', 'LineWidth', 1.5);

            % FOM line
            yl = ylim(ax);
            plot(ax, [FOM FOM], [0 yl(2)*0.95], 'k--', 'LineWidth', 1.2);
            hold(ax, 'off');

            xlabel(ax, 'TL (dB)', 'FontSize', 7);
            title(ax, sprintf('r=%.0fkm z=%.0fm\nLN3=%.2f Emp=%.2f', ...
                  r_km(ri), z_m(zi), prob_ln3, prob_emp), 'FontSize', 7);
            set(ax, 'FontSize', 7);
        end
    end
    sgtitle(sprintf('LN3 Grid Fits  |  %s  |  N=%d  FOM=%ddB', subset_labels{s}, N, FOM), ...
            'FontSize', 10);

    out_png = fullfile('figures', sprintf('LN3_grid_fits_%s.png', sn));
    exportgraphics(fig, out_png, 'Resolution', 300);
    close(fig);
    fprintf('  Saved %s\n', out_png);

    % Store LN3_prob map for composite
    LN3_map = zeros(Nz, Nr);
    for zi = 1:Nz
        for ri = 1:Nr
            LN3_map(zi,ri) = ln3prob(squeeze(TL_all(zi,ri,:)), FOM);
        end
        if mod(zi,50)==0
            fprintf('  LN3 map: %d/%d rows\n', zi, Nz);
        end
    end
    all_ln3.(sn) = LN3_map;
    fprintf('  LN3 map done: %s\n', sn);
end

%% ── COMPOSITE P(shadow) MAP ─────────────────────────────────────────────────
fn = fieldnames(all_ln3);
if ~isempty(fn)
    N_s = numel(subset_names);
    fw  = max(1400, 200*N_s);
    fig_c = figure('Name','LN3 P(shadow) All Subsets','Visible','off', ...
                   'Position',[50 50 fw 380]);
    % Need r_km/z_m — load from last successful subset
    last = load(fullfile('..','Monte_Carlo','results', ...
                         sprintf('MC_%s.mat', fn{end})), 'r_km','z_m','FOM');
    r_km = last.r_km;  z_m = last.z_m;  FOM = last.FOM;
    thresholds = [0.70 0.80 0.90 0.95];
    for s = 1:N_s
        sn = subset_names{s};
        if ~isfield(all_ln3, sn), continue; end
        ax = subplot(1, N_s, s);
        shadowThresholdMaps(ax, r_km, z_m, all_ln3.(sn), subset_labels{s}, thresholds);
    end
    sgtitle(sprintf('LN3 P(TL>FOM=%ddB)  Thresholds: 70/80/90/95%%', FOM));
    saveFig(fig_c, fullfile('figures','LN3_probability_maps'));  close(fig_c);
    fprintf('\nSaved composite: figures/LN3_probability_maps\n');
end

fprintf('\n=== LN3 complete.  Figures in LN3/figures/ ===\n');
