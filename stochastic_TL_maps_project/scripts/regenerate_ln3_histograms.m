%% regenerate_ln3_histograms.m
% Regenerate LN3 histogram figures (all 36 combos) with y-axis clipped to
% histogram range so MC bars are always visible.
%
% Root cause of empty plots: when LN3 PDF has a very narrow peak (small sigma),
% MATLAB auto-scales y-axis to the PDF peak (e.g. 50 dB^-1) while histogram bars
% are only ~0.5 dB^-1 — making them invisible. Fix: clip ylim after plotting.
%
% Saves to Methods/LN3/<scen>/<dist>/figures/LN3_<param>.png (overwrites).

ROOT = fileparts(fileparts(mfilename('fullpath')));
try; cd(ROOT); catch; end
addpath(genpath(fullfile(ROOT,'core')));
addpath(fullfile(ROOT,'binaries'));

cfg    = loadConfig();
PARAMS = {'freq','zS','svp'};
N_BINS = 30;
GRID_R = 4;
GRID_Z = 4;

total = 0;

for si = 1:numel(cfg.scenarios)
    sc = cfg.scenarios(si);
    for di = 1:numel(cfg.distributions)
        dist = cfg.distributions(di);
        FOM  = getFOM(cfg, sc.name);

        fig_dir = fullfile('Methods','LN3', sc.name, dist.name, 'figures');
        if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

        mc_res = fullfile('Methods','MC', sc.name, dist.name, 'results');

        for ki = 1:numel(PARAMS)
            sn      = PARAMS{ki};
            mc_file = fullfile(mc_res, sprintf('MC_%s.mat', sn));
            if ~isfile(mc_file)
                fprintf('  SKIP %s/%s/%s — MC file missing\n', sc.name, dist.name, sn);
                continue;
            end

            tl_cache = fullfile('Cache', sc.name, dist.name, sprintf('TL_%s.mat', sn));
            if ~isfile(tl_cache)
                fprintf('  SKIP %s/%s/%s — TL cache missing\n', sc.name, dist.name, sn);
                continue;
            end

            fprintf('[%s/%s/%s] Loading data...\n', sc.name, dist.name, sn);

            D_tl = load(tl_cache, 'TL_save');
            D_mc = load(mc_file, 'r_km','z_m','MC_P_detect');
            TL_all   = double(D_tl.TL_save);
            r_km     = D_mc.r_km;
            z_m      = D_mc.z_m;
            MC_P_det = D_mc.MC_P_detect;
            clear D_tl D_mc;
            [Nz, Nr, N] = size(TL_all);

            ri_vec = round(linspace(Nr*0.10, Nr*0.90, GRID_R));
            zi_vec = round(linspace(Nz*0.10, Nz*0.85, GRID_Z));
            ri_vec = max(1, min(Nr, ri_vec));
            zi_vec = max(1, min(Nz, zi_vec));

            fig = figure('Position', [50 50 1600 1000], 'Visible','off');
            t = tiledlayout(GRID_Z, GRID_R, 'TileSpacing','compact','Padding','compact');
            title(t, sprintf('LN3 Fit Quality | %s | %s | %s | N=%d | FOM=%ddB', ...
                             sc.name, dist.name, sn, N, FOM), ...
                  'FontSize', 11, 'Interpreter','none');

            for zi_i = 1:GRID_Z
                for ri_i = 1:GRID_R
                    zi = zi_vec(zi_i);
                    ri = ri_vec(ri_i);

                    samps = squeeze(TL_all(zi, ri, :));
                    p_emp = MC_P_det(zi, ri);
                    [p_ln3, gam, mu_ln, sig_ln] = ln3Fit(samps, FOM);
                    ks    = ksStatLN3(samps, gam, mu_ln, sig_ln);

                    ax = nexttile;

                    % Histogram first — with more bins for visibility
                    h_obj = histogram(ax, samps, N_BINS, ...
                                      'Normalization','pdf', ...
                                      'FaceColor',[0.55 0.75 0.95], ...
                                      'EdgeColor',[0.3 0.5 0.8], ...
                                      'EdgeAlpha',0.4);
                    hold(ax,'on');

                    % PDF curve
                    x_lo  = min(samps) - 0.5*(max(samps)-min(samps));
                    x_hi  = max(samps) + 0.5*(max(samps)-min(samps));
                    x_lo  = max(x_lo, gam + 1e-3);
                    x_pdf = linspace(x_lo, x_hi, 400);
                    sh    = x_pdf - gam;
                    sh(sh<=0) = NaN;
                    y_pdf = lognpdf(sh, mu_ln, sig_ln);
                    plot(ax, x_pdf, y_pdf, 'b-', 'LineWidth', 1.8);

                    xline(ax, FOM, 'r--', 'LineWidth', 1.2);
                    hold(ax,'off');

                    %% KEY FIX: clip y-axis to histogram range so bars are visible
                    hist_vals = h_obj.Values;
                    if ~isempty(hist_vals) && max(hist_vals) > 0
                        ymax_hist = max(hist_vals) * 1.30;
                        ylim(ax, [0, ymax_hist]);
                    end
                    % Also clip x-axis to data range + small margin
                    xmargin = max(0.5*(max(samps)-min(samps)), 1);
                    xlim(ax, [min(samps)-xmargin*0.3, max(samps)+xmargin*0.3]);

                    xlabel(ax,'TL (dB)','FontSize',7);
                    ylabel(ax,'pdf','FontSize',7);

                    if ks < 0.10,     tc = [0.0 0.5 0.0];
                    elseif ks < 0.15, tc = [0.8 0.5 0.0];
                    else,             tc = [0.8 0.0 0.0]; end

                    title(ax, sprintf('r=%.1fkm z=%.0fm\nP_{det}(LN3)=%.2f P_{det}(MC)=%.2f KS=%.3f', ...
                                      r_km(ri), z_m(zi), p_ln3, p_emp, ks), ...
                          'Color', tc, 'FontSize', 7, 'Interpreter','none');
                end
            end

            out_path = fullfile(fig_dir, sprintf('LN3_%s', sn));
            saveFigPNG(fig, out_path);
            close(fig);
            total = total + 1;
            fprintf('  Saved: LN3_%s.png  (%d/36)\n', sn, total);

            clear TL_all;
        end
    end
end

fprintf('\nDone. %d LN3 histogram figures regenerated.\n', total);

%% ── KS statistic helper ──────────────────────────────────────────────────────
function ks = ksStatLN3(samps, gam, mu_ln, sig_ln)
samps = sort(double(samps(:)));
n     = numel(samps);
F_emp = (1:n)' / n;
sh    = samps - gam;
sh(sh<=0) = 1e-9;
F_ln3 = logncdf(sh, mu_ln, sig_ln);
ks    = max(abs(F_emp - F_ln3));
end
