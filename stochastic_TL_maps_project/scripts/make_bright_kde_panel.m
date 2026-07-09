%% make_bright_kde_panel.m
% Generates a 4×4 grid of empirical KDE curves for TL_zS
% at 16 sampled (range, depth) locations — deep_water / Normal_5pct.
%
% No histogram, no parametric fits.  Four vivid colours, one per depth row.
% Output: for_presentation/LN3_empirical_bright_zS.png

ROOT = fileparts(fileparts(mfilename('fullpath')));
try; cd(ROOT); catch; end
addpath(genpath(fullfile(ROOT,'core')));

%% ── Config ──────────────────────────────────────────────────────────────────
SC    = 'deep_water';
DIST  = 'Normal_5pct';
PARAM = 'zS';
FOM   = getFOM(loadConfig(), SC);   % dB

TL_FILE = fullfile('Cache', SC, DIST, sprintf('TL_%s.mat', PARAM));
MC_FILE = fullfile('Methods','MC', SC, DIST,'results', sprintf('MC_%s.mat', PARAM));

OUT_DIR = fullfile('for_presentation');
if ~exist(OUT_DIR,'dir'), mkdir(OUT_DIR); end
OUT_PNG = fullfile(OUT_DIR, sprintf('LN3_empirical_bright_%s.png', PARAM));

%% ── Load data ────────────────────────────────────────────────────────────────
D    = load(TL_FILE, 'TL_save');
DMC  = load(MC_FILE, 'r_km','z_m');
TL   = double(D.TL_save);          % [Nz × Nr × N]
r_km = DMC.r_km;
z_m  = DMC.z_m;
[Nz, Nr, ~] = size(TL);

%% ── Spatial grid: 4 ranges × 4 depths ───────────────────────────────────────
ri_vec = round(linspace(Nr*0.15, Nr*0.85, 4));
zi_vec = round(linspace(Nz*0.25, Nz*0.75, 4));
ri_vec = max(1, min(Nr, ri_vec));
zi_vec = max(1, min(Nz, zi_vec));

%% ── Four vivid colours (one per depth row) ───────────────────────────────────
COLORS = [
    0.98  0.20  0.20;   % vivid red     (shallowest)
    1.00  0.55  0.00;   % vivid orange
    0.05  0.75  0.20;   % vivid green
    0.10  0.45  1.00;   % vivid blue    (deepest)
];
ALPHA_FILL = 0.22;      % filled area transparency

%% ── Figure ───────────────────────────────────────────────────────────────────
fig = figure('Position',[50 50 1400 900],'Color','white');
tl  = tiledlayout(4, 4, 'TileSpacing','tight','Padding','compact');
title(tl, ...
    sprintf('MC empirical TL | %s | %s | %s | FOM = %d dB', ...
            SC, DIST, PARAM, FOM), ...
    'FontSize', 12, 'FontWeight','bold', 'Interpreter','none');

for zi_i = 1:4
    c  = COLORS(zi_i, :);
    zi = zi_vec(zi_i);

    for ri_i = 1:4
        ri    = ri_vec(ri_i);
        samps = squeeze(TL(zi, ri, :));

        % KDE on a fine grid
        bw    = 1.06 * std(samps) * numel(samps)^(-0.2);   % Silverman's rule
        bw    = max(bw, 0.05);
        [f, xi] = ksdensity(samps, 'Bandwidth', bw, 'NumPoints', 512);

        ax = nexttile;

        % Filled area under KDE
        patch(ax, [xi, xi(end), xi(1)], [f, 0, 0], c, ...
              'FaceAlpha', ALPHA_FILL, 'EdgeColor','none');
        hold(ax,'on');

        % Bright KDE curve
        plot(ax, xi, f, '-', 'Color', c, 'LineWidth', 2.2);

        % Shade P(TL > FOM) in a stronger fill
        mask = xi >= FOM;
        if any(mask)
            xi_s = [FOM, xi(mask), xi(end)];
            f_s  = [interp1(xi, f, FOM,'linear',0), f(mask), 0];
            patch(ax, xi_s, f_s, c, 'FaceAlpha', 0.55, 'EdgeColor','none');
        end

        % FOM vertical line
        xline(ax, FOM, ':', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.4);

        hold(ax,'off');
        box(ax,'off');
        ax.TickDir = 'out';
        ax.FontSize = 7;

        % P(shadow) from empirical samples
        p_emp = mean(samps > FOM);
        xlabel(ax, 'TL (dB)', 'FontSize', 7);
        if ri_i == 1, ylabel(ax, 'pdf', 'FontSize', 7); end
        title(ax, ...
            sprintf('r=%.1f km   z=%d m\n\\itP\\rm(TL>FOM) = %.2f', ...
                    r_km(ri), round(z_m(zi)), p_emp), ...
            'FontSize', 6.5, 'FontWeight','normal');
    end
end

% Colour legend (one entry per depth row)
dummy_ax = axes(fig,'Visible','off','Position',[0 0 0 0]);
h = gobjects(4,1);
for k = 1:4
    h(k) = patch(dummy_ax, NaN, NaN, COLORS(k,:), ...
                 'FaceAlpha',0.55,'EdgeColor',COLORS(k,:),'LineWidth',2);
end
legend(dummy_ax, h, ...
    arrayfun(@(z) sprintf('z = %d m', round(z_m(zi_vec(z)))), 1:4, ...
             'UniformOutput',false), ...
    'Location','southoutside', 'Orientation','horizontal', ...
    'FontSize', 9, 'Box','off', ...
    'Position',[0.25 0.00 0.50 0.025]);

%% ── Save ─────────────────────────────────────────────────────────────────────
saveFigPNG(fig, strrep(OUT_PNG,'.png',''));
close(fig);
fprintf('Saved: %s\n', OUT_PNG);
