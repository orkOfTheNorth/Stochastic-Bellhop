%% make_ex_annotated_circles.m
% Recreate the MC E[TL] map for downslope | Normal_5pct | freq
% and overlay white circles at the 16 LN3 fit-quality sample positions.
% Saves to docs/pptx_figures/downslope_EX_annotated/image1.png

close all;
ROOT = fileparts(fileparts(mfilename('fullpath')));
try; cd(ROOT); catch, end
addpath(genpath(fullfile(ROOT,'core')));
set(0,'DefaultFigureVisible','off');

SC   = 'downslope';
DIST = 'Normal_5pct';
PARAM = 'freq';

% ── load TL cube ────────────────────────────────────────────────────────
cache_file = fullfile('Cache', SC, DIST, sprintf('TL_%s.mat', PARAM));
mc_file    = fullfile('Methods','MC', SC, DIST, 'results', ...
                      sprintf('MC_%s.mat', PARAM));

fprintf('Loading MC results...\n');
Dm   = load(mc_file, 'r_km','z_m','MC_EX');
r_km = Dm.r_km;
z_m  = Dm.z_m;
EX   = Dm.MC_EX;    % [Nz x Nr]
[Nz, Nr] = size(EX);

% ── LN3 fit-quality sample grid (same as run_ln3.m) ──────────────────
GRID_R = 4;
GRID_Z = 4;
ri_vec = round(linspace(Nr*0.10, Nr*0.90, GRID_R));
zi_vec = round(linspace(Nz*0.10, Nz*0.85, GRID_Z));
ri_vec = max(1, min(Nr, ri_vec));
zi_vec = max(1, min(Nz, zi_vec));

% physical coordinates of the 16 sample pixels
sample_r = r_km(ri_vec);   % km
sample_z = z_m(zi_vec);    % m

% ── figure ──────────────────────────────────────────────────────────────
fig = figure('Position',[50 50 900 620]);
ax  = axes(fig);

imagesc(ax, r_km, z_m, EX);
set(ax,'YDir','reverse');
colormap(ax, jet);
cb = colorbar(ax);
cb.Label.String = 'E[TL] (dB)';
clim(ax, [50 150]);

xlabel(ax,'Range (km)');
ylabel(ax,'Depth (m)');
title(ax, sprintf('MC E[TL] | %s | %s | %s | N=1000', SC, DIST, PARAM), ...
      'Interpreter','none');

% draw seafloor (downslope: depth increases linearly)
hold(ax,'on');
bathy_r = [r_km(1), r_km(end)];
bathy_d = [50, 550];   % 50m -> 550m downslope
fill(ax, [bathy_r fliplr(bathy_r)], ...
         [bathy_d, max(z_m)*[1 1]], ...
         [0.45 0.27 0.07], 'EdgeColor','none');

% ── white circles at LN3 sample positions ──────────────────────────────
for zi_i = 1:GRID_Z
    for ri_i = 1:GRID_R
        plot(ax, sample_r(ri_i), sample_z(zi_i), 'o', ...
             'MarkerSize', 14, ...
             'MarkerEdgeColor', 'white', ...
             'MarkerFaceColor', 'none', ...
             'LineWidth', 2.5);
    end
end
hold(ax,'off');

% ── save ────────────────────────────────────────────────────────────────
out_dir = fullfile('docs','pptx_figures','downslope_EX_annotated');
if ~exist(out_dir,'dir'), mkdir(out_dir); end
out_png = fullfile(out_dir, 'image1.png');

print(fig, out_png, '-dpng', '-r200');
close(fig);
fprintf('Saved: %s\n', out_png);
