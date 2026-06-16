function shadowCategoryMap(ax, r_km, z_m, prob_map, label, thresholds)
% Plot a categorical shadow-probability map with N thresholds → N+1 color bands.
%
% Band 1 (white)    — P < thr(1)          detection / below lowest threshold
% Band 2..N+1       — pale-yellow → red   ascending confidence of shadow
%
% Inputs:
%   ax         — axes handle
%   r_km       — range vector (km), 1 × Nr
%   z_m        — depth vector (m), Nz × 1
%   prob_map   — Nz × Nr probability matrix in [0,1]
%   label      — title string
%   thresholds — 1×N ascending vector, e.g. [0.50 0.60 0.70 0.80 0.90 0.95]

if nargin < 6 || isempty(thresholds)
    thresholds = [0.50, 0.60, 0.70, 0.80, 0.90, 0.95];
end
thr   = sort(thresholds(:)');   % ascending row vector
n_thr = numel(thr);

% Build integer index image: 1 = below lowest thr, 2..n_thr+1 = successive bands
idx = ones(size(prob_map));
for k = 1:n_thr
    idx(prob_map >= thr(k)) = k + 1;
end

% Colormap: white + gradient pale-yellow → deep red (n_thr+1 entries total)
% Pre-defined anchor colours; interpolated for any N ≤ 8
anchors = [
    1.00  1.00  0.80;   % pale yellow
    1.00  0.90  0.00;   % yellow
    1.00  0.65  0.00;   % yellow-orange
    1.00  0.38  0.00;   % orange
    0.90  0.18  0.00;   % deep orange
    0.75  0.05  0.05;   % red
    0.55  0.00  0.00;   % dark red
    0.35  0.00  0.00;   % very dark red
];
n_anchors = size(anchors, 1);
if n_thr <= n_anchors
    band_colors = anchors(1:n_thr, :);
else
    % Interpolate if more thresholds than anchors
    t_src = linspace(0, 1, n_anchors);
    t_dst = linspace(0, 1, n_thr);
    band_colors = interp1(t_src, anchors, t_dst);
end
cmap = [1 1 1; band_colors];   % white + N band colours

axes(ax); %#ok<MAXES>
imagesc(ax, r_km, z_m, idx);
colormap(ax, cmap);
clim(ax, [1, n_thr + 1]);
set(ax, 'YDir', 'reverse');
xlabel(ax, 'Range (km)');
ylabel(ax, 'Depth (m)');
title(ax, label, 'Interpreter', 'none');

% Legend — one patch per band (skip the white/below-threshold band)
hold(ax, 'on');
for k = 1:n_thr
    if k < n_thr
        lbl = sprintf('%.0f–%.0f%%', thr(k)*100, thr(k+1)*100);
    else
        lbl = sprintf('\\geq %.0f%%', thr(k)*100);
    end
    patch(ax, 'XData', NaN, 'YData', NaN, ...
          'FaceColor', cmap(k+1,:), 'EdgeColor', 'none', 'DisplayName', lbl);
end
legend(ax, 'Location', 'best', 'FontSize', 6);
hold(ax, 'off');
end
