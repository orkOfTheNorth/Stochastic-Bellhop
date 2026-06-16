function detectionCategoryMap(ax, r_km, z_m, prob_map, label, thresholds)
% Plot a categorical detection-probability map with N thresholds → N+1 color bands.
%
% Band 1 (white)  — P_detect < thr(1)         (no reliable detection)
% Band 2..N+1     — pale-green → dark green    (increasing detection confidence)
%
% Inputs:
%   prob_map — P(TL < FOM) = P(detection), in [0,1]

if nargin < 6 || isempty(thresholds)
    thresholds = [0.50, 0.60, 0.70, 0.80, 0.90, 0.95];
end
thr   = sort(thresholds(:)');
n_thr = numel(thr);

idx = ones(size(prob_map));
for k = 1:n_thr
    idx(prob_map >= thr(k)) = k + 1;
end

% Colormap: white + pale-green → deep green (N+1 entries)
anchors = [
    0.85  1.00  0.85;   % pale green
    0.60  0.94  0.60;   % light green
    0.30  0.85  0.30;   % medium green
    0.10  0.70  0.10;   % green
    0.00  0.55  0.00;   % dark green
    0.00  0.39  0.00;   % very dark green
    0.00  0.27  0.00;   % forest green
    0.00  0.18  0.00;   % deep forest
];
n_anchors = size(anchors, 1);
if n_thr <= n_anchors
    band_colors = anchors(1:n_thr, :);
else
    t_src = linspace(0,1,n_anchors);
    t_dst = linspace(0,1,n_thr);
    band_colors = interp1(t_src, anchors, t_dst);
end
cmap = [1 1 1; band_colors];

axes(ax); %#ok<MAXES>
imagesc(ax, r_km, z_m, idx);
colormap(ax, cmap);
clim(ax, [1, n_thr + 1]);
set(ax, 'YDir', 'reverse');
xlabel(ax, 'Range (km)');
ylabel(ax, 'Depth (m)');
title(ax, label, 'Interpreter', 'none');

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
