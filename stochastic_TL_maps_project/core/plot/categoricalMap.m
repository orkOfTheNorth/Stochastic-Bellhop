function categoricalMap(ax, r_km, z_m, prob_map, label, thresholds, color_style)
% categoricalMap — plot a categorical probability map with N thresholds → N+1 bands.
%
% Unified replacement for shadowCategoryMap and detectionCategoryMap.
%
% Band 1 (white)  — P < thresholds(1)
% Band 2..N+1     — ascending colour confidence bands
%
% Inputs:
%   ax          — axes handle
%   r_km        — range vector (km), 1 × Nr
%   z_m         — depth vector (m), Nz × 1
%   prob_map    — [Nz × Nr] probability matrix in [0, 1]
%   label       — title string
%   thresholds  — ascending vector, default [0.50 0.60 0.70 0.80 0.90 0.95]
%   color_style — 'shadow' (yellow→red) | 'detect' (pale-green→dark-green)
%                 or [N × 3] custom RGB matrix for N threshold bands

if nargin < 6 || isempty(thresholds)
    thresholds = [0.50, 0.60, 0.70, 0.80, 0.90, 0.95];
end
if nargin < 7 || isempty(color_style)
    color_style = 'shadow';
end

thr   = sort(thresholds(:)');
n_thr = numel(thr);

% Integer index image: 1 = below lowest threshold, 2..n_thr+1 = successive bands
idx = ones(size(prob_map));
for k = 1:n_thr
    idx(prob_map >= thr(k)) = k + 1;
end

% Colour anchors per style
shadow_anchors = [
    1.00  1.00  0.80;   % pale yellow
    1.00  0.90  0.00;   % yellow
    1.00  0.65  0.00;   % yellow-orange
    1.00  0.38  0.00;   % orange
    0.90  0.18  0.00;   % deep orange
    0.75  0.05  0.05;   % red
    0.55  0.00  0.00;   % dark red
    0.35  0.00  0.00;   % very dark red
];
detect_anchors = [
    1.00  0.97  0.75;   % pale yellow
    1.00  0.85  0.30;   % yellow-orange
    1.00  0.60  0.00;   % orange
    0.95  0.38  0.00;   % deep orange
    0.85  0.18  0.00;   % orange-red
    0.70  0.05  0.00;   % red
    0.50  0.00  0.00;   % dark red
    0.30  0.00  0.00;   % very dark red
];

if isnumeric(color_style)
    % Caller-supplied [N × 3] band colours
    band_colors = color_style;
else
    switch lower(color_style)
        case 'detect'
            anchors = detect_anchors;
        otherwise  % 'shadow'
            anchors = shadow_anchors;
    end
    n_anchors = size(anchors, 1);
    if n_thr <= n_anchors
        band_colors = anchors(1:n_thr, :);
    else
        t_src = linspace(0, 1, n_anchors);
        t_dst = linspace(0, 1, n_thr);
        band_colors = interp1(t_src, anchors, t_dst);
    end
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
