function shadowThresholdMaps(ax, r_km, z_m, prob_map, label, thresholds)
% Plot a continuous P(shadow) map with threshold contour overlays.
%
% ax         : axes handle
% r_km       : range vector (km)
% z_m        : depth vector (m)
% prob_map   : Nz×Nr probability matrix [0,1]
% label      : title string
% thresholds : probability thresholds to mark (default [0.70 0.80 0.90 0.95])

    if nargin < 6
        thresholds = [0.70 0.80 0.90 0.95];
    end

    pcolor(ax, r_km, z_m, prob_map); shading(ax, 'interp');
    set(ax, 'YDir', 'reverse');
    colormap(ax, parula(256)); clim(ax, [0 1]);
    cb = colorbar(ax); cb.Label.String = 'P(TL>FOM)'; cb.FontSize = 7;

    % Threshold contour colors: grey, orange, red, blue
    thr_colors = [0.5 0.5 0.5; 0.85 0.55 0.0; 0.85 0.15 0.15; 0.10 0.20 0.80];
    hold(ax, 'on');
    for k = 1:numel(thresholds)
        c = thr_colors(min(k, size(thr_colors,1)), :);
        try
            [~, hc] = contour(ax, r_km, z_m, prob_map, ...
                              [thresholds(k) thresholds(k)], ...
                              'Color', c, 'LineWidth', 1.4);
            hc.DisplayName = sprintf('%.0f%%', thresholds(k)*100);
        catch; end
    end
    hold(ax, 'off');
    legend(ax, 'show', 'Location', 'best', 'FontSize', 7);
    xlabel(ax, 'Range (km)'); ylabel(ax, 'Depth (m)');
    title(ax, label, 'FontSize', 9);
end
