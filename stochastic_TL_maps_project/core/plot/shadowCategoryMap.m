function shadowCategoryMap(ax, r_km, z_m, prob_map, label, thresholds)
% Alias for categoricalMap with color_style='shadow' (renamed during restructure).
    if nargin < 6, thresholds = []; end
    categoricalMap(ax, r_km, z_m, prob_map, label, thresholds, 'shadow');
end
