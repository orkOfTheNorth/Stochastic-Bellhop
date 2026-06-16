function overlayBathymetry(ax, bathy_m, max_depth_m)
% overlayBathymetry  Draw seafloor line + grey fill on an existing axes.
%
%   ax          — target axes (from gca or subplot)
%   bathy_m     — [N×2] matrix: [range_m, depth_m]  from bathymetryMaker
%   max_depth_m — plot depth limit (for the fill patch)

if isempty(bathy_m), return; end

r_km  = bathy_m(:,1) / 1000;
z_bty = bathy_m(:,2);

hold(ax, 'on');

% Grey filled patch below seafloor
patch(ax, [r_km; flipud(r_km)], ...
          [z_bty; ones(size(z_bty))*max_depth_m*1.05], ...
      [0.45 0.35 0.25], 'EdgeColor','none', 'FaceAlpha', 0.75, ...
      'DisplayName', 'Sediment');

% Bold seafloor line
plot(ax, r_km, z_bty, 'k-', 'LineWidth', 2, 'DisplayName', 'Seafloor');

hold(ax, 'off');
end
