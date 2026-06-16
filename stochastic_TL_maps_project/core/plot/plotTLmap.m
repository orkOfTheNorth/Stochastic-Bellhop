function plotTLmap(ax, r_km, z_m, TL, FOM)
% Standard TL heatmap on the given axes.
%
%   pcolor(jet), clim [50 150], YDir reverse, FOM contour overlay.
%
% Inputs:
%   ax   — axes handle
%   r_km — range vector (km)
%   z_m  — depth vector (m)
%   TL   — Nz x Nr TL matrix (dB)
%   FOM  — (optional) scalar dB threshold; if provided, overlays white contour

axes(ax);
pcolor(r_km, z_m, TL);
shading interp;
set(ax, 'YDir', 'reverse');
colormap(ax, jet);
colorbar(ax);
clim([50 150]);
xlabel(ax, 'Range (km)');
ylabel(ax, 'Depth (m)');

if nargin >= 5 && ~isempty(FOM)
    hold(ax, 'on');
    TL_sm = movmean(movmean(TL, 20, 2), 20, 1);
    contour(r_km, z_m, TL_sm, [FOM FOM], 'w-', 'LineWidth', 1.2);
    hold(ax, 'off');
end
end
