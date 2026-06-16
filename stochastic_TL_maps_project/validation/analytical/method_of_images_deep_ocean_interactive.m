function method_of_images_deep_ocean_interactive()
% Method of Images (Lloyd's Mirror) — deep ocean interactive
% Scale: depth=2500m, range=10km, source depth=5m
%
% Physics: exact Green's function solution (wave equation, not ray tracing)
%   P = exp(ikR1)/(R1) - exp(ikR2)/(R2)    (3D spherical spreading)
%   TL = -20*log10(|P|)
%
% Controls:
%   Sound speed slider: 1400-1600 m/s

    f = 2000;
    c = 1500;
    zs = 5;
    max_range_x = 10000;
    max_depth_z = 2500;

    fig = figure('Name', 'Method of Images — Deep Ocean Interactive', ...
                 'Position', [100, 100, 1100, 700], 'Color', 'w');
    ax = axes('Parent', fig, 'Position', [0.1, 0.15, 0.75, 0.75]);

    uicontrol('Style', 'text',   'Position', [50, 20, 150, 20], 'String', 'Sound Speed (c)');
    uicontrol('Style', 'slider', 'Min', 1400, 'Max', 1600, 'Value', c, ...
              'Position', [200, 20, 200, 20], 'Callback', @update_speed);

    render_plot();

    function render_plot()
        if ~isvalid(ax), return; end
        k = 2 * pi * f / c;
        x_vec = linspace(1, max_range_x, 800);
        z_vec = linspace(0, max_depth_z, 800);
        [X, Z] = meshgrid(x_vec, z_vec);
        R1 = sqrt(X.^2 + (Z - zs).^2);
        R2 = sqrt(X.^2 + (Z + zs).^2);
        P  = (exp(1i*k*R1)./(R1+eps)) - (exp(1i*k*R2)./(R2+eps));
        TL = -20 * log10(abs(P) + eps);
        imagesc(ax, x_vec/1000, z_vec, TL);
        colormap(ax, jet);  caxis(ax, [60 110]);
        set(ax, 'YDir', 'reverse');  axis(ax, 'tight');
        cb = colorbar(ax);  ylabel(cb, 'Transmission Loss (dB)');
        xlabel(ax, 'Range (km)');  ylabel(ax, 'Depth (m)');
        title(ax, sprintf('Method of Images (Deep Ocean)  f=%dHz  c=%.0fm/s  lambda=%.1fm', ...
              f, c, c/f));
        drawnow;
    end

    function update_speed(src, ~),  c = src.Value;  render_plot();  end
end
