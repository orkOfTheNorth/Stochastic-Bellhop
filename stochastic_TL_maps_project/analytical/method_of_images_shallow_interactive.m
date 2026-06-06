function method_of_images_shallow_interactive()
% Method of Images (Lloyd's Mirror) — shallow near-field interactive
% Scale: depth=10m, range=10m, source depth=1m
%
% Physics: P = G(r,rs) - G(r,rs_image)
%   Point:  G ~ exp(ikR)/R   (3D spherical spreading)
%   Line:   G ~ H0(kR)        (2D cylindrical, Hankel function)
%
% Controls:
%   Source type popup: Point / Line
%   Sound speed slider: 1400-1600 m/s
%   Far field checkbox: extends range to 100m

    f = 100;
    c = 1500;
    zs = 1;
    source_type = 'Point';
    view_range = 10;

    fig = figure('Name', 'Method of Images — Shallow Interactive', ...
                 'Position', [100, 100, 1100, 700]);
    ax = axes('Parent', fig, 'Position', [0.1, 0.25, 0.75, 0.65]);

    uicontrol('Style', 'popup', 'String', {'3D Point Source', 'Line Source (Y-axis)'}, ...
              'Position', [50, 50, 160, 30], 'Callback', @update_source);
    uicontrol('Style', 'slider', 'Min', 1400, 'Max', 1600, 'Value', c, ...
              'Position', [250, 55, 150, 20], 'Callback', @update_speed);
    uicontrol('Style', 'text', 'Position', [250, 30, 150, 20], 'String', 'Sound Speed (c)');
    uicontrol('Style', 'checkbox', 'String', 'Far Field (100m)', 'Value', 0, ...
              'Position', [450, 55, 120, 20], 'Callback', @update_zoom);

    render_plot();

    function render_plot()
        k = 2 * pi * f / c;
        res = 500;
        x_vec = linspace(0.01, view_range, res);
        z_vec = linspace(0.01, view_range, res);
        [X, Z] = meshgrid(x_vec, z_vec);
        R1 = sqrt(X.^2 + (Z - zs).^2);
        R2 = sqrt(X.^2 + (Z + zs).^2);
        if strcmp(source_type, 'Point')
            P = (exp(1i*k*R1)./R1) - (exp(1i*k*R2)./R2);
        else
            P = besselh(0, 1, k*R1) - besselh(0, 1, k*R2);
        end
        TL = -10 * log10(abs(P).^2 + eps);
        num_levels = 60;
        contourf(ax, x_vec, z_vec, TL, num_levels, 'LineStyle', 'none');
        cb = colorbar(ax);  ylabel(cb, 'Transmission Loss (dB)');
        set(ax, 'YDir', 'reverse');
        xlabel(ax, 'Range X (m)');  ylabel(ax, 'Depth Z (m)');
        title(ax, sprintf('%s Source  f=%dHz  c=%.0fm/s  zs=%.1fm  lambda=%.2fm', ...
              source_type, f, c, zs, c/f));
        colormap(ax, jet(num_levels));  caxis(ax, [0 50]);
        grid(ax, 'on');
    end

    function update_source(src, ~)
        types = {'Point', 'Line'};  source_type = types{src.Value};  render_plot();
    end
    function update_speed(src, ~),  c = src.Value;  render_plot();  end
    function update_zoom(src, ~)
        if src.Value == 1, view_range = 100; else, view_range = 10; end
        render_plot();
    end
end
