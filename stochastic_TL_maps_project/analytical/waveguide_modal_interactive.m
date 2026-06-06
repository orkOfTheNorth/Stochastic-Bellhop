function waveguide_modal_interactive()
% Isovelocity waveguide — normal mode solution, interactive
% D=2500m depth, range=10km, source depth=5m
%
% Modal sum:  p(x,z) = sum_m  [sin(k_zm * zs) * sin(k_zm * z) * exp(-i*k_xm*x)]
%                              / (i * k_xm * D)
%   k_zm = (m+0.5)*pi/D        (Dirichlet surface, Neumann bottom)
%   k_xm = sqrt(k^2 - k_zm^2)  (propagating if real, evanescent if imaginary)
%
% M_limit=2000 captures all propagating modes up to ~500 Hz
% Frequency slider: 10–500 Hz

    f = 15;
    c = 1500;
    zs = 5;
    D = 2500;
    max_range_x = 10000;

    fig = figure('Name', 'Waveguide Modal Solution — Interactive', ...
                 'Position', [100, 100, 1100, 700], 'Color', 'w');
    ax = axes('Parent', fig, 'Position', [0.1, 0.15, 0.75, 0.75]);

    uicontrol('Style', 'text',   'Position', [50, 20, 150, 20], 'String', 'Frequency (Hz)');
    uicontrol('Style', 'slider', 'Min', 10, 'Max', 500, 'Value', f, ...
              'Position', [200, 20, 200, 20], 'Callback', @update_params);

    render_plot();

    function render_plot()
        x_vec = linspace(0.1, max_range_x, 500);
        z_vec = linspace(0, D, 300);
        [X, Z] = meshgrid(x_vec, z_vec);
        k = 2*pi*f/c;
        P = zeros(size(X));
        n_prop = 0;
        for m = 0:2000
            k_zm    = (m + 0.5)*pi/D;
            inside  = k^2 - k_zm^2;
            if inside > 0
                k_xm  = sqrt(inside);
                n_prop = n_prop + 1;
            else
                k_xm = -1i*sqrt(abs(inside));
            end
            term = (1/(1i*k_xm*D)) .* sin(k_zm*zs) .* sin(k_zm*Z) .* exp(-1i*k_xm*X);
            P = P + term;
        end
        TL = -20*log10(abs(P));
        imagesc(ax, x_vec, z_vec, TL);
        colormap(ax, jet);  caxis(ax, [20 100]);
        set(ax, 'YDir', 'reverse');
        cb = colorbar(ax);  ylabel(cb, 'Transmission Loss (dB)');
        xlabel(ax, 'Range (m)');  ylabel(ax, 'Depth (m)');
        title(ax, sprintf('Waveguide Modal  f=%dHz  c=%.0fm/s  D=%dm  %d propagating modes', ...
              floor(f), c, D, n_prop));
        drawnow;
    end

    function update_params(src, ~),  f = src.Value;  render_plot();  end
end
