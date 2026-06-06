function sea_wave_simulation_3D()
    % --- Initial Parameters ---
    f = 100;            % Frequency (Hz)
    c = 1500;           % Speed of sound (m/s)
    zs = 1;             % Source depth (meters)
    source_type = 'Point'; 
    view_range = 10;    % Default Zoom
    
    % Create Figure
    fig = figure('Name', '3D Physics - 2D Slice Simulation', 'Position', [100, 100, 1100, 700]);
    ax = axes('Parent', fig, 'Position', [0.1, 0.25, 0.75, 0.65]);
    
    % --- UI Controls ---
    uicontrol('Style', 'popup', 'String', {'3D Point Source', 'Line Source (Y-axis)'}, ...
              'Position', [50, 50, 160, 30], 'Callback', @update_source);
    
    sld_c = uicontrol('Style', 'slider', 'Min', 1400, 'Max', 1600, 'Value', 1500, ...
              'Position', [250, 55, 150, 20], 'Callback', @update_speed);
    uicontrol('Style', 'text', 'Position', [250, 30, 150, 20], 'String', 'Sound Speed (c)');
    
    uicontrol('Style', 'checkbox', 'String', 'Far Field (100m)', 'Value', 0, ...
              'Position', [450, 55, 120, 20], 'Callback', @update_zoom);

    render_plot();

    function render_plot()
        k = 2 * pi * f / c;
        res = 500;
        
        % Slice at Y = 0
        x_vec = linspace(0.01, view_range, res);
        z_vec = linspace(0.01, view_range, res);
        [X, Z] = meshgrid(x_vec, z_vec);

       
        % Radius calculation for the Y=0 plane
        % R1 = Real Source, R2 = Image Source (Method of Images)
        R1 = sqrt(X.^2 + (Z - zs).^2); 
        R2 = sqrt(X.^2 + (Z + zs).^2);
        
        if strcmp(source_type, 'Point')
            % --- 3D Spherical Spreading Logic ---
            % P ~ e^(ikR)/R
            P = (exp(1i*k*R1)./R1) - (exp(1i*k*R2)./R2);
        else
            % --- 2D Cylindrical Spreading Logic ---
            % P ~ Hankel function (represents infinite line source)
            P = besselh(0, 1, k*R1) - besselh(0, 1, k*R2);
        end
        
        % Transmission Loss RL = -10 * log10(P^2)
        % Note: Intensity is proportional to |P|^2
        TL = -10 * log10(abs(P).^2 + eps);
        
        % Visualization
        % We use contourf (filled contours) to group the TL values into bands
        num_levels = 60; % This defines how many 'steps' or wavefronts you see
        contourf(ax, x_vec, z_vec, TL, num_levels, 'LineStyle', 'none');
        cb = colorbar(ax);
        ylabel(cb, 'Transmission Loss (dB)');
        
        set(ax, 'YDir', 'reverse'); % Z direction is DOWN
        xlabel('Horizontal Range X (m)');
        ylabel('Depth Z (m)');
        title(sprintf('%s Source: X-Z Slice at Y=0 (c=%.0f m/s)', source_type, c));
        % This limits the colormap to only have the same number of colors as your bands
        colormap(ax, jet(num_levels));
        caxis(ax, [0 50]);
        grid(ax, 'on');
        
    end

    % --- Callbacks ---
    function update_source(src, ~), types = {'Point', 'Line'}; source_type = types{src.Value}; render_plot(); end
    function update_speed(src, ~), c = src.Value; render_plot(); end
    function update_zoom(src, ~), if src.Value == 1, view_range = 100; else, view_range = 10; end; render_plot(); end
end