function sea_wave_simulation_modal()
    % --- פרמטרים התחלתיים (מבוסס על התמונה והקוד שלך) ---
    f = 15;            % תדר נמוך יותר כדי לראות את המודים בבירור
    c = 1500;          % מהירות הקול במים
    zs = 5;           % עומק המקור (z')
    D = 2500;           % עומק המים (D בתמונה)
    
    % גבולות תצוגה
    max_range_x = 5000; % 1 ק"מ
    max_depth_z = D;    
    
    % יצירת ממשק
    fig = figure('Name', 'Waveguide Modal Solution', 'Position', [100, 100, 1100, 700], 'Color', 'w');
    ax = axes('Parent', fig, 'Position', [0.1, 0.15, 0.75, 0.75]);
    
    % סליידר לשינוי תדר (משפיע ישירות על כמות המודים המתקדמים)
    uicontrol('Style', 'text', 'Position', [50, 20, 150, 20], 'String', 'Frequency (f)');
    h_freq = uicontrol('Style', 'slider', 'Min', 10, 'Max', 500, 'Value', f, ...
              'Position', [200, 20, 200, 20], 'Callback', @update_params);

    render_plot();

    function render_plot()
        % 1. הגדרת רשת
        res_x = 500; 
        res_z = 300;
        x_vec = linspace(0.1, max_range_x, res_x); % מתחיל ב-0.1 למניעת סינגולריות
        z_vec = linspace(0, D, res_z);
        [X, Z] = meshgrid(x_vec, z_vec);

        k = 2 * pi * f / c; % מספר גל כולל
        P = zeros(size(X)); % אתחול שדה הלחץ
        
        % 2. חישוב סכום המודים (טור סופי)
        % מספר המודים המקסימלי לחישוב (כולל מודים דועכים)
        M_limit = 50; 
        
        propagating_count = 0;
        
        for m = 0:M_limit
            % k_zm לפי הנוסחה בתמונה
            k_zm = (m + 0.5) * pi / D;
            
            % חישוב k_xm (מספר הגל בכיוון ההתקדמות)
            % אם k^2 > k_zm^2 המוד מתקדם. אחרת הוא דועך (Cutoff)
            inside_sqrt = k^2 - k_zm^2;
            
            if inside_sqrt > 0
                k_xm = sqrt(inside_sqrt);
                propagating_count = propagating_count + 1;
            else
                % מוד דועך (Evanescent) - k_xm הופך למדומה טהור
                % אנחנו משתמשים במינוס j כדי לקבל דעיכה e^(-|k_xm|*x)
                k_xm = -1i * sqrt(abs(inside_sqrt));
            end
            
            % הנוסחה מהתמונה: (1 / (j * kxm * D)) * sin(kzm * zs) * sin(kzm * z) * exp(-j * kxm * x)
            % שים לב: המשתנה x בקוד הוא X (מרחק מהמקור)
            mode_z = sin(k_zm * zs) * sin(k_zm * Z);
            mode_x = exp(-1i * k_xm * X);
            
            term = (1 / (1i * k_xm * D)) .* mode_z .* mode_x;
            P = P + term;
        end
        
        % 3. חישוב Transmission Loss
        TL = -20 * log10(abs(P) );
        
        % 4. ויזואליזציה
        imagesc(ax, x_vec, z_vec, TL); 
        colormap(ax, jet);
        caxis(ax, [20 70]); % התאמת ניגודיות
        
        set(ax, 'YDir', 'reverse');
        cb = colorbar(ax);
        ylabel(cb, 'Transmission Loss (dB)');
        xlabel(ax, 'Range (m)');
        ylabel(ax, 'Depth (m)');
        title(ax, sprintf('Modal Solution: f=%dHz, Propagating Modes: %d', floor(f), propagating_count));
        drawnow;
    end

    function update_params(src, ~)
        f = src.Value; 
        render_plot(); 
    end
end