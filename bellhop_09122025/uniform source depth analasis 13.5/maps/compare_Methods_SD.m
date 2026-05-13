%% compare_Methods_SD.m
clear; close all; clc;

% 1. Load Data
load('Delta_SD_Results.mat');
load('MC_SD_Results.mat');

r_km = r_grid / 1000;
z_m = z_grid;

fig = figure('Name', 'Interactive Comparison: MC vs Delta', 'Position', [100 100 900 800]);

%% 2. Expected Value Difference Map (Top Plot)
ax1 = subplot(2,1,1);
diff_Expected = MC_Expected_TL - Delta_Expected_TL;
imagesc(r_km(1,:), z_m(:,1), diff_Expected);
set(gca, 'YDir', 'reverse');
colormap(ax1, jet); colorbar;
title('Difference in Expected TL (MC - Delta) [dB]');
xlabel('Range (km)'); ylabel('Depth (m)');

%% 3. Variance Difference Map (Bottom Plot)
ax2 = subplot(2,1,2);
diff_Var = MC_Var_TL - Delta_Var_TL;
imagesc(r_km(1,:), z_m(:,1), diff_Var);
set(gca, 'YDir', 'reverse');
colormap(ax2, hot); colorbar;
title('Difference in Variance (MC - Delta)');
xlabel('Range (km)'); ylabel('Depth (m)');

%% 4. Set up the Custom Click Function
dcm_obj = datacursormode(fig);
dcm_obj.Enable = 'on';
dcm_obj.DisplayStyle = 'window';

% Pass both axes (ax1, ax2) and all four data matrices to the callback
set(dcm_obj, 'UpdateFcn', @(obj, event_obj) customDatatip(obj, event_obj, ...
    ax1, ax2, r_km, z_m, MC_Expected_TL, Delta_Expected_TL, MC_Var_TL, Delta_Var_TL));

disp('Click on either map to see the underlying values.');

%% 5. The Callback Function (Keep at the bottom of your file)
function txt = customDatatip(~, event_obj, ax1, ax2, r_km, z_m, MC_Exp, Delta_Exp, MC_Var, Delta_Var)
    % Get coordinates of the click
    pos = event_obj.Position;
    
    % Determine which subplot was clicked
    clicked_ax = event_obj.Target.Parent; 
    
    % Find the closest indices in the grid
    [~, r_idx] = min(abs(r_km(1,:) - pos(1)));
    [~, z_idx] = min(abs(z_m(:,1) - pos(2)));
    
    % Format the pop-up text based on the axis clicked
    if clicked_ax == ax1
        % User clicked the Expected TL map
        val_mc    = MC_Exp(z_idx, r_idx);
        val_delta = Delta_Exp(z_idx, r_idx);
        val_diff  = val_mc - val_delta;
        
        txt = {['Range: ', num2str(pos(1), '%.2f'), ' km'], ...
               ['Depth: ', num2str(pos(2), '%.1f'), ' m'], ...
               '--- Expected TL ---', ...
               ['MC:     ', num2str(val_mc, '%.2f'), ' dB'], ...
               ['Delta:  ', num2str(val_delta, '%.2f'), ' dB'], ...
               ['Diff:   ', num2str(val_diff, '%.2f'), ' dB']};
               
    elseif clicked_ax == ax2
        % User clicked the Variance map
        val_mc    = MC_Var(z_idx, r_idx);
        val_delta = Delta_Var(z_idx, r_idx);
        val_diff  = val_mc - val_delta;
        
        txt = {['Range: ', num2str(pos(1), '%.2f'), ' km'], ...
               ['Depth: ', num2str(pos(2), '%.1f'), ' m'], ...
               '--- Variance ---', ...
               ['MC:     ', num2str(val_mc, '%.4f')], ...
               ['Delta:  ', num2str(val_delta, '%.4f')], ...
               ['Diff:   ', num2str(val_diff, '%.4f')]};
    else
        txt = 'Click registered outside data area';
    end
end