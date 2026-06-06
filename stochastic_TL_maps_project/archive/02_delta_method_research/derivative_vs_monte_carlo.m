%% 3) Comparison: Derivative (Delta) Method vs. Monte Carlo
clear; close all; clc;

%% Load Results
try
    mc = load('monte_carlo_results.mat');
    delta = load('delta_method_results.mat');
catch
    error('Could not find the .mat files. Ensure scripts 1 and 2 have been run.');
end

r_grid = mc.r_grid;
z_grid = mc.z_grid;

%% Calculate Differences
% Difference in Expected Value
Diff_E_TL = mc.E_TL_mc - delta.E_TL_delta;

% Difference in Variance
Diff_Var_TL = mc.Var_TL_mc - delta.Var_TL_delta;

%% Visualization
fig = figure('Name', 'Method Comparison: MC vs Delta', 'Position', [200, 200, 1200, 500]);

% Plot Expected Value Difference
subplot(1,2,1);
imagesc(r_grid/1000, z_grid, Diff_E_TL);
colorbar; set(gca, 'YDir', 'reverse');
title('Difference in E[TL] (MC - Delta)');
xlabel('Range (km)'); ylabel('Depth (m)');
% Center colormap around 0
max_val_E = max(abs(Diff_E_TL(:)));
if max_val_E > 0, clim([-max_val_E, max_val_E]); end
colormap(gca, 'parula'); 

% Plot Variance Difference
subplot(1,2,2);
imagesc(r_grid/1000, z_grid, Diff_Var_TL);
colorbar; set(gca, 'YDir', 'reverse');
title('Difference in Var(TL) (MC - Delta)');
xlabel('Range (km)'); ylabel('Depth (m)');
% Center colormap around 0
max_val_V = max(abs(Diff_Var_TL(:)));
if max_val_V > 0, clim([-max_val_V, max_val_V]); end
colormap(gca, 'parula');

% Interactive Data Cursor
dcm_obj = datacursormode(fig);
set(dcm_obj, 'UpdateFcn', {@customDatatipDiff, r_grid, z_grid, Diff_E_TL, Diff_Var_TL});

%% Custom Datatip Function
function txt = customDatatipDiff(~, event_obj, r_grid, z_grid, Diff_E, Diff_Var)
    pos = get(event_obj, 'Position');
    clicked_r = pos(1) * 1000; 
    clicked_z = pos(2);
    
    [~, r_idx] = min(abs(r_grid - clicked_r));
    [~, z_idx] = min(abs(z_grid - clicked_z));
    
    txt = {sprintf('Range: %.2f km', pos(1)), ...
           sprintf('Depth: %.1f m', pos(2)), ...
           sprintf('Δ E[TL]: %.2f dB', Diff_E(z_idx, r_idx)), ...
           sprintf('Δ Var(TL): %.2f dB^2', Diff_Var(z_idx, r_idx))};
end