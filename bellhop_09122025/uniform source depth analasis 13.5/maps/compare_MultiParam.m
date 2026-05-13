%% compare_MultiParam.m
clear; close all; clc;

%% 1. Select the file to compare
% This script looks for the multi-parameter results file
multi_files = dir('Multi_UQ_Results*.mat');

if isempty(multi_files)
    error('No Multi-Parameter result files found. Run run_MultiParam_UQ.m first.');
end

% List files for the user to choose or just take the latest
[~, idx] = max([multi_files.datenum]);
target_file = multi_files(idx).name;
fprintf('Loading and comparing: %s\n', target_file);
load(target_file);

%% 2. Calculate Differences
% Difference in Expected Value (Bias)
% Note: Delta uses the nominal run as the expected value
diff_Expected = Expected_LHS - TL_nom;

% Difference in Variance (The "Linearity Gap")
diff_Var = Var_LHS - Var_Delta_Total;

%% 3. Plotting the Comparison
r_km = r_grid / 1000;
fig = figure('Name', ['Multi-Var Comparison: ', target_file], 'Position', [100 100 1000 800]);

% --- Top Plot: Expected Value Difference ---
subplot(2,1,1);
pcolor(r_km, z_grid, diff_Expected);
shading interp; set(gca, 'YDir', 'reverse');
colormap(gca, jet); colorbar;
title(['Expected Value Difference (LHS Mean - Nominal) [', target_file, ']']);
xlabel('Range (km)'); ylabel('Depth (m)');

% --- Bottom Plot: Variance Difference ---
subplot(2,1,2);
pcolor(r_km, z_grid, diff_Var);
shading interp; set(gca, 'YDir', 'reverse');
colormap(gca, hot); colorbar;
title('Variance Difference (LHS Variance - Delta Variance)');
xlabel('Range (km)'); ylabel('Depth (m)');

%% 4. Save the Comparison Data
save_name = strrep(target_file, 'Results', 'Comparison');
save(save_name, 'diff_Expected', 'diff_Var', 'r_grid', 'z_grid');
saveas(fig, strrep(save_name, '.mat', '.png'));
fprintf('Comparison saved as %s and .png\n', save_name);
% Enable the data cursor and set a custom update function
dcm = datacursormode(fig);
set(dcm, 'UpdateFcn', @customDataTip);
%%clicking
% This function goes at the very end of your file
function txt = customDataTip(~, event_obj)
    pos = get(event_obj, 'Position');
    % Get the color value (TL difference) at the clicked coordinate
    cdata = get(event_obj.Target, 'CData');
    xdata = get(event_obj.Target, 'XData');
    ydata = get(event_obj.Target, 'YData');
    
    % Find the closest index in the grid
    [~, ix] = min(abs(xdata - pos(1)));
    [~, iy] = min(abs(ydata - pos(2)));
    val = cdata(iy, ix);
    
    txt = {['Range: ', num2str(pos(1), 4), ' km'], ...
           ['Depth: ', num2str(pos(2), 4), ' m'], ...
           ['Diff: ', num2str(val, 4), ' dB']};
end