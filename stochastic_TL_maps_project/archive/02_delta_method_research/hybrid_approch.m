%% Initialization
clear; close all; clc; warning('off');
addpath('Code/Functions');

%% Base Parameters (Nominal Values)
sourceFrequency = 10000; % Hz
maxRange        = 10000; % m
sourceDepth     = 5;     % m
sourceHalfBeam  = 0;     % deg
sourceTilt      = 0;     % deg
svpType         = "const_1500"; 
bathyBase       = "const_2500";
geoAcousticsBase= [1500, 1.63, 0.07]; 
fom             = 100;

%% Define Uncertainties (5% Standard Deviation)
sig_f = sourceFrequency * 0.05;
sig_z = sourceDepth * 0.05;
sig_c = geoAcousticsBase(1) * 0.05;

%% 1. Baseline Run (To get dimensions)
sim_pars_base = {sourceFrequency, maxRange, sourceDepth, sourceHalfBeam, ...
                 sourceTilt, svpType, bathyBase, geoAcousticsBase, fom};
TL_base = simpleBellhopHazat(sim_pars_base);

[num_depths, num_ranges] = size(TL_base);
ranges = linspace(0, maxRange, num_ranges);
depths = linspace(0, 2500, num_depths);

%% 2. The Linear Component (C_bot)
fprintf('Calculating Linear Component for C_bot...\n');
% Run once with a 1-sigma shift (+5%)
geo_nudge = geoAcousticsBase; 
geo_nudge(1) = geoAcousticsBase(1) + sig_c;
sim_pars_c = sim_pars_base; 
sim_pars_c{8} = geo_nudge;

TL_c_nudge = simpleBellhopHazat(sim_pars_c);
dTL_c = TL_c_nudge - TL_base;

% Variance is the square of the difference (Delta^2)
Var_Linear = dTL_c .^ 2;

%% 3. The Stochastic Component (Reduced Monte Carlo for Freq & Depth)
N_runs = 40; % מספר ריצות מצומצם מספיק לסטטיסטיקה של 2 משתנים
fprintf('Running Reduced Monte Carlo (%d iterations) for Freq & Depth...\n', N_runs);

TL_MC_runs = zeros(num_depths, num_ranges, N_runs);

% Generate Gaussian distributions for Freq and Depth
f_dist = sourceFrequency + sig_f * randn(N_runs, 1);
z_dist = sourceDepth + sig_z * randn(N_runs, 1);

for i = 1:N_runs
    % Note: C_bot remains strictly at its nominal baseline value here!
    sim_pars_mc = {f_dist(i), maxRange, z_dist(i), sourceHalfBeam, ...
                   sourceTilt, svpType, bathyBase, geoAcousticsBase, fom};
               
    TL_MC_runs(:,:,i) = simpleBellhopHazat(sim_pars_mc);
    
    if mod(i, 10) == 0
        fprintf('Completed %d / %d runs...\n', i, N_runs);
    end
end

% Variance from the stochastic runs (Std^2)
Std_Stochastic = std(TL_MC_runs, 0, 3);
Var_Stochastic = Std_Stochastic .^ 2;

%% 4. The Hybrid Combination (Root Sum Square)
fprintf('Calculating Hybrid Total Uncertainty...\n');
Total_Variance = Var_Stochastic + Var_Linear;
Total_Std = sqrt(Total_Variance);

%% 5. Visualization
fig = figure('Name', 'Hybrid Uncertainty Analysis', 'Color', 'w', 'Position', [100, 100, 1500, 450]);

% Determine uniform color limits for comparison
max_val = std(Total_Std(:), 'omitnan') * 3;
clim_vals = [0, max_val]; % Uncertainty is absolute (positive), so starts at 0

% --- Plot 1: Linear Uncertainty (from C_bot) ---
ax1 = subplot(1,3,1);
imagesc(ranges, depths, abs(dTL_c)); % Absolute difference represents 1 sigma
set(gca, 'YDir', 'reverse');
colormap(ax1, hot); clim(ax1, clim_vals);
title('Linear Uncertainty: 1\sigma (C_{bot})');
xlabel('Range (m)'); ylabel('Depth (m)');
colorbar;

% --- Plot 2: Stochastic Uncertainty (from Freq + Depth) ---
ax2 = subplot(1,3,2);
imagesc(ranges, depths, Std_Stochastic);
set(gca, 'YDir', 'reverse');
colormap(ax2, hot); clim(ax2, clim_vals);
title('Stochastic Uncertainty: 1\sigma (Freq + Depth)');
xlabel('Range (m)'); 
colorbar;

% --- Plot 3: Hybrid Total Uncertainty ---
ax3 = subplot(1,3,3);
imagesc(ranges, depths, Total_Std);
set(gca, 'YDir', 'reverse');
colormap(ax3, hot); clim(ax3, clim_vals);
title('Total Hybrid Uncertainty: \sigma_{Total} = \surd(\sigma_{MC}^2 + \sigma_{Lin}^2)');
xlabel('Range (m)'); 
colorbar;

linkaxes([ax1, ax2, ax3], 'xy');

% Add Interactive Tooltip
dcm = datacursormode(fig);
dcm.Enable = 'on';
dcm.DisplayStyle = 'datatip';
dcm.UpdateFcn = @(obj, event_obj) tooltip_hybrid(event_obj, ranges, depths, TL_base, abs(dTL_c), Std_Stochastic, Total_Std);

fprintf('החישוב הושלם בהצלחה.\n');

%% --- Custom Tooltip Function ---
function txt = tooltip_hybrid(event_obj, ranges, depths, TL_base, std_lin, std_mc, std_tot)
    pos = event_obj.Position;
    
    [~, r_idx] = min(abs(ranges - pos(1)));
    [~, d_idx] = min(abs(depths - pos(2)));
    
    base = TL_base(d_idx, r_idx);
    s_lin = std_lin(d_idx, r_idx);
    s_mc = std_mc(d_idx, r_idx);
    s_tot = std_tot(d_idx, r_idx);
    
    txt = {
        sprintf('--- Location ---'), ...
        sprintf('Range: %.0f m | Depth: %.0f m', pos(1), pos(2)), ...
        sprintf('Nominal TL: %.1f dB', base), ...
        sprintf(''), ...
        sprintf('--- 1\\sigma Uncertainty ---'), ...
        sprintf('Linear (C_{bot}): \\pm%.2f dB', s_lin), ...
        sprintf('Non-Linear (Freq+Z): \\pm%.2f dB', s_mc), ...
        sprintf('-------------------'), ...
        sprintf('TOTAL UNCERTAINTY: \\pm%.2f dB', s_tot)
    };
end