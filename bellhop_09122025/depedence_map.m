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

%% Define the 5% Nudges
f_nudge = sourceFrequency * 1.05;
z_nudge = sourceDepth * 1.05;
c_nudge = geoAcousticsBase(1) * 1.05;

%% 1. Baseline Run
sim_pars_base = {sourceFrequency, maxRange, sourceDepth, sourceHalfBeam, ...
                 sourceTilt, svpType, bathyBase, geoAcousticsBase, fom};
TL_base = simpleBellhopHazat(sim_pars_base);

% Extract grid sizes for plotting
[num_depths, num_ranges] = size(TL_base);
ranges = linspace(0, maxRange, num_ranges);
depths = linspace(0, 2500, num_depths);

%% 2. Isolated Runs (Single Variable Changes)
% --- Freq Only ---
sim_pars_f = sim_pars_base; sim_pars_f{1} = f_nudge;
dTL_f = simpleBellhopHazat(sim_pars_f) - TL_base;

% --- Source Depth Only ---
sim_pars_z = sim_pars_base; sim_pars_z{3} = z_nudge;
dTL_z = simpleBellhopHazat(sim_pars_z) - TL_base;

% --- C_bot Only ---
geo_nudge = geoAcousticsBase; geo_nudge(1) = c_nudge;
sim_pars_c = sim_pars_base; sim_pars_c{8} = geo_nudge;
dTL_c = simpleBellhopHazat(sim_pars_c) - TL_base;

%% 3. Actual Joint Runs (Two Variables Changing Simultaneously)
% --- Freq + Source Depth ---
sim_pars_fz = sim_pars_base; sim_pars_fz{1} = f_nudge; sim_pars_fz{3} = z_nudge;
dTL_fz = simpleBellhopHazat(sim_pars_fz) - TL_base;

% --- Freq + C_bot ---
sim_pars_fc = sim_pars_base; sim_pars_fc{1} = f_nudge; sim_pars_fc{8} = geo_nudge;
dTL_fc = simpleBellhopHazat(sim_pars_fc) - TL_base;

% --- Source Depth + C_bot ---
sim_pars_zc = sim_pars_base; sim_pars_zc{3} = z_nudge; sim_pars_zc{8} = geo_nudge;
dTL_zc = simpleBellhopHazat(sim_pars_zc) - TL_base;

%% 4. Calculate Linearity Error (Error = Joint - Linear Sum)
Err_fz = dTL_fz - (dTL_f + dTL_z);
Err_fc = dTL_fc - (dTL_f + dTL_c);
Err_zc = dTL_zc - (dTL_z + dTL_c);

%% 5. Visualization - Pairwise Linearity Error Comparison
fig = figure('Name', 'Pairwise Linearity Error', 'Color', 'w', 'Position', [100, 100, 1500, 450]);

% Determine uniform color limits based on the standard deviation of all errors combined
all_errors = [Err_fz(:); Err_fc(:); Err_zc(:)];
max_err = std(all_errors, 'omitnan') * 2.5; 
clim_vals = [-max_err, max_err];

% --- Plot 1: Freq + Source Depth ---
ax1 = subplot(1,3,1);
imagesc(ranges, depths, Err_fz);
set(gca, 'YDir', 'reverse');
colormap(ax1, parula); clim(ax1, clim_vals);
title('Error: Freq + Source Depth');
xlabel('Range (m)'); ylabel('Depth (m)');
colorbar;

% --- Plot 2: Freq + C_bot ---
ax2 = subplot(1,3,2);
imagesc(ranges, depths, Err_fc);
set(gca, 'YDir', 'reverse');
colormap(ax2, parula); clim(ax2, clim_vals);
title('Error: Freq + C_{bot}');
xlabel('Range (m)'); 
colorbar;

% --- Plot 3: Source Depth + C_bot ---
ax3 = subplot(1,3,3);
imagesc(ranges, depths, Err_zc);
set(gca, 'YDir', 'reverse');
colormap(ax3, parula); clim(ax3, clim_vals);
title('Error: Source Depth + C_{bot}');
xlabel('Range (m)'); 
colorbar;

% Link the axes so zooming in one plot zooms all of them
linkaxes([ax1, ax2, ax3], 'xy');

fprintf('החישוב הושלם. תוכל כעת להשוות את חומרת האינטראקציה בין הזוגות.\n');