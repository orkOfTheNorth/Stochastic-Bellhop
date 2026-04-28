%% 4) Coupling Check: First-Order Taylor Expansion vs. Combined Perturbation
clear; close all; clc; warning('off');
addpath('Code/Functions');

%% Base Parameters
sourceFrequency = 10000;    
maxRange = 50000;           
sourceDepth = 4.5;          
sourceHalfBeam = 0;         
sourceTilt = 0;             
fom = 0;                    

% Fixed SVP profile for this test
svpType = "const_1500";     

%% Nominal Values and Perturbations
z_nom = 52.5;         % Nominal depth
c_nom = 1500;         % Nominal bottom sound speed
rho = 1.63;           % Bottom density
alpha = 0.07;         % Bottom attenuation

dz = 0.5;             % Delta z
dc = 5.0;             % Delta c_bottom (m/s)

%% 1. Nominal Run
fprintf('Running Nominal Model (z=%.1f, c=%.1f)...\n', z_nom, c_nom);
[TL_nom, r_grid, z_grid] = run_bellhop_instance(z_nom, c_nom, ...
    sourceFrequency, maxRange, sourceDepth, sourceHalfBeam, sourceTilt, svpType, rho, alpha, fom);

%% 2. Partial Derivative w.r.t Depth (z)
fprintf('Calculating dTL/dz...\n');
[TL_z_plus, ~, ~]  = run_bellhop_instance(z_nom + dz, c_nom, ...
    sourceFrequency, maxRange, sourceDepth, sourceHalfBeam, sourceTilt, svpType, rho, alpha, fom);
[TL_z_minus, ~, ~] = run_bellhop_instance(z_nom - dz, c_nom, ...
    sourceFrequency, maxRange, sourceDepth, sourceHalfBeam, sourceTilt, svpType, rho, alpha, fom);

dTL_dz = (TL_z_plus - TL_z_minus) ./ (2 * dz);

%% 3. Partial Derivative w.r.t Bottom Sound Speed (c)
fprintf('Calculating dTL/dc...\n');
[TL_c_plus, ~, ~]  = run_bellhop_instance(z_nom, c_nom + dc, ...
    sourceFrequency, maxRange, sourceDepth, sourceHalfBeam, sourceTilt, svpType, rho, alpha, fom);
[TL_c_minus, ~, ~] = run_bellhop_instance(z_nom, c_nom - dc, ...
    sourceFrequency, maxRange, sourceDepth, sourceHalfBeam, sourceTilt, svpType, rho, alpha, fom);

dTL_dc = (TL_c_plus - TL_c_minus) ./ (2 * dc);

%% 4. First-Order Taylor Approximation
% Calculate expected delta TL if both change by +dz and +dc
delta_TL_taylor = dTL_dz .* dz + dTL_dc .* dc;

%% 5. Actual Combined Perturbation
fprintf('Running Combined Perturbation (z+dz, c+dc)...\n');
[TL_combined, ~, ~] = run_bellhop_instance(z_nom + dz, c_nom + dc, ...
    sourceFrequency, maxRange, sourceDepth, sourceHalfBeam, sourceTilt, svpType, rho, alpha, fom);

% Actual change from nominal
delta_TL_actual = TL_combined - TL_nom;

%% 6. Difference / Coupling Error
coupling_error = delta_TL_actual - delta_TL_taylor;

%% Visualization
fig = figure('Name', 'Coupling Check: Taylor vs Combined', 'Position', [100, 100, 1500, 450]);

subplot(1,3,1);
imagesc(r_grid/1000, z_grid, delta_TL_taylor);
colorbar; set(gca, 'YDir', 'reverse');
title('\Delta TL (Taylor 1st Order)');
xlabel('Range (km)'); ylabel('Depth (m)');
colormap('jet');

subplot(1,3,2);
imagesc(r_grid/1000, z_grid, delta_TL_actual);
colorbar; set(gca, 'YDir', 'reverse');
title('\Delta TL (Actual Combined)');
xlabel('Range (km)'); ylabel('Depth (m)');
colormap('jet');

subplot(1,3,3);
imagesc(r_grid/1000, z_grid, coupling_error);
colorbar; set(gca, 'YDir', 'reverse');
title('Coupling Error (Actual - Taylor)');
xlabel('Range (km)'); ylabel('Depth (m)');
% Center error colormap around 0
max_err = max(abs(coupling_error(:)));
if max_err > 0, clim([-max_err, max_err]); end
colormap(gca, 'parula');

%% Helper Function for Cleaner Calls
function [TL, r, z] = run_bellhop_instance(z_bot, c_bot, freq, rmax, zsrc, beam, tilt, svp, rho, alpha, fom)
    bathy_str = sprintf('const_%.2f', z_bot);
    geo = [c_bot, rho, alpha];
    sim_pars = {freq, rmax, zsrc, beam, tilt, svp, bathy_str, geo, fom};
    [TL, r, z] = simpleBellhopHazat(sim_pars);
end