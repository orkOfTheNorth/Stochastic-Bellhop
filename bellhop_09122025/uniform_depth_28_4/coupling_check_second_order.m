%% 4) Advanced Coupling Check: Full Mixed-Derivative Calculation
clear; close all; clc; warning('off');
addpath('../Code/Functions');

%% Base Parameters
sourceFrequency = 10000;    
maxRange = 50000;           
sourceDepth = 4.5;          
sourceHalfBeam = 0;         
sourceTilt = 0;             
fom = 0;                    
svpType = "const_1500";     

%% Nominal Values and Perturbations
z_nom = 52.5;         
c_nom = 1500;         
rho = 1.63;           
alpha = 0.07;         
dz = 0.5;             
dc = 5.0;             

%% Define Master Grid for Interpolation
z_master = linspace(0, 100, 500); 
r_master = linspace(0, maxRange, 1000); 

%% Run the 4 Corners
fprintf('Running 4 Corner Perturbations for Cross-Derivative...\n');

fprintf('1/4: (+dz, +dc)\n');
TL_pp = run_bellhop_instance(z_nom + dz, c_nom + dc, ...
    sourceFrequency, maxRange, sourceDepth, sourceHalfBeam, sourceTilt, svpType, rho, alpha, fom, r_master, z_master);

fprintf('2/4: (+dz, -dc)\n');
TL_pm = run_bellhop_instance(z_nom + dz, c_nom - dc, ...
    sourceFrequency, maxRange, sourceDepth, sourceHalfBeam, sourceTilt, svpType, rho, alpha, fom, r_master, z_master);

fprintf('3/4: (-dz, +dc)\n');
TL_mp = run_bellhop_instance(z_nom - dz, c_nom + dc, ...
    sourceFrequency, maxRange, sourceDepth, sourceHalfBeam, sourceTilt, svpType, rho, alpha, fom, r_master, z_master);

fprintf('4/4: (-dz, -dc)\n');
TL_mm = run_bellhop_instance(z_nom - dz, c_nom - dc, ...
    sourceFrequency, maxRange, sourceDepth, sourceHalfBeam, sourceTilt, svpType, rho, alpha, fom, r_master, z_master);

%% Calculate the Mixed Partial Derivative (Coupling Term)
fprintf('Calculating Cross-Derivative Matrix...\n');
d2TL_dzdc = (TL_pp - TL_pm - TL_mp + TL_mm) ./ (4 * dz * dc);

% כדי לראות את ההשפעה הפיזית האמיתית (בדציבלים) של הצימוד 
% נכפיל את הנגזרת בחזרה בגדלי ההפרעות כדי לקבל את "השגיאה ב-dB"
coupling_impact_dB = d2TL_dzdc * dz * dc;

%% Visualization
fig = figure('Name', 'Advanced Coupling Check (Mixed Derivative)', 'Position', [200, 200, 800, 500]);

imagesc(r_master/1000, z_master, coupling_impact_dB);
colorbar; set(gca, 'YDir', 'reverse'); 
title('Impact of Coupling Term (\partial^2 TL / \partial z \partial c \cdot \Delta z \Delta c)');
xlabel('Range (km)'); ylabel('Depth (m)');
colormap(gca, 'parula');

% Center error colormap around 0
max_err = prctile(abs(coupling_impact_dB(:)), 98); 
if max_err > 0 && ~isnan(max_err)
    clim([-max_err, max_err]); 
end

%% Helper Function
function [TL_out] = run_bellhop_instance(z_bot, c_bot, freq, rmax, zsrc, beam, tilt, svp, rho, alpha, fom, r_mast, z_mast)
    bathy_str = sprintf('const_%.2f', z_bot);
    geo = [c_bot, rho, alpha];
    sim_pars = {freq, rmax, zsrc, beam, tilt, svp, bathy_str, geo, fom};
    
    [TL_raw, r_raw, z_raw] = simpleBellhopHazat(sim_pars);
    
    % Clean data
    TL_raw(TL_raw > 150) = NaN;
    TL_raw(isinf(TL_raw)) = NaN;
    TL_raw(TL_raw == 0) = NaN;
    
    % Interpolate
    [R_raw, Z_raw] = meshgrid(r_raw, z_raw);
    [R_mast, Z_mast] = meshgrid(r_mast, z_mast);
    TL_out = interp2(R_raw, Z_raw, TL_raw, R_mast, Z_mast, 'linear', NaN);
end