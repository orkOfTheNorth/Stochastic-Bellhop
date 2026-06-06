%% Initialization
clear; close all; clc; warning('off');
addpath('Code/Functions');

%% 1. Baseline Parameters
sourceFrequency = 15000; % Hz
maxRange = 5000; % m
sourceDepth = 5; % m
sourceHalfBeam = 0; % deg
sourceTilt = 0; % deg
svpType = "const_1500"; % string
geoAcoustics = [1600 1.63 0.07]; % [c_bottom rho alpha]
fom = 100; % dB

% Baseline Bathymetry (Ocean Depth)
bottom_depth_base = 2500; % Baseline depth in meters
bathymetryType_base = "const_" + num2str(bottom_depth_base);

% Pack parameters into the cell array
sim_pars_base = {sourceFrequency, maxRange, sourceDepth, sourceHalfBeam, ...
    sourceTilt, svpType, bathymetryType_base, geoAcoustics, fom};

%% 2. Run Baseline Simulation
disp('Running Baseline Simulation...');
[TL_base, r_grid, z_grid] = simpleBellhopHazat(sim_pars_base);

% Convert Transmission Loss (dB) to linear Incoherent Intensity
I_base = 10.^(-TL_base / 10);

%% 3. Perturb the Bottom Depth (delta H)
delta_H = 1.0; % Increase the ocean depth by 1 meter
bottom_depth_pert = bottom_depth_base + delta_H;
bathymetryType_pert = "const_" + num2str(bottom_depth_pert);

% Pack perturbed parameters
sim_pars_pert = {sourceFrequency, maxRange, sourceDepth, sourceHalfBeam, ...
    sourceTilt, svpType, bathymetryType_pert, geoAcoustics, fom};

%% 4. Run Perturbed Simulation
disp('Running Perturbed Simulation...');
% We only need the TL matrix this time, so we ignore the grid outputs
[TL_pert, ~, ~] = simpleBellhopHazat(sim_pars_pert);

% Convert perturbed TL to linear intensity
I_pert = 10.^(-TL_pert / 10);

%% 5. Calculate the Derivative Map (dI / dH)
% Finite difference formula: f'(x) = (f(x + dx) - f(x)) / dx
dI_dH = (I_pert - I_base) / delta_H;

%% 6. Plot the Sensitivity Map
% Convert the absolute derivative back to a dB scale for visualization
dI_dH_dB = 10 * log10(abs(dI_dH) + 1e-20); % 1e-20 prevents log(0) errors

figure('Name', 'Sensitivity Analysis', 'Position', [100, 100, 800, 500]);
imagesc(r_grid, z_grid, dI_dH_dB);
colorbar;
colormap('jet');
set(gca, 'YDir', 'reverse'); % Ensure depth increases downwards
xlabel('Range (km)'); 
ylabel('Depth (m)');
title('Sensitivity of Acoustic Intensity to Ocean Bottom Depth (dI/dH)');