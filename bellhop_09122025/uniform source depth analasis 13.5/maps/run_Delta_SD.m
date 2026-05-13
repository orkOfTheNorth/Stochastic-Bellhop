%% run_Delta_SD.m
clear; close all; clc; warning('off');
addpath('Code/Functions'); % Adjust if needed

%% 1. Parameters
sim_pars = {10000, 50000, 5, 0, 0, "summer", "const_35", [0.989*1500 1.63 0.07], 100};
% Uniform bounds for SD
sd_min = 4.75;
sd_max = 5.25;
var_SD = ((sd_max - sd_min)^2) / 12; % Uniform Variance

h = 0.1; % Perturbation step for derivative

%% 2. Run Nominal (Center)
disp('Running Delta Method: Nominal...');
[TL_nom, r_grid, z_grid] = simpleBellhopHazat(sim_pars);

%% 3. Run Up / Down Perturbations
disp('Running Delta Method: SD + h...');
sim_pars{3} = 5 + h;
[TL_up, ~, ~] = simpleBellhopHazat(sim_pars);

disp('Running Delta Method: SD - h...');
sim_pars{3} = 5 - h;
[TL_down, ~, ~] = simpleBellhopHazat(sim_pars);

%% 4. Calculate Variance via Delta Method
% Central finite difference derivative
dTL_dSD = (TL_up - TL_down) / (2 * h);

% Apply Delta formula
Delta_Var_TL = (dTL_dSD.^2) .* var_SD;
Delta_Expected_TL = TL_nom;

%% 5. Save Data
save('Delta_SD_Results.mat', 'Delta_Expected_TL', 'Delta_Var_TL', 'r_grid', 'z_grid');
disp('Delta Method complete and saved.');