%% run_MC_SD.m
%% Initialization
clear; close all; clc; warning('off');

% Go up one level from 'maps' and then into 'code/Functions'
addpath('../code/Functions'); 

%% 1. Parameters
sim_pars = {10000, 50000, 5, 0, 0, "summer", "const_100", [0.989*1500 1.63 0.07], 100};
N = 100;
sd_min = 4.75;
sd_max = 5.25;

% Generate N Latin Hypercube samples for Source Depth
% This ensures perfectly stratified sampling across the [4.75, 5.25] range
sd_samples = sd_min + (sd_max - sd_min) * lhsdesign(N, 1);

%% 2. Monte Carlo Loop
disp(['Starting LHS Monte Carlo with N=', num2str(N)]);
% Pre-allocate based on a dummy run to get dimensions
[TL_dummy, r_grid, z_grid] = simpleBellhopHazat(sim_pars);
TL_all = zeros(size(TL_dummy,1), size(TL_dummy,2), N);

for i = 1:N
    disp(['LHS MC Run ', num2str(i), ' / ', num2str(N)]);
    sim_pars{3} = sd_samples(i);
    [TL, ~, ~] = simpleBellhopHazat(sim_pars);
    TL_all(:,:,i) = TL;
end

%% 3. Calculate Statistics
MC_Expected_TL = mean(TL_all, 3);
MC_Var_TL = var(TL_all, 0, 3);

%% 4. Save Data
save('MC_SD_Results.mat', 'MC_Expected_TL', 'MC_Var_TL', 'r_grid', 'z_grid');
disp('LHS Monte Carlo complete and saved.');