%% Initialization
clear; close all; clc; warning('off');
addpath('Code/Functions');

%% User Parameters
sourceFrequency = 10000; % Hz, > 3000 Hz
maxRange = 50000; % m
sourceDepth = 5; % m
sourceHalfBeam = 0; % deg
sourceTilt = 0; % deg
svpType = "summer"; % string, has to be 
        % "summer", "winter", "const_<sv[m/s]>" or "custom"
bathymetryType = "const_35"; % string, has to be "const_<maxDepth[m]>", 
        % "slope_<startDepth[m]>_<endDepth[m]>" or "custom"
geoAcoustics = [0.989*1500 1.63 0.07]; % [c_bottom rho[g/cm^3] alpha[dB/lambda]]
fom = 100; % dB

%% Run
sim_pars = {sourceFrequency, maxRange, sourceDepth, sourceHalfBeam, ...
    sourceTilt, svpType, bathymetryType, geoAcoustics, fom};
simpleBellhopHazat(sim_pars)
% simpleBellhopHazat(sim_pars, "CustomSVP", svp, "CustomBathythermy", bm)