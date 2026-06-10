%% run_MC.m  –  Monte Carlo UQ  (2026-05-26 | Baseline scenario)
%
% This script sets up the configuration for the baseline scenario and
% calls the core analysis function from the library.

clearvars -except N; close all; clc;
try, cd(fileparts(mfilename('fullpath'))); catch; end

base_dir = fullfile('..');
addpath(genpath(fullfile(base_dir, 'Functions')));
addpath(base_dir);
addpath(fullfile(base_dir, 'lib')); % Add library path

%% --- Configuration for Baseline Scenario ---
config.label = 'Baseline (35m flat, 50km)';
config.bathy_type = 'const_35';
config.maxR = 50000;
config.output_base_dir = '.'; % Save results in current folder (Monte_Carlo/)
config.base_dir = base_dir;
if ~exist('N','var'), config.N = 50; else, config.N = N; end

fprintf('=== Monte Carlo | %s | N=%d ===\n', config.label, config.N);

%% --- Run Core Analysis ---
run_mc_analysis(config);
