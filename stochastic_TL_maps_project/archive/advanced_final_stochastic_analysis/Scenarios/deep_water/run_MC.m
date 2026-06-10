﻿%% run_MC.m  –  Monte Carlo UQ | Scenario: Deep Water
% This script sets up the configuration for the Deep Water scenario and
% calls the core analysis function from the library.

clearvars -except N; close all; clc; warning('off');
try, cd(fileparts(mfilename('fullpath'))); catch; end

base_dir = fullfile('..', '..');
addpath(genpath(fullfile(base_dir, 'Functions')));
addpath(base_dir);
addpath(fullfile(base_dir, 'lib')); % Add library path

%% ── SCENARIO ──────────────────────────────────────────────────────────────
config.label = 'Deep Water (2500m flat, 50km)';
config.bathy_type = 'const_2500';
config.maxR = 50000;
config.output_base_dir = '.'; % Save results in current folder
config.base_dir = base_dir;
if ~exist('N','var'), config.N = 50; else, config.N = N; end

fprintf('=== MC | %s | N=%d ===\n', config.label, config.N);

%% --- Run Core Analysis ---
run_mc_analysis(config);
