%% compare_Methods_SD.m
clear; close all; clc;

% Load Data
load('Delta_SD_Results.mat');
load('MC_SD_Results.mat');

% Convert range/depth to km for plotting
r_km = r_grid / 1000;
z_m = z_grid;

figure('Name', 'Delta vs Monte Carlo (SD Only)', 'Position', [100 100 900 800]);

% 1. Expected Value Difference
subplot(2,1,1);
diff_Expected = MC_Expected_TL - Delta_Expected_TL;
pcolor(r_km, z_m, diff_Expected);
shading interp; set(gca, 'YDir', 'reverse');
colormap(gca, jet); colorbar;
title('Difference in Expected TL (MC - Delta) [dB]');
xlabel('Range (km)'); ylabel('Depth (m)');

% 2. Variance Difference
subplot(2,1,2);
diff_Var = MC_Var_TL - Delta_Var_TL;
pcolor(r_km, z_m, diff_Var);
shading interp; set(gca, 'YDir', 'reverse');
colormap(gca, hot); colorbar;
title('Difference in Variance (MC - Delta)');
xlabel('Range (km)'); ylabel('Depth (m)');