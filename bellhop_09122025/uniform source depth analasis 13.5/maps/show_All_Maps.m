%% show_All_Maps.m
% dashboard to load and display all saved uncertainty quantification results.
clear; close all; clc;

disp('Scanning directory for saved acoustic maps...');

%% 1. DISPLAY SINGLE-VARIABLE (SD) RESULTS
if isfile('MC_SD_Results.mat') && isfile('Delta_SD_Results.mat')
    load('MC_SD_Results.mat');
    load('Delta_SD_Results.mat');
    r_km = r_grid/1000;
    
    figure('Name', 'Dashboard: Source Depth (SD) Uncertainty', 'Position', [50 50 1200 800]);
    
    % Expected Values
    subplot(2,2,1); pcolor(r_km, z_grid, MC_Expected_TL); shading interp; set(gca,'YDir','reverse');
    title('MC Expected TL (SD Only)'); colormap(gca, jet); colorbar; caxis([50 150]);
    
    subplot(2,2,2); pcolor(r_km, z_grid, Delta_Expected_TL); shading interp; set(gca,'YDir','reverse');
    title('Delta Expected TL (SD Only)'); colormap(gca, jet); colorbar; caxis([50 150]);
    
    % Variances
    subplot(2,2,3); pcolor(r_km, z_grid, MC_Var_TL); shading interp; set(gca,'YDir','reverse');
    title('MC Variance (SD Only)'); colormap(gca, hot); colorbar;
    
    subplot(2,2,4); pcolor(r_km, z_grid, Delta_Var_TL); shading interp; set(gca,'YDir','reverse');
    title('Delta Variance (SD Only)'); colormap(gca, hot); colorbar;
else
    disp('No Single-Variable SD results found.');
end

%% 2. DISPLAY MULTI-PARAMETER RESULTS
% Finds all files starting with 'Multi_UQ_Results'
multi_files = dir('Multi_UQ_Results*.mat');

for i = 1:length(multi_files)
    fname = multi_files(i).name;
    load(fname);
    r_km = r_grid/1000;
    
    figure('Name', ['Dashboard: Multi-Variable (', fname, ')'], 'Position', [100 100 1200 800]);
    
    subplot(2,2,1); pcolor(r_km, z_grid, Expected_LHS); shading interp; set(gca,'YDir','reverse');
    title(['LHS Expected (', fname, ')']); colormap(gca, jet); colorbar; caxis([50 150]);
    
    subplot(2,2,2); pcolor(r_km, z_grid, TL_nom); shading interp; set(gca,'YDir','reverse');
    title(['Delta Expected (', fname, ')']); colormap(gca, jet); colorbar; caxis([50 150]);
    
    subplot(2,2,3); pcolor(r_km, z_grid, Var_LHS); shading interp; set(gca,'YDir','reverse');
    title('LHS Variance'); colormap(gca, hot); colorbar;
    
    subplot(2,2,4); pcolor(r_km, z_grid, Var_Delta_Total); shading interp; set(gca,'YDir','reverse');
    title('Delta Total Variance'); colormap(gca, hot); colorbar;
end

%% 3. DISPLAY COMPARISON / DIFFERENCE MAPS
% Finds all files starting with 'Multi_UQ_Comparison' or 'Compare_SD'
comp_files = [dir('Multi_UQ_Comparison*.mat'); dir('Compare_SD_Results.mat')];

for j = 1:length(comp_files)
    cname = comp_files(j).name;
    load(cname);
    
    % Handle different variable naming in SD vs Multi comparison files
    if contains(cname, 'SD')
        % r_km and z_m are already in the SD compare file
    else
        r_km = r_grid/1000;
        z_m = z_grid;
    end
    
    figure('Name', ['Dashboard: Difference Maps (', cname, ')'], 'Position', [150 150 1000 700]);
    
    subplot(2,1,1); pcolor(r_km, z_m, diff_Expected); shading interp; set(gca, 'YDir', 'reverse');
    title(['Expected Value Difference: ', cname]); colormap(gca, jet); colorbar;
    
    subplot(2,1,2); pcolor(r_km, z_m, diff_Var); shading interp; set(gca, 'YDir', 'reverse');
    title(['Variance Difference: ', cname]); colormap(gca, hot); colorbar;
end

disp('All available maps have been loaded and displayed.');