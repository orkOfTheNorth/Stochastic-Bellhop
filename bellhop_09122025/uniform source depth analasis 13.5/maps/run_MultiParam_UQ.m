%% Initialization
clear; close all; clc; warning('off');

% Go up one level from 'maps' and then into 'code/Functions'
addpath('../code/Functions'); 


% --- Create UI for Selection ---
fig = uifigure('Name', 'Select Variables to Perturb', 'Position', [500 500 300 250]);
lbl = uilabel(fig, 'Position', [20 200 250 22], 'Text', 'Select parameters for UQ Analysis:');
cb_SD = uicheckbox(fig, 'Position', [40 160 200 22], 'Text', 'Source Depth (5m ± 5%)', 'Value', 1);
cb_Bathy = uicheckbox(fig, 'Position', [40 120 200 22], 'Text', 'Bathy Slope (0° ± 0.1°)', 'Value', 1);
cb_SVP = uicheckbox(fig, 'Position', [40 80 200 22], 'Text', 'SVP Summer (± 1°C)', 'Value', 1);
btn = uibutton(fig, 'push', 'Position', [100 20 100 30], 'Text', 'Run Suite', ...
    'ButtonPushedFcn', @(btn,event) runAnalysis(cb_SD.Value, cb_Bathy.Value, cb_SVP.Value, fig));
uiwait(fig);

function runAnalysis(do_SD, do_Bathy, do_SVP, fig)
    delete(fig); % Close UI
    
    %% Base Parameters
    N = 100;
    freq = 10000; maxR = 50000; fom = 100;
    geo = [0.989*1500 1.63 0.07]; % Untouched
    
    % Bounds arrays: [min, max]
    b_SD = [4.75, 5.25];
    b_SVP = [-1, 1];     % Temperature shift
    b_Bathy = [-0.1, 0.1]; % Slope angle in degrees
    
    numVars = sum([do_SD, do_Bathy, do_SVP]);
    if numVars == 0
        disp('No variables selected!'); return;
    end
    
    %% --- PART 1: Latin Hypercube Monte Carlo ---
    disp('--- Running LHS Monte Carlo ---');
    lhs_matrix = lhsdesign(N, numVars);
    
    % Map LHS [0,1] to physical bounds
    col = 1;
    if do_SD
        lhs_SD = b_SD(1) + lhs_matrix(:,col)*(b_SD(2)-b_SD(1)); col = col+1;
    else
        lhs_SD = 5 * ones(N,1);
    end
    if do_Bathy
        lhs_Bathy = b_Bathy(1) + lhs_matrix(:,col)*(b_Bathy(2)-b_Bathy(1)); col = col+1;
    else
        lhs_Bathy = 0 * ones(N,1);
    end
    if do_SVP
        lhs_SVP = b_SVP(1) + lhs_matrix(:,col)*(b_SVP(2)-b_SVP(1));
    else
        lhs_SVP = 0 * ones(N,1);
    end
    
    % Pre-allocate
    sim_pars = {freq, maxR, 5, 0, 0, "summer", "const_35", geo, fom};
    [TL_dummy, r_grid, z_grid] = simpleBellhopHazat(sim_pars);
    TL_MC = zeros(size(TL_dummy,1), size(TL_dummy,2), N);
    
    for i = 1:N
        fprintf('LHS Run %d/%d\n', i, N);
        % Update SD
        sim_pars{3} = lhs_SD(i);
        
        % Update Bathy Custom
        if lhs_Bathy(i) == 0
            sim_pars{7} = "const_35";
            custom_Bathy = [];
        else
            sim_pars{7} = "custom";
            end_depth = 35 + tan(deg2rad(lhs_Bathy(i))) * maxR;
            custom_Bathy = [linspace(0, maxR, 100).', linspace(35, end_depth, 100).'];
        end
        
        % Update SVP Custom
        if lhs_SVP(i) == 0
            sim_pars{6} = "summer";
            custom_SVP = [];
        else
            sim_pars{6} = "custom";
            temps = [0 25; 30 25; 180 17; 400 13.6; 5000 13.6];
            temps(:,2) = temps(:,2) + lhs_SVP(i); % Apply temp shift
            depths = linspace(0, 5000, 1000);     % Create custom SV profile
            t_interp = interp1(temps(:,1), temps(:,2), depths);
            S = 37;
            sv = 1499.2 + 4.6.*t_interp - 0.055.*t_interp.^2 + 0.00029*t_interp.^3 + ... 
                (1.34 - 0.01.*t_interp).*(S-35) + 0.016.*depths;
            custom_SVP = [depths.' sv.'];
        end
        
        % Run Bellhop
        [TL, ~, ~] = simpleBellhopHazat(sim_pars, 'CustomBathymetry', custom_Bathy, 'CustomSVP', custom_SVP);
        TL_MC(:,:,i) = TL;
    end
    
    Expected_LHS = mean(TL_MC, 3);
    Var_LHS = var(TL_MC, 0, 3);
    
    %% --- PART 2: MULTI-VARIABLE DELTA METHOD ---
    disp('--- Running Multi-Variable Delta Method ---');
    sim_pars = {freq, maxR, 5, 0, 0, "summer", "const_35", geo, fom};
    [TL_nom, ~, ~] = simpleBellhopHazat(sim_pars);
    
    Var_Delta_Total = zeros(size(TL_nom));
    
    if do_SD
        h = 0.1; 
        sim_pars{3} = 5 + h; [TL_up, ~, ~] = simpleBellhopHazat(sim_pars);
        sim_pars{3} = 5 - h; [TL_dn, ~, ~] = simpleBellhopHazat(sim_pars);
        sim_pars{3} = 5; % Reset
        
        d_SD = (TL_up - TL_dn) / (2*h);
        var_input = ((b_SD(2)-b_SD(1))^2)/12;
        Var_Delta_Total = Var_Delta_Total + (d_SD.^2 .* var_input);
    end
    
    if do_Bathy
        h = 0.02; % Degrees perturbation
        % Run Up (+h)
        sim_pars{7} = "custom";
        end_depth = 35 + tan(deg2rad(h)) * maxR;
        cB_up = [linspace(0, maxR, 100).', linspace(35, end_depth, 100).'];
        [TL_up, ~, ~] = simpleBellhopHazat(sim_pars, 'CustomBathymetry', cB_up);
        
        % Run Down (-h)
        end_depth = 35 + tan(deg2rad(-h)) * maxR;
        cB_dn = [linspace(0, maxR, 100).', linspace(35, end_depth, 100).'];
        [TL_dn, ~, ~] = simpleBellhopHazat(sim_pars, 'CustomBathymetry', cB_dn);
        sim_pars{7} = "const_35"; % Reset
        
        d_Bathy = (TL_up - TL_dn) / (2*h);
        var_input = ((b_Bathy(2)-b_Bathy(1))^2)/12;
        Var_Delta_Total = Var_Delta_Total + (d_Bathy.^2 .* var_input);
    end
    
    if do_SVP
        h = 0.2; % Temp shift perturbation
        depths = linspace(0, 5000, 1000);
        base_temps = [0 25; 30 25; 180 17; 400 13.6; 5000 13.6];
        
        % Function handle for SV inline
        getSVP = @(t_shift) [depths.' (1499.2 + 4.6.*interp1(base_temps(:,1), base_temps(:,2)+t_shift, depths) - ...
            0.055.*interp1(base_temps(:,1), base_temps(:,2)+t_shift, depths).^2 + ...
            0.00029*interp1(base_temps(:,1), base_temps(:,2)+t_shift, depths).^3 + ...
            (1.34 - 0.01.*interp1(base_temps(:,1), base_temps(:,2)+t_shift, depths)).*(37-35) + 0.016.*depths).'];
            
        sim_pars{6} = "custom";
        [TL_up, ~, ~] = simpleBellhopHazat(sim_pars, 'CustomSVP', getSVP(h));
        [TL_dn, ~, ~] = simpleBellhopHazat(sim_pars, 'CustomSVP', getSVP(-h));
        sim_pars{6} = "summer"; % Reset
        
        d_SVP = (TL_up - TL_dn) / (2*h);
        var_input = ((b_SVP(2)-b_SVP(1))^2)/12;
        Var_Delta_Total = Var_Delta_Total + (d_SVP.^2 .* var_input);
    end
    
    %% --- PART 3: PLOTTING ---
    figure('Name', 'LHS vs Multi-Delta Analysis', 'Position', [100 100 1200 800]);
    r_km = r_grid/1000;
    
    subplot(2,2,1); pcolor(r_km, z_grid, Expected_LHS); shading interp; set(gca,'YDir','reverse');
    title('LHS Monte Carlo: Expected TL'); colormap(gca, jet); colorbar; caxis([50 150]);
    
    subplot(2,2,2); pcolor(r_km, z_grid, TL_nom); shading interp; set(gca,'YDir','reverse');
    title('Delta Method: Expected TL (Nominal)'); colormap(gca, jet); colorbar; caxis([50 150]);
    
    subplot(2,2,3); pcolor(r_km, z_grid, Var_LHS); shading interp; set(gca,'YDir','reverse');
    title('LHS Monte Carlo: Variance'); colormap(gca, hot); colorbar;
    
    subplot(2,2,4); pcolor(r_km, z_grid, Var_Delta_Total); shading interp; set(gca,'YDir','reverse');
    title('Delta Method: Total Variance'); colormap(gca, hot); colorbar;
    %% --- PART 4: SAVE RESULTS ---
    % Create a dynamic filename based on what was checked
    varNames = "";
    if do_SD, varNames = varNames + "_SD"; end
    if do_Bathy, varNames = varNames + "_Bathy"; end
    if do_SVP, varNames = varNames + "_SVP"; end
    
    filename = strcat('Multi_UQ_Results', varNames, '.mat');
    
    disp(['Saving results to ', char(filename), '...']);
    save(filename, 'Expected_LHS', 'Var_LHS', 'TL_nom', 'Var_Delta_Total', 'r_grid', 'z_grid');
    disp('Suite Finished Successfully!');
end