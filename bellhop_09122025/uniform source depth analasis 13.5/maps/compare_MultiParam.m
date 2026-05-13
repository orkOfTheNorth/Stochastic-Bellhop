%% run_MultiParam_UQ.m (100m Base Depth & Live SVP Plotting)
clear; close all; clc; warning('off');

% Go up one level from 'maps' and then into 'code/Functions'
addpath('../code/Functions'); 

% --- Create UI for Selection ---
fig = uifigure('Name', 'Select Variables to Perturb', 'Position', [500 500 350 250]);
lbl = uilabel(fig, 'Position', [20 200 300 22], 'Text', 'Select parameters for UQ Analysis:');
cb_SD = uicheckbox(fig, 'Position', [40 160 250 22], 'Text', 'Source Depth (5m ± 5%)', 'Value', 1);
cb_Bathy = uicheckbox(fig, 'Position', [40 120 250 22], 'Text', 'Bathy Slope (0° ± 0.1°)', 'Value', 1);
cb_SVP = uicheckbox(fig, 'Position', [40 80 300 22], 'Text', 'SVP Nodes (Point-by-point ± 1°C)', 'Value', 1);
btn = uibutton(fig, 'push', 'Position', [125 20 100 30], 'Text', 'Run Suite', ...
    'ButtonPushedFcn', @(btn,event) runAnalysis(cb_SD.Value, cb_Bathy.Value, cb_SVP.Value, fig));
uiwait(fig);

function runAnalysis(do_SD, do_Bathy, do_SVP, fig)
    delete(fig); % Close UI
    
    %% Base Parameters
    N = 100;
    freq = 10000; maxR = 50000; fom = 100;
    geo = [0.989*1500 1.63 0.07]; % Untouched bottom properties
    
    % Bounds arrays: [min, max]
    b_SD = [4.75, 5.25];
    b_Bathy = [-0.1, 0.1]; 
    b_SVP = [-1, 1];       

    % SVP Base Profile Definition (5 Nodes)
    base_depths_SVP = [0; 30; 180; 400; 5000];
    base_temps_SVP  = [25; 25; 17; 13.6; 13.6];
    num_SVP_nodes   = length(base_depths_SVP);
    
    % Calculate total LHS variables needed
    numVars = do_SD + do_Bathy;
    if do_SVP
        numVars = numVars + num_SVP_nodes; 
    end
    
    if numVars == 0
        disp('No variables selected!'); return;
    end
    
    %% --- PART 1: Latin Hypercube Monte Carlo ---
    disp('--- Running LHS Monte Carlo ---');
    lhs_matrix = lhsdesign(N, numVars);
    
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
        lhs_SVP_shifts = b_SVP(1) + lhs_matrix(:, col:col+num_SVP_nodes-1)*(b_SVP(2)-b_SVP(1));
    else
        lhs_SVP_shifts = zeros(N, num_SVP_nodes);
    end
    
    % Pre-allocate based on 100m depth
    sim_pars = {freq, maxR, 5, 0, 0, "summer", "const_100", geo, fom};
    [TL_dummy, r_grid, z_grid] = simpleBellhopHazat(sim_pars);
    TL_MC = zeros(size(TL_dummy,1), size(TL_dummy,2), N);
    
    % --- LIVE PLOT SETUP ---
    high_res_depths = linspace(0, 5000, 1000)';
    if do_SVP
        fig_svp_live = figure('Name', 'Live SVP Perturbations', 'Position', [1300 200 400 700]);
        ax_svp = axes(fig_svp_live);
        hold(ax_svp, 'on');
        
        % Plot the nominal (base) profile first
        t_base = interp1(base_depths_SVP, base_temps_SVP, high_res_depths, 'linear', 'extrap');
        sv_base = 1499.2 + 4.6.*t_base - 0.055.*t_base.^2 + 0.00029*t_base.^3 + ... 
                (1.34 - 0.01.*t_base).*(37-35) + 0.016.*high_res_depths;
        plot(ax_svp, sv_base, high_res_depths, 'r', 'LineWidth', 2.5);
        
        set(ax_svp, 'YDir', 'reverse');
        ylim(ax_svp, [0 500]); % Zooming into the top 500m
        xlabel(ax_svp, 'Sound Speed (m/s)');
        ylabel(ax_svp, 'Depth (m)');
        title(ax_svp, 'Live SVP "Lightning" Perturbations');
        grid(ax_svp, 'on');
        box(ax_svp, 'on');
    end
    
    % --- THE MONTE CARLO LOOP ---
    for i = 1:N
        fprintf('LHS Run %d/%d\n', i, N);
        
        sim_pars{3} = lhs_SD(i);
        
        if lhs_Bathy(i) == 0
            sim_pars{7} = "const_100";
            custom_Bathy = [];
        else
            sim_pars{7} = "custom";
            % Base depth is now 100m
            end_depth = 100 + tan(deg2rad(lhs_Bathy(i))) * maxR;
            custom_Bathy = [linspace(0, maxR, 100).', linspace(100, end_depth, 100).'];
        end
        
        if do_SVP
            sim_pars{6} = "custom";
            perturbed_temps = base_temps_SVP + lhs_SVP_shifts(i, :)';
            
            t_interp = interp1(base_depths_SVP, perturbed_temps, high_res_depths, 'linear', 'extrap');
            S = 37;
            sv = 1499.2 + 4.6.*t_interp - 0.055.*t_interp.^2 + 0.00029*t_interp.^3 + ... 
                (1.34 - 0.01.*t_interp).*(S-35) + 0.016.*high_res_depths;
            custom_SVP = [high_res_depths sv];
            
            % Plot this specific branch live
            plot(ax_svp, sv, high_res_depths, 'Color', [0.2 0.5 0.8 0.3], 'LineWidth', 0.5);
            drawnow; % Forces MATLAB to update the figure window immediately
            
        else
            sim_pars{6} = "summer";
            custom_SVP = [];
        end
        
        [TL, ~, ~] = simpleBellhopHazat(sim_pars, 'CustomBathymetry', custom_Bathy, 'CustomSVP', custom_SVP);
        TL_MC(:,:,i) = TL;
    end
    
    Expected_LHS = mean(TL_MC, 3);
    Var_LHS = var(TL_MC, 0, 3);
    
    %% --- PART 2: MULTI-VARIABLE DELTA METHOD ---
    disp('--- Running Multi-Variable Delta Method ---');
    sim_pars = {freq, maxR, 5, 0, 0, "summer", "const_100", geo, fom};
    [TL_nom, ~, ~] = simpleBellhopHazat(sim_pars);
    
    Var_Delta_Total = zeros(size(TL_nom));
    
    if do_SD
        disp('Calculating Partial Derivative: Source Depth');
        h = 0.1; 
        sim_pars{3} = 5 + h; [TL_up, ~, ~] = simpleBellhopHazat(sim_pars);
        sim_pars{3} = 5 - h; [TL_dn, ~, ~] = simpleBellhopHazat(sim_pars);
        sim_pars{3} = 5; 
        
        d_SD = (TL_up - TL_dn) / (2*h);
        var_input = ((b_SD(2)-b_SD(1))^2)/12;
        Var_Delta_Total = Var_Delta_Total + (d_SD.^2 .* var_input);
    end
    
    if do_Bathy
        disp('Calculating Partial Derivative: Bathymetry');
        h = 0.02; 
        sim_pars{7} = "custom";
        
        % Base depth 100m for perturbations
        end_depth = 100 + tan(deg2rad(h)) * maxR;
        cB_up = [linspace(0, maxR, 100).', linspace(100, end_depth, 100).'];
        [TL_up, ~, ~] = simpleBellhopHazat(sim_pars, 'CustomBathymetry', cB_up);
        
        end_depth = 100 + tan(deg2rad(-h)) * maxR;
        cB_dn = [linspace(0, maxR, 100).', linspace(100, end_depth, 100).'];
        [TL_dn, ~, ~] = simpleBellhopHazat(sim_pars, 'CustomBathymetry', cB_dn);
        
        sim_pars{7} = "const_100"; 
        d_Bathy = (TL_up - TL_dn) / (2*h);
        var_input = ((b_Bathy(2)-b_Bathy(1))^2)/12;
        Var_Delta_Total = Var_Delta_Total + (d_Bathy.^2 .* var_input);
    end
    
    if do_SVP
        disp('Calculating Partial Derivatives: SVP Nodes (1 to 5)');
        h = 0.2; 
        sim_pars{6} = "custom";
        var_input = ((b_SVP(2)-b_SVP(1))^2)/12;
        
        for pt = 1:num_SVP_nodes
            fprintf('   -> Perturbing Node %d at depth %d m\n', pt, base_depths_SVP(pt));
            
            temps_up = base_temps_SVP; temps_up(pt) = temps_up(pt) + h;
            t_interp_up = interp1(base_depths_SVP, temps_up, high_res_depths, 'linear', 'extrap');
            sv_up = 1499.2 + 4.6.*t_interp_up - 0.055.*t_interp_up.^2 + 0.00029*t_interp_up.^3 + ... 
                (1.34 - 0.01.*t_interp_up).*(37-35) + 0.016.*high_res_depths;
            [TL_up, ~, ~] = simpleBellhopHazat(sim_pars, 'CustomSVP', [high_res_depths sv_up]);
            
            temps_dn = base_temps_SVP; temps_dn(pt) = temps_dn(pt) - h;
            t_interp_dn = interp1(base_depths_SVP, temps_dn, high_res_depths, 'linear', 'extrap');
            sv_dn = 1499.2 + 4.6.*t_interp_dn - 0.055.*t_interp_dn.^2 + 0.00029*t_interp_dn.^3 + ... 
                (1.34 - 0.01.*t_interp_dn).*(37-35) + 0.016.*high_res_depths;
            [TL_dn, ~, ~] = simpleBellhopHazat(sim_pars, 'CustomSVP', [high_res_depths sv_dn]);
            
            d_SVP_pt = (TL_up - TL_dn) / (2*h);
            Var_Delta_Total = Var_Delta_Total + (d_SVP_pt.^2 .* var_input);
        end
        sim_pars{6} = "summer"; 
    end
    
    %% --- PART 3: PLOTTING & SAVING ---
    figure('Name', 'LHS vs Multi-Delta Analysis', 'Position', [50 100 1200 800]);
    r_km = r_grid/1000;
    
    subplot(2,2,1); pcolor(r_km, z_grid, Expected_LHS); shading interp; set(gca,'YDir','reverse');
    title('LHS Monte Carlo: Expected TL'); colormap(gca, jet); colorbar; caxis([50 150]);
    
    subplot(2,2,2); pcolor(r_km, z_grid, TL_nom); shading interp; set(gca,'YDir','reverse');
    title('Delta Method: Expected TL (Nominal)'); colormap(gca, jet); colorbar; caxis([50 150]);
    
    subplot(2,2,3); pcolor(r_km, z_grid, Var_LHS); shading interp; set(gca,'YDir','reverse');
    title('LHS Monte Carlo: Variance'); colormap(gca, hot); colorbar;
    
    subplot(2,2,4); pcolor(r_km, z_grid, Var_Delta_Total); shading interp; set(gca,'YDir','reverse');
    title('Delta Method: Total Variance'); colormap(gca, hot); colorbar;
    
    % Save the Live SVP plot if it was generated
    if do_SVP
        disp('Saving Live SVP plot...');
        saveas(fig_svp_live, 'Live_SVP_Lightning_Plot.png');
        savefig(fig_svp_live, 'Live_SVP_Lightning_Plot.fig');
    end

    varNames = "";
    if do_SD, varNames = varNames + "_SD"; end
    if do_Bathy, varNames = varNames + "_Bathy"; end
    if do_SVP, varNames = varNames + "_SVP"; end
    
    filename = strcat('Multi_UQ_Results_100m', varNames, '.mat');
    
    disp(['Saving main results to ', char(filename), '...']);
    save(filename, 'Expected_LHS', 'Var_LHS', 'TL_nom', 'Var_Delta_Total', 'r_grid', 'z_grid');
    disp('Suite Finished Successfully!');
end