function [TL, r_grid, z_grid] = simpleBellhopHazat(sim_pars, options)
    % This function runs a range predicition calculation using BELLHOP.
    %% Define arguments
    arguments
        sim_pars (1, 9)
        options.CustomSVP = [];
        options.CustomBathymetry = [];
    end
    
    %% Read simulation parameters needed
    sourceFrequency = sim_pars{1}; % Hz
    maxRange = sim_pars{2}; % m
    sourceDepth = sim_pars{3}; % m
    sourceHalfBeam = sim_pars{4}; % deg
    sourceTilt = sim_pars{5}; % deg
    svpType = sim_pars{6}; % string
    bathymetryType = sim_pars{7}; % string
    geoAcoustics = sim_pars{8}; % [c_bottom/c_water rho alpha]
    FOM = sim_pars{9}; % dB
    
    %% Prepare parameters for files creation
    rd = false; % Don't use range-dependent SVP
    
    if strcmp(bathymetryType, "custom")
        bathymetry = options.CustomBathymetry;
    else
        bathymetry = bathymetryMaker(bathymetryType, maxRange);
    end
    
    maxDepth = max(bathymetry(:, 2));
    if strcmp(svpType, "custom")
        svp = options.CustomSVP;
    else
        svp = svpMaker(svpType, maxDepth);
    end
    
    sourceBeam = [sourceTilt-sourceHalfBeam sourceTilt+sourceHalfBeam];
    ttl = join([ % Title for plot
          'f=' num2str(sourceFrequency) '_' ...
          'sD=' num2str(sourceDepth) '_' ...
          'beamhalfwidth=' num2str(sourceHalfBeam) '_' ...
          'bathymetryType=' bathymetryType '_' ...
          'svpType=' svpType ...
          ]);
    envInputs = {ttl, sourceFrequency, rd, svp, geoAcoustics, ...
        sourceDepth, maxDepth, maxRange, sourceBeam};
    %% Create relevant BELLHOP files
    createBellhopEnv(envInputs) % Create bellhop.env
    createBellhopBty(bathymetry) % Create bellhop.bty
    
    %% Miscelleanous
    disp('start')
    global units
    units = 'km';
    tic % start measure time
    
    %% Clean previous shd and prt files
    filename = {'bellhop.shd', 'bellhop.prt'};
    for i = 1:length(filename)
        if isfile(filename{i})
            delete(filename{i})
        end
    end
    
    %% Run BELLHOP and plot
    bellhop('bellhop');
    plotshd('bellhop.shd', FOM);
    clim([50 150]);
    title(ttl)
    colormap(colors('bar.png'));
    hold on
    plot(bathymetry(:,1)/1000, bathymetry(:,2), 'LineStyle','--', 'Color', 'w');
    hold off
    
    %% Extract Data to Return
    % This is the new section! We read the binary file Bellhop just made 
    % and format it into standard arrays to send back to the main script.
    [~, ~, ~, ~, ~, Pos, pressure] = read_shd('bellhop.shd');
    r_grid = Pos.r.r;     % Range vector
    z_grid = Pos.r.z;     % Depth vector
    
    % Convert the complex pressure field into Transmission Loss (dB)
    % We add 'eps' to avoid log(0) errors in complete shadow zones
    TL = squeeze(-20 * log10(abs(pressure) + eps)); 

    %% Miscelleanous
    global units jkpsflag
    units = 'km';
    toc % finish measuring runtime
    
    %% Clean bellhop files
    % filename = {'bellhop.env', 'bellhop.bty', 'bellhop.ssp'};
    % for i = 1:length(filename)
    %     if isfile(filename{i})
    %         delete(filename{i})
    %     end
    % end
end