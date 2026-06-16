function [TL, r_grid, z_grid] = simpleBellhopHazat(sim_pars, options)
    % Run a range-prediction calculation using BELLHOP and return the TL field.
    %% Define arguments
    arguments
        sim_pars (1, 9)
        options.CustomSVP = [];
        options.CustomBathymetry = [];
    end

    %% Read simulation parameters
    sourceFrequency = sim_pars{1}; % Hz
    maxRange        = sim_pars{2}; % m
    sourceDepth     = sim_pars{3}; % m
    sourceHalfBeam  = sim_pars{4}; % deg
    sourceTilt      = sim_pars{5}; % deg
    svpType         = sim_pars{6}; % string
    bathymetryType  = sim_pars{7}; % string
    geoAcoustics    = sim_pars{8}; % [c_bottom/c_water rho alpha]
    FOM             = sim_pars{9}; % dB

    %% Build bathymetry and SVP
    rd = false; % range-independent SVP

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

    sourceBeam = [sourceTilt-sourceHalfBeam, sourceTilt+sourceHalfBeam];
    ttl = join([ ...
          'f=' num2str(sourceFrequency) '_' ...
          'sD=' num2str(sourceDepth) '_' ...
          'beamhalfwidth=' num2str(sourceHalfBeam) '_' ...
          'bathymetryType=' bathymetryType '_' ...
          'svpType=' svpType ...
          ]);
    envInputs = {ttl, sourceFrequency, rd, svp, geoAcoustics, ...
        sourceDepth, maxDepth, maxRange, sourceBeam};

    %% Create BELLHOP input files
    createBellhopEnv(envInputs);
    createBellhopBty(bathymetry);

    %% Remove stale output files from any previous run
    for f = {'bellhop.shd', 'bellhop.prt'}
        if isfile(f{1}), delete(f{1}); end
    end

    %% Run BELLHOP
    bellhop('bellhop');

    %% Extract TL from the shadow file
    [~, ~, ~, ~, ~, Pos, pressure] = read_shd('bellhop.shd');
    r_grid = Pos.r.r;  % range vector (m)
    z_grid = Pos.r.z;  % depth vector (m)

    % eps avoids log(0) at complete shadow-zone pixels
    TL = squeeze(-20 * log10(abs(pressure) + eps));

    %% Optional interactive plot (skipped in batch / headless runs)
    if strcmp(get(0, 'DefaultFigureVisible'), 'on')
        plotshd('bellhop.shd', FOM);
        clim([50 150]);
        title(ttl);
        colormap(colors('bar.png'));
        hold on;
        plot(bathymetry(:,1)/1000, bathymetry(:,2), 'LineStyle', '--', 'Color', 'w');
        hold off;
    end

    %% Clean up temporary Bellhop files
    for f = {'bellhop.env', 'bellhop.bty', 'bellhop.ssp', 'bellhop.shd', 'bellhop.prt'}
        if isfile(f{1}), delete(f{1}); end
    end
end
