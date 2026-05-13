function createBellhopEnv(envInputs) % Create bellhop.env file
    arguments
        envInputs (1, 9)
    end

    %% Extract Parameters from INPUTS
    ttl = envInputs{1};
    frequency = envInputs{2};
    rd = envInputs{3};
    svp = envInputs{4};
    geoAcoustics = envInputs{5}; % [c_bottom/c_water rho alpha]
    sourceDepth = envInputs{6};
    maxDepth = envInputs{7};
    range = envInputs{8};
    beam = envInputs{9};


    %% Open File
    fileID = fopen('bellhop.env', 'w');
    
    %% Title
    fprintf(fileID, ["'"]);
    fprintf(fileID, '%s', ttl);
    fprintf(fileID, ["'"]);
    fprintf(fileID, '\n');

    %% Frequency
    fprintf(fileID, [num2str(frequency) '\n']);

    %% nmedia (dummy integer < 20)
    fprintf(fileID, ['1' '\n']);

    %% Type of SVP
    if rd
        fprintf(fileID, ["'QVWT'"]); % 'CVWT' = Single, 'QVWT' = Range Dependent
    else
        fprintf(fileID, ["'CVWT'"]);
    end
    fprintf(fileID, ['\n']);

    %% SVP formatted accordingly
    fprintf(fileID, [num2str(length(svp)) ' 0.0 ' num2str(svp(end, 1)) ,'\n']);
    for i = 1:length(svp(:, 1))
        fprintf(fileID, [num2str([svp(i, 1), svp(i, 2)]) '  /\n']);
    end

    %% Filler
    fprintf(fileID, ["'A*' 0.0"]);
    fprintf(fileID, '\n');

    %% Deepest depth of SVP | c_bottom | 0.0 | rho | alpha
    svp_z0 = interp1(svp(:, 1), svp(:, 2), sourceDepth);
    fprintf(fileID, [num2str(svp(end, 1)) ' ' ...
        num2str(geoAcoustics(1)) ' 0.0 ' ...
        num2str(geoAcoustics(2)) ' ' num2str(geoAcoustics(3)) ' /\n']);

    %% Number of sources
    fprintf(fileID, ['1\n']);

    %% Source depth
    fprintf(fileID, [num2str(sourceDepth) '   /\n']);

    %% nrd
    fprintf(fileID, ['574      /\n']);

    %% 0.0 | Bottom max depth
    fprintf(fileID, ['0.0   ' num2str(maxDepth) '   /		! RD(1:NRD) (m)' '\n']);

    %% nrr
    fprintf(fileID, ['766  /   		! NR' '\n']);

    %% 0.0 | Range (km)
    fprintf(fileID, ['0.0 ' num2str(range/1000) '   /		! R(1:NR ) (km)' '\n']);

    %% 'SB' or 'CB'
    fprintf(fileID, ["'SB'"]);
    fprintf(fileID, '\n');

    %% nbeams
    fprintf(fileID, ['1500 /   ! NBEAMS' '\n']);

    %% Aperture
    if beam(1) == 0
        fprintf(fileID, [num2str([-90, 90]) '   /    ! ALPHA1, 2 (degrees)' '\n']);
    else
        fprintf(fileID, [num2str(beam) '   /    ! ALPHA1, 2 (degrees)' '\n']);
    end

    %% 0.0 | Bottom max depth | Range (km)
    fprintf(fileID, ['0.0 ' num2str(1.1*maxDepth) ' ' num2str(range/1000) '  / ! STEP (m), ZBOX (m), RBOX (km)   \n']);
    fclose(fileID);
    
end