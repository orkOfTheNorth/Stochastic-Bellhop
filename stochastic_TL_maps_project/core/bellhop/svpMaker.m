function [svp] = svpMaker(svpType, maxSVPDepth)
    % This function creates a SVP according a given SVP type.

    %% Arguments
    arguments
        svpType (1, 1) string
        maxSVPDepth (1, 1) double {mustBePositive}
    end

    %% Extract parameters
    splitStr = split(svpType, '_');
    typeStr = splitStr(1);
    
    %% Create depths column
    depths = linspace(0, maxSVPDepth, 1000);

    %% Create velocities column using svpType
    if strcmp(typeStr, "winter")
        temperatureArray = [0 17; 180 17; 400 13.6; 5000 13.6];
        sv = wilsonSV(depths, temperatureArray);
    elseif strcmp(typeStr, "summer")
        temperatureArray = [0 25; 30 25; 180 17; 400 13.6; 5000 13.6];
        sv = wilsonSV(depths, temperatureArray);
    elseif strcmp(typeStr, "const")
        svValue = str2double(splitStr(2));
        sv = svValue*ones(size(depths));
    else
        error('svpType has to be either summer, winter or const.');
    end
    svp = [depths.' sv.'];

end

%% Compute SV using Wilson's formula
function sv = wilsonSV(depths, temperatureArray)
    t = interp1(temperatureArray(:,1), temperatureArray(:,2), depths);
    S = 37; % salinity [promils]
    sv = 1499.2 + 4.6.*t - 0.055.*t.^2 + 0.00029*t.^3 + ... 
        (1.34 - 0.01.*t).*(S-35) + 0.016.*depths;
end