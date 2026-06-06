function [bathymetry] = bathymetryMaker(bathymetryType, maxRange)
    % This function creates a SVP according a given SVP type.

    %% Arguments
    arguments
        bathymetryType (1, 1) string
        maxRange (1, 1) double {mustBePositive}
    end

    %% Extract parameters
    splitStr = split(bathymetryType, '_');
    typeStr = splitStr(1);

    if strcmp(typeStr, "const")
        maxDepth = str2double(splitStr(2));
        bathymetry = [linspace(0, maxRange, 100).', maxDepth*ones(100, 1)];
    else if strcmp(typeStr, "slope")
        startDepth = str2double(splitStr(2));
        endDepth = str2double(splitStr(3));
        bathymetry = [linspace(0, maxRange, 100).', ...
            linspace(startDepth, endDepth, 100).'];
    else
        error("bathymeryType has to be const or slope.")
    end
end