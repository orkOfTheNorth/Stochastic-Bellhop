% Script to get array shapes from TL cube .mat files
% Output: shapes.txt with pipe-delimited lines: scenario|dist|varname|dims

cacheRoot = 'C:\Users\orind\Stochastic-Bellhop\stochastic_TL_maps_project\Cache\';
outFile = 'C:\Users\orind\Stochastic-Bellhop\shapes.txt';

fid = fopen(outFile, 'w');

scenarios = {'baseline','deep_water','downslope','shallow_water','upslope'};
tlNames   = {'TL_freq','TL_zS','TL_svp'};

for si = 1:numel(scenarios)
    sc = scenarios{si};
    scPath = fullfile(cacheRoot, sc);

    % Get distribution dirs (not bellhop_raw, not pce, not svp_profiles)
    d = dir(scPath);
    distDirs = {};
    for di = 1:numel(d)
        if d(di).isdir && ~startsWith(d(di).name,'.') && ...
           ~strcmp(d(di).name,'bellhop_raw') && ~strcmp(d(di).name,'pce')
            distDirs{end+1} = d(di).name; %#ok<AGROW>
        end
    end

    for di = 1:numel(distDirs)
        dist = distDirs{di};
        for ti = 1:numel(tlNames)
            tlName = tlNames{ti};
            matPath = fullfile(scPath, dist, [tlName '.mat']);
            if exist(matPath, 'file')
                try
                    info = whos('-file', matPath);
                    % Find TL variable (usually same name as file stem)
                    found = false;
                    for vi = 1:numel(info)
                        v = info(vi);
                        if numel(v.size) >= 2
                            dims = sprintf('%d', v.size(1));
                            for di2 = 2:numel(v.size)
                                dims = [dims ',' sprintf('%d', v.size(di2))]; %#ok<AGROW>
                            end
                            fprintf(fid, '%s|%s|%s|%s|%s\n', sc, dist, tlName, v.name, dims);
                            found = true;
                        end
                    end
                    if ~found
                        fprintf(fid, '%s|%s|%s|NO_ARRAY_VAR|N/A\n', sc, dist, tlName);
                    end
                catch e
                    fprintf(fid, '%s|%s|%s|ERROR|%s\n', sc, dist, tlName, e.message);
                end
            end
        end
    end
end

fclose(fid);
fprintf('Done. Written to %s\n', outFile);
