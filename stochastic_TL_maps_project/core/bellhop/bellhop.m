function bellhop( filename )
% Run the BELLHOP acoustic propagation model.
%
% Usage: bellhop( filename )
%   filename — base name of the .env file (no extension)
%
% OS-aware binary selection:
%   Windows : bellhopcuda.exe / bellhop.exe
%   Linux   : bellhopcuda     / bellhop_cuda / bellhop   (no .exe)
%
% GPU acceleration: CUDA build is preferred when a compatible GPU
% (compute capability >= 3.5) is detected; CPU build is the fallback.

runbellhop = '';

if ispc
    cuda_names = {'bellhopcuda.exe'};
    cpu_names  = {'bellhop.exe'};
else
    % Linux / macOS — AT toolbox compiles without .exe suffix
    cuda_names = {'bellhopcuda', 'bellhop_cuda'};
    cpu_names  = {'bellhop'};
end

% Prefer CUDA build when a compatible GPU is available
if gpuAvailable()
    for k = 1:numel(cuda_names)
        p = which(cuda_names{k});
        if ~isempty(p)
            runbellhop = p;
            fprintf('[bellhop] Using CUDA build: %s\n', p);
            break;
        end
    end
end

% CPU fallback
if isempty(runbellhop)
    for k = 1:numel(cpu_names)
        p = which(cpu_names{k});
        if ~isempty(p)
            runbellhop = p;
            break;
        end
    end
end

if isempty(runbellhop)
    if ispc
        hint = 'bellhop.exe / bellhopcuda.exe';
    else
        hint = 'bellhop / bellhopcuda (Linux build from atoolbox)';
    end
    error('bellhop:notFound', ...
        'Bellhop binary not found on the MATLAB path.\nExpected: %s\n%s', ...
        hint, 'Add the Bellhop/ directory: addpath(genpath(''Bellhop''))');
end

% Quote path to handle spaces; on Linux single-quotes are safer but
% MATLAB system() always uses /bin/sh so double-quotes work on both.
[status, cmdout] = system(['"' runbellhop '" ' filename]);
if status ~= 0
    error('bellhop:runFailed', ...
        'Bellhop exited with code %d for file "%s".\nOutput:\n%s', ...
        status, filename, cmdout);
end
end


function result = gpuAvailable()
% Returns true if a CUDA GPU with compute capability >= 3.5 is accessible.
try
    g = gpuDevice();
    result = str2double(g.ComputeCapability) >= 3.5;
catch
    result = false;
end
end
