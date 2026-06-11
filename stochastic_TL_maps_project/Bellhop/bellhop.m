function bellhop( filename )
% Run the BELLHOP acoustic propagation model.
%
% Usage: bellhop( filename )
%   filename — base name of the .env file (no extension)
%
% GPU acceleration: if bellhopcuda.exe is on the MATLAB path AND a
% compatible CUDA GPU (compute capability >= 3.5) is detected, the CUDA
% build is used automatically.  Falls back to bellhop.exe otherwise.

runbellhop = '';

% Prefer CUDA build when a compatible GPU is available
cuda_exe = which('bellhopcuda.exe');
if ~isempty(cuda_exe) && gpuAvailable()
    runbellhop = cuda_exe;
    fprintf('[bellhop] Using CUDA build: %s\n', cuda_exe);
end

% CPU fallback
if isempty(runbellhop)
    runbellhop = which('bellhop.exe');
    if isempty(runbellhop)
        error('bellhop:notFound', ...
            'Neither bellhopcuda.exe nor bellhop.exe found on the MATLAB path.\n%s', ...
            'Add the Bellhop/ directory to the path: addpath(genpath(''Bellhop''))');
    end
end

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
