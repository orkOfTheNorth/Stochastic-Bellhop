function [TL_all, Nz, Nr, N] = loadTLCube(sc_name, dist_name, param_name, cache_root)
% loadTLCube — load a TL field cube from the Bellhop cache.
%
%   [TL_all, Nz, Nr, N] = loadTLCube(sc_name, dist_name, param_name)
%   [TL_all, Nz, Nr, N] = loadTLCube(sc_name, dist_name, param_name, cache_root)
%
% Returns TL_all [Nz × Nr × N] as double.
% Throws an error (caller should catch and skip) if the file is missing.
%
% Default cache_root: 'Cache'

if nargin < 4 || isempty(cache_root)
    cache_root = 'Cache';
end

tl_cache = fullfile(cache_root, sc_name, dist_name, sprintf('TL_%s.mat', param_name));

D = load(tl_cache, 'TL_save');   % throws MException if file absent
TL_all = double(D.TL_save);
[Nz, Nr, N] = size(TL_all);
end
