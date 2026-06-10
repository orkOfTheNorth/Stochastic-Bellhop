function [TL, r_grid, z_grid] = bellhopCached(sim_pars, cache_dir, options)
% bellhopCached  –  caching wrapper around simpleBellhopHazat.
%
% On the first call with a given set of parameters, runs Bellhop and saves
% [TL, r_grid, z_grid] to  <cache_dir>/<key>.mat.
% On every subsequent call with identical parameters, loads from cache and
% skips Bellhop entirely.
%
% Usage (mirrors simpleBellhopHazat):
%   [TL, r, z] = bellhopCached(sim_pars, cache_dir)
%   [TL, r, z] = bellhopCached(sim_pars, cache_dir, 'CustomSVP', svp)
%   [TL, r, z] = bellhopCached(sim_pars, cache_dir, 'CustomBathymetry', bty)

    arguments
        sim_pars  (1,9)
        cache_dir string
        options.CustomSVP        = []
        options.CustomBathymetry = []
    end

    if ~exist(cache_dir, 'dir'), mkdir(cache_dir); end

    key        = makeKey(sim_pars, options.CustomSVP, options.CustomBathymetry);
    cache_file = fullfile(cache_dir, [key '.mat']);

    if isfile(cache_file)
        load(cache_file, 'TL', 'r_grid', 'z_grid');
        fprintf('  [CACHE HIT]  %s\n', key);
        return;
    end

    fprintf('  [BELLHOP RUN] %s\n', key);
    [TL, r_grid, z_grid] = simpleBellhopHazat(sim_pars, ...
        'CustomSVP',        options.CustomSVP, ...
        'CustomBathymetry', options.CustomBathymetry);

    save(cache_file, 'TL', 'r_grid', 'z_grid', '-v7.3');
end

%% ── helpers ──────────────────────────────────────────────────────────────────
function key = makeKey(sim_pars, custom_svp, custom_bathy)
    freq  = sim_pars{1};
    maxR  = sim_pars{2};
    zS    = sim_pars{3};
    geo   = sim_pars{8};   % [c_ratio  rho  alpha]
    svpT  = char(sim_pars{6});
    bthyT = char(sim_pars{7});

    key = sprintf('f%.4g_zS%.6g_R%.4g_g%.4g-%.4g-%.4g', ...
                  freq, zS, maxR, geo(1), geo(2), geo(3));

    if strcmp(svpT, 'custom') && ~isempty(custom_svp)
        key = [key '_svp' arrayHash(custom_svp)];
    else
        key = [key '_' svpT];
    end

    if strcmp(bthyT, 'custom') && ~isempty(custom_bathy)
        key = [key '_bty' arrayHash(custom_bathy)];
    else
        key = [key '_' bthyT];
    end

    % Sanitise for filesystem
    key = regexprep(key, '[^A-Za-z0-9_\-\.]', '_');
end

function h = arrayHash(A)
% Fast, toolbox-free hash of a numeric array.
% Samples up to 32 evenly-spaced elements, then does a polynomial roll.
    v   = double(A(:));
    n   = numel(v);
    idx = round(linspace(1, n, min(n, 32)));
    v   = v(idx);
    acc = 0;
    for k = 1:numel(v)
        acc = mod(acc * 1000003 + round(v(k) * 1e6), 2^31 - 1);
    end
    h = sprintf('%08x', abs(acc));
end
