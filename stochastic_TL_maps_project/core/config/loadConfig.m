function cfg = loadConfig()
% Read config.json from the project root and expose convenience helpers.
    here     = fileparts(mfilename('fullpath'));  % → core/config
    root_dir = fileparts(fileparts(here));        % → core → project root
    cfg_path = fullfile(root_dir, 'config.json');
    if ~isfile(cfg_path)
        error('config.json not found at %s', cfg_path);
    end
    cfg = jsondecode(fileread(cfg_path));

    % Convenience: N_MC from MC.N
    cfg.N_MC = cfg.MC.N;

    % Convenience: PCE max order and LOO weight
    if isfield(cfg, 'PCE')
        cfg.max_pce_order = cfg.PCE.max_order;
        cfg.LOO_alpha     = cfg.PCE.LOO_alpha;
    else
        cfg.max_pce_order = 10;
        cfg.LOO_alpha     = 0.5;
    end
end
