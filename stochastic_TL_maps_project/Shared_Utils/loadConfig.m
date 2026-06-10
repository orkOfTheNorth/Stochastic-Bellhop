function cfg = loadConfig()
% Read config.json from the project root (two levels up from Shared_Utils).
    here = fileparts(mfilename('fullpath'));
    cfg_path = fullfile(here, '..', 'config.json');
    if ~isfile(cfg_path)
        error('config.json not found at %s', cfg_path);
    end
    cfg = jsondecode(fileread(cfg_path));
end
