function fom = getFOM(cfg, scen_name)
% Return the FOM (dB) for a scenario, falling back to the global default.
    fom = cfg.nominal.FOM_dB;
    for k = 1:numel(cfg.scenarios)
        if strcmp(cfg.scenarios(k).name, scen_name)
            if isfield(cfg.scenarios(k), 'FOM_dB') && ~isempty(cfg.scenarios(k).FOM_dB)
                fom = cfg.scenarios(k).FOM_dB;
            end
            return;
        end
    end
end
