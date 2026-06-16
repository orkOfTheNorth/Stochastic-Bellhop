function run_analysis_only()
%% run_analysis_only — PCE + LN3 + Comparison on already-cached MC combos.
% Does NOT run new Bellhop/MC simulations.
% Safe to rerun: each step checks its output cache and skips if done.
close all; clc; warning('off');
try, cd(fileparts(mfilename('fullpath'))); catch; end

addpath(genpath('Shared_Utils'));
addpath(genpath('Bellhop'));

cfg = loadConfig();

res_root = fullfile('Methods', 'MC');

for si = 1:numel(cfg.scenarios)
    sc = cfg.scenarios(si);
    for di = 1:numel(cfg.distributions)
        dist = cfg.distributions(di);

        % Only process combos where ALL THREE MC result files exist
        res_dir = fullfile(res_root, sc.name, dist.name, 'results');
        if ~isfile(fullfile(res_dir,'MC_freq.mat')) || ...
           ~isfile(fullfile(res_dir,'MC_zS.mat'))   || ...
           ~isfile(fullfile(res_dir,'MC_svp.mat'))
            fprintf('  SKIP (MC missing): %s / %s\n', sc.name, dist.name);
            continue;
        end

        fprintf('\n====== ANALYSIS: %s / %s ======\n', sc.name, dist.name);
        run_pce(sc.name, dist.name);
        run_LN3(sc.name, dist.name);
        run_Comparison(sc.name, dist.name);
    end
end

fprintf('\n====== run_analysis_only complete. ======\n');
end
