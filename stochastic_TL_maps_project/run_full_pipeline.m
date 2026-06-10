function run_full_pipeline()
%% run_full_pipeline — Combo-by-combo: MC → PCE → LN3 → Comparison.
% Each (scenario, distribution) pair is fully analysed before moving on.
% Safe to restart mid-run — each step caches its own output and skips if done.
close all; clc; warning('off');
try, cd(fileparts(mfilename('fullpath'))); catch; end

addpath(genpath('Shared_Utils'));
addpath(genpath('Bellhop'));

cfg = loadConfig();

if isempty(gcp('nocreate'))
    parpool('local', 10);
end

for si = 1:numel(cfg.scenarios)
    sc = cfg.scenarios(si);
    for di = 1:numel(cfg.distributions)
        dist = cfg.distributions(di);
        fprintf('\n====== PIPELINE: %s / %s ======\n', sc.name, dist.name);

        run_MC(sc.name, dist.name);
        run_pce(sc.name, dist.name);
        run_LN3(sc.name, dist.name);
        run_Comparison(sc.name, dist.name);
    end
end

fprintf('\n====== run_full_pipeline complete. ======\n');
end
