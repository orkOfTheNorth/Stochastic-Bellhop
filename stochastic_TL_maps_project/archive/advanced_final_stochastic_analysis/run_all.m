%% run_all.m  –  Master runner  (advanced_final_stochastic_analysis)
%
% EXECUTION ORDER AND ESTIMATED TIMES
% ────────────────────────────────────
% Phase 1 – Baseline scenario (const_35, 50km):
%   A. Delta   :  7 runs           ~1 min
%   B. MC      :  7×50 = 350 runs  ~20-30 min  (N=50 per subset)
%   C. LN3     :  grid fits + map  ~5-10 min
%   D. Compare :  (no Bellhop)     ~2 min
%
% Phase 2 – 4 new scenarios (run sequentially):
%   Deep water  : 7+350 runs → ~30-50 min  (50km, deep → slower)
%   Shallow water: 7+350 runs → ~15-25 min  (5km, fast)
%   Upslope     : 7+350 runs → ~20-35 min
%   Downslope   : 7+350 runs → ~20-35 min
%
% Phase 3 – Distribution study:
%   5 configs × 4 subsets × 50 runs = 1000 runs → ~50-90 min
%
% Phase 4 – Decision Tree / Random Forest:
%   RF training on full dataset  → ~5-15 min
%
% TOTAL ESTIMATE: ~4-6 hours (most spent in MC loops)
% Progress shown per-run inside each script (% complete + ETA).
%
% SELECTIVE EXECUTION: comment out phases you don't need.
% Calls exit at the end so MATLAB closes automatically.

try, cd(fileparts(mfilename('fullpath'))); catch; end

fprintf('\n════════════════════════════════════════\n');
fprintf('  Stochastic Bellhop UQ Pipeline\n');
fprintf('  2026-05-26\n');
fprintf('════════════════════════════════════════\n\n');

%% ── PHASE 1: SCENARIO ANALYSIS ───────────────────────────────────────────
fprintf('--- PHASE 1: Scenario Analysis ---\n');

scenario_dirs = {'Scenarios/baseline', 'Scenarios/deep_water', ...
                 'Scenarios/shallow_water', 'Scenarios/upslope', 'Scenarios/downslope'};

for k = 1:numel(scenario_dirs)
    sd = scenario_dirs{k};
    fprintf('\n\n===== Running Scenario: %s =====\n', sd);
    
    fprintf('\n[1.%d.A] Delta Method...\n', k); tic;
    run(fullfile(sd,'Delta_Method','run_delta.m'));
    fprintf('  Done in %.0fs\n', toc);

    fprintf('\n[1.%d.B] Monte Carlo...\n', k); tic;
    run(fullfile(sd,'Monte_Carlo','run_MC.m'));
    fprintf('  Done in %.0fs\n', toc);

    fprintf('\n[1.%d.C] Comparison...\n', k); tic;
    run(fullfile(sd,'Comparison','run_comparison.m'));
    fprintf('  Done in %.0fs\n', toc);
end

%% ── PHASE 2: DISTRIBUTION STUDY ─────────────────────────────────────────
fprintf('\n\n--- PHASE 2: Distribution Study (~50-90 min) ---\n');
tic;
run('Distribution_Study/run_dist_study.m');
fprintf('  Done in %.0fs\n\n', toc);

%% ── PHASE 3: DECISION TREE / RANDOM FOREST ──────────────────────────────
fprintf('\n--- PHASE 3: Decision Tree / Random Forest ---\n');

fprintf('[3A] RF training...\n'); tic;
run('Decision_Tree/run_DT.m');
fprintf('  Done in %.0fs\n', toc);

fprintf('[3B] RF OOD analysis...\n'); tic;
run('Decision_Tree/run_DT_OOD.m');
fprintf('  Done in %.0fs\n', toc);

fprintf('[3C] RF explainability...\n'); tic;
run('Decision_Tree/run_DT_explain.m');
fprintf('  Done in %.0fs\n', toc);

fprintf('[3D] RF variance decomposition...\n'); tic;
run('Decision_Tree/run_DT_vardecomp.m');
fprintf('  Done in %.0fs\n\n', toc);

fprintf('\n════════════════════════════════════════\n');
fprintf('  Pipeline complete.\n');
fprintf('  All figures saved in figures/ subfolders.\n');
fprintf('════════════════════════════════════════\n');

exit;
