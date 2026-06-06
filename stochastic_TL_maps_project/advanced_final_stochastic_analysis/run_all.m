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

%% ── PHASE 1: BASELINE SCENARIO ──────────────────────────────────────────
fprintf('--- PHASE 1: Baseline (const_35, 50km) ---\n');

fprintf('[1A] Delta Method (~1 min)...\n'); tic;
run('Delta_Method/run_delta.m');
fprintf('  Done in %.0fs\n\n', toc);

fprintf('[1B] Monte Carlo N=50 (~20-30 min)...\n'); tic;
run('Monte_Carlo/run_MC.m');
fprintf('  Done in %.0fs\n\n', toc);

fprintf('[1C] LN3 grid fits (~5-10 min)...\n'); tic;
run('LN3/run_LN3.m');
fprintf('  Done in %.0fs\n\n', toc);

fprintf('[1D] Comparison...\n'); tic;
run('Comparison/run_comparison.m');
fprintf('  Done in %.0fs\n\n', toc);

%% ── PHASE 2: NEW SCENARIOS ───────────────────────────────────────────────
fprintf('--- PHASE 2: New Scenarios ---\n');

scenario_dirs = {'Scenarios/deep_water', 'Scenarios/shallow_water', ...
                 'Scenarios/upslope',    'Scenarios/downslope'};
scenario_names = {'Deep Water','Shallow Water','Upslope','Downslope'};

for k = 1:numel(scenario_dirs)
    sd = scenario_dirs{k};
    fprintf('\n[2%s] %s – Delta...\n', char('A'+k-1), scenario_names{k}); tic;
    run(fullfile(sd,'run_delta.m'));
    fprintf('  Done in %.0fs\n', toc);

    fprintf('[2%s] %s – MC (~15-50 min)...\n', char('A'+k-1), scenario_names{k}); tic;
    run(fullfile(sd,'run_MC.m'));
    fprintf('  Done in %.0fs\n', toc);

    fprintf('[2%s] %s – Comparison...\n', char('A'+k-1), scenario_names{k}); tic;
    run(fullfile(sd,'run_comparison.m'));
    fprintf('  Done in %.0fs\n', toc);
end

%% ── PHASE 3: DISTRIBUTION STUDY ─────────────────────────────────────────
fprintf('\n--- PHASE 3: Distribution Study (~50-90 min) ---\n');
tic;
run('Distribution_Study/run_dist_study.m');
fprintf('  Done in %.0fs\n\n', toc);

%% ── PHASE 4: DECISION TREE / RANDOM FOREST ──────────────────────────────
fprintf('\n--- PHASE 4: Decision Tree / Random Forest ---\n');

fprintf('[4A] RF training...\n'); tic;
run('Decision_Tree/run_DT.m');
fprintf('  Done in %.0fs\n', toc);

fprintf('[4B] RF OOD analysis...\n'); tic;
run('Decision_Tree/run_DT_OOD.m');
fprintf('  Done in %.0fs\n', toc);

fprintf('[4C] RF explainability...\n'); tic;
run('Decision_Tree/run_DT_explain.m');
fprintf('  Done in %.0fs\n', toc);

fprintf('[4D] RF variance decomposition...\n'); tic;
run('Decision_Tree/run_DT_vardecomp.m');
fprintf('  Done in %.0fs\n\n', toc);

fprintf('\n════════════════════════════════════════\n');
fprintf('  Pipeline complete.\n');
fprintf('  All figures saved in figures/ subfolders.\n');
fprintf('════════════════════════════════════════\n');

exit;
