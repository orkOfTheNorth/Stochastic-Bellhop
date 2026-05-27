%% run_all.m  –  Master runner  (2026-05-26)
% Executes the full UQ pipeline in order:
%   1. Delta Method   (7 Bellhop runs)
%   2. Monte Carlo    (1400 Bellhop runs – slow, grab a coffee)
%   3. Comparison     (loads saved results, no Bellhop)
%   4. Opens viewer   (interactive dashboard)

try, cd(fileparts(mfilename('fullpath'))); catch; end

fprintf('\n========================================\n');
fprintf('  UQ Pipeline  –  2026-05-26\n');
fprintf('========================================\n\n');

%% Step 1 – Delta Method
fprintf('[STEP 1/4]  Delta Method...\n');
tic;
run('Delta_Method/run_delta.m');
fprintf('  Done in %.1f s\n\n', toc);

%% Step 2 – Monte Carlo
fprintf('[STEP 2/4]  Monte Carlo (N=200 x 7 subsets)...\n');
tic;
run('Monte_Carlo/run_MC.m');
fprintf('  Done in %.1f s\n\n', toc);

%% Step 3 – Comparison
fprintf('[STEP 3/4]  Comparison...\n');
tic;
run('Comparison/run_comparison.m');
fprintf('  Done in %.1f s\n\n', toc);

%% Step 4 – Viewer
fprintf('[STEP 4/4]  Opening viewer...\n');
run('viewer.m');

fprintf('\n========================================\n');
fprintf('  Pipeline complete.\n');
fprintf('========================================\n');
