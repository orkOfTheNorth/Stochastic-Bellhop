%% run_all.m — Master pipeline runner
%
% Each sub-script starts with 'clear', so timing uses the global tic/toc
% (no output arg) which is unaffected by 'clear'.  Per-step elapsed is
% stashed in root appdata which also survives 'clear'.

PROJECT_ROOT = fileparts(mfilename('fullpath'));
if isempty(PROJECT_ROOT)
    PROJECT_ROOT = pwd;
end
cd(PROJECT_ROOT);

% Store root in app data so it survives 'clear' inside sub-scripts
setappdata(0, 'project_root',  PROJECT_ROOT);
setappdata(0, 'pipe_elapsed',  zeros(1, 5));

fprintf('\n');
fprintf('====================================================\n');
fprintf('   Stochastic Bellhop UQ Pipeline -- Full Run\n');
fprintf('====================================================\n\n');

% ── Step 1: Analytical sweep ──────────────────────────────────────────────
fprintf('[1/5] Analytical sweep...\n');
R = getappdata(0, 'project_root');
tic;  run(fullfile(R, 'analytical', 'run_analytical_sweep.m'));
R = getappdata(0, 'project_root');  cd(R);
e = getappdata(0, 'pipe_elapsed');  e(1) = toc;
setappdata(0, 'pipe_elapsed', e);
fprintf('[1/5] Done in %.0fs\n\n', e(1));

% ── Step 2: Delta Method ──────────────────────────────────────────────────
fprintf('[2/5] Delta Method (all scenarios x distributions)...\n');
R = getappdata(0, 'project_root');
tic;  run(fullfile(R, 'run_Delta.m'));
R = getappdata(0, 'project_root');  cd(R);
e = getappdata(0, 'pipe_elapsed');  e(2) = toc;
setappdata(0, 'pipe_elapsed', e);
fprintf('[2/5] Done in %.0fs\n\n', e(2));

% ── Step 3: Monte Carlo ───────────────────────────────────────────────────
fprintf('[3/5] Monte Carlo (N=50 per subset -- hours for new scenarios)...\n');
R = getappdata(0, 'project_root');
tic;  run(fullfile(R, 'run_MC.m'));
R = getappdata(0, 'project_root');  cd(R);
e = getappdata(0, 'pipe_elapsed');  e(3) = toc;
setappdata(0, 'pipe_elapsed', e);
fprintf('[3/5] Done in %.0fs\n\n', e(3));

% ── Step 4: LN3 diagnostic ───────────────────────────────────────────────
fprintf('[4/5] LN3 diagnostics...\n');
R = getappdata(0, 'project_root');
tic;  run(fullfile(R, 'run_LN3.m'));
R = getappdata(0, 'project_root');  cd(R);
e = getappdata(0, 'pipe_elapsed');  e(4) = toc;
setappdata(0, 'pipe_elapsed', e);
fprintf('[4/5] Done in %.0fs\n\n', e(4));

% ── Step 5: Comparison ────────────────────────────────────────────────────
fprintf('[5/5] Comparison...\n');
R = getappdata(0, 'project_root');
tic;  run(fullfile(R, 'run_Comparison.m'));
R = getappdata(0, 'project_root');  cd(R);
e = getappdata(0, 'pipe_elapsed');  e(5) = toc;
setappdata(0, 'pipe_elapsed', e);
fprintf('[5/5] Done in %.0fs\n\n', e(5));

fprintf('====================================================\n');
fprintf('   Pipeline complete.  Total: %.1f min\n', sum(e)/60);
fprintf('====================================================\n');
