%% run_all_v2.m — Full pipeline orchestrator (300-sample IS redesign)
%
% Execution order:
%   run_MC            — 300-sample Normal_10pct base + IS reweighting
%   run_Delta         — GH orders 1–10 convergence
%   run_pce           — PCE up to order 20, combined LOO-CV
%   run_LN3           — LN3 histograms + pixel-wise CDF RMSE maps
%   run_Comparison    — Inverted detection-probability comparison maps
%   run_MC_convergence— EX/Var convergence scatter, mean, derivative
%   fill_findings     — Auto-patch findings.tex [FILL] markers

STEPS = { ...
    'run_MC',             'MC (300-sample IS engine)'; ...
    'run_Delta',          'Delta GH orders 1-10'; ...
    'run_pce',            'PCE order up to 20 + LOO-CV'; ...
    'run_LN3',            'LN3 histograms + RMSE maps'; ...
    'run_Comparison',     'Comparison detection maps'; ...
    'run_MC_convergence', 'MC convergence diagnostics'; ...
    'fill_findings',      'Patch findings.tex'; ...
};

PROJECT_ROOT = fileparts(mfilename('fullpath'));
if isempty(PROJECT_ROOT), PROJECT_ROOT = pwd; end
cd(PROJECT_ROOT);
setappdata(0, 'project_root', PROJECT_ROOT);

fprintf('\n');
fprintf('====================================================\n');
fprintf('  Stochastic Bellhop UQ Pipeline v2 (IS/300)\n');
fprintf('====================================================\n\n');

n_steps = size(STEPS, 1);
elapsed = zeros(1, n_steps);

for si = 1:n_steps
    fn   = STEPS{si, 1};
    desc = STEPS{si, 2};
    fprintf('[%d/%d] %s...\n', si, n_steps, desc);
    R = getappdata(0, 'project_root');
    fpath = fullfile(R, [fn '.m']);
    t0 = tic;
    try
        run(fpath);
    catch ME
        fprintf('  !! ERROR in %s: %s\n', fn, ME.message);
    end
    elapsed(si) = toc(t0);
    R = getappdata(0, 'project_root');
    cd(R);
    fprintf('[%d/%d] Done in %.0fs\n\n', si, n_steps, elapsed(si));
end

fprintf('====================================================\n');
fprintf('  Pipeline complete.  Total: %.1f min\n', sum(elapsed)/60);
fprintf('  Step breakdown:\n');
for si = 1:n_steps
    fprintf('    %-28s %6.0f s\n', STEPS{si,1}, elapsed(si));
end
fprintf('====================================================\n');
