% run_overnight — full pipeline runner for server / headless execution.
%
% Steps (in order):
%   1.  run_MC                    — Monte Carlo (N=1000, LHS+MoM + IID+MLE, SVP ensemble)
%   2.  run_Delta                 — Delta method, 1st + 2nd order, GH orders 1-10
%   3.  run_pce                   — PCE up to order 20, LOO-CV
%   4.  run_LN3                   — LN3 histograms + pixel-wise CDF RMSE maps
%   5.  run_Comparison            — Detection-probability comparison maps
%   6.  run_MC_convergence        — MC estimator convergence diagnostics
%   7.  run_lhs_iid_convergence   — LHS+MoM vs IID+MLE convergence study
%   8.  plot_delta_order_comparison — 1st vs 2nd order correction visualisation
%   9.  plot_taylor_order_study   — PCE/Taylor order study plots
%  10.  fill_findings             — Auto-patch findings.tex [FILL] markers
%
% Errors in any step are logged but do not abort subsequent steps.
% Summary printed to console and overnight_run.log.

ROOT = fileparts(fileparts(mfilename('fullpath')));
cd(ROOT);
addpath(genpath(fullfile(ROOT, 'core')));
addpath(fullfile(ROOT, 'pipeline'));
addpath(fullfile(ROOT, 'scripts'));
addpath(fullfile(ROOT, 'binaries'));
set(0, 'DefaultFigureVisible', 'off');

STEPS = {
    'run_MC',                   'function',  'MC (N=1000, LHS+IID+SVP)';
    'run_Delta',                'function',  'Delta 1st+2nd order, GH sweep';
    'run_pce',                  'function',  'PCE order 1-20 + LOO-CV';
    'run_LN3',                  'function',  'LN3 histograms + RMSE maps';
    'run_Comparison',           'function',  'Comparison detection maps';
    'run_MC_convergence',       'function',  'MC estimator convergence';
    'run_lhs_iid_convergence',  'function',  'LHS+MoM vs IID+MLE convergence';
    'plot_delta_order_comparison', 'script', '1st vs 2nd order correction plots';
    'plot_taylor_order_study',  'script',    'PCE/Taylor order study plots';
    'fill_findings',            'script',    'Patch findings.tex';
};

n = size(STEPS, 1);
elapsed = zeros(1, n);
failed  = false(1, n);

log_path = fullfile(ROOT, 'overnight_run.log');
flog = fopen(log_path, 'w');
logboth = @(varargin) deal(fprintf(varargin{:}), fprintf(flog, varargin{:}));

logboth('\n=== Stochastic Bellhop Full Pipeline === %s\n', datestr(now));

for si = 1:n
    name  = STEPS{si, 1};
    kind  = STEPS{si, 2};
    desc  = STEPS{si, 3};
    logboth('\n[%d/%d] %s\n', si, n, desc);
    t0 = tic;
    try
        if strcmp(kind, 'function')
            feval(name);
        else
            % script — locate in pipeline/ or scripts/
            spath = fullfile(ROOT, 'pipeline', [name '.m']);
            if ~isfile(spath)
                spath = fullfile(ROOT, 'scripts', [name '.m']);
            end
            run(spath);
        end
    catch ME
        logboth('  !! ERROR: %s\n', ME.message);
        fprintf(flog, '%s\n', getReport(ME));
        failed(si) = true;
    end
    elapsed(si) = toc(t0);
    cd(ROOT);
    logboth('[%d/%d] Done in %.0f s\n', si, n, elapsed(si));
end

logboth('\n=== Pipeline complete. Total: %.1f min ===\n', sum(elapsed)/60);
for si = 1:n
    status = 'ok  ';
    if failed(si), status = 'FAIL'; end
    logboth('  %s  %-44s %6.0f s\n', status, STEPS{si,3}, elapsed(si));
end
fclose(flog);
fprintf('Log saved: %s\n', log_path);
