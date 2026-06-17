% run_overnight — batch runner for server / headless execution.
% Runs the full pipeline and logs any error to overnight_error.log.
try
    ROOT = fileparts(fileparts(mfilename('fullpath')));
    cd(ROOT);
    addpath(genpath(fullfile(ROOT, 'core')));
    addpath(fullfile(ROOT, 'pipeline'));
    addpath(fullfile(ROOT, 'binaries'));
    set(0, 'DefaultFigureVisible', 'off');

    run_mc;
    run_delta;
    run_pce;
    run_ln3;
    run_comparison;
    fill_findings;
catch ME
    fid = fopen(fullfile(ROOT, 'overnight_error.log'), 'w');
    fprintf(fid, '%s\n%s\n', ME.message, getReport(ME));
    fclose(fid);
    rethrow(ME);
end
