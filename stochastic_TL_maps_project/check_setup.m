%% check_setup.m — Server environment verification.
%
% Run this first after cloning on a new machine (SSH server or local).
% Prints a pass/fail summary for every prerequisite.
%
% Usage:
%   matlab -batch "check_setup"
%   — or interactively: check_setup

try
    cd(fileparts(mfilename('fullpath')));
catch
end

addpath(genpath('Shared_Utils'));
addpath(genpath('Bellhop'));

PASS = true;
fprintf('\n====================================================\n');
fprintf('  check_setup — Stochastic Bellhop environment check\n');
fprintf('====================================================\n\n');

% ── MATLAB version ───────────────────────────────────────────────────────────
v    = version('-release');
year = str2double(v(1:4));
PASS = chk(PASS, 'MATLAB version >= R2021a', year >= 2021, ...
    sprintf('Found R%s — upgrade if parpool/toolbox features are missing.', v));

% ── Parallel Computing Toolbox ────────────────────────────────────────────────
PASS = chk(PASS, 'Parallel Computing Toolbox', ~isempty(ver('parallel')), ...
    'Not found — parfor will run sequentially (slower but correct).');

% ── config.json ───────────────────────────────────────────────────────────────
cfg_ok = isfile('config.json');
PASS   = chk(PASS, 'config.json present', cfg_ok, ...
    'Missing — check working directory.');
if cfg_ok
    try
        loadConfig();
        PASS = chk(PASS, 'config.json parseable', true, '');
    catch me
        PASS = chk(PASS, 'config.json parseable', false, me.message);
    end
end

% ── Shared_Utils on path ─────────────────────────────────────────────────────
PASS = chk(PASS, 'Shared_Utils on MATLAB path', ~isempty(which('loadConfig')), ...
    'Run: addpath(genpath(''Shared_Utils''))');

% ── Bellhop CPU binary ───────────────────────────────────────────────────────
bhp  = which('bellhop.exe');
PASS = chk(PASS, 'bellhop.exe on MATLAB path', ~isempty(bhp), ...
    sprintf('Not found.  Expected at: %s', fullfile(pwd, 'Bellhop', 'bellhop.exe')));
if ~isempty(bhp)
    fprintf('    Path: %s\n', bhp);
end

% ── Bellhop CUDA binary ──────────────────────────────────────────────────────
bhp_cuda      = which('bellhopcuda.exe');
has_cuda_bin  = ~isempty(bhp_cuda);
PASS = chk(PASS, 'bellhopcuda.exe on MATLAB path', has_cuda_bin, ...
    'CUDA build not found — GPU acceleration will be skipped (CPU fallback active).');
if has_cuda_bin
    fprintf('    Path: %s\n', bhp_cuda);
end

% ── GPU ───────────────────────────────────────────────────────────────────────
try
    g      = gpuDevice();
    cc     = str2double(g.ComputeCapability);
    gpu_ok = cc >= 3.5;
    PASS   = chk(PASS, sprintf('GPU compute capability >= 3.5 (found %.1f)', cc), gpu_ok, ...
        'GPU present but compute capability too low for CUDA bellhop.');
    if gpu_ok
        fprintf('    GPU: %s  (%.1f GB free / %.1f GB total)\n', ...
            g.Name, g.AvailableMemory/1e9, g.TotalMemory/1e9);
    end
catch
    PASS = chk(PASS, 'GPU available', false, ...
        'No GPU detected — CPU-only mode (bellhop.exe will be used).');
end

% ── Cache and Methods directories ────────────────────────────────────────────
PASS = chk(PASS, 'Cache/ directory writable',   canWriteDir('Cache'), ...
    'Create Cache/ or fix permissions.');
PASS = chk(PASS, 'Methods/ directory writable', canWriteDir('Methods'), ...
    'Create Methods/ or fix permissions.');

% ── LaTeX for findings.tex ───────────────────────────────────────────────────
[s, ~] = system('pdflatex --version');
PASS   = chk(PASS, 'pdflatex available (for findings.tex)', s == 0, ...
    'Install MiKTeX or TeX Live to compile findings.tex.');

% ── Summary ──────────────────────────────────────────────────────────────────
fprintf('\n====================================================\n');
if PASS
    fprintf('  ALL CHECKS PASSED — ready to run.\n');
else
    fprintf('  SOME CHECKS FAILED — see warnings above.\n');
end
fprintf('====================================================\n\n');


% ─────── helpers ──────────────────────────────────────────────────────────────
function pass_out = chk(pass_in, label, ok, hint)
    if ok
        fprintf('  [OK]   %s\n', label);
        pass_out = pass_in;
    else
        fprintf('  [FAIL] %s\n', label);
        if ~isempty(hint)
            fprintf('         Hint: %s\n', hint);
        end
        pass_out = false;
    end
end

function ok = canWriteDir(d)
    if ~exist(d, 'dir'), mkdir(d); end
    tmp = fullfile(d, ['.write_test_' num2str(randi(1e6))]);
    try
        fid = fopen(tmp, 'w'); fclose(fid); delete(tmp);
        ok = true;
    catch
        ok = false;
    end
end
