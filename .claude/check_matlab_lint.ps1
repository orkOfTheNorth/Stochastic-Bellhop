# check_matlab_lint.ps1 — PostToolUse hook: run MATLAB checkcode on edited .m files
# Called by Claude Code after every Edit or Write tool use.
# Receives JSON on stdin with tool_input.file_path.

$raw = $input | Out-String
try {
    $data = $raw | ConvertFrom-Json
    $fpath = $data.tool_input.file_path
} catch {
    exit 0
}

if (-not $fpath) { exit 0 }
if ($fpath -notmatch '\.m$') { exit 0 }

# Normalize path separators for MATLAB
$mpath = $fpath -replace '\\', '/'

Write-Host "`n[mlint] $fpath"

$result = & matlab -batch "try; msgs = checkcode('$mpath','-id'); for k=1:numel(msgs); fprintf('[%s] line %d: %s\n', msgs(k).id, msgs(k).line, msgs(k).message); end; catch ME; fprintf('checkcode error: %s\n', ME.message); end" 2>&1

if ($result) {
    Write-Host ($result | Out-String)
} else {
    Write-Host "[mlint] No issues found."
}
