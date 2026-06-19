# auto_push.ps1 — PostToolUse hook: auto-commit and push edited project files

$raw = $input | Out-String
try {
    $data = $raw | ConvertFrom-Json
    $fpath = $data.tool_input.file_path
} catch {
    exit 0
}

if (-not $fpath) { exit 0 }

# Skip .claude/ directory (memory, settings, hooks)
if ($fpath -match '[\\/]\.claude[\\/]') { exit 0 }

# Only commit recognised source files
if ($fpath -notmatch '\.(m|json|tex|ps1|md)$') { exit 0 }

# Must be inside a git repo
$null = git rev-parse --show-toplevel 2>$null
if (-not $?) { exit 0 }

# Stage the file
git add "$fpath" 2>$null

# Nothing staged — file unchanged or untracked outside repo
$staged = git diff --cached --name-only 2>$null
if (-not $staged) { exit 0 }

$fname = Split-Path -Leaf $fpath
git commit -m "Auto-save: $fname" 2>$null
git push 2>$null
