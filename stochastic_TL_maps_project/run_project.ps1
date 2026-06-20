# run_project.ps1 — One-command full pipeline launcher.
#
# Usage (from stochastic_TL_maps_project\):
#   powershell -ExecutionPolicy Bypass -File run_project.ps1
#
# What it does:
#   1. Runs the full MATLAB pipeline  (run_overnight)
#   2. Compiles findings.tex → findings.pdf  (pdflatex x2)
#   3. Sends a completion email to orindig0106@gmail.com
#
# One-time setup for email:
#   - Enable Gmail "App Password" at https://myaccount.google.com/apppasswords
#   - Paste the 16-char app password into $GmailAppPassword below

param(
    [string]$GmailAppPassword = "PASTE_APP_PASSWORD_HERE"
)

$ROOT    = Split-Path -Parent $MyInvocation.MyCommand.Path
$LOG     = Join-Path $ROOT "overnight_run.log"
$TEX_DIR = Join-Path $ROOT "docs"
$NOTIFY  = "orindig0106@gmail.com"

Set-Location $ROOT

# ── Step 1: Run MATLAB pipeline ───────────────────────────────────────────────
Write-Host "`n=== Launching MATLAB pipeline ===" -ForegroundColor Cyan
Write-Host "Log: $LOG"

$matlabCmd = "cd('$ROOT'); addpath(genpath(fullfile('$ROOT','core'))); addpath(fullfile('$ROOT','pipeline')); addpath(fullfile('$ROOT','binaries')); run_overnight; exit"

$proc = Start-Process -FilePath "matlab" `
    -ArgumentList "-batch", $matlabCmd `
    -WorkingDirectory $ROOT `
    -PassThru -Wait

$matlabExit = $proc.ExitCode
if ($matlabExit -ne 0) {
    Write-Host "MATLAB exited with code $matlabExit — check $LOG" -ForegroundColor Yellow
} else {
    Write-Host "MATLAB pipeline complete." -ForegroundColor Green
}

# ── Step 2: Compile findings.tex ─────────────────────────────────────────────
Write-Host "`n=== Compiling findings.tex ===" -ForegroundColor Cyan

$pdfOK = $false
if (Get-Command pdflatex -ErrorAction SilentlyContinue) {
    # Run from project root so \graphicspath{} entries resolve relative to ROOT, not docs/
    # -output-directory sends findings.pdf / .aux / .log into docs/
    pdflatex -interaction=nonstopmode -output-directory docs docs/findings.tex | Out-Null
    pdflatex -interaction=nonstopmode -output-directory docs docs/findings.tex | Out-Null
    if (Test-Path (Join-Path $TEX_DIR "findings.pdf")) {
        Write-Host "findings.pdf compiled successfully." -ForegroundColor Green
        $pdfOK = $true
    } else {
        Write-Host "PDF compilation failed — check docs\findings.log" -ForegroundColor Yellow
    }
} else {
    Write-Host "pdflatex not found — skipping PDF compilation." -ForegroundColor Yellow
}

Set-Location $ROOT

# ── Step 3: Send email notification ──────────────────────────────────────────
Write-Host "`n=== Sending completion email ===" -ForegroundColor Cyan

$subject = "Stochastic Bellhop pipeline complete"
$body    = @"
Pipeline finished.

MATLAB exit code : $matlabExit
PDF compiled     : $pdfOK
Log file         : $LOG

Results are in:
  $ROOT\Methods\
  $ROOT\docs\findings.pdf
"@

try {
    $cred = New-Object System.Management.Automation.PSCredential(
        $NOTIFY,
        (ConvertTo-SecureString $GmailAppPassword -AsPlainText -Force)
    )
    Send-MailMessage `
        -From    $NOTIFY `
        -To      $NOTIFY `
        -Subject $subject `
        -Body    $body `
        -SmtpServer "smtp.gmail.com" `
        -Port 587 `
        -UseSsl `
        -Credential $cred
    Write-Host "Email sent to $NOTIFY" -ForegroundColor Green
} catch {
    Write-Host "Email failed: $_" -ForegroundColor Yellow
    Write-Host "Pipeline still complete — check $LOG manually." -ForegroundColor White
}

Write-Host "`n=== All done ===" -ForegroundColor Green
