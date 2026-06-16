$proj = "c:\Users\orind\Stochastic-Bellhop\stochastic_TL_maps_project"
$stage1Pid = 41416
Write-Host "$(Get-Date) Watching for MATLAB PID $stage1Pid to finish..."
# Wait for stage 1 MATLAB to exit (check every 60s)
while ($true) {
    $proc = Get-Process -Id $stage1Pid -ErrorAction SilentlyContinue
    if (-not $proc) { break }
    Start-Sleep -Seconds 60
}
Write-Host "$(Get-Date) Stage 1 done. Starting Stage 2..."
Start-Process -FilePath "matlab" `
    -ArgumentList "-batch", "cd('$proj'); run_overnight_stage2" `
    -WorkingDirectory $proj `
    -RedirectStandardOutput "$proj\run_stage2_overnight.log" `
    -RedirectStandardError  "$proj\run_stage2_overnight_err.log" `
    -WindowStyle Hidden
Write-Host "$(Get-Date) Stage 2 MATLAB launched."
