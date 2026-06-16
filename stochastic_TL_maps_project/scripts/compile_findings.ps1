# Compile findings.tex → findings.pdf
# Run from the project directory or right-click → "Run with PowerShell" in Explorer

Set-Location $PSScriptRoot

Write-Host "Pass 1/2: compiling findings.tex ..." -ForegroundColor Cyan
pdflatex -interaction=nonstopmode findings.tex

if ($LASTEXITCODE -eq 0) {
    Write-Host "Pass 2/2: resolving cross-references ..." -ForegroundColor Cyan
    pdflatex -interaction=nonstopmode findings.tex
}

if (Test-Path "$PSScriptRoot\findings.pdf") {
    Write-Host "Success: findings.pdf created/updated" -ForegroundColor Green
    # Uncomment the next line to open the PDF automatically after compiling:
    # Start-Process "$PSScriptRoot\findings.pdf"
} else {
    Write-Host "Compilation failed — check findings.log for errors" -ForegroundColor Red
}
