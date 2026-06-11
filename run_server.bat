@echo off
REM ============================================================
REM  run_server.bat — Launch Stochastic-Bellhop on Windows server
REM
REM  Usage (from repo root, via SSH or locally):
REM    run_server.bat                      -> run_analysis_only (default)
REM    run_server.bat run_analysis_only    -> PCE + LN3 + Comparison on cached MC
REM    run_server.bat run_full_pipeline    -> MC + PCE + LN3 + Comparison (long)
REM    run_server.bat check_setup          -> verify environment before running
REM    run_server.bat run_Delta            -> Delta UQ pass only
REM    run_server.bat run_MC               -> Monte Carlo pass only
REM    run_server.bat run_pce              -> PCE pass only
REM    run_server.bat run_LN3              -> LN3 pass only
REM    run_server.bat run_Comparison       -> Comparison figures only
REM    run_server.bat fill_findings        -> print LaTeX values for findings.tex
REM
REM  Output log: stochastic_TL_maps_project\run_<MODE>.log
REM ============================================================

setlocal EnableDelayedExpansion

set MODE=%~1
if "%MODE%"=="" set MODE=run_analysis_only

set PROJDIR=%~dp0stochastic_TL_maps_project
set LOGFILE=%PROJDIR%\run_%MODE%.log

echo [run_server] mode    = %MODE%
echo [run_server] projdir = %PROJDIR%
echo [run_server] log     = %LOGFILE%
echo.

REM -- Validate mode -------------------------------------------------------
if "%MODE%"=="run_analysis_only" goto :RUN
if "%MODE%"=="run_full_pipeline" goto :RUN
if "%MODE%"=="check_setup"       goto :RUN
if "%MODE%"=="run_Delta"         goto :RUN
if "%MODE%"=="run_MC"            goto :RUN
if "%MODE%"=="run_pce"           goto :RUN
if "%MODE%"=="run_LN3"           goto :RUN
if "%MODE%"=="run_Comparison"    goto :RUN
if "%MODE%"=="fill_findings"     goto :RUN

echo [run_server] ERROR: unknown mode "%MODE%"
echo   Valid modes: run_analysis_only  run_full_pipeline  check_setup
echo               run_Delta  run_MC  run_pce  run_LN3  run_Comparison  fill_findings
exit /b 1

:RUN
matlab -nosplash -nodesktop ^
  -logfile "%LOGFILE%" ^
  -batch "cd('%PROJDIR%'); addpath(genpath('Shared_Utils')); addpath(genpath('Bellhop')); %MODE%"

if %ERRORLEVEL% neq 0 (
    echo.
    echo [run_server] MATLAB exited with error %ERRORLEVEL% — check %LOGFILE%
    exit /b %ERRORLEVEL%
)

echo.
echo [run_server] Done.  Log: %LOGFILE%
endlocal
