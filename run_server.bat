@echo off
REM ============================================================
REM  run_server.bat — Launch Stochastic-Bellhop on Windows server
REM
REM  Usage (from repo root, via SSH or locally):
REM    run_server.bat                      -> run_all_v2 (full v2 pipeline, default)
REM    run_server.bat run_all_v2           -> MC (IS/300) + Delta + PCE + LN3 + Comparison
REM                                          + MC_convergence + fill_findings
REM    run_server.bat check_setup          -> verify environment before running
REM    run_server.bat run_MC               -> Monte Carlo IS pass only
REM    run_server.bat run_Delta            -> Delta GH orders 1-10 only
REM    run_server.bat run_pce              -> PCE (order 20, LOO-CV) only
REM    run_server.bat run_LN3              -> LN3 histograms + RMSE maps only
REM    run_server.bat run_Comparison       -> Comparison detection maps only
REM    run_server.bat run_MC_convergence   -> MC convergence diagnostics only
REM    run_server.bat fill_findings        -> patch findings.tex [FILL] markers
REM
REM  Legacy modes (old N=50 pipeline, still callable):
REM    run_server.bat run_analysis_only    -> PCE + LN3 + Comparison on cached MC
REM    run_server.bat run_full_pipeline    -> legacy full pipeline
REM
REM  Output log: stochastic_TL_maps_project\run_<MODE>.log
REM ============================================================

setlocal EnableDelayedExpansion

set MODE=%~1
if "%MODE%"=="" set MODE=run_all_v2

set PROJDIR=%~dp0stochastic_TL_maps_project
set LOGFILE=%PROJDIR%\run_%MODE%.log

echo [run_server] mode    = %MODE%
echo [run_server] projdir = %PROJDIR%
echo [run_server] log     = %LOGFILE%
echo.

REM -- Validate mode -------------------------------------------------------
if "%MODE%"=="run_all_v2"          goto :RUN
if "%MODE%"=="check_setup"         goto :RUN
if "%MODE%"=="run_MC"              goto :RUN
if "%MODE%"=="run_Delta"           goto :RUN
if "%MODE%"=="run_pce"             goto :RUN
if "%MODE%"=="run_LN3"             goto :RUN
if "%MODE%"=="run_Comparison"      goto :RUN
if "%MODE%"=="run_MC_convergence"  goto :RUN
if "%MODE%"=="fill_findings"       goto :RUN
if "%MODE%"=="run_analysis_only"   goto :RUN
if "%MODE%"=="run_full_pipeline"   goto :RUN

echo [run_server] ERROR: unknown mode "%MODE%"
echo   Valid modes: run_all_v2  check_setup
echo               run_MC  run_Delta  run_pce  run_LN3  run_Comparison
echo               run_MC_convergence  fill_findings
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
