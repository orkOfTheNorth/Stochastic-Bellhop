#!/usr/bin/env bash
# ============================================================
#  run_server.sh — Launch Stochastic-Bellhop on Linux / macOS
#
#  Usage (from repo root):
#    ./run_server.sh                      # full v2 pipeline (default)
#    ./run_server.sh run_all_v2           # same as above
#    ./run_server.sh check_setup          # verify environment
#    ./run_server.sh run_MC               # Monte Carlo IS pass only
#    ./run_server.sh run_Delta            # Delta GH sweep only
#    ./run_server.sh run_pce              # PCE only
#    ./run_server.sh run_LN3              # LN3 RMSE maps only
#    ./run_server.sh run_Comparison       # Comparison figures only
#    ./run_server.sh run_MC_convergence   # MC convergence diagnostics only
#    ./run_server.sh fill_findings        # patch findings.tex [FILL] markers
#
#  Log: stochastic_TL_maps_project/run_<MODE>.log
# ============================================================

set -euo pipefail

MODE="${1:-run_all_v2}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJDIR="$SCRIPT_DIR/stochastic_TL_maps_project"
LOGFILE="$PROJDIR/run_${MODE}.log"

# ── Validate mode ─────────────────────────────────────────────────────────────
valid_modes=(
    run_all_v2
    check_setup
    run_MC
    run_Delta
    run_pce
    run_LN3
    run_Comparison
    run_MC_convergence
    fill_findings
    run_analysis_only
    run_full_pipeline
)

valid=0
for m in "${valid_modes[@]}"; do
    if [[ "$MODE" == "$m" ]]; then valid=1; break; fi
done

if [[ $valid -eq 0 ]]; then
    echo "[run_server] ERROR: unknown mode '$MODE'"
    echo "  Valid modes: ${valid_modes[*]}"
    exit 1
fi

# ── Locate MATLAB ─────────────────────────────────────────────────────────────
if ! command -v matlab &>/dev/null; then
    echo "[run_server] ERROR: 'matlab' not found on PATH."
    echo "  Add MATLAB to PATH, e.g.:"
    echo "    export PATH=\$PATH:/usr/local/MATLAB/R2024a/bin"
    echo "  Or call this script as:"
    echo "    MATLAB_BIN=/path/to/matlab ./run_server.sh"
    exit 1
fi
MATLAB_BIN="${MATLAB_BIN:-matlab}"

echo "[run_server] mode    = $MODE"
echo "[run_server] projdir = $PROJDIR"
echo "[run_server] log     = $LOGFILE"
echo "[run_server] matlab  = $(command -v "$MATLAB_BIN")"
echo ""

# ── Run MATLAB ────────────────────────────────────────────────────────────────
"$MATLAB_BIN" -nosplash -nodesktop \
    -logfile "$LOGFILE" \
    -batch "cd('$PROJDIR'); addpath(genpath('Shared_Utils')); addpath(genpath('Bellhop')); $MODE"

EXIT_CODE=$?
echo ""
if [[ $EXIT_CODE -ne 0 ]]; then
    echo "[run_server] MATLAB exited with error $EXIT_CODE — check $LOGFILE"
    exit $EXIT_CODE
fi

echo "[run_server] Done.  Log: $LOGFILE"
