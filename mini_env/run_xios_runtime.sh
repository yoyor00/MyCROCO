#!/bin/bash
# Experimental XIOS runtime smoke test.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$HERE/.." && pwd)"
# shellcheck source=env.sh
source "$HERE/env.sh"

RUN_DIR="${1:-$HERE/work/xios}"
[[ -d "$RUN_DIR" ]] || { echo "missing run directory: $RUN_DIR" >&2; exit 1; }
RUN_DIR="$(cd "$RUN_DIR" && pwd)"
[[ -x "$RUN_DIR/croco" ]] || { echo "missing $RUN_DIR/croco; run ./run_checks.sh first" >&2; exit 1; }
[[ -f "$RUN_DIR/croco.nml" ]] || { echo "missing $RUN_DIR/croco.nml" >&2; exit 1; }

# Build a run-local copy of process_xios_xml.sh with absolute paths.
sed \
    -e "s|^XIOS_NAM_DIR=.*|XIOS_NAM_DIR=$RUN_DIR|" \
    -e "s|^ROOT_DIR=.*|ROOT_DIR=$REPO_ROOT|" \
    "$REPO_ROOT/XIOS/process_xios_xml.sh" > "$RUN_DIR/process_xios_xml.sh"
chmod +x "$RUN_DIR/process_xios_xml.sh"

(
  cd "$RUN_DIR"
  OCE_XIOS=TRUE ATM_XIOS=FALSE USE_OASIS=FALSE ./process_xios_xml.sh
) > "$HERE/logs/xios-xml.log" 2>&1

cat > "$RUN_DIR/xios_launch_mini.file" <<'LAUNCH'
-np 4 ./croco
-np 1 ./xios_server.exe
LAUNCH

echo "Starting XIOS runtime smoke test in $RUN_DIR"
(
  cd "$RUN_DIR"
  timeout "${XIOS_RUN_TIMEOUT:-120}" \
      mpirun --allow-run-as-root --app xios_launch_mini.file
) > "$HERE/logs/xios-runtime.log" 2>&1 || {
    code=$?
    echo "XIOS runtime FAILED (exit $code); see $HERE/logs/xios-runtime.log" >&2
    tail -40 "$HERE/logs/xios-runtime.log" >&2 || true
    exit "$code"
}

if compgen -G "$RUN_DIR/*.nc" > /dev/null; then
    echo "XIOS runtime PASS"
else
    echo "XIOS runtime finished but produced no NetCDF files" >&2
    exit 1
fi
