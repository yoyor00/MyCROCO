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

if [[ -n "${XIOS_NTIMES:-}" ]]; then
    sed -i "s/^  ntimes = .*/  ntimes = ${XIOS_NTIMES}/" "$RUN_DIR/croco.nml"
fi

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

# Use the local BASIN workaround without modifying XIOS/ or TEST_CASES/.
PATCHED_FILE_DEF="$HERE/patches/XIOS/file_def_croco.xml_full_withcpp"
if [[ -f "$PATCHED_FILE_DEF" ]]; then
    echo "Using patched XIOS file definition: $PATCHED_FILE_DEF"
    cpp -P -traditional -imacros "$RUN_DIR/cppdefs.h" \
        "$PATCHED_FILE_DEF" > "$RUN_DIR/file_def_croco.xml"
fi

USING_SERVER="${XIOS_USING_SERVER:-false}"
OUTPUT_FREQ="${XIOS_OUTPUT_FREQ:-}"
sed -i "s|<variable id=\"using_server\" type=\"bool\">[^<]*</variable>|<variable id=\"using_server\" type=\"bool\">${USING_SERVER}</variable>|" \
    "$RUN_DIR/iodef.xml"

# The patched file uses 9600s for BASIN.  Keep this override for experiments.
if [[ -n "$OUTPUT_FREQ" ]]; then
    sed -i "s/output_freq=\"9600s\"/output_freq=\"${OUTPUT_FREQ}\"/g" \
        "$RUN_DIR/file_def_croco.xml"
fi

if [[ "$USING_SERVER" == "false" ]]; then
    printf -- '-np 4 ./croco\n' > "$RUN_DIR/xios_launch_mini.file"
else
    cat > "$RUN_DIR/xios_launch_mini.file" <<'LAUNCH'
-np 4 ./croco
-np 1 ./xios_server.exe
LAUNCH
fi

echo "Starting XIOS runtime smoke test in $RUN_DIR"
ulimit -s unlimited 2>/dev/null || true
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
