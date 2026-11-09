#!/bin/bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$HERE/.." && pwd)"
# shellcheck source=env.sh
source "$HERE/env.sh"

WORK="$HERE/work"
LOGS="$HERE/logs"
JOBS="${JOBS:-$(nproc)}"
BACKENDS="${MINI_ENV_BACKENDS:-make scons}"
mkdir -p "$WORK" "$LOGS"

prepare_case() {
    local name="$1" key="$2"
    local run="$WORK/$name"
    rm -rf "$run"
    mkdir -p "$run"
    cp "$REPO_ROOT/TEST_CASES/BASIN/cppdefs_BASIN.h" "$run/cppdefs_BASIN.h"
    cp "$REPO_ROOT/TEST_CASES/BASIN/param_BASIN.h" "$run/param.h"
    cp "$REPO_ROOT/TEST_CASES/BASIN/croco_BASIN.nml" "$run/croco.nml"
    sed -i 's/^# undef  MPI$/# define MPI/' "$run/cppdefs_BASIN.h"
    if [[ "$key" == "XIOS" ]]; then
        sed -i 's|^#include "cppdefs_dev.h"|# define XIOS\n# define XIOS2\n# define SPONGE\n# define NO_SPONGE_GRID\n# define SALINITY\n# define ANA_BSFLUX\n# define CURVGRID\n# define SPHERICAL\n# define MASKING\n# define PSOURCE\n#include "cppdefs_dev.h"|' "$run/cppdefs_BASIN.h"
    elif [[ "$key" == "OASIS" ]]; then
        sed -i 's|^#include "cppdefs_dev.h"|# define OA_COUPLING\n# define CURVGRID\n# define SPHERICAL\n# define MASKING\n# define PSOURCE\n#include "cppdefs_dev.h"|' "$run/cppdefs_BASIN.h"
    fi
    printf '#include "cppdefs_BASIN.h"\n' > "$run/cppdefs.h"
    echo "$run"
}

run_backend() {
    local label="$1" run="$2" backend="$3"
    shift 3
    local log="$LOGS/${label}-${backend}.log"
    local cmd=("$REPO_ROOT/OCEAN/jobcomp" --src "$REPO_ROOT/OCEAN" --jobs "$JOBS" "$@")
    if [[ "$backend" == "scons" ]]; then
        cmd+=(--scons)
    fi
    echo "=== $label / $backend ==="
    if (cd "$run" && "${cmd[@]}") > "$log" 2>&1; then
        if [[ -x "$run/croco" ]]; then
            echo "PASS: $label / $backend"
            return 0
        fi
        echo "FAIL: $label / $backend (no croco executable)" >&2
        return 1
    fi
    echo "FAIL: $label / $backend (see $log)" >&2
    tail -40 "$log" >&2
    return 1
}

status_xios=0
status_oasis=0

xios_run="$(prepare_case xios XIOS)"
for backend in $BACKENDS; do
    run_backend xios "$xios_run" "$backend" --xios-dir "$XIOS_ROOT_DIR" || status_xios=1
done
[[ -L "$xios_run/xios_server.exe" ]] && echo "XIOS server link: OK" || { echo "XIOS server link: MISSING" >&2; status_xios=1; }

oasis_run="$(prepare_case oasis OASIS)"
for backend in $BACKENDS; do
    run_backend oasis "$oasis_run" "$backend" --prism-dir "$CROCO_PRISM_ROOT_DIR" --mpif90 mpif90 || status_oasis=1
done

echo
echo "Summary:"
[[ "$status_xios" == "0" ]] && echo "  XIOS build/link: PASS" || echo "  XIOS build/link: FAIL"
[[ "$status_oasis" == "0" ]] && echo "  OASIS build/link: PASS" || echo "  OASIS build/link: FAIL"
echo "  Full coupled run: not attempted by this experimental script"
echo "  Logs: $LOGS"
exit $(( status_xios || status_oasis ))
