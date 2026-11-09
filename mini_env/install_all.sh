#!/bin/bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DO_XIOS=0
DO_OASIS=0
FORCE=0
for arg in "$@"; do
    case "$arg" in
        --xios) DO_XIOS=1 ;;
        --oasis) DO_OASIS=1 ;;
        --all) DO_XIOS=1; DO_OASIS=1 ;;
        --force) FORCE=1 ;;
        -h|--help)
            echo "Usage: $0 [--xios] [--oasis] [--all] [--force]"
            exit 0 ;;
        *) echo "Unknown option: $arg" >&2; exit 1 ;;
    esac
done
if [[ "$DO_XIOS" == "0" && "$DO_OASIS" == "0" ]]; then
    DO_XIOS=1
    DO_OASIS=1
fi
force_args=()
[[ "$FORCE" == "1" ]] && force_args=(--force)

if [[ "$DO_XIOS" == "1" ]]; then
    "$HERE/download_xios.sh" "${force_args[@]}"
    "$HERE/install_xios.sh"
fi
if [[ "$DO_OASIS" == "1" ]]; then
    "$HERE/download_oasis.sh" "${force_args[@]}"
    "$HERE/install_oasis.sh"
fi
echo "Mini environment installed under $HERE/install"
