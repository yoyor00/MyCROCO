#!/bin/bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=versions.env
source "$HERE/versions.env"

FORCE=0
[[ "${1:-}" == "--force" ]] && FORCE=1
SRC="${XIOS_SRC_DIR:-$HERE/src/xios}"

if [[ -d "$SRC/.git" && "$FORCE" == "0" ]]; then
    echo "XIOS source already present: $SRC"
    exit 0
fi
rm -rf "$SRC"
mkdir -p "$HERE/src"
echo "Cloning XIOS from $XIOS_REPO ($XIOS_REF)"
git clone --depth 1 --branch "$XIOS_REF" "$XIOS_REPO" "$SRC"
echo "XIOS source ready: $SRC"
