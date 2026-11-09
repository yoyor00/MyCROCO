#!/bin/bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=versions.env
source "$HERE/versions.env"

FORCE=0
[[ "${1:-}" == "--force" ]] && FORCE=1
SRC="$HERE/src/oasis3-mct"

if [[ -d "$SRC/.git" && "$FORCE" == "0" ]]; then
    echo "OASIS source already present: $SRC"
    exit 0
fi
rm -rf "$SRC"
mkdir -p "$HERE/src"
echo "Cloning OASIS3-MCT from $OASIS_REPO ($OASIS_REF)"
git clone --depth 1 --branch "$OASIS_REF" "$OASIS_REPO" "$SRC"
echo "OASIS source ready: $SRC"
