#!/bin/bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$HERE/.." && pwd)"
VENV="$HERE/.venv"

python3 -m venv "$VENV"
"$VENV/bin/python" -m pip install --upgrade pip
"$VENV/bin/python" -m pip install -r "$REPO_ROOT/requirements.txt"
"$VENV/bin/python" -m pip install -r "$REPO_ROOT/OCEAN/python/scons/requirements.txt"
echo "Python environment ready: $VENV"
