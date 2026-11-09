#!/bin/bash
# Remove local experimental artifacts.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

MODE=temp
DRY_RUN=0
for arg in "$@"; do
    case "$arg" in
        --all) MODE=all ;;
        --env) MODE=env ;;
        --build) MODE=build ;;
        --source) MODE=source ;;
        --work) MODE=work ;;
        --logs) MODE=logs ;;
        --dry-run) DRY_RUN=1 ;;
        -h|--help)
            cat <<HELP
Usage: $0 [--all|--env|--build|--source|--work|--logs] [--dry-run]

Default removes temporary work, logs and downloads.
  --work    remove mini_env/work
  --logs    remove mini_env/logs
  --env     remove mini_env/.venv
  --build   remove mini_env/build and mini_env/install
  --source  remove mini_env/src
  --all     remove all generated content: src, build, install, downloads,
            work, logs and .venv
  --dry-run print what would be removed
HELP
            exit 0 ;;
        *) echo "Unknown option: $arg" >&2; exit 1 ;;
    esac
done

targets=()
case "$MODE" in
    temp)   targets=("$HERE/work" "$HERE/logs" "$HERE/downloads") ;;
    work)   targets=("$HERE/work") ;;
    logs)   targets=("$HERE/logs") ;;
    env)    targets=("$HERE/.venv") ;;
    build)  targets=("$HERE/build" "$HERE/install") ;;
    source) targets=("$HERE/src") ;;
    all)    targets=("$HERE/src" "$HERE/build" "$HERE/install" "$HERE/downloads" "$HERE/work" "$HERE/logs" "$HERE/.venv") ;;
esac

for target in "${targets[@]}"; do
    if [[ ! -e "$target" ]]; then
        continue
    fi
    if [[ "$DRY_RUN" == "1" ]]; then
        echo "would remove $target"
    else
        echo "removing $target"
        rm -rf "$target"
    fi
done

# Core dumps from failed experimental runs.
for core in "$HERE"/core "$HERE"/core.*; do
    [[ -e "$core" ]] || continue
    if [[ "$DRY_RUN" == "1" ]]; then
        echo "would remove $core"
    else
        echo "removing $core"
        rm -f "$core"
    fi
done
