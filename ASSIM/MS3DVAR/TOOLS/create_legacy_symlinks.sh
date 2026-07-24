#!/bin/bash
# create_legacy_symlinks.sh -- create old-name -> new-name symlinks for every
# file renamed by the ms3dvar_ prefix rename (see docs/ms3dvar/RENAME_PLAN.md).
#
# Purpose: lets you diff a file against the original (pre-rename) MS-3DVAR
# repo by its old name, e.g.:
#   diff ASSIM/MS3DVAR/COMMON/das_get_corr.F \
#        /path/to/original/croco_ms3dvar/da_l0_lr/das_get_corr.F
# even though the real file is now ms3dvar_get_corr.F.
#
# Purely a dev convenience -- none of these symlinks are required by the
# build. (cppdefs.h/param.h are NOT in this manifest: they keep their literal
# CROCO names as real files, since CROCO's untouched OCEAN/ source does
# #include "cppdefs.h" / #include "param.h" and jobcomp's plain glob-copy
# depends on that literal name existing.)
#
# Usage:
#   ./create_legacy_symlinks.sh            # create every missing symlink
#   ./create_legacy_symlinks.sh --remove   # remove them all
#   ./create_legacy_symlinks.sh --check    # verify only, exit 1 if anything's missing/wrong
#
# Driven entirely by TOOLS/rename_map.tsv (old_relpath<TAB>new_relpath, one
# row per renamed file, paths relative to ASSIM/MS3DVAR/).

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MS3DVAR_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
MANIFEST="${SCRIPT_DIR}/rename_map.tsv"

MODE="create"
case "$1" in
    --remove) MODE="remove" ;;
    --check)  MODE="check" ;;
    "") ;;
    *) echo "Usage: $0 [--remove|--check]"; exit 1 ;;
esac

[ -f "$MANIFEST" ] || { echo "ERROR: manifest not found: $MANIFEST"; exit 1; }

n_created=0
n_ok=0
n_removed=0
n_errors=0

while IFS=$'\t' read -r old new; do
    [ -z "$old" ] && continue
    old_path="${MS3DVAR_DIR}/${old}"
    new_path="${MS3DVAR_DIR}/${new}"
    new_base="$(basename "$new")"

    case "$MODE" in
        remove)
            if [ -L "$old_path" ]; then
                rm "$old_path"
                n_removed=$((n_removed + 1))
            fi
            ;;

        check)
            if [ -L "$old_path" ] && [ "$(readlink "$old_path")" = "$new_base" ]; then
                n_ok=$((n_ok + 1))
            else
                echo "MISSING/WRONG: $old -> $new_base"
                n_errors=$((n_errors + 1))
            fi
            ;;

        create)
            if [ ! -e "$new_path" ] && [ ! -L "$new_path" ]; then
                echo "ERROR: target does not exist, cannot link: $new"
                n_errors=$((n_errors + 1))
                continue
            fi
            if [ -L "$old_path" ]; then
                if [ "$(readlink "$old_path")" = "$new_base" ]; then
                    n_ok=$((n_ok + 1))
                    continue
                else
                    echo "ERROR: $old already a symlink to something else ($(readlink "$old_path")), not touching"
                    n_errors=$((n_errors + 1))
                    continue
                fi
            fi
            if [ -e "$old_path" ]; then
                echo "ERROR: $old exists as a real file, not overwriting"
                n_errors=$((n_errors + 1))
                continue
            fi
            ln -s "$new_base" "$old_path"
            n_created=$((n_created + 1))
            ;;
    esac
done < "$MANIFEST"

case "$MODE" in
    create) echo "Created: $n_created, already OK: $n_ok, errors: $n_errors" ;;
    remove) echo "Removed: $n_removed" ;;
    check)  echo "OK: $n_ok, missing/wrong: $n_errors" ;;
esac

[ "$n_errors" -eq 0 ]
