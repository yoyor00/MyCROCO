#!/bin/bash
# ms3dvar_compile_scale.sh — Build one (or all) MS3DVAR scale(s) in the current directory
#
# This script is designed to be called by WOOM (Workflow manager for Ocean Models).
# It is invoked by the compile_filter/lr/mr/ms/all tasks defined in WOOM/tasks.cfg.
# The compiler environment (FC, flags, NetCDF paths) and source paths are set by
# the ms3dvar_compile environment declared in WOOM/hosts.cfg.
#
# Usage:
#   ms3dvar_compile_scale.sh <SCALE> <MPC_PATH> [--user-files DIR]
#
#   SCALE       : FILTER | LR | MR | MS | ALL
#                 ALL builds every scale, each in its own subdirectory of
#                 the current directory (./filter/, ./lr/, ./mr/, ./ms/),
#                 mirroring jobcomp's per-scale directory layout. Stops at
#                 the first failure (no partial-success summary, unlike
#                 jobcomp -- a WOOM task is meant to have a clear pass/fail).
#   MPC_PATH    : path to the pre-built mpc binary (artifact from compile_mpc task)
#   --user-files DIR : directory containing per-scale user file overrides,
#                       in DIR/<scale-lowercase>/ (e.g. DIR/lr/foo.F). Files
#                       there are copied in last, after every other source,
#                       so they take precedence. Optional; if omitted, no
#                       override is applied. Independent of MS3DVAR_SRC:
#                       these are not part of the versioned MS3DVAR sources,
#                       they can live anywhere -- WOOM points this at
#                       WOOM/user_files/ (see WOOM/tasks.cfg).
#                       May appear before, after, or between the positional
#                       arguments (mirrors jobcomp, which accepts positionals
#                       and options in either order).
#
# Environment variables (set by WOOM hosts.cfg ms3dvar_compile env):
#   MS3DVAR_SRC, CROCO_SRC, FC, CPP1, CFT1, FFLAGS1, CPPFLAGS1, LDFLAGS1, NPROCS
#
# WOOM runs this script from the task run_dir (already set as cwd).
# All source files are copied in; the built executable(s) stay in run_dir
# (or run_dir/<scale>/ for ALL) and are declared as WOOM artifacts.

set -e

POSITIONAL=()
USER_FILES_DIR=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --user-files)
            USER_FILES_DIR="${2:?--user-files requires a directory argument}"
            shift 2
            ;;
        --user-files=*)
            USER_FILES_DIR="${1#*=}"
            shift
            ;;
        -*)
            echo "Unknown option: $1" >&2
            exit 1
            ;;
        *)
            POSITIONAL+=("$1")
            shift
            ;;
    esac
done

SCALE="${POSITIONAL[0]:?Usage: $0 SCALE MPC_PATH [--user-files DIR]}"
MPC_PATH="${POSITIONAL[1]:?MPC_PATH required}"

: "${MS3DVAR_SRC:?MS3DVAR_SRC env var is not set}"
: "${CROCO_SRC:?CROCO_SRC env var is not set}"

# ---- Build logic for a single scale, run in the current directory --------
compile_one_scale() {
    local SCALE="$1"
    local MPC_PATH="$2"
    local USER_FILES_DIR="$3"

    echo "===== MS3DVAR compile: ${SCALE} ====="
    echo "MS3DVAR_SRC : ${MS3DVAR_SRC}"
    echo "CROCO_SRC   : ${CROCO_SRC}"
    echo "Compiler    : ${FC}"
    echo "NPROCS      : ${NPROCS}"
    echo ""

    # ---- Copy CROCO ocean sources ------------------------------------
    ls "${CROCO_SRC}"/*.F   >/dev/null 2>&1 && \cp "${CROCO_SRC}"/*.F   .
    ls "${CROCO_SRC}"/*.F90 >/dev/null 2>&1 && \cp "${CROCO_SRC}"/*.F90 . || true
    ls "${CROCO_SRC}"/*.h   >/dev/null 2>&1 && \cp "${CROCO_SRC}"/*.h   .

    # Save a modified CROCO cppdefs.h as croco_cppdefs.h before scale files
    # overwrite it (mirrors jobcomp: the scale's cppdefs.h does
    # #include "croco_cppdefs.h" and expects the two final includes below
    # stripped, since it re-includes them itself at the end).
    \cp "${CROCO_SRC}/cppdefs.h" croco_cppdefs.h
    sed -i '/^#.*include .*cppdefs_dev/{d;}' croco_cppdefs.h 2>/dev/null || true
    sed -i '/^#.*include .*set_global_definitions/{d;}' croco_cppdefs.h 2>/dev/null || true

    # Save CROCO param.h before scale files overwrite it (included as croco_param.h)
    \cp "${CROCO_SRC}/param.h" croco_param.h

    # ---- Copy MS3DVAR common sources ---------------------------------
    ls "${MS3DVAR_SRC}/COMMON/"*.F   >/dev/null 2>&1 && \cp "${MS3DVAR_SRC}/COMMON/"*.F   .
    ls "${MS3DVAR_SRC}/COMMON/"*.F90 >/dev/null 2>&1 && \cp "${MS3DVAR_SRC}/COMMON/"*.F90 . || true
    ls "${MS3DVAR_SRC}/COMMON/"*.h   >/dev/null 2>&1 && \cp "${MS3DVAR_SRC}/COMMON/"*.h   .

    # ---- Copy scale-specific sources and Makefile --------------------
    \cp "${MS3DVAR_SRC}/${SCALE}/Makefile" .
    ls "${MS3DVAR_SRC}/${SCALE}/"*.h   >/dev/null 2>&1 && \cp "${MS3DVAR_SRC}/${SCALE}/"*.h   . || true
    ls "${MS3DVAR_SRC}/${SCALE}/"*.F   >/dev/null 2>&1 && \cp "${MS3DVAR_SRC}/${SCALE}/"*.F   . || true
    ls "${MS3DVAR_SRC}/${SCALE}/"*.F90 >/dev/null 2>&1 && \cp "${MS3DVAR_SRC}/${SCALE}/"*.F90 . || true

    # ---- Apply user files for this scale, if any ---------------------
    # Copied last so they take precedence over everything above.
    if [ -n "$USER_FILES_DIR" ]; then
        local SCALE_LOWER
        SCALE_LOWER=$(echo "${SCALE}" | tr '[:upper:]' '[:lower:]')
        local SCALE_USER_FILES_DIR="${USER_FILES_DIR}/${SCALE_LOWER}"
        if [ -d "$SCALE_USER_FILES_DIR" ]; then
            ls "$SCALE_USER_FILES_DIR"/* >/dev/null 2>&1 && \cp -f "$SCALE_USER_FILES_DIR"/* .
            echo "User files applied from: ${SCALE_USER_FILES_DIR}"
        fi
    fi

    # ---- Configure Makefile.inc to use current directory -------------
    \cp "${MS3DVAR_SRC}/Makefile.inc" Makefile.inc
    sed -i 's|^COMMON_DIR = .*|COMMON_DIR = .|' Makefile.inc
    sed -i 's|^CROCO_SRC = .*|CROCO_SRC = .|'  Makefile.inc
    sed -i 's|^VPATH = .*|VPATH = .|'          Makefile.inc

    # ---- Install mpc preprocessor ------------------------------------
    \cp "${MPC_PATH}" mpc

    # ---- Generate Makedefs from the generic template -----------------
    # Substitutes $(CPP1), $(CFT1), $(FFLAGS1), $(CPPFLAGS1), $(LDFLAGS1)
    # with the values exported by the ms3dvar_compile env in hosts.cfg.
    rm -f Makedefs flags.tmp
    printf 's?$(FFLAGS1)?%s?g\n'   "${FFLAGS1}"   >> flags.tmp
    printf 's?$(LDFLAGS1)?%s?g\n'  "${LDFLAGS1}"  >> flags.tmp
    printf 's?$(CPP1)?%s?g\n'      "${CPP1}"       >> flags.tmp
    printf 's?$(CFT1)?%s?g\n'      "${CFT1}"       >> flags.tmp
    printf 's?$(CPPFLAGS1)?%s?g\n' "${CPPFLAGS1}"  >> flags.tmp
    sed -f flags.tmp "${MS3DVAR_SRC}/Makedefs.generic" > Makedefs
    rm -f flags.tmp

    # ---- Build ---------------------------------------------------------
    # depend must run sequentially before the parallel build (mirrors
    # jobcomp): it generates Make.depend (module-use ordering, e.g.
    # buffer.mod before any file that USEs it), which -j can otherwise race
    # past since $(OBJS) targets don't themselves depend on it.
    make depend
    make -j${NPROCS} all

    echo ""
    echo "===== ${SCALE} built successfully ====="
}

SCALE_UPPER=$(echo "${SCALE}" | tr '[:lower:]' '[:upper:]')

if [ "$SCALE_UPPER" = "ALL" ]; then
    # Resolve to absolute paths before cd'ing into per-scale subdirectories.
    MPC_PATH="$(cd "$(dirname "$MPC_PATH")" && pwd)/$(basename "$MPC_PATH")"
    if [ -n "$USER_FILES_DIR" ]; then
        USER_FILES_DIR="$(cd "$USER_FILES_DIR" && pwd)"
    fi
    for s in FILTER LR MR MS; do
        s_lower=$(echo "$s" | tr '[:upper:]' '[:lower:]')
        mkdir -p "$s_lower"
        ( cd "$s_lower" && compile_one_scale "$s" "$MPC_PATH" "$USER_FILES_DIR" )
    done
else
    compile_one_scale "$SCALE_UPPER" "$MPC_PATH" "$USER_FILES_DIR"
fi
