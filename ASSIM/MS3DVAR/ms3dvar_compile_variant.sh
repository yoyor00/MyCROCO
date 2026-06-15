#!/bin/bash
# ms3dvar_compile_variant.sh — Build one MS3DVAR variant
#
# Called by WOOM compile_filter/lr/mr/ms tasks (see WOOM/tasks.cfg).
# The compiler environment is set by the ms3dvar_compile env in hosts.cfg.
#
# Usage:
#   ms3dvar_compile_variant.sh <VARIANT> <MPC_PATH> [HEADERS_DIR]
#
#   VARIANT     : FILTER | LR | MR | MS
#   MPC_PATH    : path to the pre-built mpc binary (from compile_mpc artifact)
#   HEADERS_DIR : optional directory of pycroco-generated .h files (from
#                 generate_headers artifact).  Staged into CWD with cp -n so
#                 any file already present (user override) is not overwritten.
#
# Source priority (lowest → highest, all gathered into compile/ before make):
#   1. CROCO OCEAN sources  ($CROCO_SRC)
#   2. MS3DVAR COMMON sources
#   3. Variant-specific sources
#   4. pycroco-generated headers  (HEADERS_DIR, via CWD — cp -n)
#   5. User's own files in CWD    (highest priority — .h .h90 .F .F90).
#      Placed manually or copied by generate_headers without overwriting user files.
#
# The built binary is placed in CWD (the WOOM run_dir) for artifact pickup.
# All intermediate files live in compile/ and do not pollute CWD.
#
# Environment variables (set by WOOM hosts.cfg ms3dvar_compile env):
#   CROCO_SRC, FC, CPP1, CFT1, FFLAGS1, CPPFLAGS1, LDFLAGS1, NPROCS

set -eo pipefail

VARIANT="${1:?Usage: $0 VARIANT MPC_PATH [HEADERS_DIR]}"
MPC_PATH="${2:?MPC_PATH required}"
HEADERS_DIR="${3:-}"

: "${CROCO_SRC:?CROCO_SRC env var is not set}"

MS3DVAR_SRC="${CROCO_SRC}/../ASSIM/MS3DVAR"

echo "===== MS3DVAR compile: ${VARIANT} ====="
echo "MS3DVAR_SRC : ${MS3DVAR_SRC}"
echo "CROCO_SRC   : ${CROCO_SRC}"
echo "Compiler    : ${FC}"
echo "NPROCS      : ${NPROCS}"
echo ""

# ---- Stage generated headers into CWD (user files take priority) ----
# cp -n never overwrites, so a file the user already placed in CWD wins.
if [[ -n "${HEADERS_DIR}" ]]; then
    shopt -s nullglob
    for f in "${HEADERS_DIR}"/*.h "${HEADERS_DIR}"/*.h90; do
        cp -n "$f" . && echo "  staged $(basename "$f") from generate_headers" || true
    done
    shopt -u nullglob
fi

# ---- Create build directory -----------------------------------------
mkdir -p compile

# ---- Gather sources into compile/ (priority order 1→3) --------------

# 1. CROCO OCEAN
shopt -s nullglob
for ext in F F90 h; do
    files=("${CROCO_SRC}"/*."${ext}")
    (( ${#files[@]} > 0 )) && \cp "${files[@]}" compile/
done
shopt -u nullglob
# Save param.h separately: param_ms3dvar.h includes it as croco_ocean_param.h
\cp "${CROCO_SRC}/param.h" compile/croco_ocean_param.h

# 2. MS3DVAR COMMON
shopt -s nullglob
for ext in F F90 h; do
    files=("${MS3DVAR_SRC}/COMMON/"*."${ext}")
    (( ${#files[@]} > 0 )) && \cp "${files[@]}" compile/
done
shopt -u nullglob

# 3. Variant-specific
\cp "${MS3DVAR_SRC}/${VARIANT}/Makefile" compile/
shopt -s nullglob
for ext in F F90 h; do
    files=("${MS3DVAR_SRC}/${VARIANT}/"*."${ext}")
    (( ${#files[@]} > 0 )) && \cp "${files[@]}" compile/
done
shopt -u nullglob

# ---- Apply CWD overrides last (priority 4+5: generated + user) ------
# Source files (.h, .h90, .F, .F90) in CWD override CROCO/MS3DVAR defaults.
# Mirrors jobcomp user-file handling.  Place custom sources here before running.
shopt -s nullglob
cwd_files=(*.h *.h90 *.F *.F90)
shopt -u nullglob
if (( ${#cwd_files[@]} > 0 )); then
    \cp "${cwd_files[@]}" compile/
    echo "Applied ${#cwd_files[@]} override file(s) from working directory"
fi

# ---- Configure Makefile.inc -----------------------------------------
\cp "${MS3DVAR_SRC}/Makefile.inc" compile/Makefile.inc
sed -i 's|^COMMON_DIR = .*|COMMON_DIR = .|' compile/Makefile.inc
sed -i 's|^CROCO_SRC = .*|CROCO_SRC = .|'  compile/Makefile.inc
sed -i 's|^VPATH = .*|VPATH = .|'          compile/Makefile.inc

# ---- Install mpc preprocessor ---------------------------------------
\cp "${MPC_PATH}" compile/mpc

# ---- Generate Makedefs ----------------------------------------------
{
    printf 's?$(FFLAGS1)?%s?g\n'   "${FFLAGS1}"
    printf 's?$(LDFLAGS1)?%s?g\n'  "${LDFLAGS1}"
    printf 's?$(CPP1)?%s?g\n'      "${CPP1}"
    printf 's?$(CFT1)?%s?g\n'      "${CFT1}"
    printf 's?$(CPPFLAGS1)?%s?g\n' "${CPPFLAGS1}"
} > compile/flags.tmp
sed -f compile/flags.tmp "${MS3DVAR_SRC}/Makedefs.generic" > compile/Makedefs
rm -f compile/flags.tmp

# ---- Build ----------------------------------------------------------
# Stale .mod files cause spurious "m2c: No such file or directory" errors.
(
    cd compile
    rm -f *.mod
    make -j"${NPROCS}" all
)

# ---- Move binary to CWD for WOOM artifact ---------------------------
shopt -s nullglob
binaries=(compile/ms3dvar-*)
shopt -u nullglob
if (( ${#binaries[@]} == 0 )); then
    echo "ERROR: no ms3dvar-* binary found in compile/" >&2
    exit 1
fi
mv "${binaries[@]}" .

echo ""
echo "===== ${VARIANT} built successfully ====="
