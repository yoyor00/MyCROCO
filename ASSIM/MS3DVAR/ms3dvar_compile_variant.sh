#!/bin/bash
# ms3dvar_compile_variant.sh — Build one MS3DVAR variant in the current directory
#
# This script is designed to be called by WOOM (Workflow manager for Ocean Models).
# It is invoked by the compile_filter/lr/mr/ms tasks defined in WOOM/tasks.cfg.
# The compiler environment (FC, flags, NetCDF paths) and source paths are set by
# the ms3dvar_compile environment declared in WOOM/hosts.cfg.
#
# Usage:
#   ms3dvar_compile_variant.sh <VARIANT> <MPC_PATH>
#
#   VARIANT   : FILTER | LR | MR | MS
#   MPC_PATH  : path to the pre-built mpc binary (artifact from compile_mpc task)
#
# Environment variables (set by WOOM hosts.cfg ms3dvar_compile env):
#   MS3DVAR_SRC, CROCO_SRC, FC, CPP1, CFT1, FFLAGS1, CPPFLAGS1, LDFLAGS1, NPROCS
#
# WOOM runs this script from the task run_dir (already set as cwd).
# All source files are copied in; the built executable stays in run_dir
# and is declared as a WOOM artifact.

set -e

VARIANT="${1:?Usage: $0 VARIANT MPC_PATH}"
MPC_PATH="${2:?MPC_PATH required}"

: "${MS3DVAR_SRC:?MS3DVAR_SRC env var is not set}"
: "${CROCO_SRC:?CROCO_SRC env var is not set}"

echo "===== MS3DVAR compile: ${VARIANT} ====="
echo "MS3DVAR_SRC : ${MS3DVAR_SRC}"
echo "CROCO_SRC   : ${CROCO_SRC}"
echo "Compiler    : ${FC}"
echo "NPROCS      : ${NPROCS}"
echo ""

# ---- Copy CROCO ocean sources ----------------------------------------
ls "${CROCO_SRC}"/*.F   >/dev/null 2>&1 && \cp "${CROCO_SRC}"/*.F   .
ls "${CROCO_SRC}"/*.F90 >/dev/null 2>&1 && \cp "${CROCO_SRC}"/*.F90 . || true
ls "${CROCO_SRC}"/*.h   >/dev/null 2>&1 && \cp "${CROCO_SRC}"/*.h   .

# Save CROCO param.h before variant files overwrite it (included as croco_ocean_param.h)
\cp "${CROCO_SRC}/param.h" croco_ocean_param.h

# ---- Copy MS3DVAR common sources -------------------------------------
ls "${MS3DVAR_SRC}/COMMON/"*.F   >/dev/null 2>&1 && \cp "${MS3DVAR_SRC}/COMMON/"*.F   .
ls "${MS3DVAR_SRC}/COMMON/"*.F90 >/dev/null 2>&1 && \cp "${MS3DVAR_SRC}/COMMON/"*.F90 . || true
ls "${MS3DVAR_SRC}/COMMON/"*.h   >/dev/null 2>&1 && \cp "${MS3DVAR_SRC}/COMMON/"*.h   .

# ---- Copy variant-specific sources and Makefile ----------------------
\cp "${MS3DVAR_SRC}/${VARIANT}/Makefile" .
ls "${MS3DVAR_SRC}/${VARIANT}/"*.h   >/dev/null 2>&1 && \cp "${MS3DVAR_SRC}/${VARIANT}/"*.h   . || true
ls "${MS3DVAR_SRC}/${VARIANT}/"*.F   >/dev/null 2>&1 && \cp "${MS3DVAR_SRC}/${VARIANT}/"*.F   . || true
ls "${MS3DVAR_SRC}/${VARIANT}/"*.F90 >/dev/null 2>&1 && \cp "${MS3DVAR_SRC}/${VARIANT}/"*.F90 . || true

# ---- Configure Makefile.inc to use current directory -----------------
\cp "${MS3DVAR_SRC}/Makefile.inc" Makefile.inc
sed -i 's|^COMMON_DIR = .*|COMMON_DIR = .|' Makefile.inc
sed -i 's|^CROCO_SRC = .*|CROCO_SRC = .|'  Makefile.inc
sed -i 's|^VPATH = .*|VPATH = .|'          Makefile.inc

# ---- Install mpc preprocessor ----------------------------------------
\cp "${MPC_PATH}" mpc

# ---- Generate Makedefs from the generic template ---------------------
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

# ---- Build -----------------------------------------------------------
make -j${NPROCS} all

echo ""
echo "===== ${VARIANT} built successfully ====="
