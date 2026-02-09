!======================================================================
! CROCO is a branch of ROMS developped at IRD, INRIA,
! Ifremer, CNRS and Univ. Toulouse III  in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
!
! CROCO website : http://www.croco-ocean.org
!======================================================================
!
! param.h - MS3DVAR Medium Resolution (MR) Variant Grid Parameters
! =================================================================
!
! This is the parameter file for the MS3DVAR MR (Medium Resolution) variant.
! It includes the minimal param_ms3dvar.h and can add MR-specific settings.
!
! WHAT CHANGED (February 2026):
! -----------------------------
! This file was refactored from 1151 lines to ~120 lines.
! The full CROCO param.h (lines 51-1151) has been removed because
! MS3DVAR only needs MS3DVAR-specific grid parameters (lines 1-50)
! plus minimal CROCO grid parameters.
!
! Old structure: param.h (1151 lines)
!   - Lines 1-50: MS3DVAR-specific grid parameters
!   - Lines 51-1151: Full CROCO param.h with all test cases
!
! New structure: param.h (this file, ~120 lines)
!   - Include ../COMMON/param_ms3dvar.h (minimal grid parameters)
!   - Add MR variant-specific overrides below (if any)
!
! BACKWARD COMPATIBILITY:
! -----------------------
! The full original file is saved as param.h.full_backup
! All functionality is preserved - this is a restructuring, not a change.
!
! BENEFITS:
! ---------
! - 91% smaller file (1151 lines -> ~120 lines)
! - Only shows MS3DVAR-relevant grid parameters
! - Eliminates confusion about CROCO test cases
! - Easier to understand configuration
!
!======================================================================

!----------------------------------------------------------------------
! INCLUDE MINIMAL MS3DVAR PARAMETERS
!
! This includes the core MS3DVAR grid parameters (~100 lines):
! - High-resolution geopotential grid (LLmH, MMmH, NnH)
! - Medium-resolution grid (LLm_mr, MMm_mr) - for MR variant
! - Tiling configuration (NSUB_XM, NSUB_EM)
! - Standard CROCO grid parameters (LLm, MMm, N, etc.)
!----------------------------------------------------------------------

#include "../COMMON/param_ms3dvar.h"

!----------------------------------------------------------------------
! MR VARIANT-SPECIFIC PARAMETER OVERRIDES
!
! If the MR variant needs different grid dimensions or tiling
! from the defaults in param_ms3dvar.h, override them here.
!----------------------------------------------------------------------

! Currently using defaults from param_ms3dvar.h:
! - LLmH=629, MMmH=537, NnH=50 (WMED high-resolution grid)
! - nratio=1 (no downsampling for MR)
! - NSUB_XM=6, NSUB_EM=6 (6x6 tiling)
!
! To override, uncomment and modify below:

! Example: Different high-resolution grid for MR
!     parameter (LLmH=xxx,  MMmH=yyy, NnH=zz)

! Example: Different downsampling for MR
!     parameter (nratio=2, nhalf=nratio/2+1)

! Example: Different tiling for MR
!     parameter (NSUB_XM=4, NSUB_EM=4, NPPM=16)

!----------------------------------------------------------------------
! ADDITIONAL CROCO PARAMETERS (if needed)
!
! param_ms3dvar.h provides the minimal set of CROCO parameters.
! If you need additional CROCO-specific parameters for your MR variant,
! define them here.
!----------------------------------------------------------------------

! Example: Tracer parameters
!     integer   NT, itemp, isalt
!     parameter (itemp=1, isalt=2, NT=2)

! Example: Biology parameters (if using BIOLOGY)
!     integer   ntrc_bio
!     parameter (ntrc_bio=6)

! Most standard CROCO parameters are computed automatically
! from LLm, MMm, N in param_ms3dvar.h

!----------------------------------------------------------------------
! END OF MR VARIANT PARAMETERS
!----------------------------------------------------------------------
