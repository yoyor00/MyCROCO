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
/*
 * cppdefs.h - MS3DVAR Low Resolution (LR) Variant Configuration
 * ==============================================================
 *
 * This is the CPP configuration file for the MS3DVAR LR (Low Resolution)
 * variant. It includes the minimal cppdefs_ms3dvar.h and adds LR-specific
 * settings.
 *
 * WHAT CHANGED (February 2026):
 * -----------------------------
 * This file was refactored from 2203 lines to ~200 lines.
 * The full CROCO cppdefs.h (lines 1-2139) has been removed because
 * MS3DVAR only needs ~5 CROCO defines + ~60 DAS_* defines.
 *
 * Old structure: cppdefs.h (2203 lines)
 *   - Lines 1-2139: Full CROCO config (test cases, physics, biology, etc.)
 *   - Lines 2142-2203: MS3DVAR-specific (DAS_*) defines
 *
 * New structure: cppdefs.h (this file, ~200 lines)
 *   - Include ../COMMON/cppdefs_ms3dvar.h (minimal core + DAS_* defines)
 *   - Add LR variant-specific overrides below
 *
 * BACKWARD COMPATIBILITY:
 * -----------------------
 * The full original file is saved as cppdefs.h.full_backup
 * All functionality is preserved - this is a restructuring, not a change.
 *
 * BENEFITS:
 * ---------
 * - 93% smaller file (2203 lines -> 150 lines + overrides)
 * - Only shows relevant MS3DVAR options
 * - Eliminates confusion about CROCO test cases
 * - Easier to understand and maintain
 * - Reduces risk of enabling incompatible options
 */

/*=====================================================================
 * CROCO STANDARD CONFIGURATION FILE
 *
 * The file included below is the standard croco cppdefs.h
 * to which has been removed the final include, namely
 * cppdefs_dev.h and set_global_definitions.h
 *
 * Also a factorization of the standard croco cppdefs.h file
 * is performed when the user either selects the REGIONAL or the
 * COASTAL configuration type.
 * In this case, the jobcomp will modify the standard croco cppdefs
 * replacing all config definition lines in between 
 * the REGIONAL (or COASTAL) markups by including
 * the user defined configuration file cppdefs_config.h
 * 
 * See path ASSIM/MS3DVAR/COMMON/cppdefs_config.h
 * Copied to  ASSIM/MS3DVAR/FILTER/Compile/cppdefs_config.h
 *
 * See path OCEAN/cppdefs.h
 * Copied to ASSIM/MS3DVAR/LR/Compile/cppdefs_croco.h
 *=====================================================================*/
#include "cppdefs_croco.h"

/*=====================================================================
 * MANDATORY CONFIGURATION CHANGES FOR ALL MS3DVAR SCALE EXECUTABLES
 *
 *=====================================================================*/

/*=====================================================================
 * Parallelization options
 *
 * OpenMP (shared memory) is mandatory
 * MPI (distributed memory) doesn't work yet.
 * Note: DAS is incompatible with CROCO's AUTOTILING + OpenMP combination.
 *=====================================================================*/
#define OPENMP          /* OpenMP parallelization */
#undef  MPI             /* MPI parallelization */

/*=====================================================================
 * Grid and coordinate options
 *
 * Expected grid configuration for MS3DVAR (TOCHECK)
 *=====================================================================*/
#define MASKING         /* Land/sea masking */
#define NEW_S_COORD     /* New s-coordinate system (Vtransform=2) */
#define SPHERICAL       /* Spherical (geographic) coordinate grid */
#define CURVGRID        /* Curvilinear (non-orthogonal) grid */

/*=====================================================================
 * Time management
 *=====================================================================*/
#undef USE_CALENDAR    /* Use calendar for time management */

/*=====================================================================
 * XIOS I/O server
 *=====================================================================*/
#undef  XIOS

/*=====================================================================
 * LR VARIANT CPP-switch : DO NOT CHANGE
 *=====================================================================*/
#define DAS_LR

/*=====================================================================
 * FULL STANDARD MS3DVAR CONFIGURATION
 *
 * This includes the core MS3DVAR CPP defines (~150 lines)
 * Settings below the include will override defaults in cppdefs_ms3dvar.h
 *
 * See path ASSIM/MS3DVAR/COMMON/cppdefs_ms3dvar.
 * Copied to  ASSIM/MS3DVAR/LR/Compile/cppdefs_ms3dvar.h
 *=====================================================================*/
#include "cppdefs_ms3dvar.h"

/*=====================================================================
 * LR VARIANT-SPECIFIC OVERRIDES
 *
 * Customize MS3DVAR DAS_* options specifically for LR variant.
 * These override the defaults in cppdefs_ms3dvar.h
 *=====================================================================*/

/* No LR-specific overrides currently - using defaults from cppdefs_ms3dvar.h */

/*=====================================================================
 * END OF CONFIGURATION
 *=====================================================================*/

/*=====================================================================
 * DO NOT TOUCH THE BELOW INCLUDED FILES
 *
 * These files contain additional preprocessor logic and global
 * definitions. They must be included last.
 *=====================================================================*/
#include "cppdefs_dev.h"
#include "set_global_definitions.h"
#include "das_set_global_def.h"
