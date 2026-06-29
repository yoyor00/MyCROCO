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
 * cppdefs.h - MS3DVAR Medium Resolution (MR) Variant Configuration
 * ================================================================
 *
 * This is the CPP configuration file for the MS3DVAR MR (Medium Resolution)
 * variant. It includes :
 *
 * 1. A slightly modified version of the standard Croco header
 *    "OCEAN/croco/cppdefs.h" which is copied to 
 *    "ASSIM/MS3DVAR/${SCALE}/cppdefs_croco.h" with on the fly 
 *    modifications introduced at compilation stage, namely : 
 *    - the final standard includes of "OCEAN/croco/cppdefs_dev.h" and 
 *      "set_global_definitions.h" are removed; because these instructions
 *      must be brought at the very end of the current final header 
 *      "ASSIM/MS3DVAR/${SCALE}/cppdefs.h". This change is mandatory.
 *    - The "ASSIM/MS3DVAR/${SCALE}/cppdefs_croco.h" header may be itself 
 *      a factorized header file, in the case the user run a REGIONAL 
 *      or COASTAL configuration. 
 *      This is a proposal to follow @sraynaud's work. 
 *      But this latter change is not mandatory and in any case 
 *      it has to be decided at the Croco-dev level.
 * 2. Very few redefinition of cpp switches, such as OPENMP,
 *    in order to fullfill a set of mandatory MS3DVAR cpp keys.
 * 3. An MS3DVAR cpp switch saying which scale-variant is currently active : 
 *    DAS_FILTER, DAS_LR , DAS_MR or DAS_MS
 * 4. The main MS3DVAR header file "ASSIM/MS3DVAR/COMMON/cppdefs_ms3dvar.h" 
 *    with MS3DVAR CPP switches shared among all scales 
 * 5. MS3DVAR Scale-variant specific cpp switches such as :
 *    DAS_READ_INC,...
 * 6. The final standard croco includes of "OCEAN/croco/cppdefs_dev.h" 
 *    and "set_global_definitions.h" followed by the MS3DVAR final
 *    include "ASSIM/MS3DVAR/COMMON/das_set_global_def.h"
 *
 * SCALE = FILTER, LR, MR or MS
 *
 * WHAT CHANGED (SPRING 2026):
 * -----------------------------
 * This file was refactored from 2203 lines to ~200 lines.
 * The full CROCO cppdefs.h (lines 1-2139) has been removed
 *
 * BACKWARD COMPATIBILITY:
 * -----------------------
 * All functionality is preserved - this is a restructuring, not a change.
 *
 * BENEFITS:
 * ---------
 * - Only shows relevant MS3DVAR options
 * - Easier to understand and maintain
 * - Reduces risk of enabling incompatible options
 */

/*=====================================================================
 * CROCO STANDARD CONFIGURATION FILE INCLUDED
 *
 * A slightly modified version of the standard Croco header
 * "OCEAN/croco/cppdefs.h" is copied to 
 * "ASSIM/MS3DVAR/${SCALE}/cppdefs_croco.h" with on the fly 
 * modifications introduced at compilation stage, namely : 
 *  - the final standard includes of "OCEAN/croco/cppdefs_dev.h" and 
 *    "set_global_definitions.h" are removed; because these instructions
 *    must be brought at the very end of the current final header 
 *    "ASSIM/MS3DVAR/$SCALE/cppdefs.h". This change is mandatory.
 *  - Also a factorization of  "ASSIM/MS3DVAR/$SCALE/cppdefs_croco.h" 
 *    may itself be performed when the user either selects the REGIONAL 
 *    or the COASTAL configuration type.
 *    In this case, the jobcomp will modify this header file
 *    eliminating all configuration definition lines in between 
 *    the REGIONAL (or COASTAL) markups and replacing them by including
 *    the user defined configuration file "ASSIM/MS3DVAR/COMMON/cppdefs_config.h"
 *    This is ONLY A PROPOSAL to follow @sraynaud's work. 
 *    But this latter change is NOT MANDATORY and in any case 
 *    it has TO BE DECIDED AT the Croco-dev level.
 *
 * SCALE = FILTER, LR, MR or MS
 *
 * See path "OCEAN/cppdefs.h"
 * Copied to "ASSIM/MS3DVAR/${SCALE}/Compile/cppdefs_croco.h"
 * 
 * See path "ASSIM/MS3DVAR/COMMON/cppdefs_config.h"
 * Copied to "ASSIM/MS3DVAR/${SCALE}/Compile/cppdefs_config.h"
 *
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
 * MR VARIANT CPP-switch : DO NOT CHANGE
 *=====================================================================*/
#define DAS_MR

/*=====================================================================
 * FULL STANDARD MS3DVAR CONFIGURATION
 *
 * This includes the core MS3DVAR CPP defines (~150 lines)
 * Settings below the include will override defaults in cppdefs_ms3dvar.h
 *
 * See path ASSIM/MS3DVAR/COMMON/cppdefs_ms3dvar.
 * Copied to  ASSIM/MS3DVAR/MR/Compile/cppdefs_ms3dvar.h
 *=====================================================================*/
#include "cppdefs_ms3dvar.h"

/*=====================================================================
 * MR VARIANT-SPECIFIC OVERRIDES
 *
 * Customize MS3DVAR DAS_* options specifically for MR variant.
 * These override the defaults in cppdefs_ms3dvar.h
 *=====================================================================*/

/* MR-specific overrides currently - using defaults from cppdefs_ms3dvar.h */

#define  DAS_JASONSSH          /* Jason altimeter series */

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
