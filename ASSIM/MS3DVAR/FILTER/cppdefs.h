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
 * cppdefs.h - MS3DVAR FILTER Scale Configuration
 * ==============================================================
 *
 * This is the CPP configuration file for the MS3DVAR FILTER
 * scale. It includes :
 *
 * 1. A slightly modified version of the standard Croco header
 *    "OCEAN/croco/cppdefs.h" which is copied to 
 *    "ASSIM/MS3DVAR/${SCALE}/croco_cppdefs.h" with on the fly 
 *    modifications introduced at the compilation stage, namely : 
 *    - the final standard includes of "OCEAN/croco/cppdefs_dev.h" and 
 *      "set_global_definitions.h" are removed; because to be effective
 *      it is mandatory for those two instructions to be brought 
 *      at the very end of the current cppdefs header file
 *      "ASSIM/MS3DVAR/${SCALE}/cppdefs.h". 
 * 2. The definition of the MS3DVAR CPP switch saying which 
 *    scale is currently active : DAS_FILTER or DAS_LR or 
 *    DAS_MR or DAS_MS.
 * 3. The core MS3DVAR header file "ASSIM/MS3DVAR/COMMON/ms3dvar_cppdefs.h" 
 *    which performs :
 *    - the redefinition of a very small number of Croco CPP keys in order
 *      to fullfill the current set of mandatory MS3DVAR CPP keys, e.g.,
 *      defining OPENMP, undefing MPI.
 *    - the definition of the core MS3DVAR CPP switches basically 
 *      shared among all scales. 
 * 4. The Specific MS3DVAR scale CPP switches such as :
 *    DAS_READ_INC,...
 * 5. At the very end, the mandatory final standard Croco cppdefs includes 
 *    of "OCEAN/croco/cppdefs_dev.h" and "set_global_definitions.h" 
 *    followed by the MS3DVAR mandatory final include 
 *    "ASSIM/MS3DVAR/COMMON/ms3dvar_set_global_def.h".
 *
 * where the $SCALE parameter refer to either FILTER, LR, MR or MS scale.
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
 * "ASSIM/MS3DVAR/${SCALE}/croco_cppdefs.h" with on the fly 
 * modifications introduced at compilation stage, namely : 
 *  - the final standard includes of "OCEAN/croco/cppdefs_dev.h" and 
 *    "set_global_definitions.h" are removed; because these instructions
 *    must be brought at the very end of the current final header 
 *    "ASSIM/MS3DVAR/$SCALE/cppdefs.h". This change is mandatory.
 *
 * SCALE = FILTER, LR, MR or MS
 *
 * See path "OCEAN/cppdefs.h"
 * Copied to "ASSIM/MS3DVAR/${SCALE}/Compile/croco_cppdefs.h"
 * 
 *=====================================================================*/
#include "croco_cppdefs.h"

/*=====================================================================
 * FILTER SCALE CPP-switch :
 *
 * The following MS3DVAR CPP key is positioning the scale 
 * for further use.
 *
 * DO NOT CHANGE or REMOVE, 
 * DO NOT MOVE THIS INSTRUCTION ACCROSS THIS FILE.
 *
 * The scale CPP key has to be set to either DAS_FILTER or DAS_LR, 
 * or DAS_MR or DAS_MS, according to the MS3DVAR SCALE currently
 * defined (here FILTER scale => DAS_FILTER ccp key)
 *
 * All the MS3DVAR CPP keys all start with pattern DAS_*
 *=====================================================================*/
#define DAS_FILTER

/*=====================================================================
 * CORE STANDARD MS3DVAR CONFIGURATION
 *
 * The next included file contains the core MS3DVAR CPP 
 * switches (~150 lines) which are shared among all scales
 * (FILTER, LR, MR and MS) enabling the definition of
 * part of the Multi-Scale 3D Variational DA scheme.
 *
 * All the MS3DVAR CPP keys all start with pattern DAS_*
 *
 * See path "ASSIM/MS3DVAR/COMMON/cppdefs_ms3dvar"
 * Copied to "ASSIM/MS3DVAR/${SCALE}/Compile/ms3dvar_cppdefs.h"
 *
 * Note that DAS cpp keys settings below this include 
 * will override the core defaults.
 *
 *=====================================================================*/
#include "ms3dvar_cppdefs.h"

/*=====================================================================
 * FILTER SCALE-SPECIFIC OVERRIDES
 *
 * Customized MS3DVAR DAS_* options associated to the MS scale.
 * The following MS3DVAR CPP keys override the defaults previously defined 
 * in through the included file "ASSIM/MS3DVAR/COMMON/ms3dvar_cppdefs.h"
 *
 *=====================================================================*/

/* No FILTER-specific overrides currently - using defaults from ms3dvar_cppdefs.h */

/*=====================================================================
 * END OF CONFIGURATION
 *=====================================================================*/

/*=====================================================================
 * DO NOT TOUCH THE BELOW INCLUDED FILES
 *
 * These files contain additional preprocessor logic and global
 * definitions. They must be included at the very end.
 *=====================================================================*/
#include "cppdefs_dev.h"
#include "set_global_definitions.h"
#include "ms3dvar_set_global_def.h"
