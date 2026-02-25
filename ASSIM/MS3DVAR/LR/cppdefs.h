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
 * LR VARIANT: CONFIGURATION NAME
 *
 * Define your regional configuration. This must match the
 * configuration in param.h
 *=====================================================================*/
#undef  COASTAL         /* COASTAL Applications */
#define REGIONAL        /* REGIONAL Applications */

#if defined REGIONAL
/*=====================================================================
 * REGIONAL CONFIGURATION: WMED (West Mediterranean)
 *=====================================================================*/
# define WMED

/*---------------------------------------------------------------------
 * PARALLELIZATION
 *---------------------------------------------------------------------*/
# define OPENMP         /* OpenMP parallelization */
# undef  MPI            /* MPI parallelization */

/*---------------------------------------------------------------------
 * NON-HYDROSTATIC AND NESTING
 *---------------------------------------------------------------------*/
# undef  NBQ            /* Non-Boussinesq/non-hydrostatic */
# undef  CROCO_QH       /* Quasi-hydrostatic */
# undef  AGRIF          /* Adaptive grid refinement */
# undef  AGRIF_2WAY     /* Two-way nesting */

/*---------------------------------------------------------------------
 * COUPLING
 *---------------------------------------------------------------------*/
# undef  OA_COUPLING    /* Ocean-Atmosphere coupling via OASIS */
# undef  OW_COUPLING    /* Ocean-Wave coupling via OASIS */

/*---------------------------------------------------------------------
 * OPEN BOUNDARY CONDITIONS
 *---------------------------------------------------------------------*/
# undef  TIDES          /* Tidal forcing */
# define OBC_EAST       /* Open boundary east */
# define OBC_WEST       /* Open boundary west */
# undef  OBC_NORTH      /* Open boundary north */
# undef  OBC_SOUTH      /* Open boundary south */

/*---------------------------------------------------------------------
 * APPLICATIONS
 *---------------------------------------------------------------------*/
# undef  BIOLOGY        /* Biology module */
# undef  FLOATS         /* Lagrangian floats */
# undef  STATIONS       /* Station output */
# undef  SEDIMENT       /* Sediment transport */
# undef  MUSTANG        /* Wave-sediment interactions */
# undef  BBL            /* Bottom boundary layer */

/*---------------------------------------------------------------------
 * TRACERS
 *---------------------------------------------------------------------*/
# define TEMPERATURE    /* Temperature tracer */
# define SALINITY       /* Salinity tracer (REQUIRED for MS3DVAR) */
# undef  PASSIVE_TRACER /* Passive tracers */

/*---------------------------------------------------------------------
 * EQUATION OF STATE
 *---------------------------------------------------------------------*/
# define NONLIN_EOS     /* Nonlinear equation of state */
# undef  SPLIT_EOS      /* Split equation of state */

/*---------------------------------------------------------------------
 * ADVECTION
 *---------------------------------------------------------------------*/
# define UV_ADV         /* Advection of momentum */
# undef  UV_COR         /* Coriolis term */

/*---------------------------------------------------------------------
 * LATERAL BOUNDARY CONDITIONS (if OBC defined)
 *---------------------------------------------------------------------*/
# ifdef OBC_EAST
#  define FRC_BRY       /* Forcing at boundaries */
#  ifdef FRC_BRY
#   define ANA_BRY      /* Analytical boundary conditions */
#   define Z_FRC_BRY    /* SSH forcing at boundaries */
#   define OBC_M2CHARACT /* M2 characteristic BCs */
#   define OBC_REDUCED_PHYSICS /* Reduced physics at boundaries */
#   define M2_FRC_BRY   /* 2D momentum forcing */
#   undef  M3_FRC_BRY   /* 3D momentum forcing */
#   define T_FRC_BRY    /* Tracer forcing */
#  endif
# endif

/*---------------------------------------------------------------------
 * VERTICAL MIXING
 *---------------------------------------------------------------------*/
# define GLS_MIXING     /* Generic Length Scale mixing */

/*---------------------------------------------------------------------
 * SOURCES/SINKS
 *---------------------------------------------------------------------*/
# define PSOURCE        /* Point sources/sinks */
# undef  PSOURCE_MASS   /* Include mass in point sources */
# define ANA_PSOURCE    /* Analytical point sources */

/*---------------------------------------------------------------------
 * GRID AND MASKING
 *---------------------------------------------------------------------*/
# define ANA_GRID       /* Analytical grid */
# define MASKING        /* Land/sea masking */

/*---------------------------------------------------------------------
 * INITIAL CONDITIONS
 *---------------------------------------------------------------------*/
# define ANA_INITIAL    /* Analytical initial conditions */
# define ZCLIMATOLOGY   /* Enable clmname for LR initial conditions file */

/*---------------------------------------------------------------------
 * SURFACE FORCING (Analytical for standalone testing)
 *---------------------------------------------------------------------*/
# define ANA_SMFLUX     /* Analytical surface momentum flux */
# define ANA_STFLUX     /* Analytical surface tracer flux */
# define ANA_SSFLUX     /* Analytical surface salinity flux */
# define ANA_SRFLUX     /* Analytical surface radiation flux */
# define ANA_BTFLUX     /* Analytical bottom temperature flux */
# define ANA_BSFLUX     /* Analytical bottom salinity flux */

/*---------------------------------------------------------------------
 * PERIODICITY
 *---------------------------------------------------------------------*/
# define NS_PERIODIC    /* North-South periodic */
# undef  EW_PERIODIC    /* East-West periodic */

/*---------------------------------------------------------------------
 * FILE OPTIONS
 *---------------------------------------------------------------------*/
# undef  FLOATS         /* Lagrangian floats */
# define NO_FRCFILE     /* No forcing file */

/*---------------------------------------------------------------------
 * CALENDAR AND TIME
 *---------------------------------------------------------------------*/
# define USE_CALENDAR   /* Use calendar for time management */

/*---------------------------------------------------------------------
 * COORDINATE SYSTEM
 *---------------------------------------------------------------------*/
# define NEW_S_COORD    /* New s-coordinate (Vtransform=2) */

#endif /* REGIONAL */

/*=====================================================================
 * INCLUDE MINIMAL MS3DVAR CONFIGURATION
 *
 * This includes the core MS3DVAR CPP defines (~150 lines)
 * Settings above will override defaults in cppdefs_ms3dvar.h
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
