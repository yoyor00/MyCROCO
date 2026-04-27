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
 * cppdefs_ms3dvar.h - MS3DVAR CPP Configuration
 * ==============================================
 *
 * This file contains ONLY the CPP defines needed by MS3DVAR.
 * It is independent of CROCO's full cppdefs.h to avoid confusion
 * and maintenance issues.
 *
 * The full CROCO cppdefs.h has 2203 lines with configurations for:
 * - Academic test cases (BASIN, CANYON, UPWELLING, etc.)
 * - Physics modules (BIOLOGY, PISCES, SEDIMENT, etc.)
 * - Coupling options (WRF, WW3, OASIS, etc.)
 *
 * MS3DVAR only needs ~5 core CROCO defines plus ~60 DAS_* specific defines.
 * This minimal file contains exactly what's needed.
 *
 * USAGE:
 * ------
 * Variant-specific files (LR/cppdefs.h, MR/cppdefs.h, etc.) should
 * include this file and override specific settings as needed.
 */

/*=====================================================================
 * CONFIGURATION SELECTION
 *
 * Define your configuration name here. This should match your
 * regional configuration in param_ms3dvar.h
 *=====================================================================*/
#undef  BASIN           /* Basin Example */
#undef  CANYON          /* Canyon Example */
#undef  UPWELLING       /* Upwelling Example */
#undef  REGIONAL        /* REGIONAL Applications */

/* Select your configuration - default to REGIONAL */
/* This will be overridden by variant-specific cppdefs.h */
#define REGIONAL

/*=====================================================================
 * REQUIRED CROCO CORE OPTIONS
 *
 * These options are always required for MS3DVAR to function.
 * DO NOT modify these unless you know what you're doing.
 *=====================================================================*/
#define SOLVE3D         /* 3D solution (REQUIRED for MS3DVAR) */
#define SALINITY        /* Salinity tracer (REQUIRED for MS3DVAR) */

/*=====================================================================
 * PARALLELIZATION OPTIONS
 *
 * Choose between OpenMP (shared memory) or MPI (distributed memory).
 * Note: DAS is incompatible with CROCO's AUTOTILING + OpenMP combination.
 *=====================================================================*/
#define OPENMP          /* OpenMP parallelization */
#undef  MPI             /* MPI parallelization */

/*=====================================================================
 * GRID AND COORDINATE OPTIONS
 *
 * Basic grid configuration for MS3DVAR.
 *=====================================================================*/
#define MASKING         /* Land/sea masking */
#define NEW_S_COORD     /* New s-coordinate system (Vtransform=2) */
#define SPHERICAL       /* Spherical (geographic) coordinate grid */
#define CURVGRID        /* Curvilinear (non-orthogonal) grid */

/*=====================================================================
 * BOUNDARY CONDITIONS
 *
 * Configure your domain boundaries.
 * These settings are typically overridden in variant-specific files.
 *=====================================================================*/
#undef  OBC_EAST        /* Open boundary east */
#undef  OBC_WEST        /* Open boundary west */
#undef  OBC_NORTH       /* Open boundary north */
#undef  OBC_SOUTH       /* Open boundary south */
#undef  EW_PERIODIC     /* East-West periodic boundaries */
#undef  NS_PERIODIC     /* North-South periodic boundaries */

/*=====================================================================
 * ANALYTICAL FUNCTIONS
 *
 * For standalone MS3DVAR testing without external forcing files.
 * Typically disabled when assimilating into real ocean model runs.
 *=====================================================================*/
#undef  ANA_GRID        /* Analytical grid */
#undef  ANA_INITIAL     /* Analytical initial conditions */
#undef  ANA_SMFLUX      /* Analytical surface momentum flux */
#undef  ANA_STFLUX      /* Analytical surface tracer flux */
#undef  ANA_SSFLUX      /* Analytical surface salinity flux */
#undef  ANA_SRFLUX      /* Analytical surface radiation flux */
#undef  ANA_BTFLUX      /* Analytical bottom temperature flux */
#undef  ANA_BSFLUX      /* Analytical bottom salinity flux */

/*=====================================================================
 * MIXING OPTIONS
 *
 * Vertical mixing scheme (if needed by MS3DVAR).
 *=====================================================================*/
#undef  GLS_MIXING      /* Generic Length Scale mixing */

/*=====================================================================
 * TIME MANAGEMENT
 *=====================================================================*/
#define USE_CALENDAR    /* Use calendar for time management */

/*=====================================================================
 * MS3DVAR DATA ASSIMILATION OPTIONS
 *
 * This section contains all DAS_* specific CPP defines for configuring
 * the Multi-Scale 3D Variational data assimilation system.
 *=====================================================================*/

/*---------------------------------------------------------------------
 * CORE DATA ASSIMILATION SETTINGS
 *---------------------------------------------------------------------*/
#define DAS                    /* Enable data assimilation (REQUIRED) */
#undef  DAS_READ_INC          /* Read increments from file */
#undef  DAS_DBLE_BKG          /* Use double background for single-scale */

/*---------------------------------------------------------------------
 * DYNAMICAL CONSTRAINTS
 *
 * Control how strongly the analysis respects dynamical balance.
 * Choose ONE geostrophic option and ONE hydrostatic option.
 *---------------------------------------------------------------------*/

/* Geostrophic Balance */
#define DAS_GEOS_STRONG       /* Strong geostrophic constraint (recommended) */
#undef  DAS_GEOS_WEAK         /* Weak geostrophic constraint */

/* Hydrostatic Balance */
#undef  DAS_HYDRO_STRONG      /* Strong hydrostatic constraint */
#undef  DAS_HYDRO_WEAK        /* Weak hydrostatic constraint */

/*---------------------------------------------------------------------
 * BACKGROUND ERROR COVARIANCE
 *
 * Configuration of background error statistics.
 *---------------------------------------------------------------------*/
#define DAS_BVAR_CORR         /* Background error correlations (recommended) */
#undef  DAS_ANA_BVAR          /* Analytical background variance (for testing) */

/*---------------------------------------------------------------------
 * OBSERVATION QUALITY CONTROL
 *
 * Options for filtering or weighting observations.
 *---------------------------------------------------------------------*/
#define DAS_DISCOAST          /* Discount coastal observations (recommended) */

/*---------------------------------------------------------------------
 * SATELLITE SEA SURFACE TEMPERATURE (SST) OBSERVATIONS
 *
 * Enable the satellite SST products you want to assimilate.
 *---------------------------------------------------------------------*/
#define DAS_MCSST             /* Multi-channel SST (AVHRR, MODIS, VIIRS) */
#define DAS_SSTMCMSK          /* SST multi-channel mask for quality control */
#undef  DAS_GOES_SST          /* GOES geostationary satellite SST */
#undef  DAS_TMISST            /* TMI (TRMM Microwave Imager) SST */
#undef  DAS_FDN_MCSST         /* Foundation temperature from multi-channel SST */

/*---------------------------------------------------------------------
 * SATELLITE SEA SURFACE HEIGHT (SSH) OBSERVATIONS
 *
 * Enable the altimetry products you want to assimilate.
 *---------------------------------------------------------------------*/
#define DAS_SWOTSSH           /* SWOT (Surface Water and Ocean Topography) */
#define DAS_SWOTSSHMSK        /* SWOT SSH mask for quality control */
#undef  DAS_JASONSSH          /* Jason altimeter series */
#undef  DAS_TPSSH             /* TOPEX/Poseidon altimeter */
#undef  DAS_ERS2SSH           /* ERS-2 altimeter */
#undef  DAS_GFOSSH            /* GFO altimeter */

/*---------------------------------------------------------------------
 * SATELLITE SEA SURFACE SALINITY (SSS) OBSERVATIONS
 *
 * Enable for SMOS, SMAP, or Aquarius SSS products.
 *---------------------------------------------------------------------*/
#undef  DAS_SATSSSS           /* Satellite sea surface salinity */

/*---------------------------------------------------------------------
 * IN-SITU OBSERVATIONS
 *
 * Enable in-situ observation platforms (ARGO, gliders, CTD, etc.)
 *---------------------------------------------------------------------*/
#define DAS_INSITU            /* Enable in-situ observations (master switch) */

/* Gliders */
#undef  DAS_WHOIGLIDER        /* WHOI (Woods Hole) gliders */
#undef  DAS_SIOGLIDER         /* SIO (Scripps) gliders */

/* CTD (Conductivity-Temperature-Depth) */
#undef  DAS_PTSURCDT          /* Point surface CTD */
#undef  DAS_MARTINCDT         /* Martin CTD */
#undef  DAS_MAPPEDCDT         /* Gridded/mapped CTD */
#undef  DAS_PRFCDT            /* Profile CTD */

/* Ships and Floats */
#undef  DAS_SHIPSST           /* Ship-based SST measurements */
#undef  DAS_NPSFLIGHT         /* NPS flight/aircraft observations */

/* Autonomous Underwater Vehicles (AUV) */
#undef  DAS_CALPOLYAUV        /* Cal Poly AUV */
#undef  DAS_DORADOAUV         /* Dorado AUV */

/*---------------------------------------------------------------------
 * CURRENT OBSERVATIONS
 *
 * Direct velocity measurements from HF radar, ADCP, etc.
 *---------------------------------------------------------------------*/
#undef  DAS_HFRADAR           /* HF radar surface currents */
#undef  DAS_CURRENT_UV        /* Direct UV current observations */
#undef  DAS_SYNTHRADAR        /* Synthetic radar for testing */
#undef  DAS_HUGEHRADAR        /* Support for very large radar datasets */
#undef  DAS_HFHEADER          /* HF radar header processing */

/*---------------------------------------------------------------------
 * MAPPED/GRIDDED OBSERVATIONS
 *
 * For pre-gridded observation products.
 *---------------------------------------------------------------------*/
#undef  DAS_MAPPED            /* Mapped/gridded observations */

/*---------------------------------------------------------------------
 * PROCESSING AND ANALYSIS OPTIONS
 *
 * Advanced options for tuning the analysis procedure.
 *---------------------------------------------------------------------*/
#define DAS_PWTHZETATOT       /* Weight total sea surface height */
#define DAS_ZETAHADJCORR      /* SSH adjustment with correlations */
#define DAS_REDUCED_ZCORR     /* Reduced vertical correlation lengths */
#undef  DAS_BKGZETA           /* Use background SSH in analysis */
#undef  DAS_PSICHI_COSTREG    /* Streamfunction cost regularization */
#undef  DAS_ADIMPSICHI        /* Adimensional psi/chi formulation */
#undef  DAS_COUNT_PSICHI      /* Count psi/chi observations for diagnostics */
#undef  DAS_UVAGEO_DIAG       /* Diagnose ageostrophic velocity components */
#undef  DAS_D2CUV             /* Second-order current analysis */

/*---------------------------------------------------------------------
 * DIAGNOSTIC OUTPUT
 *
 * Options for debugging and detailed output.
 *---------------------------------------------------------------------*/
#undef  DAS_VERBOSE_INSITU    /* Verbose output for in-situ observations */
#undef  DAS_WRT_ADJ_INI       /* Write adjoint initial state */
#undef  DAS_WRT_ADJ_END       /* Write adjoint final state */

/*=====================================================================
 * INCLUDE ADDITIONAL DEFINITIONS
 *
 * These files contain additional preprocessor logic and global
 * definitions. They must be included last.
 *=====================================================================*/
#include "cppdefs_dev.h"
#include "set_global_definitions.h"
#include "das_set_global_def.h"
