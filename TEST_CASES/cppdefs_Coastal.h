# define COASTAL
# define VILAINE
                      /* Parallelization */
# undef  OPENMP
# undef  MPI
                      /* Open Boundary Conditions */
# define TIDES
# undef  OBC_EAST
# define OBC_WEST
# undef  OBC_NORTH
# define OBC_SOUTH
                      /* Applications */
# undef  BIOLOGY
# undef  STATIONS
# undef  PASSIVE_TRACER
# undef  SEDIMENT
# undef  BBL
# define MUSTANG
                      /* I/O server */
# undef  XIOS
                     /* Custion IO */
# undef  FILLVAL
                      /* Calendar */

# define USE_CALENDAR
                      /* dedicated croco.log file */
# undef  LOGFILE
/*!
!-------------------------------------------------
! PRE-SELECTED OPTIONS
!
! ADVANCED OPTIONS ARE IN CPPDEFS_DEV.H
!-------------------------------------------------
*/
                      /* Parallelization */
# ifdef MPI
#  define NC4PAR
#  undef  MPI_NOLAND
#  undef  MPI_TIME
# endif
                      /* Non-hydrostatic options */
# ifdef NBQ
#  define W_HADV_WENO5
#  define W_VADV_WENO5
# endif
                      /* Grid configuration */
# define ANA_INITIAL
# define CURVGRID
# define SPHERICAL
# define MASKING
# define WET_DRY
# define NEW_S_COORD
                      /* Model dynamics */
# define SOLVE3D
# define UV_COR
# define UV_ADV
                      /* Equation of State */
# define SALINITY
# define NONLIN_EOS
                      /* Lateral Momentum Advection (default UP3) */
# undef  UV_HADV_UP3
# define UV_HADV_WENO5
                      /* Lateral Explicit Momentum Mixing */
# define UV_VIS2
# ifdef UV_VIS2
#  define UV_VIS_SMAGO
# endif
                      /* Vertical Momentum Advection  */
# undef  UV_VADV_SPLINES
# define UV_VADV_WENO5
                      /* Lateral Tracer Advection (default UP3) */
# undef  TS_HADV_UP3
# define TS_HADV_WENO5
                      /* Lateral Explicit Tracer Mixing */
# define TS_DIF2
# define TS_MIX_S
                      /* Vertical Tracer Advection  */
# undef  TS_VADV_SPLINES
# define TS_VADV_WENO5
                      /* Sponge layers for UV and TS */
# define SPONGE
                      /* Semi-implicit Vertical Tracer/Mom Advection */
# undef  VADV_ADAPT_IMP
                      /* Bottom friction in fast 3D step */
# define LIMIT_BSTRESS
# undef  BSTRESS_FAST
                      /* Vertical Mixing */
# define GLS_MIXING
                      /* Surface Forcing */
# define BULK_FLUX
# ifdef BULK_FLUX
#  undef  BULK_ECUMEV0
#  undef  BULK_ECUMEV6
#  undef  BULK_WASP
#  define BULK_GUSTINESS
#  undef  BULK_LW
#  undef  SST_SKIN
#  undef  ANA_DIURNAL_SW
#  define ONLINE
#  ifdef ONLINE
#   define AROME
#   undef  ERA_ECMWF
#  endif
#  define READ_PATM
#  ifdef READ_PATM
#   define OBC_PATM
#  endif
# else
#  undef  QCORRECTION
#  undef  SFLX_CORR
#  undef  SFLX_CORR_COEF
#  undef  ANA_DIURNAL_SW
# endif
# undef  ANA_SSFLUX
# define ANA_STFLUX
                      /* Lateral Forcing */
# undef  ANA_BRY
# define FRC_BRY
# ifdef FRC_BRY
#  define Z_FRC_BRY
#  define M2_FRC_BRY
#  undef  M3_FRC_BRY
#  define T_FRC_BRY
# endif
                      /* Bottom Forcing */
# define ANA_BSFLUX
# define ANA_BTFLUX
                      /* Point Sources - Rivers */
# define PSOURCE
# undef  PSOURCE_MASS
# define PSOURCE_NCFILE
# ifdef PSOURCE_NCFILE
#  define PSOURCE_NCFILE_TS
# endif
                      /* Open Boundary Conditions */
# ifdef TIDES
#  define M2FILTER_NONE
#  define SSH_TIDES
#  define UV_TIDES
#  undef  POT_TIDES
#  define TIDES_MAS
#  ifndef UV_TIDES
#   define OBC_REDUCED_PHYSICS
#  endif
#  define TIDERAMP
# endif
# define OBC_M2CHARACT
# define OBC_M3ORLANSKI
# define OBC_TORLANSKI
                      /* Input/Output */
# undef  AVERAGES
# undef  AVERAGES_K
# undef  OUTPUTS_SURFACE
# undef  HOURLY_VELOCITIES
/*
!           Applications:
!---------------------------------
! Biology, Stations,
! Passive tracer, Sediments, BBL
!---------------------------------
!
   Quasi-monotone lateral advection scheme (WENO5)
   for passive/biology/sediment tracers
*/
# if defined PASSIVE_TRACER || defined BIOLOGY || defined SEDIMENT \
                            || defined SUBSTANCE || defined MUSTANG
#  define BIO_HADV_WENO5
# endif
                      /*     USGS Sediment model     */
# ifdef SEDIMENT
#  define SUSPLOAD
#  define BEDLOAD
#  define MORPHODYN
# endif
                      /*   MUSTANG Sediment model     */
# ifdef MUSTANG
#  undef  key_MUSTANG_V2
#  undef  key_MUSTANG_bedload
#  undef  MORPHODYN
#  undef  key_tauskin_c_upwind
#  define WAVE_OFFLINE
# endif
/*
!
!==========================================================
!              IDEALIZED CONFIGURATIONS
!==========================================================
!
*/

#include "cppdefs_dev.h"
#include "set_global_definitions.h"
