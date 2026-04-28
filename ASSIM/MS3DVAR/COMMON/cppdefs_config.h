/*=====================================================================
 * REGIONAL CONFIGURATION: WMED (West Mediterranean)
 *=====================================================================*/

                      /* Configuration Name */
# define WMED
                      /* Parallelization */
# undef  OPENMP
# undef  MPI
                      /* Non-hydrostatic option */
# undef  NBQ
# undef  CROCO_QH
                      /* Nesting */
# undef  AGRIF
# undef  AGRIF_2WAY
                      /* OA and OW Coupling via OASIS (MPI) */
# undef  OA_COUPLING
# ifdef OA_COUPLING
#  define READ_PATM
#  define OBC_PATM
#  undef  OA_GRID_UV
# endif
# undef  OW_COUPLING
# ifdef OW_COUPLING
#  undef OW_COUPLING_FULL
#  undef WAVE_SMFLUX
# endif
                      /* Wave-current interactions */
# undef  MRL_WCI
                      /* Open Boundary Conditions */
# undef  TIDES
# define OBC_EAST
# define OBC_WEST
# undef OBC_NORTH
# undef OBC_SOUTH
                      /* Applications */
# undef  BIOLOGY
# undef  STATIONS
# undef  PASSIVE_TRACER
# undef  SEDIMENT
# undef  MUSTANG
# undef  BBL
                      /* I/O server */
# undef  XIOS
                      /* Calendar */
# undef  USE_CALENDAR
                      /* dedicated croco.log file */
# define LOGFILE
/*
!-------------------------------------------------
! PRE-SELECTED OPTIONS
!
! ADVANCED OPTIONS ARE IN CPPDEFS_DEV.H
!-------------------------------------------------
*/
                      /* Parallelization */
# ifdef MPI
#  undef  PARALLEL_FILES
#  undef  NC4PAR
#  undef  MPI_NOLAND
#  undef  MPI_TIME
# endif
                      /* Grid configuration */
# define CURVGRID
# define SPHERICAL
# define MASKING
# undef  WET_DRY
# define NEW_S_COORD
                      /* Model dynamics */
# define SOLVE3D
# define UV_COR
# define UV_ADV
                      /* Equation of State */
# define SALINITY
# define NONLIN_EOS
                      /* Lateral Forcing */
# undef CLIMATOLOGY
# ifdef CLIMATOLOGY
#  define ZCLIMATOLOGY
#  define M2CLIMATOLOGY
#  define M3CLIMATOLOGY
#  define TCLIMATOLOGY

#  define ZNUDGING
#  define M2NUDGING
#  define M3NUDGING
#  define TNUDGING
#  undef  ROBUST_DIAG
# endif

# define FRC_BRY
# ifdef FRC_BRY
#  define Z_FRC_BRY
#  define M2_FRC_BRY
#  define M3_FRC_BRY
#  define T_FRC_BRY
# endif
                      /* Surface Forcing */
/*
! Bulk flux algorithms (options)
! by default : COARE3p0 paramet with GUSTINESS effects
!
! To change bulk param, define one the following keys (exclusive) :
! - define BULK_ECUMEV0 : ECUME_v0 param
! - define BULK_ECUMEV6 : ECUME_v6 param
! - define BULK_WASP    : WASP param
! Note : gustiness effects can be added for all params
!        by defining BULK_GUSTINESS
*/
# define BULK_FLUX
# ifdef BULK_FLUX
#  undef  BULK_ECUMEV0
#  undef  BULK_ECUMEV6
#  undef  BULK_WASP
#  define BULK_GUSTINESS
#  define BULK_LW
#  undef  SST_SKIN
#  undef  ANA_DIURNAL_SW
#  undef  ONLINE
#  ifdef ONLINE
#   undef  AROME
#   undef  ERA_ECMWF
#  endif
#  undef READ_PATM
#  ifdef READ_PATM
#   define OBC_PATM
#  endif
#  undef  ABL1D
#  ifdef  ABL1D
#   undef  ANA_ABL_LSDATA
#   undef  ANA_ABL_VGRID
#   define STRESS_AT_RHO_POINTS
#   define ABL_NUDGING
#   define ABL_NUDGING_DYN
#   define ABL_NUDGING_TRA
#   undef  ABL_DYN_RESTORE_EQ
#   undef  SFLUX_CFB
#  endif
# else
#  define QCORRECTION
#  define SFLX_CORR
#  undef  SFLX_CORR_COEF
#  define ANA_DIURNAL_SW
# endif
# define  SFLUX_CFB
# undef  SEA_ICE_NOFLUX
                      /* Lateral Momentum Advection (default UP3) */
# define UV_HADV_UP3
# undef  UV_HADV_UP5
# undef  UV_HADV_WENO5
                      /* Lateral Explicit Momentum Mixing */
# undef  UV_VIS2
# ifdef UV_VIS2
#  define UV_VIS_SMAGO
# endif
                      /* Vertical Momentum Advection */
# define UV_VADV_SPLINES
# undef  UV_VADV_WENO5
                      /* Lateral Tracer Advection (default UP3) */
# undef  TS_HADV_UP3
# define TS_HADV_RSUP3
# undef  TS_HADV_UP5
# undef  TS_HADV_WENO5
                      /* Lateral Explicit Tracer Mixing */
# undef  TS_DIF2
# undef  TS_DIF4
# undef  TS_MIX_S
                      /* Vertical Tracer Advection  */
# define TS_VADV_SPLINES
# undef  TS_VADV_AKIMA
# undef  TS_VADV_WENO5
                      /* Sponge layers for UV and TS */
# define SPONGE
                      /* Semi-implicit Vertical Tracer/Mom Advection */
# define  VADV_ADAPT_IMP
                      /* Bottom friction in fast 3D step */
# define LIMIT_BSTRESS
# undef  BSTRESS_FAST
                      /* Vertical Mixing */
# undef  BODYFORCE
# define LMD_MIXING
# undef  GLS_MIXING
# ifdef LMD_MIXING
#  define LMD_SKPP
#  define LMD_BKPP
#  define LMD_RIMIX
#  define LMD_CONVEC
#  define LMD_NONLOCAL
#  undef  LMD_DDMIX
#  undef  LMD_LANGMUIR
# endif
                      /* Wave-current interactions */
# ifdef OW_COUPLING
#  define MRL_WCI
#  define BBL
# endif
# ifdef MRL_WCI
#  ifndef OW_COUPLING
#   undef  WAVE_OFFLINE
#   define ANA_WWAVE
#   undef  WKB_WWAVE
#  endif
#  undef  WAVE_ROLLER
#  define WAVE_STREAMING
#  define WAVE_FRICTION
#  define WAVE_RAMP
#  ifdef WKB_WWAVE
#   undef  WKB_OBC_NORTH
#   undef  WKB_OBC_SOUTH
#   define WKB_OBC_WEST
#   undef  WKB_OBC_EAST
#  endif
# endif
                      /* Bottom Forcing */
# define ANA_BSFLUX
# define ANA_BTFLUX
                      /* Point Sources - Rivers */
# define PSOURCE
# define PSOURCE_NCFILE
# ifdef PSOURCE_NCFILE
#  undef PSOURCE_NCFILE_TS
# endif
                      /* Open Boundary Conditions */
# ifdef TIDES
#  define SSH_TIDES
#  define UV_TIDES
#  define POT_TIDES
#  undef  TIDES_MAS
#  define TIDERAMP
# endif
# define OBC_M2CHARACT
# undef  OBC_M2ORLANSKI
# define OBC_M3ORLANSKI
# define OBC_TORLANSKI
# undef  OBC_M2SPECIFIED
# undef  OBC_M3SPECIFIED
# undef  OBC_TSPECIFIED
                      /* Input/Output */
# define AVERAGES
# define AVERAGES_K
# undef  OUTPUTS_SURFACE
# undef  HOURLY_VELOCITIES
                     /* Exact restart */
# undef EXACT_RESTART
/*
!                        Diagnostics
!--------------------------------------------
! 3D Tracer & momentum balance
! 2D Mixing layer balance
! Depth-mean vorticity and energy balance
! Eddy terms
!--------------------------------------------
!
*/
# undef DO_NOT_OVERWRITE
# undef RESTART_DIAGS

# undef DIAGNOSTICS_TS
# undef DIAGNOSTICS_UV
# ifdef DIAGNOSTICS_TS
#  undef  DIAGNOSTICS_TS_ADV
#  undef  DIAGNOSTICS_TS_MLD
#  ifdef DIAGNOSTICS_TS_MLD
#   define DIAGNOSTICS_TS_MLD_CRIT
#  endif

# endif

# undef DIAGNOSTICS_TSVAR
# ifdef DIAGNOSTICS_TSVAR
#  define DIAGNOSTICS_TS
#  define DIAGNOSTICS_TS_ADV
# endif

# undef  DIAGNOSTICS_VRT
# undef  DIAGNOSTICS_EK
# ifdef DIAGNOSTICS_EK
#  undef DIAGNOSTICS_EK_FULL
#  undef DIAGNOSTICS_EK_MLD
# endif

# undef DIAGNOSTICS_BARO
# undef DIAGNOSTICS_PV
# undef DIAGNOSTICS_DISS
# ifdef DIAGNOSTICS_DISS
#  define DIAGNOSTICS_PV
# endif

# undef DIAGNOSTICS_EDDY

# undef TENDENCY
# ifdef TENDENCY
#  define DIAGNOSTICS_UV
# endif
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
                                               || defined MUSTANG
#  define BIO_HADV_WENO5
# endif
                      /*   Choice of Biology models   */
# ifdef BIOLOGY
#  define PISCES
#  undef  BIO_NChlPZD
#  undef  BIO_N2ChlPZD2
#  undef  BIO_BioEBUS
                      /*   Biology options    */
#  ifdef PISCES
#   undef  DIURNAL_INPUT_SRFLX
#   define key_pisces
#   define key_ligand
#   undef key_pisces_quota
#   undef key_pisces_npzd
#   undef key_sediment
#  endif
#  ifdef BIO_NChlPZD
#   define OXYGEN
#  endif
#  ifdef BIO_BioEBUS
#   define NITROUS_OXIDE
#  endif
                      /*   Biology diagnostics    */
#  define DIAGNOSTICS_BIO
#  if defined DIAGNOSTICS_BIO && defined PISCES
#   define key_trc_diaadd
#  endif
# endif
                      /*   Stations recording    */
# ifdef STATIONS
#  define ALL_SIGMA
# endif
                      /*   USGS Sediment model     */
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
#  define key_sand2D
#  define MUSTANG_CORFLUX
#  undef  key_tauskin_c_upwind
#  undef  WAVE_OFFLINE
# endif
