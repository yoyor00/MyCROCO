!======================================================================
! CROCO is derived from the ROMS-AGRIF branch of ROMS.
! ROMS-AGRIF was developed by IRD and Inria. CROCO also inherits
! from the UCLA branch (Shchepetkin et al.) and the Rutgers
! University branch (Arango et al.), both under MIT/X style license.
! Copyright (C) 2005-2026 CROCO Development Team
! License: CeCILL-2.1 - see LICENSE.txt
!
! CROCO website : https://www.croco-ocean.org
!======================================================================
!
/*
   This is "cppdefs_full.h": FULLY DOCUMENTED MODEL CONFIGURATION FILE
   ==== == ================= ====== ========== ===== ============= ====
!
!  This file lists every user-facing CPP key with a short description.
!  It serves as a reference companion to croco_full.nml.
!  Default values mirror OCEAN/cppdefs.h (REGIONAL configuration).
!
!  Advanced options (rarely changed) are in cppdefs_dev.h.
!  Keys automatically derived from others should NOT be changed here.
*/

/*
!======================================================================
!  CONFIGURATION NAME
!  ------------------
!  Exactly one configuration-case key must be defined.
!  Defining more than one is an error.
!======================================================================
*/

/* ── Realistic configurations ─────────────────────────────────────── */
# define REGIONAL        /* Realistic regional simulations              */
# undef  COASTAL         /* Coastal configuration (alternative to REGIONAL) */

/* ── Classical idealised cases ────────────────────────────────────── */
# undef  UPWELLING       /* Upwelling example                           */
# undef  BASIN           /* Basin example                               */
# undef  CANYON          /* Canyon example                              */
# ifdef  CANYON
#  undef  CANYON_STRAT   /* Stratified canyon variant                   */
# endif
# undef  EQUATOR         /* Equatorial channel example                  */
# undef  GRAV_ADJ        /* Gravitational Adjustment example (also with NBQ) */
# undef  ACOUSTIC        /* Acoustic wave example (requires NBQ)        */
# undef  INNERSHELF      /* Inner Shelf example                         */
# ifdef  INNERSHELF
#  define INNERSHELF_EKMAN /* Ekman layer variant                       */
#  define INNERSHELF_APG   /* Atmospheric pressure gradient variant     */
# endif
# undef  OVERFLOW        /* Gravitational/Overflow example              */
# undef  SEAMOUNT        /* Seamount example                            */
# undef  SHELFRONT       /* Shelf Front example                         */
# undef  SOLITON         /* Equatorial Rossby Wave example              */
# undef  INTERNAL        /* Internal tides example                      */
# ifdef  INTERNAL
#  undef  BODYTIDE       /* Tidal body force variant                    */
# endif
# undef  VORTEX          /* Baroclinic Vortex example                   */
# undef  JET             /* Jet example                                 */
# undef  THACKER         /* Thacker analytical solution example         */
# ifdef  THACKER
#  undef  THACKER_2DV    /* 2D vertical plane variant                   */
# endif
# undef  TANK            /* Tank example                                */
# undef  ISOLITON        /* Nonlinear internal solitary wave example    */
# undef  IGW             /* Internal Gravity Wave example (also with NBQ) */
# ifdef  IGW
#  undef  EXPERIMENT3    /* 3rd experiment variant (default)            */
#  undef  ONLINE_ANALYSIS /* Online spectral analysis output            */
# endif
# undef  KILPATRICK      /* SST fronts and island wakes example         */
# undef  KH_INST         /* Kelvin-Helmholtz shear instability example  */
# ifdef  KH_INST
#  undef  KH_INSTY       /* Y-direction (2D horizontal) variant         */
#  undef  KH_INST3D      /* 3D variant                                  */
# endif

/* ── Nearshore / wave cases ───────────────────────────────────────── */
# undef  RIP             /* Rip current example                         */
# ifdef  RIP
#  undef  RIP_TOPO_2D    /* Simplified 2D cross-shore topography        */
#  undef  BISCA          /* Biscarrosse beach config (no rip channel)  */
#  undef  GRANDPOPO      /* Grand Popo beach initial conditions         */
# endif
# undef  FLASH_RIP       /* Flash rip current example                   */
# undef  SHOREFACE       /* Shoreface cross-shore dynamics example      */
# undef  SWASH           /* Swash zone wave runup example               */
# ifdef  SWASH
#  define SWASH_GLOBEX_B3  /* GLOBEX B3 beach profile (default)        */
#  undef  SWASH_GLOBEX_B2  /* GLOBEX B2 beach profile                  */
#  undef  SWASH_GLOBEX_A3  /* GLOBEX A3 beach profile                  */
# endif
# undef  SANDBAR         /* Sandbar cross-shore dynamics (with MRL_WCI) */
# ifdef  SANDBAR
#  undef  SANDBAR_OFFSHORE /* Offshore bar variant                      */
#  undef  SANDBAR_ONSHORE  /* Onshore bar variant                       */
# endif

/* ── Sediment / morphodynamics cases ─────────────────────────────── */
# undef  DUNE            /* 2D dune migration (SEDIMENT or MUSTANG)     */
# ifdef  DUNE
#  undef  ANA_DUNE        /* Analytical dune initial conditions          */
#  undef  DUNE3D          /* 3D dune variant                             */
# endif
# undef  TIDAL_FLAT      /* Tidal flat morphodynamics example           */
# undef  ESTUARY         /* Estuary sediment dynamics example           */
# undef  RIVER           /* River / point-source test case              */
# undef  SEAGRASS        /* Seagrass bed wave-damping example           */
# undef  SED_TOY         /* Sediment 0D/1D toy model                   */
# ifdef  SED_TOY
#  undef  SED_TOY_ROUSE      /* Rouse equilibrium profile variant       */
#  undef  SED_TOY_RESUSP     /* Resuspension variant                    */
#  undef  SED_TOY_CONSOLID   /* Consolidation variant                   */
#  undef  SED_TOY_FLOC_0D    /* 0D flocculation variant                 */
#  undef  SED_TOY_FLOC_1D    /* 1D flocculation variant                 */
# endif
# undef  MOVING_BATHY    /* Moving bathymetry (NBQ) example             */

/* ── Single-column mixing cases ──────────────────────────────────── */
/* Use SINGLE_COLUMN + one forcing flag                               */
# undef  SINGLE_COLUMN   /* Single-column model base key                */
# ifdef  SINGLE_COLUMN
#  undef  KATO_PHILLIPS        /* Kato-Phillips mixed-layer forcing      */
#  undef  WILLIS_DEARDORFF     /* Willis-Deardorff convection forcing    */
#  undef  DIURNAL_CYCLE        /* Diurnal shortwave heating cycle        */
#  undef  FORCED_DBLEEK        /* Forced double-Ekman layer              */
#  undef  FORCED_EKBBL         /* Forced Ekman bottom boundary layer     */
#  undef  FORCED_NONROTBBL     /* Forced non-rotating bottom BBL         */
#  undef  FORCED_OSCNONROTBBL  /* Oscillating non-rotating BBL           */
# endif

/* ── Tracer advection verification cases ─────────────────────────── */
/* Use TS_HADV_TEST + one variant flag                                */
# undef  TS_HADV_TEST    /* Tracer advection test base key              */
# ifdef  TS_HADV_TEST
#  undef  SOLID_BODY_PER  /* Rotating solid-body (repeating)            */
#  undef  SOLID_BODY_ROT  /* Rotating solid-body (one rotation)         */
#  undef  DIAGONAL_ADV    /* Diagonal advection variant                  */
                          /* Without any variant: body-force case        */
# endif


/*
!======================================================================
!  PARALLELIZATION
!  ---------------
!  Only one of OPENMP / MPI may be defined at a time.
!  XIOS requires MPI.
!======================================================================
*/
# undef  OPENMP          /* Open-MP shared-memory parallelism           */
# undef  MPI             /* MPI distributed-memory parallelism          */
# undef  OPENACC
# ifdef MPI
   /* Activated automatically when MPI is defined – do not edit below  */
#  undef  PARALLEL_FILES  /* Parallel I/O: each MPI task writes its own file */
#  undef  NC4PAR          /* NetCDF-4 parallel I/O (requires HDF5/NC4)   */
#  undef  MPI_NOLAND      /* Skip land-only tiles (load balancing)       */
#  undef  MPI_TIME        /* Time MPI communications                     */
# endif

/*
!======================================================================
!  I/O SERVER
!======================================================================
*/
# undef  XIOS            /* Use the external XIOS I/O server (needs MPI) */

/*
!======================================================================
!  CALENDAR
!======================================================================
*/
# undef  USE_CALENDAR    /* Absolute date/time mode; requires start date
                            in croco.nml &croco_time_stepping            */

/*
!======================================================================
!  LOGGING
!======================================================================
*/
# undef  LOGFILE         /* Redirect STDOUT to a dedicated croco.log file */

/*
!======================================================================
!  NON-HYDROSTATIC OPTION
!======================================================================
*/
# undef  NBQ             /* Non-Boussinesq (non-hydrostatic) solver     */
# ifdef  NBQ
#  undef  NBQ_PRECISE    /* More precise NBQ time-stepping (slower)     */
# endif
# undef  CROCO_QH        /* Quasi-hydrostatic option (advanced)         */

/*
!======================================================================
!  NESTING (AGRIF)
!  ---------------
!  AGRIF requires a separate compilation with the AGRIF pre-processor.
!======================================================================
*/
# undef  AGRIF           /* Adaptive Mesh Refinement / nesting (1-way by default) */
# undef  AGRIF_2WAY      /* 2-way nesting: update parent solution from child */

/*
!======================================================================
!  COUPLING
!  --------
!  OA_COUPLING and OW_COUPLING both require OASIS-MCT and MPI.
!======================================================================
*/
# undef  OA_COUPLING     /* Ocean-Atmosphere coupling via OASIS-MCT     */
# ifdef OA_COUPLING
   /* Derived from OA_COUPLING – do not edit manually                   */
#  define READ_PATM      /* Read atmospheric pressure from coupling field */
#  define OBC_PATM       /* Apply atmospheric pressure at open boundaries */
#  undef  OA_GRID_UV     /* Interpolate coupling fields on U/V grids     */
# endif

# undef  OW_COUPLING     /* Ocean-Wave coupling via OASIS-MCT (WAVEWATCH III) */
# ifdef OW_COUPLING
#  undef  OW_COUPLING_FULL /* Full wave-current coupling (otherwise one-way) */
#  undef  WAVE_SMFLUX    /* Apply wave-induced surface momentum flux     */
# endif

/*
!======================================================================
!  OPEN BOUNDARY CONDITIONS – DIRECTION FLAGS
!  ------------------------------------------
!  Define the sides of the domain that have open (radiative) boundaries.
!======================================================================
*/
# undef  OBC_EAST        /* Open eastern  boundary                      */
# define OBC_WEST        /* Open western  boundary                      */
# define OBC_NORTH       /* Open northern boundary                      */
# define OBC_SOUTH       /* Open southern boundary                      */

/*
!======================================================================
!  TIDES
!  -----
!  TIDES activates the tidal forcing machinery; the sub-keys select
!  which tidal fields are processed.
!======================================================================
*/
# undef  TIDES           /* Activate tidal forcing at open boundaries   */
# ifdef TIDES
#  define SSH_TIDES      /* Read and apply tidal sea-level data         */
#  define UV_TIDES       /* Read and apply tidal current data           */
#  define POT_TIDES      /* Include tidal potential forcing             */
#  undef  TIDES_MAS      /* Tidal mask (advanced)                       */
#  define TIDERAMP       /* Ramp tidal forcing over 1 day at start      */
# endif

/*
!======================================================================
!  APPLICATIONS
!======================================================================
*/
# undef  BIOLOGY         /* Biogeochemical modeling (choose model below) */
# undef  STATIONS        /* High-frequency output at fixed station points */
# undef  PASSIVE_TRACER  /* Add one passive (non-reactive) tracer       */
# undef  SEDIMENT        /* USGS sediment transport model               */
# undef  MUSTANG         /* MUSTANG sediment dynamics model             */
# undef  BBL             /* Bottom Boundary Layer parameterization      */
# undef  OBSTRUCTION     /* Sub-grid flow obstruction (e.g. seagrass)   */
# undef  FILLVAL         /* Write fill value in masked land points      */

/*
!======================================================================
!  STOCHASTIC / ENSEMBLE
!======================================================================
*/
# undef  STOGEN          /* Stochastic forcing generator                */
# undef  ENSEMBLE        /* Ensemble mode                               */

/*
!======================================================================
!  GRID CONFIGURATION
!======================================================================
*/
# define CURVGRID        /* Curvilinear coordinate transformation       */
# define SPHERICAL       /* Longitude/latitude grid positioning         */
# define MASKING         /* Land masking                                */
# undef  WET_DRY         /* Wetting-and-drying scheme                   */
# define NEW_S_COORD     /* New (Song & Haidvogel 1994) vertical S-coord */
# undef  EW_PERIODIC     /* East-West periodic boundary conditions      */
# undef  NS_PERIODIC     /* North-South periodic boundary conditions    */

/*
!======================================================================
!  MODEL DYNAMICS
!======================================================================
*/
# define SOLVE3D         /* Solve 3D primitive equations (baroclinic)   */
# define UV_COR          /* Coriolis terms in momentum equations        */
# define UV_ADV          /* Advection terms in momentum equations       */
# undef  NO_FRCFILE      /* Run without any external forcing files      */
# undef  NO_TEMPERATURE  /* Suppress temperature tracer                 */
# undef  NO_TRACER       /* Suppress all tracers (barotropic only)     */
# undef  M2FILTER_NONE   /* Disable barotropic time filter
                            (set automatically with NBQ or M3FAST)     */
# undef  ZONAL_NUDGING   /* Zonal mean flow nudging (used with JET)    */

/*
!======================================================================
!  EQUATION OF STATE
!======================================================================
*/
# define SALINITY        /* Salinity as an active tracer                */
# define NONLIN_EOS      /* Nonlinear equation of state                 */
# undef  SPLIT_EOS       /* Split EOS into adiabatic + compressible parts
                            to reduce pressure-gradient errors          */

/*
!======================================================================
!  LATERAL FORCING (CLIMATOLOGY & OPEN BOUNDARY)
!  -----------------------------------------------
!  CLIMATOLOGY activates processing of 3D climatological files used
!  for interior nudging.  FRC_BRY activates open-boundary forcing files.
!======================================================================
*/
# undef  CLIMATOLOGY     /* Process 2D/3D climatological data for nudging
                            and open boundary forcing                   */
# ifdef CLIMATOLOGY
#  define ZCLIMATOLOGY   /* Process sea-level climatology               */
#  define M2CLIMATOLOGY  /* Process barotropic velocity climatology     */
#  define M3CLIMATOLOGY  /* Process baroclinic velocity climatology     */
#  define TCLIMATOLOGY   /* Process tracer (T, S, ...) climatology      */
   /* Nudging towards climatology in the interior                       */
#  define ZNUDGING       /* Nudging layer for sea level                 */
#  define M2NUDGING      /* Nudging layer for barotropic velocities     */
#  define M3NUDGING      /* Nudging layer for baroclinic velocities     */
#  define TNUDGING       /* Nudging layer for tracers                   */
#  undef  ROBUST_DIAG    /* Strong tracer nudging over the whole domain
                            (diagnostic / relaxation simulations)       */
# endif

# define FRC_BRY         /* Process 1D/2D boundary forcing files       */
# ifdef FRC_BRY
#  define Z_FRC_BRY      /* Open-boundary forcing for sea level         */
#  define M2_FRC_BRY     /* Open-boundary forcing for barotropic velocities */
#  define M3_FRC_BRY     /* Open-boundary forcing for baroclinic velocities */
#  define T_FRC_BRY      /* Open-boundary forcing for tracers           */
# endif

/*
!======================================================================
!  SURFACE FORCING
!  ---------------
!  BULK_FLUX: compute air-sea fluxes from atmospheric variables.
!  Default parameterization is COARE3p0 with gustiness.
!  Without BULK_FLUX: use prescribed flux fields + QCORRECTION.
!======================================================================
*/
# define BULK_FLUX        /* Bulk air-sea flux formulation              */
# ifdef BULK_FLUX
   /* Bulk parameterization choice (mutually exclusive; default COARE3p0) */
#  undef  BULK_ECUMEV0   /* ECUME_v0 parameterization                   */
#  undef  BULK_ECUMEV6   /* ECUME_v6 parameterization                   */
#  undef  BULK_WASP      /* WASP parameterization                       */
   /* Bulk flux options */
#  define BULK_GUSTINESS  /* Add gustiness correction to wind speed     */
#  define BULK_LW         /* Long-wave radiation feedback from model SST */
#  undef  SST_SKIN        /* Prognostic skin SST (Zeng & Beljaars 2005) */
#  undef  ANA_DIURNAL_SW  /* Analytical diurnal cycle of shortwave radiation
                             (use only if forcing has no diurnal cycle)  */
   /* Atmospheric input method */
#  undef  ONLINE          /* Read native atm files and interpolate online
                             to the CROCO grid (default: cubic)         */
#  ifdef ONLINE
#   undef  AROME          /* AROME high-resolution atmospheric model input */
#   undef  ERA_ECMWF      /* ERA / ECMWF reanalysis input               */
#  endif
   /* Atmospheric pressure */
#  undef  READ_PATM       /* Read surface atmospheric pressure          */
#  ifdef READ_PATM
#   define OBC_PATM       /* Apply atm pressure gradient at open boundaries */
#  endif
   /* Atmospheric Boundary Layer */
#  undef  ABL1D           /* 1D ABL model (column model coupled to ocean) */
#  ifdef  ABL1D
#   undef  ANA_ABL_LSDATA /* Analytical large-scale ABL data            */
#   undef  ANA_ABL_VGRID  /* Analytical ABL vertical grid               */
#   define STRESS_AT_RHO_POINTS /* Compute surface stress at rho points */
#   define ABL_NUDGING         /* Nudge ABL towards large-scale data     */
#   define ABL_NUDGING_DYN     /* Nudge dynamics (wind)                  */
#   define ABL_NUDGING_TRA     /* Nudge thermodynamics (temperature)     */
#   undef  ABL_DYN_RESTORE_EQ  /* Restore ABL dynamics towards equil.   */
#   undef  SFLUX_CFB           /* Coupled surface flux correction        */
#  endif
# else
   /* Without BULK_FLUX: use prescribed fluxes with corrections         */
#  define QCORRECTION     /* Linear bulk SST correction of heat flux    */
#  define SFLX_CORR       /* Freshwater flux correction towards model SSS */
#  undef  SFLX_CORR_COEF  /* Prescribed correction coefficient          */
#  define ANA_DIURNAL_SW  /* Analytical diurnal shortwave modulation    */
# endif
# undef  SFLUX_CFB        /* Surface flux coupled feedback (standalone) */
# undef  SEA_ICE_NOFLUX   /* Zero flux under sea ice                    */

/*
!======================================================================
!  LATERAL MOMENTUM ADVECTION
!  --------------------------
!  Choose exactly one horizontal advection scheme for momentum.
!  Default (none defined): 3rd-order upstream (UP3) built-in.
!======================================================================
*/
# define UV_HADV_UP3     /* 3rd-order upstream biased advection (default) */
# undef  UV_HADV_UP5     /* 5th-order upstream biased advection         */
# undef  UV_HADV_C2      /* 2nd-order centred  (use with explicit mixing) */
# undef  UV_HADV_C4      /* 4th-order centred  (use with explicit mixing) */
# undef  UV_HADV_C6      /* 6th-order centred  (use with explicit mixing) */
# undef  UV_HADV_WENO5   /* 5th-order WENO quasi-monotone advection     */

/*
!======================================================================
!  LATERAL MOMENTUM MIXING
!======================================================================
*/
# undef  UV_VIS2         /* Laplacian  horizontal viscosity             */
# undef  UV_VIS4         /* Bilaplacian horizontal viscosity            */
# ifdef UV_VIS2
#  define UV_VIS_SMAGO   /* Smagorinsky parameterization of viscosity (2D) */
# endif
# undef  UV_VIS_SMAGO_3D  /* 3D Smagorinsky parameterization            */
# undef  UV_MIX_GEO      /* Mix momentum on geopotential (z) surfaces   */
# undef  UV_MIX_S        /* Mix momentum on iso-sigma surfaces          */

/*
!======================================================================
!  VERTICAL MOMENTUM ADVECTION
!======================================================================
*/
# define UV_VADV_SPLINES /* Splines vertical advection for momentum     */
# undef  UV_VADV_WENO5   /* WENO5 vertical advection for momentum       */

/*
!======================================================================
!  VERTICAL VELOCITY ADVECTION (W)
!  --------------------------------
!  Only relevant with NBQ or wave-resolving simulations.
!  Set automatically by cppdefs_dev.h when UV_HADV_WENO5 is defined.
!======================================================================
*/
# undef  W_HADV_WENO5    /* 5th-order WENOZ horizontal advection for W  */
# undef  W_VADV_WENO5    /* 5th-order WENOZ vertical advection for W    */

/*
!======================================================================
!  LATERAL TRACER ADVECTION
!  ------------------------
!  Choose exactly one horizontal advection scheme for tracers.
!======================================================================
*/
# undef  TS_HADV_UP3     /* 3rd-order upstream biased advection         */
# define TS_HADV_RSUP3   /* Rotated split 3rd-order upstream biased (recommended) */
# undef  TS_HADV_RSUP5   /* Rotated split 5th-order upstream biased, reduced
                            dispersion/diffusion                        */
# undef  TS_HADV_UP5     /* 5th-order upstream biased advection         */
# undef  TS_HADV_C4      /* 4th-order centred advection                 */
# undef  TS_HADV_C6      /* 6th-order centred advection                 */
# undef  TS_HADV_WENO5   /* 5th-order WENOZ quasi-monotone for all tracers */
# undef  BIO_HADV_WENO5  /* 5th-order WENOZ for passive/biology/sediment
                            tracers only (set automatically – see below) */

/*
!======================================================================
!  LATERAL TRACER MIXING
!======================================================================
*/
# undef  TS_DIF2         /* Laplacian  horizontal diffusion of tracers  */
# undef  TS_DIF4         /* Bilaplacian horizontal diffusion of tracers */
# undef  TS_MIX_ISO      /* Mix along isopycnal (isoneutral) surfaces   */
# undef  TS_MIX_GEO      /* Mix along geopotential (z) surfaces         */
# undef  TS_MIX_S        /* Mix along iso-sigma surfaces                */
# undef  TS_MIX_IMP      /* Stabilizing correction for rotated diffusion
                            (use with TS_MIX_ISO or TS_MIX_GEO)       */

/*
!======================================================================
!  VERTICAL TRACER ADVECTION
!======================================================================
*/
# define TS_VADV_SPLINES /* Splines vertical advection for tracers      */
# undef  TS_VADV_WENO5   /* WENO5 vertical advection for tracers        */

/*
!======================================================================
!  SEMI-IMPLICIT VERTICAL ADVECTION
!======================================================================
*/
# undef  VADV_ADAPT_IMP  /* Adaptive semi-implicit vertical advection   */

/*
!======================================================================
!  SPONGE LAYERS
!  -------------
!  Enhance viscosity and diffusivity near open boundaries to damp
!  spurious reflections.
!======================================================================
*/
# define SPONGE          /* Sponge layers near lateral open boundaries  */

/*
!======================================================================
!  BOTTOM STRESS
!======================================================================
*/
# define LIMIT_BSTRESS   /* Limit bottom stress to prevent negative transport */
# undef  BSTRESS_FAST    /* Compute bottom stress at 3D fast (baroclinic) time steps */

/*
!======================================================================
!  VERTICAL MIXING
!  ---------------
!  Choose exactly one of LMD_MIXING or GLS_MIXING (or neither for ANA).
!======================================================================
*/
# undef  BODYFORCE        /* Apply surface and bottom stresses as body forces
                             instead of boundary conditions              */
# undef  ANA_VMIX         /* Analytical viscosity/diffusivity (no TKE scheme) */

# define LMD_MIXING       /* KPP: Large-McWilliams-Doney turbulence closure */
# ifdef LMD_MIXING
#  define LMD_SKPP        /* KPP surface boundary layer mixing           */
#  define LMD_BKPP        /* KPP bottom  boundary layer mixing           */
#  define LMD_RIMIX       /* Shear-instability interior mixing           */
#  define LMD_CONVEC      /* Convective adjustment interior mixing       */
#  define LMD_NONLOCAL    /* Nonlocal (counter-gradient) transport in SKPP */
#  undef  LMD_DDMIX       /* Double-diffusion interior mixing            */
#  undef  LMD_LANGMUIR    /* Langmuir-cell parameterization              */
# endif

# undef  GLS_MIXING       /* Generic Length Scale turbulence scheme     */
# ifdef  GLS_MIXING
#  undef  GLS_MIXING_3D   /* Full 3D GLS (vs. columnar 1D version)      */
#  ifdef  GLS_MIXING_3D
#   define GLS_KOMEGA     /* K-omega model (K-epsilon if not defined)   */
#  endif
# endif

/*
!======================================================================
!  WAVE-CURRENT INTERACTIONS
!  -------------------------
!  MRL_WCI activates the wave-current interaction module.
!  Wave forcing can come from coupling (OW_COUPLING), an offline file,
!  the WKB internal wave model, or analytical values.
!======================================================================
*/
# undef  MRL_WCI          /* Wave-Current Interaction model (Mellor-type) */
# ifdef OW_COUPLING
#  define MRL_WCI
#  define BBL
# endif
# ifdef MRL_WCI
#  undef  MRL_CEW          /* Wave-driven longshore current energetics   */
#  ifndef OW_COUPLING
#   undef  WAVE_OFFLINE   /* Read pre-computed wave fields from a file  */
#   define ANA_WWAVE      /* Analytical (constant) wave forcing (Hs,Tp,Dir) */
#   undef  WKB_WWAVE      /* Internal WKB wave model                    */
#  endif
#  undef  WAVE_ROLLER     /* Wave roller model                          */
#  define WAVE_STREAMING  /* Bottom wave streaming (wave-induced current) */
#  define WAVE_FRICTION   /* Bottom friction in WKB model               */
#  define WAVE_RAMP       /* Ramp wave forcing using wave_ramp parameter */
#  ifdef WKB_WWAVE
    /* WKB open boundary settings (choose offshore boundary)            */
#   undef  WKB_OBC_NORTH  /* WKB open boundary – north side             */
#   undef  WKB_OBC_SOUTH  /* WKB open boundary – south side             */
#   define WKB_OBC_WEST   /* WKB open boundary – west side              */
#   undef  WKB_OBC_EAST   /* WKB open boundary – east side              */
#  endif
   /* Wave breaking parameterization (for WKB)                         */
#  undef  WAVE_BREAK_CT93 /* Thornton & Guza (1983) wave breaking       */
#  undef  WAVE_BREAK_TG86 /* Church & Thornton (1993) wave breaking     */
#  undef  WAVE_BREAK_TG86A /* Alternate Church-Thornton formulation     */
#  undef  ANA_BRY_WKB     /* Analytical wave boundary forcing (croco.in) */
# endif

/*
!======================================================================
!  WAVE MAKER
!  ----------
!  WAVE_MAKER: wave-resolving internal wave generation (ideal cases,
!  e.g. SANDBAR, SWASH, RIP).
!======================================================================
*/
# undef  WAVE_MAKER          /* Wave-maker boundary source               */
# ifdef  WAVE_MAKER
#  undef  WAVE_MAKER_SPECTRUM /* Multi-frequency spectrum input          */
#  undef  WAVE_MAKER_DSPREAD  /* Directional spreading                  */
# endif
# undef  WAVE_MAKER_INTERNAL  /* Internal (immersed) wave maker         */

/*
!======================================================================
!  BOTTOM FORCING
!======================================================================
*/
# define ANA_BSFLUX       /* Analytical bottom salinity flux (typically 0) */
# define ANA_BTFLUX       /* Analytical bottom temperature flux (typically 0) */

/*
!======================================================================
!  ANALYTICAL FORCING
!  ------------------
!  Use these keys when running test cases without external NetCDF files.
!  Each ANA_ key replaces a file read with an inline Fortran function.
!======================================================================
*/
# undef  ANA_GRID         /* Analytical grid (no grid NetCDF file)      */
# undef  ANA_INITIAL      /* Analytical initial conditions               */
# undef  ANA_BRY          /* Analytical open boundary conditions         */
# undef  ANA_SSH          /* Analytical sea level (replaces ZCLIMATOLOGY) */
# undef  ANA_M2CLIMA      /* Analytical barotropic velocity climatology  */
# undef  ANA_M3CLIMA      /* Analytical baroclinic velocity climatology  */
# undef  ANA_TCLIMA       /* Analytical tracer climatology               */
# undef  ANA_SMFLUX       /* Analytical surface momentum (wind) flux     */
# undef  ANA_SRFLUX       /* Analytical shortwave radiation flux         */
# undef  ANA_STFLUX       /* Analytical surface heat flux                */
# undef  ANA_SSFLUX       /* Analytical surface salinity flux            */
# undef  ANA_SST          /* Analytical SST (for QCORRECTION)           */
# undef  ANA_TIDES        /* Analytical tidal forcing (RIP / TIDAL_FLAT) */
# undef  ANA_JET          /* Analytical jet initial condition (JET case) */
# undef  ANA_MORPHODYN    /* Analytical morphodynamic forcing           */
# undef  ANA_PSOURCE      /* Analytical vertical profiles for point sources */

/*
!======================================================================
!  POINT SOURCES – RIVERS
!======================================================================
*/
# undef  PSOURCE          /* River / point-source mass and tracer input  */
# undef  PSOURCE_NCFILE   /* Read time-varying river transports from NetCDF */
# ifdef PSOURCE_NCFILE
#  undef  PSOURCE_NCFILE_TS /* Read time-varying river concentrations from NetCDF */
# endif
# undef  PSOURCE_MASS      /* Inject river as vertical volume flux at rho
                              points (instead of flux through u,v faces) */

/*
!======================================================================
!  OPEN BOUNDARY CONDITIONS – SCHEME CHOICES
!  ------------------------------------------
!  Choose one scheme per variable class (M2, M3, T).
!======================================================================
*/
# define OBC_M2CHARACT    /* Characteristic methods for barotropic velocities */
# undef  OBC_M2ORLANSKI   /* Radiative OBC for barotropic velocities     */
# undef  OBC_M2SPECIFIED  /* Specified (Dirichlet) OBC for barotropic velocities */
# undef  OBC_VOLCONS      /* Enforce volume conservation at open boundaries
                             (use with OBC_M2ORLANSKI)                  */
# undef  OBC_SPECIFIED_WEST /* Specified (Dirichlet) OBC on western boundary */
# undef  OBC_REDUCED_PHYSICS /* Simplified physics at OBCs (for tidal flat etc.) */

# define OBC_M3ORLANSKI   /* Radiative OBC for baroclinic velocities     */
# undef  OBC_M3SPECIFIED  /* Specified OBC for baroclinic velocities     */

# define OBC_TORLANSKI    /* Radiative OBC for tracers                  */
# undef  OBC_TSPECIFIED   /* Specified OBC for tracers                  */
# undef  OBC_TUPWIND      /* Upwind OBC for tracers                     */

/*
!======================================================================
!  INPUT / OUTPUT
!======================================================================
*/
# define AVERAGES         /* Time-averaged output fields                 */
# define AVERAGES_K       /* Time-averaged vertical mixing fields        */
# undef  OUTPUTS_SURFACE  /* Surface-only output snapshot fields         */
# undef  HOURLY_VELOCITIES /* Additional hourly velocity output          */

/*
!======================================================================
!  RESTART
!======================================================================
*/
# undef  EXACT_RESTART    /* Write extra fields for exact bit-reproducible restart */

/*
!======================================================================
!  DIAGNOSTICS
!  -----------
!  Budget diagnostics for tracers, momentum, vorticity, energy, etc.
!======================================================================
*/
# undef  DO_NOT_OVERWRITE  /* Keep existing diagnostic files (append mode) */
# undef  RESTART_DIAGS     /* Restart diagnostic accumulation from file  */

/* Tracer budget */
# undef  DIAGNOSTICS_TS    /* Tracer equation budget terms               */
# ifdef DIAGNOSTICS_TS
#  undef  DIAGNOSTICS_TS_ADV  /* Use advective (flux) formulation        */
#  undef  DIAGNOSTICS_TS_MLD  /* Integrate budgets over the mixed layer  */
#  ifdef DIAGNOSTICS_TS_MLD
#   undef  DIAGNOSTICS_TS_MLD_DENS /* Use density criterion for MLD     */
#  endif
# endif

# undef  DIAGNOSTICS_TSVAR  /* Tracer variance budget                   */
# ifdef DIAGNOSTICS_TSVAR
#  define DIAGNOSTICS_TS
#  define DIAGNOSTICS_TS_ADV
# endif

/* Momentum budget */
# undef  DIAGNOSTICS_UV    /* Momentum equation budget terms             */
# undef  DIAGNOSTICS_VRT   /* Depth-mean vorticity budget                */
# undef  DIAGNOSTICS_KE    /* Kinetic energy budget                      */
# undef  DIAGNOSTICS_BARO  /* Barotropic kinetic energy                  */

/* Potential vorticity / dissipation */
# undef  DIAGNOSTICS_PV    /* Potential vorticity budget                 */
# undef  DIAGNOSTICS_DISS  /* Dissipation (activates DIAGNOSTICS_PV)    */
# ifdef DIAGNOSTICS_DISS
#  define DIAGNOSTICS_PV
# endif

/* Eddy terms */
# undef  DIAGNOSTICS_EDDY  /* Reynolds-stress eddy terms                 */

/*
!======================================================================
!  APPLICATIONS – DETAILS
!======================================================================
*/

/* Quasi-monotone advection for passive/biology/sediment tracers       */
/* (set automatically when PASSIVE_TRACER, BIOLOGY, SEDIMENT or MUSTANG
    is defined – do not edit unless you know what you are doing)       */
# if defined PASSIVE_TRACER || defined BIOLOGY || defined SEDIMENT \
                                               || defined MUSTANG
#  define BIO_HADV_WENO5  /* 5th-order WENOZ quasi-monotone for passive tracers */
# endif

/*
!----------------------------------------------------------------------
!  BIOGEOCHEMICAL MODELS
!  ---------------------
!  Exactly one model must be chosen when BIOLOGY is defined.
!----------------------------------------------------------------------
*/
# ifdef BIOLOGY
#  define PISCES          /* 24-component PISCES biogeochemical model    */
#  undef  BIO_NChlPZD     /* 5-component NPZD model                     */
#  undef  BIO_N2ChlPZD2   /* 7-component NPZD model                     */
#  undef  BIO_BioEBUS     /* 12-component NPZD model (EBUS-oriented)    */

   /* PISCES options */
#  ifdef PISCES
#   undef  DIURNAL_INPUT_SRFLX /* Diurnal cycle modulation of shortwave input */
#   define key_pisces           /* Activates PISCES code               */
#   define key_ligand           /* Include ligand cycling              */
#   undef  key_pisces_quota     /* Quota-based nutrient model          */
#   undef  key_pisces_npzd      /* Simplified NPZD-PISCES variant      */
#   undef  key_sediment         /* PISCES sediment module              */
#  endif
#  ifdef BIO_NChlPZD
#   define OXYGEN                /* Oxygen tracer in BIO_NChlPZD       */
#  endif
#  ifdef BIO_BioEBUS
#   define NITROUS_OXIDE         /* N2O tracer in BIO_BioEBUS          */
#  endif

   /* Biology diagnostics */
#  define DIAGNOSTICS_BIO         /* Store biological flux diagnostics */
#  if defined DIAGNOSTICS_BIO && defined PISCES
#   define key_trc_diaadd         /* Additional PISCES tracer diagnostics */
#  endif
# endif

/*
!----------------------------------------------------------------------
!  STATIONS OUTPUT
!----------------------------------------------------------------------
*/
# ifdef STATIONS
#  define ALL_SIGMA       /* Output all sigma levels at stations         */
# endif

/*
!----------------------------------------------------------------------
!  USGS SEDIMENT MODEL (SEDIMENT)
!----------------------------------------------------------------------
*/
# ifdef SEDIMENT
   /* Transport components */
#  define SUSPLOAD        /* Suspended load transport                    */
#  define BEDLOAD         /* Bedload transport                           */
#  define MORPHODYN       /* Morphodynamics (bed evolution)              */
   /* Analytical / parameter inputs */
#  undef  ANA_SEDIMENT    /* Analytical sediment ripple and bed parameters */
#  undef  ANA_BPFLUX      /* Analytical kinematic bottom flux of sediment */
   /* Bedload formulation (mutually exclusive)                          */
#  undef  BEDLOAD_MPM     /* Bedload: Meyer-Peter-Muller formulation     */
#  undef  BEDLOAD_WULIN   /* Bedload: Wu & Lin formulation               */
#  undef  BEDLOAD_MARIEU  /* Bedload: Marieu formulation                 */
#  undef  BEDLOAD_WENO5   /* Bedload: WENO 5th-order advection           */
   /* Avalanching */
#  undef  SLOPE_NEMETH    /* Avalanching: Nemeth et al. (2006)           */
#  undef  SLOPE_LESSER    /* Avalanching: Lesser et al. (2004)           */
   /* Bed composition */
#  undef  COHESIVE_BED    /* Purely cohesive (mud) bed                   */
#  undef  MIXED_BED       /* Mixed cohesive/non-cohesive bed             */
   /* Flocculation (requires COHESIVE_BED or MIXED_BED)                */
#  undef  SED_FLOCS       /* Flocculation of cohesive particles         */
#  undef  SED_DEFLOC      /* De-flocculation process                    */
#  undef  FLOC_TURB_DISS  /* Turbulent dissipation of flocs             */
#  undef  FLOC_BBL_DISS   /* BBL dissipation of flocs                   */
   /* Bed friction */
#  undef  SED_TAU_CD_CONST /* Constant drag coefficient for bed stress  */
   /* Critical shear stress */
#  undef  TAU_CRIT_WULIN  /* Wu & Lin critical shear stress formulation  */
# endif

/*
!----------------------------------------------------------------------
!  MUSTANG SEDIMENT MODEL
!----------------------------------------------------------------------
*/
# ifdef MUSTANG
#  undef  key_MUSTANG_V2       /* MUSTANG version 2 (advanced features) */
#  undef  key_MUSTANG_bedload  /* MUSTANG bedload transport             */
#  undef  key_MUSTANG_flocmod  /* MUSTANG flocculation model            */
#  undef  key_tauskin_c_upwind /* Upwind skin friction for MUSTANG      */
#  undef  MORPHODYN             /* Morphodynamics (usually off)         */
#  undef  WAVE_OFFLINE          /* Wave offline forcing with MUSTANG    */
# endif

/*
!======================================================================
!  BOTTOM BOUNDARY LAYER (BBL) OPTIONS
!  ------------------------------------
!  Sub-options when BBL is defined (wave-dependent roughness).
!======================================================================
*/
# ifdef BBL
#  undef  BBL_BREAKING_STIR /* Wave breaking stirring in BBL            */
#  undef  ANA_BSEDIM      /* Analytical bed parameters (when SEDIMENT not defined) */
#  undef  ANA_WWAVE       /* Analytical wave forcing for BBL            */
#  undef  Z0_BL           /* Compute bedload roughness for ripple predictor */
#  undef  Z0_RIP           /* Bedform (ripple) roughness for sandy beds */
#  undef  Z0_BIO           /* Biogenic bedform roughness for silty beds */
# endif

/*
!======================================================================
!  INCLUDE ADVANCED OPTIONS AND GLOBAL DEFINITIONS
!  (Do not edit the two lines below)
!======================================================================
*/
#include "cppdefs_dev.h"
#include "set_global_definitions.h"
