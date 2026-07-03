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
! This is include file "scalars.h"
!---------------------------------
!
!  The following common block contains time variables and indices
! for 2D (k-indices) and 3D (n-indices) computational engines. Since
! they are changed together, they are placed into the same cache line
! despite their mixed type, so that only one cachene is being
! invalidated and has to be propagated accross the cluster.
!
! Note that the real values are placed first into the common block
! before the integer variables. This is done to prevent the
! possibility of misallignment of the 8-byte objects in the case
! when an uneven number of 4-byte integers is placed before a 8-byte
! real (in the case when default real size is set to 8bytes).
! Thought misallignment is not formally a violation of fortran
! standard, it may cause performance degradation and/or make compiler
! issue a warning message (Sun, DEC Alpha) or even crash (Alpha).
!
! time        Time since initialization [seconds];
! time_start  Initialization time [seconds];
! tdays       Time since initialization [days];
! dtfast      Time step for 2D (barotropic) mode [seconds];
!
      real dtfast, time, time2, time_start, tdays, start_time
#ifdef USE_CALENDAR
      real time_end
      character*19 date
#endif
      integer iic, kstp, krhs, knew, next_kstp
#ifdef SOLVE3D
     &      , iif, nstp, nrhs, nnew, nbstep3d
#endif
# if defined OA_COUPLING || defined OW_COUPLING
#  ifdef AGRIF
     &     , it_inside_root
#  endif
# endif

#ifdef WKB_WWAVE
     &      , wstp, wnew
#endif
      logical PREDICTOR_2D_STEP
      common /time_indices/  dtfast, time, time2,time_start, tdays,
     &     start_time,
     &     iic, kstp, krhs, knew, next_kstp,
#ifdef SOLVE3D
     &                       iif, nstp, nrhs, nnew, nbstep3d,
#endif
# if defined OA_COUPLING || defined OW_COUPLING
#  ifdef AGRIF
     &      it_inside_root,
#  endif
# endif

#ifdef WKB_WWAVE
     &                       wstp, wnew,
#endif
     &                       PREDICTOR_2D_STEP
#ifdef USE_CALENDAR
      common /time_indices2/ time_end,
     &                       date
#endif

!
! Slowly changing variables: these are typically set in the beginning
! of the run and either remain unchanged, or are changing only in
! association with the I/0.
!
! xl, el   Physical size (m) of domain box in the XI-,ETA-directions.
!
!
! sc_r     S-coordinate independent variable, [-1 < sc < 0] at
!             vertical RHO-points
! sc_w     S-coordinate independent variable, [-1 < sc < 0] at
!             vertical W-points.
! Cs_r     Set of S-curves used to stretch the vertical coordinate
!             lines that follow the topography at vertical RHO-points.
! Cs_w     Set of S-curves used to stretch the vertical coordinate
!             lines that follow the topography at vertical W-points.
!
! rho0     Boussinesque Approximation Mean density [kg/m^3].
! R0       Background constant density anomaly [kg/m^3] used in
!                                      linear equation of state.
! T0,S0    Background temperature (Celsius) and salinity [PSU]
!                          values used in analytical fields;
! Tcoef    Thermal expansion coefficient in linear EOS;
! Scoef    Saline contraction coefficient in linear EOS;
!
! rdrg     Linear bottom drag coefficient.
! rdrg2    Quadratic bottom drag coefficient.
! Cdb_max  Maximum bottom drag coefficient allowed.
! Cdb_min  Minimum bottom drag coefficient to avoid the
!                law-of-the-wall to extend indefinitely.
! Zobt      Bottom roughness (m).
!
! gamma2   Slipperiness parameter, either 1. (free-slip)
!
! mld_crit_D     Density criterion to define the ML [kg/m^3] 
! mld_crit_T     Temperature criterion to define the ML [Celsius]
!
! ntstart  Starting timestep in evolving the 3D primitive equations;
!                              usually 1, if not a restart run.
! nsta     Number of timesteps between storage of station data.
! navg     Number of timesteps between storage of time-averaged
!                                                           fields.
! ntsavg   Starting timestep for accumulation of output time-
!                                                 averaged fields.
!
! levsfrc  Deepest level to apply surface momentum stress as
!                                                 bodyforce.
! levbfrc  Shallowest level to apply bottom momentum stress as
!                                                 bodyforce.
! got_tini Logical switch used at initialisation
!              If TRUE, the tracer is present in the initial file
!              If FALSE, the tracer needs an analytical value
!
! got_inised Logical switch used at initialisation  of sediments
!              If TRUE, the sediment var. is in the initial file
!              If FALSE, the sed. var. gets analytical value from file
!
! got_inibed Logical switch used at initialisation of ripple height, length
!              If TRUE, the ripple var. is in the initial file
!              If FALSE, the ripple var. is obtained from file (ifdef also SEDIMENT)
!                        the ripple var. is set in ana_bsedim (ifndef SEDIMENT)
!
      real time_avg, time2_avg
     &               , xl, el
#ifdef SOLVE3D
# ifndef M3FAST_SEDLAYERS
      real  sc_w(0:N), Cs_w(0:N), sc_r(N), Cs_r(N)
# else
      real  sc_w(-N_sl  :N), Cs_w(-N_sl  :N) 
     &    , sc_r(-N_sl+1:N), Cs_r(-N_sl+1:N)
# endif
      real  rx0, rx1
      real weight(6,0:NWEIGHT)

#endif
      integer numthreads,     ntstart
     &      , nfast
#ifdef EXACT_RESTART
     &     ,  forw_start
#endif

#if defined SOLVE3D && defined TRACERS
      logical got_tini(NT)
#endif
#ifdef SEDIMENT
      logical got_inised(3)
#endif
#ifdef BBL
      logical got_inibed(2)
#endif

#if defined DIAGNOSTICS_EDDY && ! defined XIOS
      integer nwrtdiags_eddy
      logical ldefdiags_eddy
#endif


#ifdef  BAND_DEBUG         
      character(len=50) :: chkbandname
      character(len=50) :: fileline
#endif
      common /scalars_main/
     &             time_avg, time2_avg
     &           , xl, el
#ifdef SOLVE3D
     &           , sc_w,      Cs_w,      sc_r,    Cs_r
     &           , rx0,       rx1
     &                      , weight
#endif
     &      , numthreads,     ntstart
     &      , nfast
#ifdef EXACT_RESTART
     &       , forw_start
#endif
#if defined SOLVE3D && defined TRACERS
     &                      , got_tini
#endif
#ifdef SEDIMENT
     &                      , got_inised
#endif
#ifdef BBL
     &                      , got_inibed
#endif
#if defined DIAGNOSTICS_EDDY && ! defined XIOS
     &                      , ldefdiags_eddy, nwrtdiags_eddy
#endif  
#ifdef  BAND_DEBUG         
       common /scalchkbandname/ chkbandname
#endif       
!
!-----------------------------------------------------------------------
! This following common block contains a set of globally accessable
! variables in order to allow information exchange between parallel
! threads working on different subdomains.
!
! Global summation variables are declared with 16 byte precision
! to avoid accumulation of roundoff errors, since roundoff error
! depends on the order of summation, which is undeterministic in
! the case of summation between the parallel threads; not doing so
! would make it impossible to pass an ETALON CHECK test if there is
! a feedback of these sums into the dynamics of the model, such as
! in the case when global mass conservation is enforced.
!
!  One sunny spring day, sometime in 1989 an american tourist, who
! happened to be an attorney, was walking along a Moscow street.
! Because it was the period of 'Perestroika' (which literally means
! 'remodelling'), so that a lot of construction was going on in
! Moscow, dozens of holes and trenches were open on the street. He
! felt into one of them, broke his leg, ended up in a hospital and
! complaining: In my country if a construction firm would not place
! little red flags around the construction zone to warn passers-by
! about the danger, I will sue em for their negligence! The doctor,
! who was performing surgery on his leg replied to him: Did not you
! see the one big red flag above the whole country in the first place?
!
! WARNING: FRAGILE ALIGNMENT SEQUENCE: In the following common block:
! since real objects are grouped in pairs and integer*4 are grouped
! in quartets, it is guaranteed that 16 Byte objects are aligned
! in 16 Byte boundaries and 8 Byte objects are aligned in 8 Byte
! boundaries. Removing or introduction of variables with violation
! of parity, as well as changing the sequence of variables in the
! common block may cause violation of alignment.
!-----------------------------------------------------------------------
!
      logical synchro_flag
      common /sync_flag/ synchro_flag

      integer may_day_flag  ! This is a shared variable among nested grids
      integer tile_count, first_time, bc_count
#ifdef BIOLOGY
     &      , bio_count
#endif
      common /communicators_i/
     &        may_day_flag, tile_count, first_time, bc_count
#ifdef BIOLOGY
     &      , bio_count
#endif

      real hmin, hmax, grdmin, grdmax, Cu_min, Cu_max
      common /communicators_r/
     &     hmin, hmax, grdmin, grdmax, Cu_min, Cu_max

#ifdef SPHERICAL
      real lonmin, lonmax, latmin, latmax
      common /communicators_lonlat/
     &     lonmin, lonmax, latmin, latmax
#endif

      real*QUAD Cu_Adv3d,  Cu_W, Cu_Nbq_X, Cu_Nbq_Y, Cu_Nbq_Z
      integer i_cx_max, j_cx_max, k_cx_max
      common /diag_vars/ Cu_Adv3d,  Cu_W,
     &        i_cx_max, j_cx_max, k_cx_max
      real*QUAD volume, avgke, avgpe, avgkp, bc_crss

#ifdef OBC_VOLCONS
     &        , bc_flux, ubar_xs
#endif
#if defined BIOLOGY && defined TRACERS
     &        , global_sum(0:2*NT+1)
#endif
#ifdef RESET_RHO0
     &        , avg_vol, avg_rho
#endif

      common /communicators_rq/
     &          volume, avgke, avgpe, avgkp, bc_crss
#ifdef OBC_VOLCONS
     &        , bc_flux,  ubar_xs
#endif
#ifdef BIOLOGY
     &        , global_sum
#endif
#ifdef RESET_RHO0
     &        , avg_vol, avg_rho
#endif

#ifdef MPI
!
! MPI related variables
! === ====== =========
!
      logical EAST_INTER2, WEST_INTER2, NORTH_INTER2, SOUTH_INTER2
      logical EAST_INTER, WEST_INTER, NORTH_INTER, SOUTH_INTER
      logical CORNER_SW,CORNER_NW,CORNER_NE,CORNER_SE
      integer mynode, mynode2, ii,jj, p_W,p_E,p_S,p_N, p_SW,p_SE,
     & p_NW,p_NE,NNODES2
      common /comm_setup/ mynode, mynode2, ii,jj, p_W,p_E,p_S,p_N,
     & p_SW,p_SE, p_NW,p_NE, EAST_INTER, WEST_INTER, NORTH_INTER,
     & SOUTH_INTER, EAST_INTER2, WEST_INTER2, NORTH_INTER2, SOUTH_INTER2,
     & CORNER_SW,CORNER_NW,CORNER_NE,CORNER_SE,NNODES2

#endif

!
! Physical constants:
! ======== ==========

      real pi, deg2rad, rad2deg
      parameter (pi=3.14159265358979323846, deg2rad=pi/180.,
     &                                      rad2deg=180./pi)
!
! Earth radius [m]; Earth rotation [rad/s]; Acceleration of gravity [m/s^2];
! duration of the day in seconds and its inverse; Julian offset day.

      real Eradius, Erotation, g, day2sec,sec2day, jul_off,
     &     year2day,day2year
      parameter (Eradius=6371315.0,  Erotation=7.292115090e-5,
     &           day2sec=86400., sec2day=1./86400.,
     &           year2day=365.25, day2year=1./365.25,
     &           jul_off=2440000.)
!
! Acceleration of gravity (nondimensional for Soliton problem)
!
#ifdef SOLITON
      parameter (g=1.)
#else
      parameter (g=9.81)
#endif
!
!  Specific heat [Joules/kg/degC] for seawater, it is approximately
!  4000, and varies only slightly (see Gill, 1982, Appendix 3).
!
      real Cp
      parameter (Cp=3985.0)

      real vonKar
      parameter (vonKar=0.41)
!
!   FillValue (Needed if the FILLVAL key is defined)
!   (See fillvalue.F subroutine)
      real spval
      parameter (spval=-999.0)
      logical mask_val
      parameter (mask_val = .true.)
