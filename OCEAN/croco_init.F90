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
#include "cppdefs.h"
!
!=======================================================================
!  MODULE croco_init
!
!  Purpose : Model initialization phase. croco_initialize() carries out
!            every step from "read the namelist" to "write the initial
!            history file" before handing control back to the 
!            time-stepping driver (croco.F90).

!=======================================================================
!
      module croco_init

      implicit none
      private
      public :: croco_initialize

      contains

      subroutine croco_initialize(ierr)
!
      use param
      use scalars
      use ncscrum
#if defined CVTK_DEBUG || defined CVTK_DEBUG_ADVANCED || \
    defined CVTK_DEBUG_PERFRST
      use debug
#endif
#ifdef XIOS
      use xios           ! XIOS module
#endif
#ifdef ENSEMBLE
      use ensmpi         ! Ensemble module
#endif
#ifdef OPENACC
      use openacc
#endif
      use buffer, only : init_buffer
      use croco_namelist_read, only : read_nml, read_nml_fname
      use croco_namelist, only : dt, ndtfast, ntimes, ninfo
      use croco_namelist, only : ldefhis, nrrec
#ifdef STATIONS
      use croco_namelist, only : ldefsta
#endif
#if defined ABL1D && !defined XIOS
      use croco_namelist, only : ldefablhis
#endif
#ifdef USE_CALENDAR
      use croco_namelist, only : end_date
#endif
#if defined OA_COUPLING || defined OW_COUPLING
      use mod_prism
#endif
#ifdef PISCES
      use pisces_ini
      use trcnam_pisces
#endif
#ifdef STOGEN
      use stomod, only : sto_mod_init
      use stoexternal, only : ocean_2_stogen
#endif
#ifdef SUBSTANCE
      use substance, only : substance_read_alloc
      use substance, only : substance_surfcell
#endif
#ifdef MUSTANG
      use plug_MUSTANG_CROCO, only : mustang_init_main
#endif
#ifdef OBSTRUCTION
      use plug_OBSTRUCTIONS, only : obst_init_main
      use plug_OBSTRUCTIONS, only : obst_update_main
#endif
#ifdef ONLINE_ANALYSIS
      use module_interface_oa, only : init_parameter_oa
      use module_interface_oa, only : if_oa, if_oa_fast_mode
#endif
!
      integer, intent(out) :: ierr
      integer :: tile
#ifdef WKB_WWAVE
      integer :: winterp
#endif
#ifdef USE_CALENDAR
      real(kind=8) :: tool_datosec
#endif
#ifdef MPI
      include 'mpif.h'
#endif
#ifdef ENSEMBLE
      integer :: mpi_comm_all
#endif
!
#include "private_scratch.h"
#include "nbq.h"
! provides ocean_grid_comm/oasis_time (the "world" communicator macro
! used below is redefined to it by cppdefs_dev.h under coupling)
#include "mpi_cpl.h"
#include "grid.h"
#include "ocean2d.h"
#ifdef WKB_WWAVE
#include "wkb_wwave.h"
#endif
#ifdef STATIONS
#include "sta.h"
#include "nc_sta.h"
#endif
!
      ierr = 0
!
#ifdef JEANZAY
! Initialise OpenACC...
# ifdef OPENACC
        call initialisation_openacc
# endif
! ... before initialise MPI
#endif
!
!----------------------------------------------------------------------
!  Initialize communicators and subgrids decomposition:
!  MPI parallelization, XIOS server, OASIS coupling, AGRIF nesting,
!----------------------------------------------------------------------
!
#ifdef MPI
# if (!defined AGRIF && !defined OA_COUPLING && !defined OW_COUPLING)
      call MPI_Init (ierr)
# endif
!
!  XIOS, OASIS and AGRIF: split MPI communicator
!  (XIOS with OASIS is not done yet)
!
# if (defined XIOS && !defined AGRIF)
      call xios_initialize( "crocox",return_comm=MPI_COMM_WORLD )
# elif (defined OA_COUPLING && !defined AGRIF)
      call cpl_prism_init  ! If AGRIF --> call cpl_prism_init in zoom.F
# elif (defined OW_COUPLING && !defined AGRIF)
      call cpl_prism_init  ! In AGRIF case, cpl_prism_init is in zoom.F
# elif defined AGRIF
      call Agrif_MPI_Init(MPI_COMM_WORLD)
# endif
!
!  ENSEMBLE: further split CROCO communicator
!
# if defined ENSEMBLE
#  if defined XIOS
      mpi_comm_all = MPI_COMM_WORLD
#  else
      mpi_comm_all = mpi_comm_world   ! true world communicator (keep lowercase to avoid cpp replacement!)
#  endif
      call ens_comm_set ( )
# endif
#endif /* MPI */
!
!  Initialize AGRIF nesting (Agrif_Init_Grids() itself stays in
!  croco.F90: it is the one call conv is only known to handle
!  correctly when it is made directly from the PROGRAM unit)
!
#ifdef AGRIF
      call declare_zoom_variables()
#endif
!
!  Setup MPI domain decomposition
!
#ifdef MPI
      call MPI_Setup (ierr)
      if (ierr /= 0) return
#endif
!
!  Initialize debug procedure
!
#if defined CVTK_DEBUG || defined CVTK_DEBUG_ADVANCED || \
    defined CVTK_DEBUG_PERFRST
      call debug_ini
#endif
!
      call init_buffer
!
!----------------------------------------------------------------------
!  Read in tunable model parameters from croco namelist file
!----------------------------------------------------------------------
!
      call read_nml_fname ()
      call read_nml (ierr)
      if (ierr /= 0) return
      call read_inp (ierr)
      if (ierr /= 0) return
!
!----------------------------------------------------------------------
!  Initialize global model parameters
!----------------------------------------------------------------------
!
!  Global scalar variables
!
      call init_scalars (ierr)
      if (ierr /= 0) return
!
#ifdef SOLVE3D
!
!  PISCES biogeochemeical model parameters
!
# if defined BIOLOGY && defined PISCES
      call trc_nam_pisces
# endif
!
!  Read sediment initial values and parameters from sediment.in file
!
# ifdef SEDIMENT
#  ifdef AGRIF
      if (Agrif_lev_sedim == 0) call init_sediment
#  else
      call init_sediment
#  endif
# endif
#endif

#ifdef SUBSTANCE
!
!  Substance var need for MUSTANG and BIOLink
!
      call substance_read_alloc()
#endif
!
! Online spectral analysis module
!
#ifdef ONLINE_ANALYSIS
# ifdef MPI
      call init_parameter_oa( io_unit_oa=stdout,                       &
                               if_print_node_oa=(mynode==0),           &
                               mynode_oa=mynode,                       &
                               comm_oa=MPI_COMM_WORLD,                 &
                               dti_oa=dt, kount0_oa=ntstart-1,         &
                               nt_max_oa=ntimes, dtf_oa=dtfast,        &
                               ntf_max_oa=ndtfast,                     &
                               ntiles=NSUB_X*NSUB_E)
# else
      call init_parameter_oa( io_unit_oa=stdout,                       &
                               if_print_node_oa=.true.,                &
                               mynode_oa=0, comm_oa=0,                 &
                               dti_oa=dt, kount0_oa=ntstart-1,         &
                               nt_max_oa=ntimes, dtf_oa=dtfast,        &
                               ntf_max_oa=ndtfast,                     &
                               ntiles=NSUB_X*NSUB_E)
# endif
#endif
!
!----------------------------------------------------------------------
!  Create parallel threads;
!  initialize (FIRST-TOUCH) model global arrays (most of them
!  are just set to to zero).
!----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call init_arrays (tile)
      enddo
!
! Copy to device(s)
!
#if defined OPENACC
!$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call copy_to_devices (tile)
      enddo
#endif
!
!----------------------------------------------------------------------
!  Set horizontal grid, model bathymetry and Land/Sea mask
!----------------------------------------------------------------------
!
#ifdef ANA_GRID
!
!  Set grid analytically
!
!$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call ana_grid (tile)
      enddo
# if defined CVTK_DEBUG || defined CVTK_DEBUG_ADVANCED
!$OMP BARRIER
!$OMP MASTER
       call check_tab2d(h(:,:),'h initialisation #1','r',              &
            ondevice=.TRUE.)
!$OMP END MASTER
# endif
#else
!
!  Read grid from GRID NetCDF file
!
      call get_grid
      if (may_day_flag /= 0) return
#endif
!
!  Compute various metric term combinations.
!
!$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call setup_grid1 (tile)
      enddo
!
!$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call setup_grid2 (tile)
      enddo
!
!----------------------------------------------------------------------
!  Setup vertical grid variables setup vertical S-coordinates
!  and fast-time averaging for coupling of
!  split-explicit baroropic mode.
!----------------------------------------------------------------------
!
#ifdef SOLVE3D
!
!  Set vertical S-coordinate functions
!
      call set_scoord
!
!  Set fast-time averaging for coupling of split-explicit baroropic mode.
!
      call set_weights
!
!  Create three-dimensional S-coordinate system,
!  which may be needed by ana_initial
!  (here it is assumed that free surface zeta=0).
!
!$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call set_depth (tile)
      enddo
!
!  Make grid diagnostics
!
!$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call grid_stiffness (tile)
      enddo
#endif

#if defined OPENACC
!$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call copy_from_devices (tile)
      enddo
#endif

#if defined OPENACC
!$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call copy_to_devices_2 (tile)
      enddo
#endif
!
!----------------------------------------------------------------------
!  Set initial conditions for momentum and tracer variables
!----------------------------------------------------------------------
!
!  Read from NetCDF file
!
#ifdef ANA_INITIAL
      if (nrrec /= 0) then            ! read from NetCDF file
#endif
#ifdef EXACT_RESTART
        call get_initial (nrrec-1, 2) ! Set initial conditions
                                       ! in case of restart
!$OMP BARRIER
# ifdef SOLVE3D
        do tile=0,NSUB_X*NSUB_E-1
          call set_depth (tile)       !<-- needed to initialize Hz_bak
        enddo
!$OMP BARRIER
# endif
#endif
        call get_initial (nrrec, 1)   ! Set initial conditions
#ifdef ANA_INITIAL
      else  ! nrrec.eq.0
# if defined OA_COUPLING || defined OW_COUPLING
        call cpl_prism_define
        oasis_time = 0
        MPI_master_only write(*,*)'CPL-CROCO: OASIS_TIME',oasis_time
# endif
      endif
#endif
                                ! Set initial model clock: at this
      time=start_time          ! moment "start_time" (global scalar)
      tdays=time*sec2day       ! is set by get_initial or analytically
                                ! --> copy it into threadprivate "time"

#ifdef USE_CALENDAR
      time_end=tool_datosec(end_date)
      ntimes=int((time_end-time)/dt)
      MPI_master_only write(stdout,*)                                  &
           'Ntimes from date_start and date_end:',ntimes
      ntimes=ntimes+ntstart
#endif
!
!  Set initial conditions analytically for ideal cases
!  or for tracer variables not present in NetCDF file
!
!$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call ana_initial (tile)
      enddo
!
!----------------------------------------------------------------------
!  Initialize specific PISCES variables
!----------------------------------------------------------------------
!
#if defined BIOLOGY && defined PISCES
!$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call  pisces_ini_tile (tile)
      enddo
#endif
      if (may_day_flag /= 0) return
!
!----------------------------------------------------------------------
!  Bottom sediment parameters for BBL or SEDIMENT model
!----------------------------------------------------------------------
!
#if (defined BBL && defined ANA_BSEDIM) || defined SEDIMENT
!
!  --- Set analytically ---
!
!$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1

# if defined BBL && defined ANA_BSEDIM
#  ifdef AGRIF
        if (Agrif_lev_sedim == 0) call ana_bsedim (tile)
#  else
        call ana_bsedim (tile)
#  endif
# endif
!
# ifdef SEDIMENT
#  ifdef AGRIF
        if (Agrif_lev_sedim == 0) call ana_sediment (tile)
#  else
        call ana_sediment (tile)
#  endif
# endif
      enddo
#endif
!
!  --- Read from NetCDF file (to do) ---
!
#if defined BBL && !defined ANA_BSEDIM && !defined SEDIMENT
# ifdef AGRIF
      if (Agrif_lev_sedim == 0) call get_bsedim
# else
      call get_bsedim
# endif
#endif
!
#if defined SEDIMENT && !defined ANA_SEDIMENT
# ifdef AGRIF
      if (Agrif_lev_sedim == 0) call get_sediment
# else
      call get_sediment
# endif
#endif
!
!----------------------------------------------------------------------
!  SUBSTANCE : computing cell surfaces need for MUSTANG and BIOLink
!----------------------------------------------------------------------
#ifdef SUBSTANCE
      call substance_surfcell
#endif
!
!----------------------------------------------------------------------
!  STOGEN: initialization
!----------------------------------------------------------------------
!
#ifdef STOGEN
      do tile=0,NSUB_X*NSUB_E-1
        call ocean_2_stogen (tile)    ! provide parameters to STOGEN module
        call sto_mod_init             ! initialize STOGEN module
      enddo
#endif
!
!----------------------------------------------------------------------
!  MUSTANG : initialization
!----------------------------------------------------------------------
!
#ifdef MUSTANG
!$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call MUSTANG_init_main (tile)
      enddo
#endif
!
!----------------------------------------------------------------------
! OBSTRUCTION : initialization
!----------------------------------------------------------------------
!
#ifdef OBSTRUCTION
!$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call obst_init_main ()
        call obst_update_main (tile)
      enddo
#endif
!
!----------------------------------------------------------------------
!  Finalize grid setup
!----------------------------------------------------------------------
!
!  Finalize vertical grid now that zeta is knowned
!  zeta is also corrected here for Wetting/Drying
!  in both 2D and 3D cases
!
#if defined SOLVE3D || defined WET_DRY
!$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call set_depth (tile)
      enddo
#endif
!
!----------------------------------------------------------------------
!  Initialize diagnostic fields: mass flux, rho, omega
!----------------------------------------------------------------------
!
#ifdef SOLVE3D
!$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call set_HUV (tile)
# ifdef RESET_RHO0
        call reset_rho0 (tile)
# endif
      enddo

!$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call omega (tile)
        call rho_eos (tile)
      enddo
#endif
!
!----------------------------------------------------------------------
!  Initialize surface wave variables
!----------------------------------------------------------------------
!
#ifdef MRL_WCI
!$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call mrl_wci (tile)
      enddo
#endif
!
!----------------------------------------------------------------------
!  Set nudging coefficients
!  for sea surface height, momentum and tracers
!----------------------------------------------------------------------
!
#if defined TNUDGING  || defined ZNUDGING  \
  || defined M2NUDGING || defined M3NUDGING \
                       || defined SPONGE
!$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call set_nudgcof (tile)
      enddo
#endif
!
!----------------------------------------------------------------------
!  Read initial climatology fields from NetCDF file
!  (or any external oceanic forcing data)
!  in 3D interior or 2D boundary arrays
!----------------------------------------------------------------------
!
#if defined TCLIMATOLOGY && !defined ANA_TCLIMA
      call get_tclima
#endif
#if defined M2CLIMATOLOGY && !defined ANA_M2CLIMA ||\
   (defined M3CLIMATOLOGY && !defined ANA_M3CLIMA)
      call get_uclima
#endif
#if defined ZCLIMATOLOGY && !defined ANA_SSH
      call get_ssh
#endif
#if defined FRC_BRY && !defined ANA_BRY
      call get_bry
# ifdef BIOLOGY
      call get_bry_bio
# endif
#endif
#if !defined ANA_BRY_WKB && defined WKB_WWAVE
      call get_bry_wkb
#endif
!
!----------------------------------------------------------------------
!  Set analytical initial climatology fields
!  (or any external oceanic forcing data)
!  for sea surface height, momentum and tracers
!----------------------------------------------------------------------
!
#if (defined ZCLIMATOLOGY  && defined ANA_SSH)     || \
    (defined M2CLIMATOLOGY && defined ANA_M2CLIMA) || \
    (defined M3CLIMATOLOGY && defined ANA_M3CLIMA) || \
     defined TCLIMATOLOGY
!$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
# ifdef TCLIMATOLOGY
        call ana_tclima (tile)
# endif
# if defined M2CLIMATOLOGY && defined ANA_M2CLIMA ||\
    (defined M3CLIMATOLOGY && defined ANA_M3CLIMA)
        call ana_uclima (tile)
# endif
# if defined ZCLIMATOLOGY && defined ANA_SSH
        call ana_ssh (tile)
# endif
      enddo
#endif
!
!----------------------------------------------------------------------
!  Read surface forcing from NetCDF file
!----------------------------------------------------------------------
!
      call get_vbc
!
!----------------------------------------------------------------------
!  Read tidal harmonics from NetCDF file
!----------------------------------------------------------------------
!
#if defined SSH_TIDES || defined UV_TIDES || defined POT_TIDES
      call get_tides
#endif
!
!----------------------------------------------------------------------
! OA "Stand Alone" module : second initialization step (spatial domain)
!----------------------------------------------------------------------
!
#ifdef ONLINE_ANALYSIS
      if ( if_oa.eqv..true. ) then
!$OMP PARALLEL DO PRIVATE(tile)
         do tile=0,NSUB_X*NSUB_E-1
           call online_spectral_diags(tile,-1)
         enddo
      endif
#endif
!
!----------------------------------------------------------------------
!  Initialize XIOS I/O server
!----------------------------------------------------------------------
!
#ifdef XIOS
!$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call init_xios(tile)
      enddo
#endif
!
      if (may_day_flag /= 0) return
!
#ifdef ABL1D
!$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call abl_ini (tile)
      enddo
#endif
!
!----------------------------------------------------------------------
!  Initialization for stations
!
!  It is done here and not in init_scalars since it must be done only
!  once (whether child levels exist or not)
!----------------------------------------------------------------------
!
#ifdef STATIONS
      nrecsta=0
      ncidsta=-1
      staspval=1.E15  ! nodata flag for float variables.
      stadeltap2c=2.5 ! distance from the boundary at which a
                      ! float is transfered from parent to child
      call init_arrays_sta
      call init_sta
# ifdef SPHERICAL
      call interp_r2d_sta_ini (lonr(START_2D_ARRAY), istalon)
      call interp_r2d_sta_ini (latr(START_2D_ARRAY), istalat)
# else
      call interp_r2d_sta_ini (  xr(START_2D_ARRAY), istalon)
      call interp_r2d_sta_ini (  yr(START_2D_ARRAY), istalat)
# endif
# ifdef SOLVE3D
      call fill_sta_ini ! fills in trackaux for ixgrd,iygrd,izgrd
                        ! and ifld (either izgrd or ifld is meaningful)
# endif
      if (ldefsta) call wrt_sta
#endif /* STATIONS */
!
!----------------------------------------------------------------------
!  WKB surface wave model:
!
!  initialization and spinup to equilibrium
!----------------------------------------------------------------------
!
#ifdef WKB_WWAVE
!$OMP BARRIER
!$OMP MASTER
        MPI_master_only write(stdout,'(/1x,A/)')                       &
             'WKB: started steady wave computation.'
!$OMP END MASTER
      iic=0
      winfo=1
      iwave=1
      thwave=1.D+10
# ifndef ANA_BRY_WKB
!$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call set_bry_wkb (tile)   ! set boundary forcing
      enddo
# endif
# ifdef MRL_CEW
!$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call wkb_cew_prep (tile)  ! prepare coupling mode
      enddo
!$OMP BARRIER
      wint=0
      do winterp=1,interp_max
        wint=wint+1
        if (wint > 2) wint=1
!$OMP PARALLEL DO PRIVATE(tile)
        do tile=0,NSUB_X*NSUB_E-1
          call wkb_uvfield (tile, winterp)
        enddo
      enddo
!$OMP BARRIER
!$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call wkb_cew_finalize (tile)
      enddo
!$OMP BARRIER
# endif /* MRL_CEW */
!
!  Spinup: intergrate wave model to equilibrium
!
      do while (iwave <= 50000 .and. thwave >= 1.D-10)
        wstp=wnew
        wnew=wstp+1
        if (wnew >= 3) wnew=1
!$OMP PARALLEL DO PRIVATE(tile)
        do tile=0,NSUB_X*NSUB_E-1
# ifdef WAVE_OFFLINE
          if (iwave == 1) call set_wwave(tile)
# endif
          call wkb_wwave (tile)
        enddo
        call wkb_diag (0)
        iwave=iwave+1
        thwave=max(av_wac,av_wkn)
      enddo
# if defined CVTK_DEBUG || defined CVTK_DEBUG_ADVANCED
!$OMP BARRIER
!$OMP MASTER
      call check_tab2d(wac(:,:,wnew),'wac initialisation #1','r')
!$OMP END MASTER
# endif
!
!  Re-initialize wave forcing terms
!
      first_time=0
!$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call mrl_wci (tile)
      enddo
!$OMP BARRIER
!$OMP MASTER
        MPI_master_only write(stdout,'(/1x,A/)')                       &
             'WKB: completed steady wave computation.'
!$OMP END MASTER
#endif /* WKB_WWAVE */
!
!----------------------------------------------------------------------
!  Set initial non-Boussinesq (or fast 3D) parameters and variables
!----------------------------------------------------------------------
!
#ifdef M3FAST
# ifdef NBQ_MASS
!
! Re-evaluate Hz and Huon,Hvom using
! previously computed density rho
!
!$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call set_depth (tile)
      enddo

!$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call set_HUV (tile)
      enddo
# endif
!
! Set initial NBQ param. & var.
!
!$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call initial_nbq(tile)
      enddo
#endif
!
!----------------------------------------------------------------------
! OA "Stand Alone" module
!----------------------------------------------------------------------
!
#ifdef ONLINE_ANALYSIS
      if_online_analysis : if ( if_oa.eqv..true. ) then
!$OMP PARALLEL DO PRIVATE(tile)
         do tile=0,NSUB_X*NSUB_E-1
! BLXD online_spectral_diags/output_oa 1st argument his the model mode 0/1 = Slow/Fast mode
           call online_spectral_diags(tile,0)
         enddo
        call output_oa(0)
! BLXD Fast mode analysis
        if(if_oa_fast_mode.eqv..true.) then
!$OMP PARALLEL DO PRIVATE(tile)
         do tile=0,NSUB_X*NSUB_E-1
! BLXD online_spectral_diags/output_oa 1st argument his the model mode 0/1 = Slow/Fast mode
           call online_spectral_diags(tile,1)
         enddo
! BLXD output_oa source now set to be only called within the slow mode
!      even when the spectral online diagnostics have been calculated in the fast mode
         call output_oa(1)
      endif
      endif if_online_analysis
#endif
!
!----------------------------------------------------------------------
!  Write initial fields into history NetCDF files
!----------------------------------------------------------------------
!
#ifdef XIOS
      if (nrrec == 0) then
!$OMP PARALLEL DO PRIVATE(tile)
        do tile=0,NSUB_X*NSUB_E-1
          call send_xios_diags(tile)
        enddo
      endif
#else
      if (ldefhis .and. wrthis(indxTime)) call wrt_his
#endif

#ifdef ABL1D
!$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call abl_ini (tile)
      enddo
# ifndef XIOS
      if (ldefablhis .and. wrtabl(indxTime)) call wrt_abl_his
# endif
#endif
!
      if (may_day_flag /= 0) return
!
      end subroutine croco_initialize

      end module croco_init
