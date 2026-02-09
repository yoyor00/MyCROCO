!      implicit none
!
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
!
! CROCO website : http://www.croco-ocean.org
!======================================================================
!
#include "cppdefs.h"
!
      program main
!
!======================================================================
!
!                     MS-3DVAR DRIVER
!
!    
!
!======================================================================
!
!#ifdef XIOS
!         USE xios           ! XIOS module
!#endif
      USE buffer, ONLY : init_buffer

!# BLXD check why implicit none has been removed    
!     implicit none
      integer tile, subs, trd, ierr
#include "param.h"
#include "das_param.h"
#include "das_param_lr.h"
!#include "private_scratch.h"
#include "scalars.h"
#include "ncscrum.h"
#include "grid.h"
#include "ocean3d.h"
#include "ocean2d.h"
#include "mpi_cpl.h"
!#ifdef AGRIF
!      Type(Agrif_pgrid),pointer  :: parcours
!#endif
#include "das_ocean.h"
#include "das_ocean_smooth.h"
#include "das_ocean_lr.h"
#include "das_ocean9.h"
#include "das_covar.h"
#include "das_lbgfs.h"
#include "das_innov.h"

      integer i,j,k,max_tile !, ios, nsyn
!
! minimization
!
      REAL EPS,GTOL
      INTEGER LP,IPRINT(2),IFLAG ,POINT,ITER,MP,ICALL
      INTEGER MAXITER
      LOGICAL DIAGCO
      EXTERNAL VA15CD
      COMMON /VA15DD/MP,LP, GTOL
      real costpp
!      LOGICAL, parameter :: flag_lr_inc = .false.
#if defined DAS_REDUCED_ZCORR
      real fredvertcorr
#endif
#ifdef MPI
      include 'mpif.h'
!     real*8 start_time2, start_time1, exe_time
#endif
#ifdef AGRIF
      integer size_XI,size_ETA,se,sse, sz,ssz
      external :: step
# include "zoom.h"
# include "dynparam.h"
#endif
#ifdef USE_CALENDAR
      character*19 :: tool_sectodat
#endif
#ifdef AUTOTILING
! Initial values of NSUB_X, NSUB_E
      NSUB_X = 1
      NSUB_E = NPP
#endif
!
! BLXD warning check if dynderivparam include with size_XI, size_ETA, N2d, N3d...
! required for DAS (warning : dynderivparam must sometimes replace
! DAS private_scratch since it no longer holds siz_XI, size_ETA,...)
#include "dynderivparam.h"
!
!----------------------------------------------------------------------
!  Initialize communicators and subgrids decomposition:
!  MPI parallelization, XIOS server, OASIS coupling, AGRIF nesting,
!----------------------------------------------------------------------
!
#ifdef MPI
# if (!defined AGRIF)
      call MPI_Init (ierr)
# endif
!
!  XIOS, OASIS and AGRIF: split MPI communicator
!  (XIOS with OASIS is not done yet)
!
!# if (defined XIOS && !defined AGRIF)
!      call xios_initialize( "crocox",return_comm=MPI_COMM_WORLD )
!# elif defined AGRIF
#if defined AGRIF
      call Agrif_MPI_Init(MPI_COMM_WORLD)
# endif
#endif /* MPI */
!
!  Initialize AGRIF nesting
!
#ifdef AGRIF
      call Agrif_Init_Grids()
      call declare_zoom_variables()
#endif
!
!  Initialize automatic tiling
!
#ifdef AUTOTILING
      call init_auto_tiling
#endif
!
!  Setup MPI domain decomposition
!
#ifdef MPI
!     start_time1=PMPI_Wtime()
      call MPI_Setup (ierr)
      if (ierr.ne.0) goto 100        !--> ERROR
c**   call MPI_Test
c**   goto 100
#endif
!
!  Initialize debug procedure
!
#if defined RVTK_DEBUG || defined RVTK_DEBUG_ADVANCED || \
    defined RVTK_DEBUG_PERFRST

      call debug_ini
#endif
!
#define CR  !
!? #define SINGLE NSUB_X*NSUB_E,NSUB_X*NSUB_E !!!!
      call init_buffer
!
!----------------------------------------------------------------------
!  Read in tunable model parameters in roms.in file
!----------------------------------------------------------------------
!
      call read_inp (ierr)
      if (ierr.ne.0) goto 100
!
!----------------------------------------------------------------------
!  Initialize global model parameters
!----------------------------------------------------------------------
!
!  Gobal scalar variables
!
      call init_scalars (ierr)
      if (ierr.ne.0) goto 100

      call das_init_scalars

!
!
!----------------------------------------------------------------------
!  Create parallel threads; start timers for each thread;
!  initialize (FIRST-TOUCH) model global arrays (most of them
!  are just set to to zero).
!----------------------------------------------------------------------
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call start_timers()
        call init_arrays (tile)
        call das_init_arrays (tile)  ! initialize (FIRST-TOUCH) DAS
      enddo
CR      write(*,*) '-11' MYID
!
!----------------------------------------------------------------------
!  Set horizontal grid, model bathymetry and Land/Sea mask
!----------------------------------------------------------------------
!
#ifdef ANA_GRID
!
!  Set grid analytically
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call ana_grid (tile)
      enddo
# if defined RVTK_DEBUG || defined RVTK_DEBUG_ADVANCED
C$OMP BARRIER
C$OMP MASTER
       call check_tab2d(h(:,:),'h initialisation #1','r')
C$OMP END MASTER
# endif
#else
!
!  Read grid from GRID NetCDF file
!
      call get_grid
      if (may_day_flag.ne.0) goto 99 !-->  EXIT
#endif
!
!  Compute various metric term combinations.
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call setup_grid1 (tile)
      enddo
CR      write(*,*) '-10' MYID
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call setup_grid2 (tile)
      enddo
CR      write(*,*) ' -9' MYID
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
!
!  Create three-dimensional S-coordinate system,
!  which may be needed by ana_initial
!  (here it is assumed that free surface zeta=0).
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call set_depth (tile)
      enddo
CR      write(*,*) ' -8' MYID
!
!  Make grid diagnostics
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call grid_stiffness (tile)
      enddo
CR      write(*,*) ' -7' MYID
#endif
!
!----------------------------------------------------------------------
!  Set initial conditions for momentum and tracer variables
!----------------------------------------------------------------------
!
!  Read from NetCDF file
!

#ifdef DAS_DBLE_BKG 
      write(*,*)'DBLE BKG : get_initial with clmname'
      call get_initial (nrrec, 1, clmname )
#else
!      call get_initial (nrrec, 1)
      call get_initial (nrrec, 1, ininame)
#endif
 
      iic=0        ! required by set_depth TODO KEEP ?
     
                                ! Set initial model clock: at this
      time=start_time           ! moment "start_time" (global scalar)
      tdays=time*sec2day        ! is set by get_initial or analytically
                                ! --> copy it into threadprivate "time"
!
!#ifdef USE_CALENDAR
!      ! compute start time and ntimes
!# ifdef ANA_INITIAL
!      ! if no initial file set start_date from croco.in
!      if (nrrec.eq.0) time=tool_datosec(run_start_date)
!# endif
!      time_end=tool_datosec(run_end_date)
!      ntimes=int((time_end-time)/dt)
!      MPI_master_only write(stdout,*)
!     &     'Ntimes from date_start and date_end:',ntimes
!      ntimes=ntimes+ntstart
!#endif
!
!  Set initial conditions analytically for ideal cases
!  or for tracer variables not present in NetCDF file
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call ana_initial (tile)
      enddo
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
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call set_depth (tile)
      enddo
CR      write(*,*)  ' -5' MYID
#endif

!
!----------------------------------------------------------------------
!  DAS geopotential grid setup
!----------------------------------------------------------------------
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call das_setup1 (tile)
      enddo
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call das_setup2 (tile)
      enddo
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1      ! Interpolation from S-coordinate system
        call das_s2z_3d (tile)       ! to z-coordinates
      enddo


#ifdef DAS_DBLE_BKG 
!
!----------------------------------------------------------------------
!  DAS getting bkg from instantaneous croco_rst : t, u, v, zeta
!----------------------------------------------------------------------
!
      write(*,*)'DBLE BKG : get_initial with ininame'
      call get_initial (nrrec, 1)
!     use time from restart in this case 
!     avoid ncap2 on croco_avg time centered at noon
!     iic=0        ! required by set_depth
!    
!                               ! Set initial model clock: at this
!     time=start_time           ! moment "start_time" (global scalar)
!     tdays=time*sec2day        ! is set by get_initial or analytically
#endif


!
!----------------------------------------------------------------------
!  DAS getting increment
!----------------------------------------------------------------------
!

!      no_lr_inc : if (flag_lr_inc) then
#if defined DAS_READ_INC

      print*, '... Bef Reading below scale increment of analysis'
      print*,flag_lr_s_da, flag_lr_da

      call das_get_initial        ! read the increment from the low resolution DA

      print*, '... Aft Reading below scale increment of analysis'
      print*,flag_lr_s_da, flag_lr_da

!
      if (flag_lr_da .eq. 0 ) then
C$OMP PARALLEL DO PRIVATE(tile)
        do tile=0,NSUB_X*NSUB_E-1
!          call das_interp_lr (tile)
          call das_add_lr (tile)
        enddo
      endif
#endif
!      endif no_lr_inc
! BLXD : sea over land
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1      ! just for calculation possible at land grid point
        call das_fill9 (tile)
      enddo

!      print*, 'T'
!      print*, (t_lr_s(i,Mm+1,ndas,1),i=0,Lm+1)
!      print*, 'S'
!      print*, (t_lr_s(i,Mm+1,ndas,2),i=0,Lm+1)
!      print*, 'u'
!      print*, (u_lr_s(i,Mm+1,ndas),i=0,Lm+1)
!      print*, 'v'
!      print*, (v_lr_s(i,Mm+1,ndas),i=0,Lm+1)
!      print*, 'Z'
!      print*, (zeta_lr_s(i,Mm+1),i=0,Lm+1)

!C$OMP PARALLEL DO PRIVATE(tile)
!      do tile=0,NSUB_X*NSUB_E-1      ! Interpolation from S-coordinate system
!        call das_z2s (tile)          ! to z-coordinates
!      enddo      

!
!----------------------------------------------------------------------
!  Initialize diagnostic fields: mass flux, rho, omega
!----------------------------------------------------------------------
!
#ifdef SOLVE3D
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call set_HUV (tile)
# ifdef RESET_RHO0
        call reset_rho0 (tile)
# endif
      enddo
CR      write(*,*)  ' -4' MYID
#endif
!!
!!----------------------------------------------------------------------
!!  Read initial climatology fields from NetCDF file
!!  (or any external oceanic forcing data)
!!  in 3D interior or 2D boundary arrays
!!----------------------------------------------------------------------
!!
!#if defined TCLIMATOLOGY && !defined ANA_TCLIMA
!      call get_tclima
!#endif
!#if defined M2CLIMATOLOGY && !defined ANA_M2CLIMA ||\
!   (defined M3CLIMATOLOGY && !defined ANA_M3CLIMA)
!      call get_uclima
!#endif
!#if defined ZCLIMATOLOGY && !defined ANA_SSH
!      call get_ssh
!#endif
!!
!!----------------------------------------------------------------------
!!  Set analytical initial climatology fields
!!  (or any external oceanic forcing data)
!!  for sea surface height, momentum and tracers
!!----------------------------------------------------------------------
!!
!#if (defined ZCLIMATOLOGY  && defined ANA_SSH)     || \
!    (defined M2CLIMATOLOGY && defined ANA_M2CLIMA) || \
!    (defined M3CLIMATOLOGY && defined ANA_M3CLIMA) || \
!     defined TCLIMATOLOGY
!C$OMP PARALLEL DO PRIVATE(tile)
!      do tile=0,NSUB_X*NSUB_E-1
!# ifdef TCLIMATOLOGY
!        call ana_tclima (tile)
!# endif
!# if defined M2CLIMATOLOGY && defined ANA_M2CLIMA ||\
!    (defined M3CLIMATOLOGY && defined ANA_M3CLIMA)
!        call ana_uclima (tile)
!# endif
!# if defined ZCLIMATOLOGY && defined ANA_SSH
!        call ana_ssh (tile)
!# endif
!      enddo
!#endif
!!
!CR      write(*,*) ' -2' MYID
!
!----------------------------------------------------------------------
!  Read surface forcing from NetCDF file
!----------------------------------------------------------------------
!
!      call get_vbc
!
!----------------------------------------------------------------------
!  Read tidal harmonics from NetCDF file
!----------------------------------------------------------------------
!
!#if defined SSH_TIDES || defined UV_TIDES || defined POT_TIDES
!      call get_tides
!#endif
!
!CR      write(*,*) ' -1' MYID
!
!
!----------------------------------------------------------------------
!  Initialize XIOS I/O server
!----------------------------------------------------------------------
!
!#ifdef XIOS
!C$OMP PARALLEL DO PRIVATE(tile)
!      do tile=0,NSUB_X*NSUB_E-1
!        call init_xios(tile)
!      enddo
!#endif
!
!      if (may_day_flag.ne.0) goto 99 !-->  EXIT
!
!
!
#undef CR
!#ifdef MPI
!      call MPI_Barrier(MPI_COMM_WORLD, ierr)
!#endif

#ifdef AGRIF
      iind = -1
      grids_at_level = -1
      sortedint = -1
      call computenbmaxtimes
#endif

!#ifdef MPI_TIME
! Start Chrono
!      start_time1=PMPI_Wtime()
!#endif

!
!... prepare some basic variables
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1   
        call das_extrapolate9 (tile)
        call das_rho_eos9 (tile)   
#if defined DAS_BKGZETA
        call das_diag_steric_sshb (tile)     ! zetab_das9
#endif
      enddo
!
!
! ... get the file names of observations and error covariances
!
      call das_read_inp
!
! ... get the forecast error covariance
!
# ifdef DAS_ANA_BVAR
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call das_ana_bvar(tile)
      enddo
# else
      call das_get_bvar
# endif

      do k=1,ndas
        do j=0,Mm+1
          do i=0,Lm+1
!good             bt_das(i,j,k,itemp)=1.0*bt_das(i,j,k,itemp)+0.03
!good             bt_das(i,j,k,isalt)=0.6*bt_das(i,j,k,isalt)+0.01
!             bt_das(i,j,k,itemp)=1.0*bt_das(i,j,k,itemp)  !v1
!             bt_das(i,j,k,isalt)=0.6*bt_das(i,j,k,isalt)  !v1
!             bt_das(i,j,k,itemp)=0.85*bt_das(i,j,k,itemp)+0.35
!             bt_das(i,j,k,isalt)=0.25*bt_das(i,j,k,isalt)+0.05
!             bt_das(i,j,k,itemp)=6.*bt_das(i,j,k,itemp)
!            bt_das(i,j,k,isalt)=0.8*bt_das(i,j,k,isalt)
!             bt_das(i,j,k,isalt)=5.*bt_das(i,j,k,isalt)


!             bt_das(i,j,k,isalt)=0.65*bt_das(i,j,k,isalt)
!             bpsi_das(i,j,k)=1.00*bpsi_das(i,j,k)
!             bchi_das(i,j,k)=0.25*bchi_das(i,j,k)
! HF2
!             bpsi_das(i,j,k)=0.5*bpsi_das(i,j,k)
!             bchi_das(i,j,k)=0.25*bchi_das(i,j,k)
!bench        bpsi_das(i,j,k)=0.35*10.*bpsi_das(i,j,k)
!             bchi_das(i,j,k)=0.35*5.*bchi_das(i,j,k)
# if defined DAS_ADIMPSICHI
             bpsi_das(i,j,k)=0.28*10.*bpsi_das(i,j,k)/adp
             bchi_das(i,j,k)=0.28*5.*bchi_das(i,j,k)/adp
# else
             bpsi_das(i,j,k)=0.28*10.*bpsi_das(i,j,k)
             bchi_das(i,j,k)=0.28*5.*bchi_das(i,j,k)
# endif
          enddo
        enddo
      enddo
      do j=0,Mm+1
        do i=0,Lm+1
!           bz_das(i,j)=1.5*bz_das(i,j)
           bz_das(i,j)=bz_das(i,j)  !v1
! BLDX after GADT.nc issue
!
!           bz_das(i,j)=0.8*(0.6*exp(-h(i,j)*h(i,j)/(120.*120.))+1.0)
!     &                 * bz_das(i,j) 
!            bz_das(i,j)=bz_das(i,j)+0.01
        enddo
      enddo
      call das_set_geo_ratio

!      write(*,*) 'geo_ratio'
!      write(*,*) (georatio(i,50),i=1,LLm)
!      write(*,*) 'geo_ratio'
!      write(*,*) (georatio(200,j),j=1,MMm)

      write(*,*) 'bt_var'
      write(*,*) (bt_das(50,50,k,itemp),k=1,ndas)
      write(*,*) 'bs_var'
      write(*,*) (bt_das(30,40,k,isalt),k=1,ndas)
      write(*,*) 'bpsi_var'
      write(*,*) (bpsi_das(30,40,k),k=1,ndas)
      write(*,*) 'bchi_var'
      write(*,*) (bchi_das(30,40,k),k=1,ndas)
      write(*,*) 'bz_var'
      write(*,*) (bz_das(30,k), k=20,40)
!
# ifdef DAS_BVAR_CORR
      call das_get_corr
# endif
!
#if defined DAS_REDUCED_ZCORR
      do j=1,NDAS
      do i=1,NDAS
        fredvertcorr = (1.+tanh(-(z_das(i)+125.)/150.))/2.
        ctZ_das(i,j,itemp)=fredvertcorr*ctZ_das(i,j,itemp)
        ctZ_das(i,j,isalt)=fredvertcorr*ctZ_das(i,j,isalt)
        cuZ_das(i,j)=fredvertcorr*cuZ_das(i,j)
        cvZ_das(i,j)=fredvertcorr*cvZ_das(i,j)
      enddo
      enddo
#endif

# if defined DAS_JASONSSH || defined DAS_TPSSH\
 || defined DAS_ERS2SSH || defined DAS_SWOTSSH
      call das_get_ssh_ref
# endif

!
! ... read in all the observational data
!
                       ! IN-SITU, whoi, sio, ptsur, martin vertical profiles
                       !          NPS flight sst
!
      call das_get_TS_obs
!
!      write (*,*), ' TS Prf Obs time'
!      write (*,*), (time_obs(i),i=1,prf_num_TS_obs)
!      write (*,*), ' TS Prf Obs I'
!      write (*,*), (I_obs(i),i=1,prf_num_TS_obs)
!      write (*,*), ' TS Prf Obs salt'
!      write (*,*), (s_obs_raw(prf_num_TS_obs,k),k=1,ndas-1)
!
      call das_get_insitu
                       ! HF radar
      call das_get_hfradar
!
!
! ... get CalPoly & DORADO auv
!
      call das_get_auv
!      print*,'calpoly auv ', num_cal
!      print*,'dorado  auv ', num_dor
!
                       ! SST from TMI and dist_coast

! ... get distance to coast at rho points 
      call das_get_d2c
!      write(*,*) 'dist_coast'
!      write(*,*) (dist_coast(i,20),i=0,Lm+1)

# ifdef DAS_TMISST
      call das_get_tmi
# endif

                       ! SST from AVHRR
# ifdef DAS_MCSST
      call das_get_mc
# endif
# ifdef DAS_GOES_SST
      call das_get_goes
# endif
                       ! SSH from JASON-1
# ifdef DAS_SWOTSSH
      call das_get_swot
# endif

# ifdef DAS_JASONSSH
      call das_get_js1
# endif
                       ! SSH from TOPEX/POSAIDON
# ifdef DAS_TPSSH
      call das_get_tp
# endif
!
                       ! SSS from Aquarius
# ifdef DAS_SATSSSS
      call das_get_sats
# endif
!
      write (*,'(9(2x,A))') 'obs num:', 'hfradar','hfradar6','mapped',
     &                       'flight','whoi','sio',
     &                       'ptsur','martin'
      write (*,'(10x,8I7)') num_hfradar,num_hfradar6,prf_num_map,
     &                      num_flight,prf_num_whoi,prf_num_sio,
     &                      prf_num_ptsur,prf_num_martn

                                          ! when no any observation
                                          ! data assimilation not executed
! BLXD Why .not. flag_sats
      if(.not. flag_tmi .and. .not. flag_mc .and. .not. flag_goes
     &   .and. .not. flag_swot
     &   .and. .not. flag_js1 .and. .not. flag_tp 
     &   .and. .not. flag_sio .and. .not. flag_whoi
     &   .and. .not. flag_ptsur .and. .not. flag_martn 
     &   .and. .not. flag_map .and. .not. flag_prof 
     &   .and. .not. flag_ship
     &   .and. .not. flag_cal .and. .not. flag_dor
     &   .and. .not. flag_flight
     &   .and. .not. flag_hfradar.and. .not. flag_hfradar6) then
        write(*,'(6x, A, 2(/6x, A))') 
     &              '!!! NO ANY OBSERVATION AVAILABLE',
     &              '!!! DATA ASSIMILATION NOT EXECUTED',
     &              '!!! NO ANALYSIS FILE'

C$OMP PARALLEL DO PRIVATE(tile)
        do tile=0,NSUB_X*NSUB_E-1
          CALL DAS_CA2ZERO_S (tile)
        enddo

        goto 98       !!! The increment from DA lr will be added 
                      !!!  to the background

      endif

!      if (flag_mc) then
!        do i=1,Lm
!        do j=1,Mm
!        if (mc_mask(i,j) .gt. 0.0 ) then
!          write (*,*) sst_mc(i,j)
!        endif
!        enddo
!        enddo
!      endif
!
! ... compute the innovation vectors
!
!
      max_tile=max(numthreads,prf_num_martn)
C$OMP PARALLEL DO PRIVATE(k), SHARED(max_tile)
      do k=1,max_tile
        call DAS_INNOV_MARTN(k)
      enddo
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call das_innov_TS_obs(tile)
      enddo
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call das_innov_sio(tile)
        call das_innov_whoi(tile)
        call das_innov_ptsur(tile)
        call das_innov_flight(tile)
        call das_innov_calauv(tile)
        call das_innov_dorauv(tile)
# if defined DAS_SATSSSS
        call das_innov_sats(tile)
# endif
!
#if defined DAS_MCSST || defined DAS_TMISST
        call das_innov_sst(tile)
# endif
#if defined DAS_HFRADAR
        call das_innov_hfradar (tile)
#endif
      enddo
!
! smooth gridded T/S
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_TS_DAS_SMOOTH (tile)     ! smooth zeta_das, same as smooth zeta_s
      enddo
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call das_innov_map(tile)      ! smoothed T/S
      enddo
!
! smooth gridded SSH
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_ZETA_DAS_SM (tile)     ! smooth zeta_das, same as smooth zeta_s
      enddo
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
# ifdef DAS_SWOTSSH

        call das_innov_swot(tile)
# endif
      enddo
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
# if defined DAS_JASONSSH
        call das_innov_js1(tile)
# endif
# if defined DAS_TPSSH
        call das_innov_tp(tile)
# endif
      enddo

!
! print out innovations
!

!
! START MINIMIZATION
!
      ICALL=0
      IFLAG=0
      MAXITER=43
      IPRINT(1)= 1
      IPRINT(2)= 0
      DIAGCO= .FALSE.
      EPS= 1.0D-3     ! minimization criterion
                      ! |grad|/max(1, |x|) < EPS, terminated
      LP = 6
      MP = 6

!
! ... initialize the first guess of increment to be zero
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_CA2ZERO (tile)
      enddo
      CALL DAS_CA2VEC    ! control variable vector
!
! ... minimization interation starts
!
  10  CONTINUE
!
! ... initialize the cost function
!
      COSTF=0.0
!
! ... compute the matrix-vector product 
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL das_cross_ts(tile)
      enddo

C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_MATVEC(tile)  
      enddo

!
! ... terms of the SSH
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_VARVEC (tile)      ! t_s
      enddo
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_COST_REG(tile)
      enddo

C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_LIN_RHO_EOS (tile) ! rho_s
        CALL DAS_DIAG_SSH_STERIC (tile)    ! zeta_h, steric SSHs and zeta_s
      enddo
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_ZETA_S_SM (tile)     ! smooth zeta_sm
      enddo
!
#if defined DAS_PWTHZETATOT
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_ZETA_H_SM (tile)
      enddo
#endif
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
# ifdef DAS_SWOTSSH
        CALL DAS_COST_SWOT(tile)
# endif
      enddo
!
! also smoothed
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
# if defined DAS_JASONSSH || defined DAS_TPSSH
        CALL DAS_COST_JS1(tile)   ! cost function SSH JS1+TP
# endif
      enddo
!
! in the pressure, not non-steric SSH zeta_h
! for a high resolution model, should be conservative to correct T/S using velocities
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_DIAG_PRESSURE (tile)     ! geostrophic p, T/S are geostrophic?
      enddo
!
! smoothed for calculation of geostrophic velocity
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_PRESSURE_SMOOTH (tile)
      enddo
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_COPY_PRESSURE_SM (tile)
      enddo

       costpp=costf
!      write(*,*) 'zeta_s+H'
!      write(*,*) (zeta_s(53,j),j=53,Mm)
#if defined DAS_HFRADAR
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_GEOS_UV_P (tile)      !u_s, v_s before psichi
        CALL DAS_GEO_RATIO (tile)
        CALL DAS_PSICHI_UV (tile)      ! add psi/chi
      enddo
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_COST_HFRADAR (tile)
      enddo
!C$OMP PARALLEL DO PRIVATE(tile)
!      do tile=0,NSUB_X*NSUB_E-1
!        CALL DAS_COST_REG_UV(tile)   ! cost function SST
!      enddo
#endif
      costpp=costf-costpp
      write(*,*) '     HF COST=', costpp
!
! ... term of the SST
!
#if defined DAS_MCSST || defined DAS_TMISST
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_COST_SST(tile)   ! cost function SST
      enddo
#endif
!
! flight
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_COST_FLIGHT(tile)   ! cost function SST
      enddo
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_COST_CALAUV(tile)   
      enddo
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_COST_DORAUV(tile)   
      enddo
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_COST_GLIDER_SIO(tile)   
      enddo
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_COST_GLIDER_WHOI(tile)
      enddo
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL das_cost_TS_obs(tile)
      enddo
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_COST_GLIDER_PTSUR(tile)
      enddo
!
! smoothed    t_s => t_sm
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_TS_S_SM(tile)
      enddo
!
!C$OMP PARALLEL DO PRIVATE(tile)
!      do tile=0,NSUB_X*NSUB_E-1
!        CALL DAS_COPY_TS_S_SMOOTH(tile)
!      enddo
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_COST_MAP(tile)
      enddo
!
! satellite
!
!
#if defined DAS_SATSSSS
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_COST_SATS(tile)   ! cost function SSS
      enddo
#endif
!
! ... gliders
!
      CALL DAS_COST_GLIDER
!
! ... term of the background
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_COSTB(tile)   ! cost function background
      enddo
!
! ... compute the gradient with respect to u, v, t, s
!
! ... clean the adjoint arrays
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_ADJ_ZERO(tile)
      enddo
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_ADJ_COSTB(tile)
      enddo
!
! glider
!
      CALL DAS_ADJ_CGLIDER
!
#if defined DAS_SATSSSS
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_ADJ_CSATS(tile)     ! SSS
      enddo
#endif
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_ADJ_MAP(tile)
      enddo
!C$OMP PARALLEL DO PRIVATE(tile)
!      do tile=0,NSUB_X*NSUB_E-1
!        CALL DAS_ADJ_COPY_TS_S_SMOOTH(tile)
!      enddo
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_ADJ_TS_S_SM(tile)
      enddo
!
! CalPoly AUV
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_ADJ_GLIDER_PTSUR(tile)
      enddo
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL das_adj_TS_obs(tile)
      enddo
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_ADJ_GLIDER_WHOI(tile)
        CALL DAS_ADJ_GLIDER_SIO(tile)
        CALL DAS_ADJ_DAUV(tile)
        CALL DAS_ADJ_CAUV(tile)
      enddo
!
! flight
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_ADJ_CFLIGHT(tile)
      enddo
!
#if defined DAS_MCSST || defined DAS_TMISST
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_ADJ_CSST(tile)
      enddo
#endif
!
#if defined DAS_HFRADAR
!C$OMP PARALLEL DO PRIVATE(tile)
!      do tile=0,NSUB_X*NSUB_E-1
!        CALL DAS_ADJ_COST_REG_UV(tile)
!      enddo
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_ADJ_CHFRADAR (tile)
      enddo

C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_ADJ_PSICHI_UV (tile)
      enddo

      !IF(ICALL.EQ.MAXITER)THEN
      !  call das_adj_wrt
      !  zeta_s_adj, etc nul
      !ENDIF

C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_ADJ_GEO_RATIO (TILE)
      enddo
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_ADJ_GEOS_UV_P (tile)
      enddo
#endif   /* DAS_HFRADAR */

C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_ADJ_COPY_PRESSURE_SM (tile)
      enddo
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_ADJ_PRESSURE_SM (tile)
      enddo

C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_ADJ_DIAG_P (tile)
      enddo

!      write(*,*) 'zeta_s_adj'
!      write(*,*) (zeta_s_adj(10,j),j=50,Mm)
!      write(*,*) 'zeta_s_adj'
!      write(*,*) (zeta_s_adj(50,j),j=50,Mm)
!
! BLXD clean this double cpp keys 
!      openMP overhead if ifdef within the PARALLEL REGION?
# if defined DAS_JASONSSH || defined DAS_TPSSH
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_ADJ_CJS1(tile)   ! cost function SSH:JS1+TP
      enddo
#endif
#if defined DAS_SWOTSSH
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_ADJ_CSWOT(tile)   ! cost function SSH
      enddo
#endif
#if defined DAS_PWTHZETATOT
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_ADJ_ZETA_H_SM (tile)
      enddo
#endif
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_ADJ_ZETA_S_SM (tile)
      enddo
!
!      write(*,*) 'zeta_s 50'
!      write(*,*) (zeta_s(53,j),j=50,Mm)
!
!      write(*,*) 'ssh_js1 53'
!      write(*,*) (ssh_js1(53,j),j=50,Mm)
!
!      write(*,*) 'zeta_s_adj cjs1'
!      write(*,*) (zeta_s_adj(53,j),j=50,Mm)

C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_ADJ_DIAG_SSH_STERIC (tile)
        CALL DAS_ADJ_RHO_EOS (tile)
      enddo
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_ADJ_COST_REG (tile)
        CALL DAS_ADJ_VARVEC (tile)
      enddo
!
! ... compute the matrix-vector product
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_ADJ_MATVEC(tile)
      enddo

C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_ADJ_CROSS_TS(tile)
      enddo

#if defined DAS_WRT_ADJ_INI
      IF(ICALL.EQ.0)THEN
      call das_theadj_wrt(1)
      ENDIF
#endif

!
! ... clean the intermediate adjoint arrays
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL DAS_ADJ_ZERO_S(tile)
      enddo
!
      CALL DAS_GA2VEC    ! gradient vector
!
      CALL mbfgs(NDIM,MSAVE,COSTF,DIAGCO,IPRINT,EPS,POINT,IFLAG)

      CALL DAS_VEC2CA        !vector to control variables

      IF(IFLAG.LE.0) GO TO 50
      ICALL=ICALL + 1
      IF(ICALL.GT.MAXITER) GO TO 50
  
      GOTO 10
 
  50  CONTINUE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   minimization done
!
      print*, '     FLAG=', IFLAG
      if (IFLAG .lt. 0) then
        write(*, '(6x, A, (/6x, A))')
     &   '!!! MINIMIZATION NOT COMPLETED',
     &   '!!! NO SS ANALYSIS INC GENERATED'

C$OMP PARALLEL DO PRIVATE(tile)
        do tile=0,NSUB_X*NSUB_E-1
          CALL DAS_CA2ZERO_S (tile)
        enddo
        goto 98

      else
        write(*, '(6x, A)')
     &   ' MINIMIZATION COMPLETED'
      endif
!
! ... compute incremental analysis
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL das_cross_ts(tile)
      enddo

C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1   !temp and salt
        call das_matvec (tile)
      enddo
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call das_varvec (tile)       ! t_s
      enddo

!      write(*,*) 'zeta_s'
!      write(*,*) (zeta_s(53,j),j=50,Mm)
!
#define DAS_SMOOTH
#ifdef DAS_SMOOTH
                    ! smooth t_s to t_sm
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call das_innov_smooth (tile)
      enddo
                    ! then t_sm to t_s
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call das_copy_smooth (tile)
      enddo
#endif
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call das_lin_rho_eos (tile)  ! rho_s
        call das_diag_ssh_steric (tile)     ! zeta_s -> zeta_s
      enddo
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call das_zeta_smooth (tile)  ! zeta_sm   use sz_rad
      enddo
#if defined DAS_PWTHZETATOT
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL das_zeta_h_sm (tile)
      enddo
#endif
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call das_copy_zeta_smooth (tile)   !zeta_sm -> zeta_s
        call das_diag_pressure (tile)
      enddo
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL das_pressure_smooth (tile)     ! use sz_rad
      enddo
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        CALL das_copy_pressure_sm (tile)       !use sz_rad
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!
#undef PRINF_GEO

C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call das_geos_uv_p (tile)      ! u_s and v_s
        call das_geo_ratio (tile)
#if defined DAS_D2CUV
        call das_d2c_uv (tile)
#endif
#if !defined PRINF_GEO
        call das_psichi_uv (tile)      ! u_s,vs -> us,vs
#endif    
      enddo
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call das_innov_uv_smooth (tile)
      enddo
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1                          
        call das_copy_uv_smooth (tile)                  
      enddo

#ifdef PRINF_GEO
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call das_diag_ssh (tile)     ! zeta
      enddo
#endif   
                          ! increment from LR if yes
  98  continue                       
      
!
! adding the increment for zeta
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
!        if(flag_lr_inc)
!     &    call das_add_zeta_lr_s (tile)  !the large scale component, bias
#if defined DAS_READ_INC
        call das_add_zeta_lr_s (tile)  !the large scale component, bias
#endif
        call das_analysis_zeta (tile)
      enddo
!
! ... After updating zeta, z_r is updated before vertical interpolations
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1      ! S-coordinate system, which
        call set_depth (tile)        ! needed in analysis
      enddo
!
! ...Vertical interp 3D variables and add the increment
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
!       if(flag_lr_inc)
!    &    call das_add_lr_s (tile)
#if defined DAS_READ_INC
        call das_add_lr_s (tile)
#endif
        call das_analysis (tile)
      enddo

C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call rho_eos (tile)
      enddo

      write(*,*) 'rho 50 '
      write(*,*) (rho(50,50,k), k=1,n)
!
! ... create new restart file
! ... with the analysis from 3DVAR
!
      iic=ntstart
      call wrt_rst
!      call das_wrt_rst   
!...
!
      write(*,*) '     iic==', iic

      write(*,*) 't_s j=20 '
      write(*,*) (t_s(i,20,ndas,1), i=0,60)

!      write(*,*) 'psi 80 '
!      write(*,*) (psi_s(15,j,ndas), j=1,Mm)
!
!      write(*,*) 'chi 80 '
!      write(*,*) (chi_s(15,j,ndas), j=1,Mm)
! ... with the analysis from 3DVAR
!
#if defined DAS_WRT_ADJ_END
      call das_theadj_wrt(0)
#endif
      call das_wrt


                                     ! JOB DONE  
  99  continue                       ! SHUTDOWN:
!#ifdef MPI_TIME
!      start_time2=PMPI_Wtime()
!      exe_time = start_time2 - start_time1
!      MPI_master_only write(stdout,*) '*******************************'
!      MPI_master_only write(stdout,*) 'time time-stepping : ',exe_time
!      MPI_master_only write(stdout,*) '*******************************'
!#endif
!C$OMP PARALLEL DO PRIVATE(tile)            ! Stop timers and DAS Commented out
!      do tile=0,NSUB_X*NSUB_E-1            ! close netCDF files.
!        call stop_timers()
!      enddo
      call closecdf
!#ifdef XIOS
!      call iom_context_finalize( "crocox")   ! needed for XIOS+AGRIF
!#endif
#ifdef AGRIF
!
!  Close the netcdf files also for the child grids
!
      parcours=>Agrif_Curgrid%child_list % first
      do while (associated(parcours))
        Call Agrif_Instance(parcours % gr)
        call closecdf
!# ifdef XIOS
!        call iom_context_finalize( "crocox") ! needed for XIOS+AGRIF
!# endif
        parcours => parcours % next
      enddo
#endif

 100  continue

#ifdef MPI
      if (ierr.ne.0) call mpi_abort (MPI_COMM_WORLD, ierr)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)  ! XIOS

      start_time2=PMPI_Wtime()
      exe_time = start_time2 - start_time1
      if (mynode.eq.0) print *,'exe_time =',exe_time

!# if defined XIOS
!                                ! case XIOS + (OASIS / no OASIS)
!                                !      > MPI finalize done by XIOS
!      call xios_finalize()      !      > if OASIS, finalize is done by XIOS
!                                !      > if AGRIF : done using iom_context_finalize
!!!#   if !defined AGRIF
!!!      call MPI_Finalize (ierr)  !       if No AGRIF : MPI_Finalize is needed
!!!#   endif
!
!#  if  !defined OA_COUPLING && !defined OW_COUPLING  && !defined AGRIF
!       call MPI_Finalize (ierr)      !  if No coupling + No AGRIF : MPI_Finalize is needed
!#  endif
!
!# elif defined OA_COUPLING || defined OW_COUPLING
!                                          ! case no XIOS + OASIS
!      call prism_terminate_proto(ierr)    !   > Finalize OASIS3 (without XIOS)
!!!!      call MPI_Finalize (ierr)         !   > Finalize CROCO (without XIOS)
!# else
      call MPI_Finalize (ierr)            !   > Finalize CROCO (without XIOS)
!# endif
#endif

      stop
      end
