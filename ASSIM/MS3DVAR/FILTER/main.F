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
!                     MS-3DVAR FILTER DRIVER
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
!#include "private_scratch.h"
#include "scalars.h"
#include "ncscrum.h"
#include "ocean3d.h"
#include "ocean2d.h"
#include "mpi_cpl.h"
!#ifdef AGRIF
!      Type(Agrif_pgrid),pointer  :: parcours
!#endif
#include "das_ocean.h"
#include "das_covar.h"

      integer i,j,k,max_tile
      integer lstr, lenstr
      character fname*120
!
! minimization
!
!      REAL EPS,GTOL
!      INTEGER LP,IPRINT(2),IFLAG ,POINT,ITER,MP,ICALL
!      INTEGER MAXITER
!      LOGICAL DIAGCO
!      EXTERNAL VA15CD
!      COMMON /VA15DD/MP,LP, GTOL

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
!  Create parallel threads;
!  initialize (FIRST-TOUCH) model global arrays (most of them
!  are just set to to zero).
!----------------------------------------------------------------------
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
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

!#ifdef DAS_DBLE_BKG 
!      write(*,*)'DBLE BKG : get_initial with clmname'
!      call get_initial (nrrec, 1, clmname )
!#else
!      call get_initial (nrrec, 1)
      call get_initial (nrrec, 1, ininame)
!#endif
 
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

!BLXD None used
!#if defined VIS_GRID || defined DIF_GRID
!C$OMP PARALLEL DO PRIVATE(tile)
!      do tile=0,NSUB_X*NSUB_E-1      ! Rescale horizontal mixing
!        call visc_rescale (tile)     ! coefficients according to
!      enddo                          ! local grid size.
!CR      write(*,*) ' -4' MYID 
!#endif
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
!#ifdef DAS_DBLE_BKG 
!
!----------------------------------------------------------------------
!  DAS getting bkg from instantaneous croco_rst : t, u, v, zeta
!----------------------------------------------------------------------
!
!      write(*,*)'DBLE BKG : get_initial with ininame'
!      call get_initial (nrrec, 1)
!     use time from restart in this case 
!     avoid ncap2 on croco_avg time centered at noon
!     iic=0        ! required by set_depth
!    
!                               ! Set initial model clock: at this
!     time=start_time           ! moment "start_time" (global scalar)
!     tdays=time*sec2day        ! is set by get_initial or analytically
!#endif

!C$OMP PARALLEL DO PRIVATE(tile)
!      do tile=0,NSUB_X*NSUB_E-1      ! Interpolation from S-coordinate system
!        call das_z2s (tile)          ! to z-coordinates
!      enddo      

C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1      
        call das_filter (tile) 
      enddo

      lstr=lenstr(rstname)
      fname=rstname(1:lstr)
      call netcdf_filtered(fname)

!
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
