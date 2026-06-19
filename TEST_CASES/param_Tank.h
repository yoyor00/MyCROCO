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
!----------------------------------------------------------------------
! Dimensions of Physical Grid and array dimensions
!----------------------------------------------------------------------
!
! LLm,MMm  Number of the internal points of the PHYSICAL grid.
!          in the XI- and ETA-directions [physical side boundary
!          points and peroodic ghost points (if any) are excluded].
! N        Number of vertical levels.
!
      integer, parameter :: LLm0 = 50
      integer, parameter :: MMm0 = 1
      integer, parameter :: N = 50

!
!----------------------------------------------------------------------
! ABL1D
!----------------------------------------------------------------------
!
# if defined ABL1D
! N_abl Number of layer use in ABL1D
      integer, parameter :: N_abl=50
# endif

!
!----------------------------------------------------------------------
! Tides
!----------------------------------------------------------------------
!
#if defined SSH_TIDES || defined UV_TIDES || defined POT_TIDES
                                 ! Number of tides
                                 ! ====== == =====
      integer, parameter :: Ntides = 8
#endif

!
!----------------------------------------------------------------------
! Point sources, Floast, Stations
!----------------------------------------------------------------------
!
#if defined PSOURCE || defined PSOURCE_MASS || defined PSOURCE_NCFILE
                                  ! Number of point sources
      integer, parameter :: Msrc = 30        ! ====== == ===== =======
#endif
#ifdef STATIONS
       integer, parameter :: NS = 5          ! Number of output stations
       integer, parameter :: Msta = 1000     ! Maximum number of stations
#endif


#ifdef SOLVE3D
!
!----------------------------------------------------------------------
! Number of tracers
!----------------------------------------------------------------------
!
!
# ifdef PASSIVE_TRACER
      integer, parameter :: ntrc_pas = 1
# else
      integer, parameter :: ntrc_pas = 0
# endif

/*! === SUBSTANCE ===*/
!
# if defined SUBSTANCE
! ntrc_subs : number of advected substances (not fixed, neither benthic)
!
      integer, parameter :: ntrc_subs = 0
      integer, parameter :: ntfix = 0
      integer, parameter :: ntrc_substot = ntrc_subs+ntfix
      integer, parameter :: itsubs1 = itemp+ntrc_salt+ntrc_pas+ntrc_bio+1
      integer, parameter :: itsubs2 = itemp+ntrc_salt+ntrc_pas+ntrc_bio+ntrc_subs
# else
      integer, parameter :: ntrc_subs = 0
      integer, parameter :: ntrc_substot = 0
# endif /* SUBSTANCE */

!
# ifdef SEDIMENT
! NSAND          Number of sand classes
! NMUD           Number of mud classes
! NGRAV          Number of gravel classes (not implemented...)
! NST            Number of sediment (tracer) size classes
!
      integer, parameter :: NSAND = 2
      integer, parameter :: NMUD = 0
      integer, parameter :: NGRAV = 0
      integer, parameter :: ntrc_sed = NSAND+NMUD+NGRAV
      integer, parameter :: NST = ntrc_sed
# else
      integer, parameter :: ntrc_sed = 0
# endif /* SEDIMENT */

!
!----------------------------------------------------------------------
! Number of sediment layer
!----------------------------------------------------------------------
!
# ifdef SEDIMENT
! NLAY           Number of layers in sediment bed
!
      integer, parameter :: NLAY = 1
# endif /* SEDIMENT */

# ifdef MUSTANG
! Vertical dimension (ksdmin:ksdmax) of variables in sediment
! (ksdmax : max number of layers)
!
      integer, parameter :: ksdmin = 1
      integer, parameter :: ksdmax = 10
# endif /* MUSTANG */

# if defined SEDIMENT && defined AGRIF
      integer, parameter :: Agrif_lev_sedim = 0
# endif

#endif /* SOLVE3D */


!
!----------------------------------------------------------------------
! Number of layers in Sediment (SL)
!----------------------------------------------------------------------
!
      integer, parameter :: N_sl = 0

!
!----------------------------------------------------------------------
! Max time increment for computing bottom stress at the 3D fast time
! steps
!----------------------------------------------------------------------
!
#ifdef BSTRESS_FAST
      integer, parameter :: inc_faststep_max = 10
#endif



!
!----------------------------------------------------------------------
! MPI related variables
!----------------------------------------------------------------------
!
      integer Lmmpi,Mmmpi,iminmpi,imaxmpi,jminmpi,jmaxmpi
      common /comm_setup_mpi1/ Lmmpi,Mmmpi
      common /comm_setup_mpi2/ iminmpi,imaxmpi,jminmpi,jmaxmpi
!
! Domain subdivision parameters
! ====== =========== ==========
!
! NPP            Maximum allowed number of parallel threads;
! NSUB_X,NSUB_E  Number of SHARED memory subdomains in XI- and
!                                                ETA-directions;
! NNODES        Total number of MPI processes (nodes);
! NP_XI,NP_ETA  Number of MPI subdomains in XI- and ETA-directions;
!
      integer NSUB_X, NSUB_E, NPP
#ifdef MPI
      integer NP_XI, NP_ETA, NNODES
#if defined(SPLITTING_X) && defined(SPLITTING_ETA)
      parameter (NP_XI=SPLITTING_X,  NP_ETA=SPLITTING_ETA,  NNODES=NP_XI*NP_ETA)
#else
      parameter (NP_XI=1,  NP_ETA=4,  NNODES=NP_XI*NP_ETA)
#endif
      parameter (NPP=1)
      parameter (NSUB_X=1, NSUB_E=1)
#ifdef OPENACC
      integer my_acc_device
      logical compute_on_device
      common/comm_my_device/my_acc_device,compute_on_device
#endif
#elif defined OPENMP
#if defined(SPLITTING_X) && defined(SPLITTING_ETA)
      parameter (NPP=SPLITTING_X*SPLITTING_ETA)
      parameter (NSUB_X=SPLITTING_X, NSUB_E=SPLITTING_ETA)
#else
      parameter (NPP=4)
      parameter (NSUB_X=1, NSUB_E=NPP)
#endif

#else
      parameter (NPP=1)
      parameter (NSUB_X=1, NSUB_E=NPP)
#ifdef OPENACC
      integer,parameter :: my_acc_device = 0
      logical,parameter :: compute_on_device = .true.
#endif
#endif



#include "param_dev.h"
