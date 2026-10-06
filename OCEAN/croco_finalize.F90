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
!  MODULE croco_finalize
!
!  Purpose : Shutdown sequence
!            driver: close netCDF files (root and, under AGRIF, every
!            child grid), finalize STOGEN/XIOS/OASIS/MPI, then stop.
!            Called once from croco.F90, whether the time-stepping
!            loop ran to completion, exited early (may_day_flag), or
!            never ran at all (ierr/=0 from croco_initialize).
!=======================================================================
!
      module croco_finalize

      implicit none
      private
      public :: croco_finalize_run

      contains

      subroutine croco_finalize_run(ierr)
!
      use scalars
#ifdef STOGEN
      use stomod, only : sto_mod_finalize
#endif
#ifdef XIOS
      use xios
#endif
#if defined OA_COUPLING || defined OW_COUPLING
      use mod_prism
#endif
!
      integer, intent(inout) :: ierr
#ifdef AGRIF
      type(Agrif_pgrid), pointer :: parcours
#endif
#ifdef MPI
      include 'mpif.h'
#endif
!
! provides ocean_grid_comm (the "world" communicator macro used below
! is redefined to it by cppdefs_dev.h under coupling)
#include "mpi_cpl.h"
!
!  ierr == 0 here means croco_initialize succeeded (whether or not
!  the time loop itself ran or exited early via may_day_flag): only
!  then is there anything open to close. On a "hard" croco_initialize
!  failure (ierr /= 0), skip straight to the MPI abort/finalize
!  epilogue below, which is what actually reacts to ierr/=0.
!
      if (ierr == 0) then
      call closecdf
#ifdef STOGEN
      call sto_mod_finalize
#endif
#ifdef XIOS
      call iom_context_finalize( "crocox")   ! needed for XIOS+AGRIF
#endif
#ifdef AGRIF
!
!  Close the netcdf files also for the child grids
!
      parcours=>Agrif_Curgrid%child_list % first
      do while (associated(parcours))
        call Agrif_Instance(parcours % gr)
        call closecdf
# ifdef XIOS
        call iom_context_finalize( "crocox") ! needed for XIOS+AGRIF
# endif
        parcours => parcours % next
      enddo
#endif
      endif ! ierr == 0

      if (may_day_flag /= 0) ierr=1

#ifdef MPI
      if (ierr /= 0) call mpi_abort (MPI_COMM_WORLD, ierr)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)  ! XIOS

# if defined XIOS
                                ! case XIOS + (OASIS / no OASIS)
                                !      > MPI finalize done by XIOS
      call xios_finalize()      !      > if OASIS, finalize is done by XIOS
                                !      > if AGRIF : done using iom_context_finalize

#  if !defined OA_COUPLING && !defined OW_COUPLING  && !defined AGRIF
       call MPI_Finalize (ierr)      !  if No coupling + No AGRIF : MPI_Finalize is needed
#  endif

# elif defined OA_COUPLING || defined OW_COUPLING
                                          ! case no XIOS + OASIS
      call prism_terminate_proto(ierr)    !   > Finalize OASIS3 (without XIOS)
# else
                                          ! case no XIOS + no OASIS
      call MPI_Finalize (ierr)            !   > Finalize CROCO (without XIOS)
# endif
#endif

      stop
      end subroutine croco_finalize_run

      end module croco_finalize
