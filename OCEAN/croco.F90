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
      program croco
!
!======================================================================
!
!                     OCEAN MODEL MAIN DRIVER
!
!    Advances forward the equations for all nested grids, if any.
!
!    This is the slim replacement for the historical main.F: it owns
!    the MPI/XIOS/OASIS/AGRIF bootstrap, the time-stepping loop and
!    the shutdown sequence. The model initialization phase itself
!    (namelist read through the first history write) lives in
!    croco_init.F90 (module croco_init, subroutine croco_initialize).
!
!======================================================================
!
      use param
      use scalars
      use ncscrum
      use croco_init, only : croco_initialize
      use croco_finalize, only : croco_finalize_run
      use croco_namelist, only : ntimes, ninfo
#ifdef USE_CALENDAR
      use croco_namelist, only : end_date
#endif
!
      implicit none
!
      integer :: ierr
      integer :: iifroot, iicroot
#ifdef USE_CALENDAR
      character(len=19) :: tool_sectodat
      real(kind=8) :: tool_datosec
#endif
#ifdef AGRIF
      integer :: size_XI, size_ETA, se, sse, sz, ssz
      external :: step
#endif
#ifdef MPI
      include 'mpif.h'
#endif
!
#include "private_scratch.h"
#include "nbq.h"
#include "ocean2d.h"
! provides ocean_grid_comm (the "world" communicator macro used below
! is redefined to it by cppdefs_dev.h under coupling)
#include "mpi_cpl.h"
#ifdef AGRIF
#include "zoom.h"
#include "dynparam.h"
#endif
!
#include "dynderivparam.h"
!
!  Must stay directly in the PROGRAM unit: 
!  conv outlines the rest of the program into Sub_Loop_croco, whose 
!  very first call indexes into Agrif_tabvars_i/Agrif_tabvars_r ; 
!  arrays that Agrif_Init_Grids() itself allocates.
#ifdef AGRIF
      call Agrif_Init_Grids()
#endif
!
!----------------------------------------------------------------------
!  Model initialization: namelist read through the first history
!  write. See croco_init.F90 for the full sequence.
!----------------------------------------------------------------------
!
      call croco_initialize (ierr)
!
!  Only run the time loop if croco_initialize fully succeeded (ierr)
!  and did not itself flag a fatal problem (may_day_flag, e.g. a bad
!  grid file) -- either way, croco_finalize_run below still runs.
!
      if (ierr == 0 .and. may_day_flag == 0) then
!
!**********************************************************************
!                                                                     *
!             *****   ********   *****   ******  ********             *
!            **   **  *  **  *  *   **   **   ** *  **  *             *
!            **          **    **   **   **   **    **                *
!             *****      **    **   **   **   *     **                *
!                 **     **    *******   *****      **                *
!            **   **     **    **   **   **  **     **                *
!             *****      **    **   **   **   **    **                *
!                                                                     *
!**********************************************************************
!
!
      MPI_master_only write(stdout,'(/1x,A27/)')                       &
                      'MAIN: started time-stepping.'
      next_kstp=kstp
      time_start=time
#ifdef USE_CALENDAR
      time_end=tool_datosec(end_date)
#endif

#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif

#ifdef SOLVE3D
      iif = -1
      nbstep3d = 0
#endif
      iic = ntstart

#ifdef AGRIF
      iind = -1
# if defined OA_COUPLING || defined OW_COUPLING
      it_inside_root = 1
# endif
      grids_at_level = -1
      sortedint = -1
      call computenbmaxtimes
#endif

      time_stepping: do iicroot=ntstart,ntimes+1

#ifdef USE_CALENDAR
        if (mod(iicroot-1,ninfo) == 0) then
          MPI_master_only write(stdout,'(a)') tool_sectodat(time)
        endif
        if (time > time_end) exit time_stepping
#endif

#ifdef SOLVE3D
# ifndef AGRIF
        do iifroot = 0,nfast+2
# else
        nbtimes = 0
        do while (nbtimes <= nbmaxtimes)
# endif
#endif

#ifdef AGRIF
          call Agrif_Step(step)
#else
          call step()
#endif

#ifdef SOLVE3D
        enddo
#endif
        if (may_day_flag /= 0) exit time_stepping
      enddo time_stepping                  !-->  end of time step

      endif ! ierr == 0 .and. may_day_flag == 0

      call croco_finalize_run (ierr)
      end program croco
