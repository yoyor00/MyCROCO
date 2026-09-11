! Copyright (C) 2024-2026 IFREMER
! License: CeCILL-C
! See LICENSES/LICENSE_LAGRANGIAN_FOIL.txt

MODULE toolmpi

#include "cppdefs.h"
#include "toolcpp.h"
   !!======================================================================
   !!                   ***  MODULE toolmpi  ***
   !!
   !! NPROCS= total no of procs using  =>MYPROCS
   !! NPROCJ= number of division on J direction =>MYPROCJ
   !! NPROCI= number of division on I direction =>MYPROCI
   !! MPI_PROC_NULL= BLANK point, i.e. no mpi process is used for that local area
   !! NBBLANK=COUNT(PROC==IGROUDPOINT). i.e. total number of non-wet points
   !! IWETPOINT= wet point i.e. COUNT(PROC==IWETPOINT)==NPROCS
   !! NPROCI * NPROCJ=  NPROCS+NBBLANK
   !!
   !!======================================================================

#if defined MPI && (defined LAGRANGIAN)
   !! * Modules used
   USE comtraj, ONLY: rsh, rlg
   USE mpi
   USE module_lagrangian !parametres mpi
   IMPLICIT NONE

   PRIVATE

   !! * Accessibility
   PUBLIC :: ex_traj

   interface MPI_loc2glob
      module procedure MPI_loc2glob_real
      module procedure MPI_loc2glob_integer
   end interface MPI_loc2glob
   interface MPI_glob2loc
      module procedure MPI_glob2loc_real
      module procedure MPI_glob2loc_integer
   end interface MPI_glob2loc
   interface MPI_gather_sort_var
      module procedure MPI_gather_sort_int_var
      module procedure MPI_gather_sort_real4_var
      module procedure MPI_gather_sort_real8_var
   end interface MPI_gather_sort_var

   PUBLIC :: ADD_ALL_MPI_INT, MPI_SETUP_LAG, MPI_loc2glob, MPI_glob2loc
   PUBLIC :: ADD_ALL_MPI_REAL
   PUBLIC :: MPI_gather_sort_counts, MPI_gather_sort_perm, MPI_gather_sort_var
   INTEGER :: nprocs
   INTEGER :: miniproc, minjproc, maxiproc, maxjproc
   INTEGER :: myi_uproc, myj_uproc
   INTEGER :: myi_dproc, myj_dproc
   INTEGER :: myiproc, myjproc
   INTEGER, ALLOCATABLE, DIMENSION(:, :) :: PROC_NUM
   INTEGER, PARAMETER                  :: IGROUNDPOINT = -100
   INTEGER, PARAMETER                  :: IWETPOINT = -2
CONTAINS

   !!======================================================================

   SUBROUTINE init_mpi_type_size
      !&E---------------------------------------------------------------------
      !&E                 ***  ROUTINE  init_mpi_type_size  ***
      !&E
      !&E ** Purpose : Defines the MPI type for 'rsh' and 'rlg' kind real
      !&E              variables.
      !&E
      !&E       !  2014-12 (M. Honnorat) Refactor MPI routines
      !&E---------------------------------------------------------------------
      !! * Modules used
      USE comtraj, ONLY: rsh, rlg, type_mpi_rsh, type_mpi_rlg

      IMPLICIT NONE

      !! * Local declarations
      INTEGER :: ierr_mpi

      !!----------------------------------------------------------------------
      !! * Executable part
      IF (rsh == 4) THEN
         type_mpi_rsh = MPI_REAL
      ELSE IF (rsh == 8) THEN
         type_mpi_rsh = MPI_DOUBLE_PRECISION
      ELSE
         WRITE (*, *) 'MPI exchange works only for rsh=4 and 8, you gave rsh=', rsh
         CALL MPI_FINALIZE(ierr_mpi); STOP
      END IF

      IF (rlg == 4) THEN
         type_mpi_rlg = MPI_REAL
      ELSE IF (rlg == 8) THEN
         type_mpi_rlg = MPI_DOUBLE_PRECISION
      ELSE
         WRITE (*, *) 'MPI exchange works only for rlg=4 and 8, you gave rlg=', rlg
         CALL MPI_FINALIZE(ierr_mpi); STOP
      END IF

   END SUBROUTINE init_mpi_type_size

   !!======================================================================

   SUBROUTINE MPI_loc2glob_real(idx, idy)
      !&E---------------------------------------------------------------------
      !&E                 ***  ROUTINE DEFINE_UDPROC  ***
      !&E
      !&E ** Purpose : compute local indexes from global indexes
      !&E
      !&E---------------------------------------------------------------------
      IMPLICIT NONE
      REAL(kind=rsh), INTENT(inout) :: idx, idy

      idx = idx + REAL(iminmpi, kind=rsh) - 1.0_rsh
      idy = idy + REAL(jminmpi, kind=rsh) - 1.0_rsh
   END SUBROUTINE

   SUBROUTINE MPI_loc2glob_integer(idx, idy)
      !&E---------------------------------------------------------------------
      !&E                 ***  ROUTINE DEFINE_UDPROC  ***
      !&E
      !&E ** Purpose : compute local indexes from global indexes
      !&E
      !&E---------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(inout)           :: idx, idy

      idx = idx + iminmpi - 1
      idy = idy + jminmpi - 1
   END SUBROUTINE

   SUBROUTINE MPI_glob2loc_real(idx, idy)
      !&E---------------------------------------------------------------------
      !&E                 ***  ROUTINE MPI_glob2loc_real  ***
      !&E
      !&E ** Purpose : compute global indexes from local indexes
      !&E
      !&E---------------------------------------------------------------------
      IMPLICIT NONE
      REAL(kind=rsh), INTENT(inout)    :: idx, idy

      idx = idx - REAL(iminmpi, kind=rsh) + 1.0_rsh
      idy = idy - REAL(jminmpi, kind=rsh) + 1.0_rsh
   END SUBROUTINE

   SUBROUTINE MPI_glob2loc_integer(idx, idy)
      !&E---------------------------------------------------------------------
      !&E                 ***  ROUTINE MPI_glob2loc_integer  ***
      !&E
      !&E ** Purpose : compute global indexes from local indexes
      !&E
      !&E---------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(inout)           :: idx, idy

      idx = idx - iminmpi + 1
      idy = idy - jminmpi + 1
   END SUBROUTINE

   SUBROUTINE MPI_SETUP_LAG(ix, je)
      !&E---------------------------------------------------------------------
      !&E                 ***  ROUTINE MPI_SETUP_LAG  ***
      !&E
      !&E ** Purpose : define MPI setup relative to lagrangian
      !&E
      !&E---------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(in) :: ix, je
      !!----------------------------------------------------------------------
      !! * Executable part
      nprocs = nnodes
      miniproc = 0; minjproc = 0
      maxiproc = np_xi - 1
      maxjproc = np_eta - 1
      ALLOCATE (PROC_NUM(miniproc:maxiproc, minjproc:maxjproc))
      CALL define_proc_num(proc_num, miniproc, maxiproc, minjproc, maxjproc)
      CALL init_mpi_type_size()
      CALL DEFINE_UDPROC()

   END SUBROUTINE MPI_SETUP_LAG

   SUBROUTINE define_proc_num(proc, mini, maxi, minj, maxj)
      !&E---------------------------------------------------------------------
      !&E                 ***  ROUTINE DEFINE_NUM_PROC  ***
      !&E
      !&E ** Purpose : define proc number even in MPINOLAND case
      !&E
      !&E---------------------------------------------------------------------
      !! * Modules used
      use netcdf
      IMPLICIT NONE
      INTEGER, DIMENSION(mini:maxi, minj:maxj), INTENT(inout) :: PROC
      INTEGER, INTENT(in) :: mini, maxi, minj, maxj
      REAL(KIND=rlg), DIMENSION(0:LLm + 1, 0:MMm + 1)  ::   zmask
      INTEGER :: i, j, nerr, inu, inu2
      INTEGER :: lbx, ubx, lby, uby
      INTEGER :: ncid, varid, ierr
      INTEGER :: Istrmpi, Iendmpi, Jstrmpi, Jendmpi, i_X, j_E
      INTEGER :: chunk_size_X, margin_X, chunk_size_E, margin_E
# ifndef MP_3PTS
      INTEGER, parameter  :: Npts = 2
# else
      INTEGER, parameter  :: Npts = 3
# endif
      !!----------------------------------------------------------------------
      !! * Executable part

      PROC(:, :) = MPI_PROC_NULL

#ifndef MPI_NOLAND
      inu = 0
      do j = minj, maxj
         do i = mini, maxi
            PROC(i, j) = inu
            inu = inu + 1
         end do
      end do
#else
      ! read mask
      nerr = nf90_open('croco_grd.nc', NF90_NOWRITE, ncid)
      nerr = nf90_inq_varid(ncid, 'mask_rho', varid)
      nerr = nf90_get_var(ncid, varid, zmask)
      if (nerr /= nf90_noerr) then
         write (*, *) 'Reading mask file failed'
         CALL MPI_FINALIZE(ierr); STOP
      end if
      nerr = nf90_close(ncid)
      inu = 0; inu2 = 0
      do j = minj, maxj
         do i = mini, maxi
            !
            j_E = inu2/NP_XI
            i_X = inu2 - j_E*NP_XI
            !
            chunk_size_X = (LLm + NP_XI - 1)/NP_XI
            margin_X = (NP_XI*chunk_size_X - LLm)/2
            chunk_size_E = (MMm + NP_ETA - 1)/NP_ETA
            margin_E = (NP_ETA*chunk_size_E - MMm)/2
            !
            istrmpi = 1 + i_X*chunk_size_X - margin_X !-Npts
            iendmpi = istrmpi + chunk_size_X - 1  ! +Npts
            istrmpi = max(istrmpi, 1)
            iendmpi = min(iendmpi, LLm)
            !
            jstrmpi = 1 + j_E*chunk_size_E - margin_E !-NPTS
            jendmpi = jstrmpi + chunk_size_E - 1 !  +Npts
            jstrmpi = max(jstrmpi, 1)
            jendmpi = min(jendmpi, Mmm)
            !
            lbx = max(istrmpi - Npts, 1)
            ubx = min(iendmpi + Npts, LLm)
            lby = max(jstrmpi - Npts, 1)
            uby = min(jendmpi + Npts, Mmm)
            if (sum(zmask(lbx:ubx, lby:uby)) > 0.) then
               PROC(i, j) = inu
               inu = inu + 1
            end if
            inu2 = inu2 + 1
         end do
      end do
#endif
      ! i.e. each process have different value!!!
      DO j = minj, maxj
         DO i = mini, maxi
            IF (proc(i, j) .EQ. mynode) THEN
               myiproc = i
               myjproc = j
            END IF
         END DO
      END DO

   END SUBROUTINE define_proc_num

   SUBROUTINE DEFINE_UDPROC

      !&E---------------------------------------------------------------------
      !&E                 ***  ROUTINE DEFINE_UDPROC  ***
      !&E
      !&E ** Purpose : define neighbors proc numbers
      !&E
      !&E---------------------------------------------------------------------
      !! * Modules used
      IMPLICIT NONE

      !!----------------------------------------------------------------------
      !! * Executable part

      IF (myjproc == maxjproc) THEN
         myj_uproc = MPI_PROC_NULL
      ELSE IF (PROC_NUM(myiproc, myjproc + 1) .EQ. MPI_PROC_NULL) THEN
         myj_uproc = MPI_PROC_NULL
      ELSE
         myj_uproc = PROC_NUM(myiproc, myjproc + 1)
      END IF

      IF (myjproc == minjproc) THEN
         myj_dproc = MPI_PROC_NULL
      ELSE IF (PROC_NUM(myiproc, myjproc - 1) .EQ. MPI_PROC_NULL) THEN
         myj_dproc = MPI_PROC_NULL
      ELSE
         myj_dproc = PROC_NUM(myiproc, myjproc - 1)
      END IF

      IF (myiproc == maxiproc) THEN
         myi_uproc = MPI_PROC_NULL
      ELSE IF (PROC_NUM(myiproc + 1, myjproc) .EQ. MPI_PROC_NULL) THEN
         myi_uproc = MPI_PROC_NULL
      ELSE
         myi_uproc = PROC_NUM(myiproc + 1, myjproc)
      END IF

      IF (myiproc == miniproc) THEN
         myi_dproc = MPI_PROC_NULL
      ELSE IF (PROC_NUM(myiproc - 1, myjproc) .EQ. MPI_PROC_NULL) THEN
         myi_dproc = MPI_PROC_NULL
      ELSE
         myi_dproc = PROC_NUM(myiproc - 1, myjproc)
      END IF

      RETURN

   END SUBROUTINE DEFINE_UDPROC

   SUBROUTINE ADD_ALL_MPI_INT(value)
      !&E---------------------------------------------------------------------
      !&E                 ***  ROUTINE  ADD_ALL_MPI_INT  ***
      !&E
      !&E ** Purpose :  Compute global sum of an integer over all MPI processes.
      !&E
      !&E ** Description :
      !&E  Performs an MPI_ALLREDUCE with MPI_SUM on 'value'.
      !&E  On entry:  local integer contribution.
      !&E  On exit :  global sum available on all ranks.
      !&E
      !&E ** History :
      !&E---------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(inout) :: value

      !! * Local declarations
      INTEGER :: ierr_mpi
      INTEGER :: value_temp

      !CALL MPI_ALLREDUCE(MPI_IN_PLACE,value,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_ALLREDUCE(value, value_temp, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr_mpi)
      value = value_temp

   END SUBROUTINE ADD_ALL_MPI_INT

   !!======================================================================

   SUBROUTINE ADD_ALL_MPI_REAL(value)
      !&E---------------------------------------------------------------------
      !&E                 ***  ROUTINE  ADD_ALL_MPI_REAL  ***
      !&E
      !&E ** Purpose :  Compute global sum of a real over all MPI processes.
      !&E
      !&E ** Description :
      !&E  Performs an MPI_ALLREDUCE with MPI_SUM on 'value'.
      !&E  The MPI datatype (type_mpi_rlg) matches REAL(kind=rlg).
      !&E  On entry:  local real contribution.
      !&E  On exit :  global sum available on all ranks.
      !&E
      !&E ** History :
      !&E       !  2026-02 (C. Menu)  Adapted from ADD_ALL_MPI_INT
      !&E---------------------------------------------------------------------
      USE comtraj, ONLY: rlg, type_mpi_rlg
      IMPLICIT NONE
      REAL(kind=rlg), INTENT(inout) :: value

      !! * Local declarations
      INTEGER :: ierr_mpi
      REAL(kind=rlg) :: value_temp

      CALL MPI_ALLREDUCE(value, value_temp, 1, type_mpi_rlg, MPI_SUM, MPI_COMM_WORLD, ierr_mpi)
      value = value_temp
   END SUBROUTINE ADD_ALL_MPI_REAL

   !!======================================================================

   SUBROUTINE MPI_gather_sort_counts(nb_part_local, counts, displs, nb_part_total)
      !&E---------------------------------------------------------------------
      !&E                 ***  ROUTINE  MPI_gather_sort_counts  ***
      !&E
      !&E ** Purpose : Compute, on every rank, the per-rank particle counts,
      !&E              the Gatherv displacements and the total particle count,
      !&E              as a preliminary step to MPI_gather_sort_perm/MPI_gather_sort_var.
      !&E
      !&E---------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(in)                            :: nb_part_local
      INTEGER, DIMENSION(0:NNODES - 1), INTENT(out)  :: counts, displs
      INTEGER, INTENT(out)                           :: nb_part_total

      INTEGER :: ierr_mpi, i

      counts = 0
      CALL MPI_ALLGATHER(nb_part_local, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr_mpi)

      displs(0) = 0
      DO i = 1, NNODES - 1
         displs(i) = displs(i - 1) + counts(i - 1)
      END DO
      nb_part_total = SUM(counts)

   END SUBROUTINE MPI_gather_sort_counts

   !!======================================================================

   SUBROUTINE MPI_gather_sort_perm(nb_part_local, num_local, counts, displs, nb_part_total, iperm)
      !&E---------------------------------------------------------------------
      !&E                 ***  ROUTINE  MPI_gather_sort_perm  ***
      !&E
      !&E ** Purpose : Gather the local NUM values onto MASTER and build the
      !&E              permutation (valid on MASTER only) that orders the
      !&E              gathered particles by ascending NUM. counts/displs/
      !&E              nb_part_total must come from a prior call to
      !&E              MPI_gather_sort_counts.
      !&E
      !&E ** Description : particle%num is assumed to be a unique key >= 1
      !&E              (never reused), so gathered particles can be placed
      !&E              directly by index instead of being comparison-sorted.
      !&E              The upper bound on NUM is computed here as a collective
      !&E              MAX over every rank's local particles, rather than taken
      !&E              from the caller, since NUM's range is not necessarily
      !&E              the same as the current active particle count (e.g. some
      !&E              NUM values may never come active, e.g. when particles are
      !&E              numbered from their position in an input file that some
      !&E              particles get filtered out of).
      !&E
      !&E---------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(in)                             :: nb_part_local
      INTEGER, DIMENSION(nb_part_local), INTENT(in)   :: num_local
      INTEGER, DIMENSION(0:NNODES - 1), INTENT(in)    :: counts, displs
      INTEGER, INTENT(in)                             :: nb_part_total
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(out) :: iperm

      INTEGER, DIMENSION(:), ALLOCATABLE :: num_glob, pos_of_num
      INTEGER :: ierr_mpi, i, k, nb_part_max, local_max

      local_max = 0
      IF (nb_part_local > 0) local_max = MAXVAL(num_local)
      CALL MPI_ALLREDUCE(local_max, nb_part_max, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr_mpi)

      IF (MASTER) THEN
         ALLOCATE (num_glob(nb_part_total))
      ELSE
         ALLOCATE (num_glob(0))
      END IF

      CALL MPI_GATHERV(num_local, nb_part_local, MPI_INTEGER, &
                        num_glob, counts, displs, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr_mpi)

      IF (MASTER) THEN
         ALLOCATE (pos_of_num(nb_part_max))
         pos_of_num = 0
         DO i = 1, nb_part_total
            pos_of_num(num_glob(i)) = i
         END DO

         ALLOCATE (iperm(nb_part_total))
         k = 0
         DO i = 1, nb_part_max
            IF (pos_of_num(i) /= 0) THEN
               k = k + 1
               iperm(k) = pos_of_num(i)
            END IF
         END DO
         DEALLOCATE (pos_of_num)
      ELSE
         ALLOCATE (iperm(0))
      END IF

      DEALLOCATE (num_glob)

   END SUBROUTINE MPI_gather_sort_perm

   !!======================================================================

   SUBROUTINE MPI_gather_sort_int_var(nb_part_local, var_local, counts, displs, nb_part_total, iperm, var_glob)
      !&E---------------------------------------------------------------------
      !&E                 ***  ROUTINE  MPI_gather_sort_int_var  ***
      !&E
      !&E ** Purpose : Gather one INTEGER particle field onto MASTER and
      !&E              reorder it with the permutation from MPI_gather_sort_perm.
      !&E              var_glob is only meaningful (size nb_part_total) on
      !&E              MASTER; it is a zero-size array on every other rank.
      !&E
      !&E---------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(in)                             :: nb_part_local
      INTEGER, DIMENSION(nb_part_local), INTENT(in)   :: var_local
      INTEGER, DIMENSION(0:NNODES - 1), INTENT(in)    :: counts, displs
      INTEGER, INTENT(in)                             :: nb_part_total
      INTEGER, DIMENSION(:), INTENT(in)               :: iperm
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(out) :: var_glob

      INTEGER, DIMENSION(:), ALLOCATABLE :: var_raw
      INTEGER :: ierr_mpi, i

      IF (MASTER) THEN
         ALLOCATE (var_raw(nb_part_total), var_glob(nb_part_total))
      ELSE
         ALLOCATE (var_raw(0), var_glob(0))
      END IF

      CALL MPI_GATHERV(var_local, nb_part_local, MPI_INTEGER, &
                        var_raw, counts, displs, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr_mpi)

      IF (MASTER) THEN
         DO i = 1, nb_part_total
            var_glob(i) = var_raw(iperm(i))
         END DO
      END IF

      DEALLOCATE (var_raw)

   END SUBROUTINE MPI_gather_sort_int_var

   SUBROUTINE MPI_gather_sort_real4_var(nb_part_local, var_local, counts, displs, nb_part_total, iperm, var_glob)
      !&E ** Purpose : Same as MPI_gather_sort_int_var, for REAL(kind=4) fields.
      IMPLICIT NONE
      INTEGER, INTENT(in)                                  :: nb_part_local
      REAL(kind=4), DIMENSION(nb_part_local), INTENT(in)   :: var_local
      INTEGER, DIMENSION(0:NNODES - 1), INTENT(in)         :: counts, displs
      INTEGER, INTENT(in)                                  :: nb_part_total
      INTEGER, DIMENSION(:), INTENT(in)                    :: iperm
      REAL(kind=4), DIMENSION(:), ALLOCATABLE, INTENT(out) :: var_glob

      REAL(kind=4), DIMENSION(:), ALLOCATABLE :: var_raw
      INTEGER :: ierr_mpi, i

      IF (MASTER) THEN
         ALLOCATE (var_raw(nb_part_total), var_glob(nb_part_total))
      ELSE
         ALLOCATE (var_raw(0), var_glob(0))
      END IF

      CALL MPI_GATHERV(var_local, nb_part_local, MPI_REAL, &
                        var_raw, counts, displs, MPI_REAL, 0, MPI_COMM_WORLD, ierr_mpi)

      IF (MASTER) THEN
         DO i = 1, nb_part_total
            var_glob(i) = var_raw(iperm(i))
         END DO
      END IF

      DEALLOCATE (var_raw)

   END SUBROUTINE MPI_gather_sort_real4_var

   SUBROUTINE MPI_gather_sort_real8_var(nb_part_local, var_local, counts, displs, nb_part_total, iperm, var_glob)
      !&E ** Purpose : Same as MPI_gather_sort_int_var, for REAL(kind=8) fields.
      IMPLICIT NONE
      INTEGER, INTENT(in)                                  :: nb_part_local
      REAL(kind=8), DIMENSION(nb_part_local), INTENT(in)   :: var_local
      INTEGER, DIMENSION(0:NNODES - 1), INTENT(in)         :: counts, displs
      INTEGER, INTENT(in)                                  :: nb_part_total
      INTEGER, DIMENSION(:), INTENT(in)                    :: iperm
      REAL(kind=8), DIMENSION(:), ALLOCATABLE, INTENT(out) :: var_glob

      REAL(kind=8), DIMENSION(:), ALLOCATABLE :: var_raw
      INTEGER :: ierr_mpi, i

      IF (MASTER) THEN
         ALLOCATE (var_raw(nb_part_total), var_glob(nb_part_total))
      ELSE
         ALLOCATE (var_raw(0), var_glob(0))
      END IF

      CALL MPI_GATHERV(var_local, nb_part_local, MPI_DOUBLE_PRECISION, &
                        var_raw, counts, displs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr_mpi)

      IF (MASTER) THEN
         DO i = 1, nb_part_total
            var_glob(i) = var_raw(iperm(i))
         END DO
      END IF

      DEALLOCATE (var_raw)

   END SUBROUTINE MPI_gather_sort_real8_var

   !!======================================================================

   SUBROUTINE ex_traj(d_give, u_give, r_give, l_give)

      !&E---------------------------------------------------------------------
      !&E                 ***  ROUTINE ex_traj  ***
      !&E
      !&E ** Purpose :  Exchange trajectories
      !&E
      !&E ** Description : Exchange trajectories
      !&E
      !&E  When a trajectory is leaving the MPI domain, its variable 'limitbye' has
      !&E  previously been set to a value between 1 and 8, according to the border
      !&E  it touches (see 'traj_3d' subroutine in TRAJ/traject3d.F90).
      !&E  [ Basically, 1=SW, 2=S, 3=SE, 4=E, 5=NE, 6=N, 7=NW, 8=W ].
      !&E
      !&E  The arguments '[durl]_give' contain the number of outgoing particles for
      !&E  each side (down, up, right, left) of the current domain.
      !&E
      !&E  The actual transfert of particles is performed by 'ex_traj_1d', in
      !&E  four movements :
      !&E
      !&E        1:            2:            3:            4:
      !&E         +--v--+       +--^--+       +-----+       +-----+
      !&E         |     |       |     |      <|-   <|-     -|>   -|>
      !&E         +--v--+       +--^--+       +-----+       +-----+
      !&E
      !&E ** History :
      !&E       !  2014-12 (M. Honnorat) Adapt for IBM upgrade
      !&E---------------------------------------------------------------------
      !! * Modules used
      IMPLICIT NONE

      !! * Arguments
      INTEGER, INTENT(inout)   :: d_give, u_give
      INTEGER, INTENT(inout)   :: r_give, l_give

      !! * Local declarations
      INTEGER  :: d_get, u_get
      INTEGER  :: r_get, l_get

      !!----------------------------------------------------------------------
      !! * Executable part
      d_get = 0; u_get = 0
      l_get = 0; r_get = 0

      ! Send number of outgoing particles to neighbors and reveive number of incoming particles
      !  . first with up & down domains
      IF (myj_uproc /= MPI_PROC_NULL) CALL sendrecv_one_int(myj_uproc, u_give, u_get)
      IF (myj_dproc /= MPI_PROC_NULL) CALL sendrecv_one_int(myj_dproc, d_give, d_get)

      ! Carry out the exchange with up & down domains...
      !  . first communicate from the top down
      !  . then  communicate from the bottom up
      CALL ex_traj_1d(d_give, myj_dproc, u_get, myj_uproc, l_give, r_give, limit_min=1, limit_max=3)
      CALL ex_traj_1d(u_give, myj_uproc, d_get, myj_dproc, l_give, r_give, limit_min=5, limit_max=7)

      ! Send number of outgoing particles to neighbors and reveive number of incoming particles
      !  . then with left & right domains
      IF (myi_uproc /= MPI_PROC_NULL) CALL sendrecv_one_int(myi_uproc, r_give, r_get)
      IF (myi_dproc /= MPI_PROC_NULL) CALL sendrecv_one_int(myi_dproc, l_give, l_get)

      ! Carry out the exchange with left & right domains...
      !  . first communicate from right to left
      !  . then  communicate from left to right
      CALL ex_traj_1d(l_give, myi_dproc, r_get, myi_uproc, l_give, r_give, limit_min=8, limit_max=8)
      CALL ex_traj_1d(r_give, myi_uproc, l_get, myi_dproc, l_give, r_give, limit_min=4, limit_max=4)

   END SUBROUTINE ex_traj

   !!======================================================================

   SUBROUTINE ex_traj_1d(i_give, give_proc, i_get, get_proc, l_give, r_give, limit_min, limit_max)

      !&E---------------------------------------------------------------------
      !&E                 ***  ROUTINE ex_traj_1d  ***
      !&E
      !&E ** Purpose : Exchange trajectories in one direction
      !&E
      !&E ** Description : Actual exchange of particles across a boundary
      !&E
      !&E  This subroutine carries out the exchange of particle data between neighbotring MPI domains,
      !&E  in one direction, either from left to right, from right to left, from top to bottom or
      !&E  from bottom to top.
      !&E
      !&E    - 'i_give' is the number of particles to send to MPI neighbor 'give_proc'
      !&E    - 'i_get'  is the number of particles to recieve from MPI neighbor 'get_proc'
      !&E    - 'l_give' and 'r_give' (arguments from 'ex_traj') are updated in corner cases
      !&E    - 'limit_min' and 'limit_max' are the bounds (from 1 to 8) of the 'limitbye' value to consider.
      !&E
      !&E ** History :
      !&E       !  2014-12 (M. Honnorat) Adapt for IBM upgrade
      !&E---------------------------------------------------------------------
      !! * Modules used
      USE comtraj, ONLY: patches, get_patch, type_patch, type_particle, type_mpi_particle, enlarge_patch
      IMPLICIT NONE

      !! * Arguments
      INTEGER, INTENT(in)      :: i_give, give_proc
      INTEGER, INTENT(in)      :: i_get, get_proc
      INTEGER, INTENT(inout)   :: l_give, r_give
      INTEGER, INTENT(in)      :: limit_min, limit_max

      !! * Local declarations
      INTEGER                                          :: ierr_mpi
      INTEGER                                          :: l_give_proc, l_get_proc
      INTEGER                                          :: index_give, index_get
      INTEGER                                          :: i, n, m
      TYPE(type_particle), ALLOCATABLE, DIMENSION(:)   :: get, give
      INTEGER, ALLOCATABLE, DIMENSION(:)   :: get_n, give_n
      TYPE(type_patch), POINTER                     :: patch
      TYPE(type_particle), POINTER                     :: particle

      !!----------------------------------------------------------------------
      !! * Executable part

      IF ((i_give == 0) .and. (i_get == 0)) RETURN   ! Nothing to exchange

      l_give_proc = give_proc
      l_get_proc = get_proc
      ! mynode,l_give_proc,i_give,l_get_proc,i_get,limit_min,limit_max
      IF (i_get == 0) l_get_proc = MPI_PROC_NULL
      IF (i_give == 0) l_give_proc = MPI_PROC_NULL

      ALLOCATE (give(i_give), give_n(i_give), get(i_get), get_n(i_get))

      IF (i_give /= 0) THEN
         ! We have some particles to send.
         ! First, we'll pack them in array 'give'.
         i = 0

         loop_patches: DO n = 1, patches%nb
            patch => get_patch(patches, n)
            DO m = 1, patch%nb_part_alloc
               particle => patch%particles(m)
               ! Consider only particles leaving by the wright side
               IF ((particle%limitbye >= limit_min) .and. &
                   (particle%limitbye <= limit_max)) THEN
                  i = i + 1
                  give_n(i) = n           ! patch number
                  give(i) = particle    ! particle data
                  particle%active = .FALSE.   ! particle is no more active for the current domain.
                  particle%limitbye = 0
                  IF (i == i_give) exit loop_patches
               END IF
            END DO
         END DO loop_patches

      END IF   ! ( i_give /= 0 )

      index_give = 0
      index_get = 0

      ! index = 1XXXYYYZ where XXX is id of sender proc, YYY is id of reciever proc, Z is limit_min
      IF (l_give_proc /= MPI_PROC_NULL) index_give = 10000000 + 10000*mynode + 10*l_give_proc + limit_min
      IF (l_get_proc /= MPI_PROC_NULL) index_get = 10000000 + 10000*l_get_proc + 10*mynode + limit_min

      ! We can now exchange the data with the neighbors
      ! . first the patch numbers...
      CALL MPI_SENDRECV(give_n, i_give, MPI_INTEGER, l_give_proc, index_give, &
                        get_n, i_get, MPI_INTEGER, l_get_proc, index_get, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr_mpi)

      ! . then the  particle data...
      CALL MPI_SENDRECV(give, i_give, type_mpi_particle, l_give_proc, index_give, &
                        get, i_get, type_mpi_particle, l_get_proc, index_get, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr_mpi)

      IF (i_get /= 0) THEN
         ! Some particles are incoming :
         !  we can then unpack the data
         DO i = 1, i_get

            ! get the patch number
            n = get_n(i)     ! patch number
            patch => get_patch(patches, n)

            ! find a place somewhere in memory to store this new particle
            m = 0
            particle => NULL()
            DO WHILE (.NOT. ASSOCIATED(particle))
               m = m + 1
               IF (m > patch%nb_part_alloc) THEN
                  PRINT_DBG*, mynode, 'max number of particles reached. Increase patch % nb_part_alloc.', m, patch%nb_part_alloc
                  CALL enlarge_patch(patch, 1)
               END IF
               IF (.NOT. patch%particles(m)%active) THEN
                  particle => patch%particles(m)
               END IF
            END DO

            ! unpack particle data
            particle = get(i)
            ! update r_give, l_give and particle % limitbye (here we handle 'corner' cases).
            SELECT CASE (particle%limitbye)
            CASE (1)
               particle%limitbye = 8
               l_give = l_give + 1
            CASE (2)
               particle%limitbye = 0
            CASE (3)
               particle%limitbye = 4
               r_give = r_give + 1
            CASE (4)
               particle%limitbye = 0
            CASE (5)
               particle%limitbye = 4
               r_give = r_give + 1
            CASE (6)
               particle%limitbye = 0
            CASE (7)
               particle%limitbye = 8
               l_give = l_give + 1
            CASE (8)
               particle%limitbye = 0
            END SELECT
         END DO
      END IF   ! ( i_get /= 0 )

      DEALLOCATE (give, give_n, get, get_n)

   END SUBROUTINE ex_traj_1d

   SUBROUTINE sendrecv_one_int(dest, give, get)

      !&E---------------------------------------------------------------------
      !&E                 ***  ROUTINE sendrecv_one_int  ***
      !&E
      !&E ** Purpose : send and receive integer between procs
      !&E
      !&E---------------------------------------------------------------------
      !! * Modules used
      IMPLICIT NONE

      !! * Arguments
      INTEGER, INTENT(in)  :: dest
      INTEGER, INTENT(in)  :: give
      INTEGER, INTENT(out) :: get

      !! * Local declarations
      INTEGER :: ierr_mpi
      INTEGER :: send_tag, recv_tag

      !!----------------------------------------------------------------------
      !! * Executable part

      send_tag = 300*mynode + dest
      recv_tag = 300*dest + mynode

      CALL MPI_SENDRECV(give, 1, MPI_INTEGER, dest, send_tag, &
                        get, 1, MPI_INTEGER, dest, recv_tag, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr_mpi)

   END SUBROUTINE sendrecv_one_int

#endif

END MODULE toolmpi
