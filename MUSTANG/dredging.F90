! MUSTANG - This software is governed by the CeCILL-C license
! see LICENSE_MUSTANG.txt

MODULE dredging
   !============================================================================
   !!                   ***  MODULE  dredging  ***
   !!
   !! This module manage dredging effect on sediment.
   !! This is done by removing sediments in areas where a sediment height is
   !! exceeded and optionnaly dumping them in other areas of the domain.
   !!
   !! This module is linked to MUSTANG and SUBSTANCE modules.
   !!
   !! This module has been adapted and generalized in 2024 by F. Grasso and
   !! S. Le Gac from the work of J.P. Lemoine in a MARS version for ARES
   !! hindcast in the SEINE
   !! estuary.
   !============================================================================

#include "cppdefs.h"

#if defined MUSTANG

   USE module_substance ! for GLOBAL_2D_ARRAY, h, t
   USE module_MUSTANG ! for z_w
   USE comsubstance, ONLY: lchain, rsh, rlg, riosh, nvp, isand1, surf_cell
   USE comMUSTANG, ONLY: dredging_location_file, dredging_settings_file, &
                         dredging_out_file, dredging_dt, &
                         hsed, ksmi, ksma, cv_sed, dzs

   IMPLICIT NONE
   PRIVATE

   PUBLIC dredging_init_param, dredging_init_hsed0 ! called by initMUSTANG
   PUBLIC dredging_main ! called by MUSTANG_update
   PUBLIC l_dredging ! boolean , TRUE to activate dredging

   ! Shared module variables
   LOGICAL :: l_dredging ! boolean , TRUE to activate dredging

   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: dredg_mass_byclass
   REAL(KIND=rsh), DIMENSION(:, :), ALLOCATABLE :: dredg_mass_byclass_byloc
   REAL(KIND=rsh), DIMENSION(:, :), ALLOCATABLE :: dump_mass_byclass_byloc
   REAL(KIND=rsh), DIMENSION(:, :), ALLOCATABLE :: dredg_mass_byclass_byloc_cum

   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: dump_surface

   INTEGER, DIMENSION(:, :, :), ALLOCATABLE :: dredg_loc
   INTEGER, DIMENSION(:, :, :), ALLOCATABLE :: dump_loc
   CHARACTER(lchain), DIMENSION(:), ALLOCATABLE :: dredg_name
   CHARACTER(lchain), DIMENSION(:), ALLOCATABLE :: dump_name
   INTEGER, DIMENSION(:, :), ALLOCATABLE :: dredg_flag
   INTEGER, DIMENSION(:), ALLOCATABLE :: dredg_dump_flag
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: dredg_depth
   REAL(KIND=rsh), DIMENSION(:, :), ALLOCATABLE :: dredg_hsed_init
   INTEGER :: n_dredg_loc
   INTEGER :: n_dump_loc

CONTAINS

   !============================================================================
   SUBROUTINE dredging_init_param(imin, imax, jmin, jmax)
      !!------------------------------------------------------------------------
      !!                 *** SUBROUTINE dredging_init_param  ***
      !!
      !! ** Purpose : Initialization of dredging parameters and variables
      !!
      !!------------------------------------------------------------------------
      INTEGER, INTENT(IN) :: imin, imax, jmin, jmax

      l_dredging = .TRUE.
      IF (TRIM(dredging_location_file) == "" .OR. &
          TRIM(dredging_settings_file) == "") THEN
         l_dredging = .FALSE.
      ELSE
         !! reading association dredging/dumping areas and dredging depths
         CALL dredging_get_dimension_from_settings_file
         IF (n_dredg_loc == 0) THEN
            l_dredging = .FALSE.
         ELSE
            !! allocate all arrays
            CALL dredging_alloc
            !! reading association dredging/dumping areas and dredging depths
            CALL dredging_read_settings_file
            !! reading dredging_location and dumping_location
            CALL dredging_read_location_file
            !! initialize dump_surface and checking compatibilities
            CALL dredging_init_var(imin, imax, jmin, jmax)
         END IF ! test if there is at least one dredging area

      END IF ! test if there is dreddging
      RETURN

   END SUBROUTINE dredging_init_param
   !============================================================================

   SUBROUTINE dredging_get_dimension_from_settings_file
      !!------------------------------------------------------------------------
      !!    *** SUBROUTINE dredging_get_dimension_from_settings_file ***
      !!
      !! ** Purpose : Get n_dredg_locand n_dump_loc values from file
      !!
      !!------------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER :: i, ios, idump
      CHARACTER(len=1000) :: line
      CHARACTER(len=lchain) :: temp_dredging, temp_dumping
      CHARACTER(len=lchain), DIMENSION(:), ALLOCATABLE :: temp_dumping_name
      REAL(kind=rsh) :: temp_depth
      LOGICAL :: new_dump_zone, line_ok

      OPEN (unit=51, file=dredging_settings_file, &
            status='old', action='read', iostat=ios)
      IF (ios /= 0) THEN
         WRITE (stdout, *) "Error opening dredging_settings_file : ", &
            dredging_settings_file
         STOP 1
      END IF

      READ (51, '(A)') ! reading header line
      ! retrieve n_dredg_loc
      n_dredg_loc = 0
      ios = 0
      line_ok = .TRUE.
      DO WHILE (line_ok)
         READ (51, '(A)', iostat=ios) line
         IF (ios == 0) THEN
            READ (line, *, iostat=ios) temp_dredging, temp_depth, temp_dumping
            IF (ios == 0) THEN
               ! If reading is successful (iostat is 0), count this line
               n_dredg_loc = n_dredg_loc + 1
            ELSE
               line_ok = .FALSE.
            END IF
         ELSE
            line_ok = .FALSE.
         END IF
      END DO

      REWIND (51)
      READ (51, '(A)') ! reading header line
      ! retrieve n_dump_loc
      n_dump_loc = 0
      ALLOCATE (temp_dumping_name(n_dredg_loc))
      temp_dumping_name(:) = "none"
      DO i = 1, n_dredg_loc
         READ (51, '(A)', iostat=ios) line
         READ (line, *, iostat=ios) temp_dredging, temp_depth, temp_dumping

         new_dump_zone = .true.
         IF (temp_dumping == "none") THEN
            new_dump_zone = .false.
         ELSE
            DO idump = 1, i - 1
               IF (TRIM(temp_dumping) == TRIM(temp_dumping_name(idump))) THEN
                  new_dump_zone = .false.
               END IF
            END DO
         END IF
         IF (new_dump_zone) THEN
            n_dump_loc = n_dump_loc + 1
         END IF
         temp_dumping_name(i) = temp_dumping
      END DO

      CLOSE (51)

      DEALLOCATE (temp_dumping_name)

   END SUBROUTINE dredging_get_dimension_from_settings_file
   !============================================================================

   SUBROUTINE dredging_read_settings_file

      INTEGER :: i, ios, iz, idump
      CHARACTER(len=1000) :: line
      CHARACTER(len=lchain) :: temp_dredging, temp_dumping
      REAL(kind=rsh) :: temp_depth
      LOGICAL :: new_dump_zone

      OPEN (unit=51, file=dredging_settings_file, &
            status='old', action='read', iostat=ios)
      IF (ios /= 0) THEN
         WRITE (stdout, *) "Error dredging_read_settings_file"
         WRITE (stdout, *) "Error opening dredging_settings_file : ", &
            dredging_settings_file
         STOP 1
      END IF

      READ (51, '(A)') ! reading header line

      idump = 0
      dump_name(:) = ""
      DO i = 1, n_dredg_loc
         READ (51, '(A)', iostat=ios) line
         READ (line, *, iostat=ios) temp_dredging, temp_depth, temp_dumping

         dredg_name(i) = TRIM(temp_dredging)
         dredg_depth(i) = temp_depth

         ! finding index of corresponding dumping zone
         IF (temp_dumping == "none") THEN
            dredg_dump_flag(i) = 0
         ELSE
            new_dump_zone = .true.
            DO iz = 1, idump
               IF (TRIM(temp_dumping) == TRIM(dump_name(iz))) THEN
                  new_dump_zone = .false.
                  dredg_dump_flag(i) = iz
               END IF
            END DO
            IF (new_dump_zone) THEN
               idump = idump + 1
               dump_name(idump) = TRIM(temp_dumping)
               dredg_dump_flag(i) = idump
            END IF

         END IF
      END DO
      CLOSE (51)

   END SUBROUTINE dredging_read_settings_file
   !============================================================================

   SUBROUTINE dredging_read_location_file
      !!------------------------------------------------------------------------
      !!                 *** SUBROUTINE dredging_read_location_file ***
      !!
      !! ** Purpose : Reading netcdf file containing dradging and dumping
      !!              locations
      !!
      !!------------------------------------------------------------------------
      USE netcdf
      IMPLICIT NONE

      INTEGER iz, ncid, varid, status
      CHARACTER(lchain) name
      REAL tmp(GLOBAL_2D_ARRAY)

      ! Open input NetCDF file.
      status = NF90_OPEN(dredging_location_file, NF90_NOWRITE, ncid)
      IF (status /= NF90_NOERR) THEN
         WRITE (stdout, *) "Error opening dredging location file"
         WRITE (stdout, *) TRIM(NF90_STRERROR(status))
         STOP 1
      END IF

      DO iz = 1, n_dredg_loc
         name = TRIM(dredg_name(iz))
         status = NF90_INQ_VARID(ncid, name, varid)
         IF (status /= NF90_NOERR) THEN
            WRITE (stdout, *) "Unable to retrieve dredg variable id: ", iz, " ", name
            WRITE (stdout, *) TRIM(NF90_STRERROR(status))
            STOP 1
         END IF
         CALL dredging_ncget2D(ncid, varid, tmp)
         dredg_loc(iz, GLOBAL_2D_ARRAY) = INT(tmp(GLOBAL_2D_ARRAY))
      END DO

      DO iz = 1, n_dump_loc
         name = TRIM(dump_name(iz))
         status = NF90_INQ_VARID(ncid, name, varid)
         IF (status /= NF90_NOERR) THEN
            WRITE (stdout, *) "Unable to retrieve dump variable id: ", iz, " ", name
            WRITE (stdout, *) TRIM(NF90_STRERROR(status))
            STOP 1
         END IF
         CALL dredging_ncget2D(ncid, varid, tmp)
         dump_loc(iz, GLOBAL_2D_ARRAY) = INT(tmp(GLOBAL_2D_ARRAY))
      END DO

      ! Close input NetCDF file.
      status = NF90_CLOSE(ncid)
      IF (status /= NF90_NOERR) THEN
         WRITE (stdout, *) "Error closing dredging location file"
         WRITE (stdout, *) TRIM(NF90_STRERROR(status))
         STOP 1
      END IF

      RETURN

   END SUBROUTINE dredging_read_location_file
   !============================================================================

   SUBROUTINE dredging_init_var(imin, imax, jmin, jmax)
      !! Initialize dredg_flag from all dredg_loc.
      !! In case of overlapping dredging zone, the deepest dredging is keeped
      INTEGER, INTENT(IN) :: imin, imax, jmin, jmax

      INTEGER :: i, j, iz, tmp_flag

      DO j = jmin, jmax
         DO i = imin, imax
            tmp_flag = 0
            DO iz = 1, n_dredg_loc
               IF (dredg_loc(iz, i, j) > 0) THEN
                  IF (tmp_flag == 0) THEN
                     tmp_flag = iz
                  ELSE IF (dredg_depth(iz) > dredg_depth(tmp_flag)) THEN
                     ! keeping the deepest
                     tmp_flag = iz
                  END IF
               END IF
            END DO
            dredg_flag(i, j) = tmp_flag

            DO iz = 1, n_dump_loc
               IF (dump_loc(iz, i, j) > 0) THEN
                  dump_surface(iz) = dump_surface(iz) + surf_cell(i, j)
               END IF
            END DO
         END DO !i
      END DO !j

   END SUBROUTINE dredging_init_var
   !============================================================================

   SUBROUTINE dredging_alloc

      ALLOCATE (dredg_mass_byclass(nvp))
      ALLOCATE (dredg_mass_byclass_byloc(n_dredg_loc, nvp))
      ALLOCATE (dump_mass_byclass_byloc(n_dump_loc, nvp))
      ALLOCATE (dredg_mass_byclass_byloc_cum(n_dredg_loc, nvp))
      ALLOCATE (dredg_loc(n_dredg_loc, GLOBAL_2D_ARRAY))
      ALLOCATE (dump_loc(n_dump_loc, GLOBAL_2D_ARRAY))
      ALLOCATE (dredg_name(n_dredg_loc))
      ALLOCATE (dump_name(n_dump_loc))
      ALLOCATE (dredg_depth(n_dredg_loc))
      ALLOCATE (dredg_dump_flag(n_dredg_loc))
      ALLOCATE (dredg_flag(GLOBAL_2D_ARRAY))
      ALLOCATE (dredg_hsed_init(GLOBAL_2D_ARRAY))
      ALLOCATE (dump_surface(n_dump_loc))

      dredg_mass_byclass(:) = 0.0_rsh
      dredg_mass_byclass_byloc(:, :) = 0.0_rsh
      dump_mass_byclass_byloc(:, :) = 0.0_rsh
      dredg_mass_byclass_byloc_cum(:, :) = 0.0_rsh

      dredg_loc(:, :, :) = 0
      dump_loc(:, :, :) = 0
      dredg_dump_flag(:) = 0
      dredg_flag(:, :) = 0
      dredg_depth(:) = 0.0_rsh
      dredg_hsed_init(:, :) = 0.0_rsh
      dump_surface(:) = 0.0_rsh

   END SUBROUTINE dredging_alloc
   !============================================================================

   SUBROUTINE dredging_init_hsed0

      dredg_hsed_init(:, :) = hsed(:, :)

   END SUBROUTINE dredging_init_hsed0
   !============================================================================

   SUBROUTINE dredging_main(imin, imax, jmin, jmax)
      INTEGER, INTENT(IN) :: imin, imax, jmin, jmax

      CALL dredging_compute_mass(imin, imax, jmin, jmax)

      CALL dumping_compute_mass

      dredg_mass_byclass_byloc_cum(:, :) = &
         dredg_mass_byclass_byloc_cum(:, :) + &
         dredg_mass_byclass_byloc(:, :)

      CALL dredging_mpi_mass

      CALL dredging_dump_mass

      CALL dredging_mpi_waterconcentration

      CALL dredging_output

   END SUBROUTINE dredging_main
   !============================================================================

   SUBROUTINE dredging_compute_mass(imin, imax, jmin, jmax)
      ! computation of dredged masses
      INTEGER, INTENT(IN) :: imin, imax, jmin, jmax

      INTEGER :: i, j, iv, iz, k

      dredg_mass_byclass_byloc(:, :) = 0.0_rsh
      DO j = jmin, jmax
         DO i = imin, imax
            hsed(i, j) = 0.0_rsh
            DO k = ksmi(i, j), ksma(i, j)
               hsed(i, j) = hsed(i, j) + dzs(k, i, j)
            END DO
            iz = dredg_flag(i, j)
            IF (iz > 0) THEN ! in dredging area
               DO WHILE ((h(i, j) - (hsed(i, j) - dredg_hsed_init(i, j))) &
                         .LE. (dredg_depth(iz)) .AND. ksma(i, j) > 0)
                  k = ksma(i, j)
                  DO iv = 1, nvp
                     dredg_mass_byclass_byloc(iz, iv) = &
                        dredg_mass_byclass_byloc(iz, iv) + &
                        dzs(k, i, j)*cv_sed(iv, k, i, j)*surf_cell(i, j)
                     cv_sed(iv, k, i, j) = 0.0_rsh
                  END DO
                  hsed(i, j) = hsed(i, j) - dzs(k, i, j)
                  dzs(k, i, j) = 0.0_rsh
                  ksma(i, j) = ksma(i, j) - 1
               END DO
            END IF
         END DO
      END DO

   END SUBROUTINE dredging_compute_mass
   !============================================================================

   SUBROUTINE dumping_compute_mass
      ! computation of dumping masses

      INTEGER :: iz, idump

      dump_mass_byclass_byloc(:, :) = 0.0_rsh
      DO iz = 1, n_dredg_loc
         DO idump = 1, n_dump_loc
            IF (dredg_dump_flag(iz) == idump) THEN
               dump_mass_byclass_byloc(idump, :) = &
                  dredg_mass_byclass_byloc(iz, :)
            END IF
         END DO
      END DO

   END SUBROUTINE dumping_compute_mass

   !============================================================================

   SUBROUTINE dredging_dump_mass
      ! transfer of dredged masses to dumping areas

      INTEGER :: i, j, idump, iv

      DO idump = 1, n_dump_loc
         IF (sum(dump_mass_byclass_byloc(:, idump)) > 0.0_rsh) THEN
            ! if there is mass to dump
            DO iv = isand1, nvp ! gravels are treated separately !! TODO add gravels in flx_ws_loc !!
               t(:, :, 1, nstp, itsubs1 - 1 + iv) = &
                  t(:, :, 1, nstp, itsubs1 - 1 + iv) + &
                  REAL(dump_loc(idump, :, :), rsh) &
                  *dump_mass_byclass_byloc(iv, idump) &
                  /dump_surface(idump)/(z_w(:, :, 1) - z_w(:, :, 0))
            END DO

         END IF
      END DO

   END SUBROUTINE dredging_dump_mass
   !============================================================================

   SUBROUTINE dredging_mpi_mass

      !! TODO

   END SUBROUTINE dredging_mpi_mass
   !============================================================================

   SUBROUTINE dredging_mpi_waterconcentration

      !! TODO

   END SUBROUTINE dredging_mpi_waterconcentration
   !============================================================================

   SUBROUTINE dredging_output

      !! TODO

   END SUBROUTINE dredging_output
   !============================================================================

   SUBROUTINE dredging_ncget2D(ncid, varid, tmp)
      !!------------------------------------------------------------------------
      !!                 *** SUBROUTINE dredging_ncget2D  ***
      !!
      !! ** Purpose : Purpose : get var from netcdf
      !!              (from nf_fread.F translate in f90)
      !!------------------------------------------------------------------------
      USE netcdf
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ncid, varid
      REAL, INTENT(INOUT) ::  tmp(GLOBAL_2D_ARRAY)

      INTEGER :: imin, imax, jmin, jmax, start(2), count(2), status

      jmin = 0
      imin = 0
      start(1) = 1
      start(2) = 1

#ifdef MPI
      if (ii .gt. 0) THEN
         start(1) = 1 - imin + iminmpi
         imin = 1
      end if
      if (ii .eq. NP_XI - 1) THEN
         imax = Lmmpi + 1
      else
         imax = Lmmpi
      end if
      if (jj .gt. 0) THEN
         start(2) = 1 - jmin + jminmpi
         jmin = 1
      end if
      if (jj .eq. NP_ETA - 1) THEN
         jmax = Mmmpi + 1
      else
         jmax = Mmmpi
      end if
#else
      imax = Lm + 1
      jmax = Mm + 1
#endif

      count(1) = imax - imin + 1
      count(2) = jmax - jmin + 1

      status = NF90_GET_VAR(ncid, varid, tmp(imin:imax, jmin:jmax), &
                            start, count)
      IF (status /= NF90_NOERR) THEN
         WRITE (stdout, *) "Error NF90_GET_VAR in dredging"
         WRITE (stdout, *) TRIM(NF90_STRERROR(status))
         STOP 1
      END IF

#ifdef MPI
# define LOCALLM Lmmpi
# define LOCALMM Mmmpi
#else
# define LOCALLM Lm
# define LOCALMM Mm
#endif
#if defined EW_PERIODIC || defined NS_PERIODIC  || defined MPI /* exchange needed */
      CALL exchange_r2d_tile(1, LOCALLM, 1, LOCALMM, tmp)
#endif /* MPI */

   END SUBROUTINE dredging_ncget2D

# endif /* ifdef MUSTANG */
   !============================================================================
END MODULE dredging
