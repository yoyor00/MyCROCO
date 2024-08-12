MODULE main_HYBIOSED

#include "cppdefs.h"

#if defined HYBIOSED
   !!===========================================================================
   !!                   ***  MODULE  HYBIOSED  ***
   !!
   !! ** Purpose : interactions and feedbacks between hydrodynamics, biology
   !!              and sediment dynamics
   !!
   !!              Strong dependancy with hydro model for mesh size,
   !!              position of u,v, name of variable, order of indexes (i,j,k),
   !!              MPI exchange...
   !!
   !! ** Description :
   !!
   !!     subroutine hbs_update  ! Main subroutine
   !!
   !!===========================================================================
   USE com_HYBIOSED
   USE HYBIOSED1DV, ONLY: hbs1dv_update_ws_coeff, hbs1dv_update_root_level

   IMPLICIT NONE

   PUBLIC hbs_update
   PRIVATE

CONTAINS
   !!===========================================================================
   SUBROUTINE hbs_update(limin, limax, ljmin, ljmax)
      !!------------------------------------------------------------------------
      !!                 *** SUBROUTINE hbs_update  ***
      !!
      !! ** Purpose : Update of ***TODO***
      !!------------------------------------------------------------------------
      USE module_HYBIOSED  ! needed for u, v, nstp, dt, time
# if defined OBSTRUCTION
      USE com_OBSTRUCTIONS, ONLY: obst_a3d, obst_s3d, obst_nbvar
# endif
      USE comMUSTANG, ONLY: ksmi, ksma, dzs
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: limin, limax, ljmin, ljmax

      REAL(KIND=rsh), DIMENSION(hbs_kmax)   :: hbs_uz
      REAL(KIND=rsh), DIMENSION(hbs_kmax)   :: hbs_vz
      REAL(KIND=rsh), DIMENSION(hbs_kmax)   :: hbs_s3d
      REAL(KIND=rsh), DIMENSION(hbs_kmax)   :: hbs_a3d
      REAL(KIND=rsh) :: hbs_hsed, hbs_dhsed
      CHARACTER(len=19) :: tool_sectodat
      CHARACTER(len=19)   :: cdate
      INTEGER :: i, j, k, ijour, imois, ian, iheure, iminu, isec, jjulien
      INTEGER  :: tool_julien

      DO j = ljmin, ljmax
         DO i = limin, limax
            ! **************************************
            ! * COMPUTES 3D VELOCITIES AT RHO POINT
            ! **************************************
            hbs_uz(:) = (u(i, j, :, nstp) + u(i + 1, j, :, nstp))/2.
            hbs_vz(:) = (v(i, j, :, nstp) + v(i, j + 1, :, nstp))/2.
# if defined OBSTRUCTION
            hbs_a3d(:) = obst_a3d(obst_nbvar + 3, :, i, j)
            hbs_s3d(:) = obst_s3d(obst_nbvar + 3, :, i, j)
# else
            hbs_a3d(:) = 0.
            hbs_s3d(:) = 0.
# endif
            ! *******************************
            ! * UPDATE CORRECTION FACTOR FOR
            ! * SETTLING VELOCITIES
            ! *******************************
            CALL hbs1dv_update_ws_coeff(hbs_uz(:), hbs_vz(:), &
                                        hbs_s3d(:), hbs_a3d(:), &
                                        hbs_ws_trapp(:, i, j), &
                                        hbs_ws_block(:, i, j))

            IF (hbs_nbvar > 0) THEN

               ! **************************************
               ! * COMPUTES sediment height variation
               ! **************************************
               hbs_hsed = 0.0_rsh
               DO k = ksmi(i, j), ksma(i, j)
                  hbs_hsed = hbs_hsed + dzs(k, i, j)
               END DO
               hbs_dhsed = hbs_hsed - hbs_hsed_prev(i, j)
               ! save hsed for next step
               hbs_hsed_prev(i, j) = hbs_hsed

               ! **************************************
               ! * COMPUTES day of year
               ! **************************************
               cdate = tool_sectodat(time)
               CALL tool_decompdate(cdate, ijour, imois, ian, iheure, iminu, isec)
               jjulien = tool_julien(ijour, imois, ian) - &
                         tool_julien(1, 1, ian) + 1

               ! ******************************************
               ! * UPDATE EROSION/DEPOSITION EFFECT ON ROOT
               ! ******************************************
               CALL hbs1dv_update_root_level(hbs_zup_root0(:, i, j), &
                                             hbs_thick_root0(:, i, j), &
                                             hbs_zup_root(:, i, j), &
                                             hbs_thick_root(:, i, j), &
                                             hbs_dz_root(:, i, j), &
                                             hbs_dthick_root(:, i, j), &
                                             hbs_dhsed, &
                                             hbs_position_wat(:, i, j), &
                                             hbs_position_bed(:, i, j), &
                                             dt, jjulien)
               ! TODO : check for debug to add
               ! save hbs_zup_root and hbs_thick_root for next step
               hbs_zup_root0(:, i, j) = hbs_zup_root(:, i, j)
               hbs_thick_root0(:, i, j) = hbs_thick_root(:, i, j)
            END IF ! hbs_nbvar > 0
         END DO ! i
      END DO ! j

   END SUBROUTINE hbs_update

#endif /* HYBIOSED */
END MODULE main_HYBIOSED
