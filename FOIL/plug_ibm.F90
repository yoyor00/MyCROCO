! Copyright (C) 2024-2026 IFREMER
! License: CeCILL-C
! See LICENSES/LICENSE_LAGRANGIAN_FOIL.txt

#include "cppdefs.h"

#ifdef FOIL

MODULE plug_ibm
   ! interface between croco and ibm module

   USE module_ibm
   USE ibm, ONLY: ibm_init, ibm_3d

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: ibm_init_main
   PUBLIC :: ibm_update_main

   ! ====================================================================
CONTAINS

   SUBROUTINE ibm_init_main(tile)

      INTEGER :: tile
# include "compute_tile_bounds.h"

      ! zeta(:,:,knew), NOT Zt_avg1: at init time (main.F) Zt_avg1 is still
      ! stale ("zeta=0" placeholder) until the next set_depth call, unlike
      ! per-timestep below where pre_step3d has just refreshed it.
      CALL ibm_init(zeta(:, :, knew), t(:, :, :, nstp, isalt), t(:, :, :, nstp, itemp), Istr, Iend, Jstr, Jend)

   END SUBROUTINE

   !======================================================================
   SUBROUTINE ibm_update_main(tile)

      INTEGER :: tile
# include "compute_tile_bounds.h"

      ! Zt_avg1 (fresh here, see ibm_init_main above) and u/v at nstp
      ! (same instant as `time`, see plug_lagrangian.F90), sliced once so
      ! ibm_3d and everything it calls only see already-time-correct
      ! 2D/3D fields, never a raw multi-slot array.
      CALL ibm_3d(Zt_avg1, u(:, :, :, nstp), v(:, :, :, nstp), &
                  t(:, :, :, nstp, isalt), t(:, :, :, nstp, itemp), Istr, Iend, Jstr, Jend)

   END SUBROUTINE

   !=========================================================================
END MODULE plug_ibm

#else

MODULE plug_ibm_empty
END MODULE plug_ibm_empty

#endif /* FOIL */
