MODULE plug_HYBIOSED

#include "cppdefs.h"

#if defined HYBIOSED
   !!===========================================================================
   !!                   ***  MODULE  plug_HYBIOSED  ***
   !!
   !! interface between CROCO and HYBIOSED module
   !!===========================================================================

   USE module_HYBIOSED

   USE init_HYBIOSED, ONLY: hbs_init
   USE main_HYBIOSED, ONLY: hbs_update

   IMPLICIT NONE

   PRIVATE

   PUBLIC hbs_update_main
   PUBLIC hbs_init_main

CONTAINS
!
!-------------------------------------------------------------------------------
   SUBROUTINE hbs_update_main(tile)

      INTEGER :: tile
# include "ocean2d.h"
# include "compute_tile_bounds.h"

      CALL hbs_update(Istr, Iend, Jstr, Jend)

   END SUBROUTINE

!-------------------------------------------------------------------------------
   SUBROUTINE hbs_init_main()

      REAL :: h0fond
# include "ocean2d.h"

# ifdef WET_DRY
      h0fond = D_wetdry
# else
      h0fond = 0.
# endif

      CALL hbs_init(h0fond)

   END SUBROUTINE

#endif /* HYBIOSED */

END MODULE plug_HYBIOSED

