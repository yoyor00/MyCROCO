#include "cppdefs.h"

#if defined OBSTRUCTION

      MODULE plug_OBSTRUCTIONS
      ! interface between croco and obstruction module

      USE module_OBSTRUCTIONS

      USE initOBSTRUCTIONS, ONLY : OBSTRUCTIONS_init
      USE OBSTRUCTIONS, ONLY : OBSTRUCTIONS_update

      IMPLICIT NONE

      PRIVATE

      PUBLIC   OBSTRUCTIONS_update_main
      PUBLIC   OBSTRUCTIONS_init_main

contains
!
!-----------------------------------------------------------------------
      SUBROUTINE OBSTRUCTIONS_update_main (tile)

      INTEGER :: tile
# include "ocean2d.h"
# include "compute_tile_bounds.h"

      CALL OBSTRUCTIONS_update(Istr, Iend, Jstr, Jend)

      END SUBROUTINE

!-----------------------------------------------------------------------
      SUBROUTINE OBSTRUCTIONS_init_main ()

      REAL :: h0fond
# include "ocean2d.h"
      
# ifdef WET_DRY
            h0fond = D_wetdry
# else
            h0fond = 0.
# endif

      CALL OBSTRUCTIONS_init(h0fond)
      END SUBROUTINE

!-----------------------------------------------------------------------
      END MODULE plug_OBSTRUCTIONS

#else

      MODULE plug_OBSTRUCTIONS_empty
      END MODULE plug_OBSTRUCTIONS_empty

#endif /* OBSTRUCTIONS */
