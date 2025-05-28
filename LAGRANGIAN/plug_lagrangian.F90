#include "cppdefs.h"

#if defined LAGRANGIAN || defined DEB_IBM

MODULE plug_LAGRANGIAN
   ! interface between croco and lagrangian module

   USE module_lagrangian
   USE trajinitsave, ONLY: LAGRANGIAN_init 
   USE traject3d,ONLY : LAGRANGIAN_update
#ifdef MPI
   USE toolmpi,ONLY : MPI_SETUP_LAG
#endif

   IMPLICIT NONE

   PRIVATE

   PUBLIC LAGRANGIAN_update_main
   PUBLIC LAGRANGIAN_init_main
#ifdef MPI
   PUBLIC LAGRANGIAN_MPI_SETUP_main
#endif

contains
!
!-----------------------------------------------------------------------
   SUBROUTINE LAGRANGIAN_update_main(tile)

      INTEGER :: tile
# include "compute_tile_bounds.h"
      CALL LAGRANGIAN_update(zeta,u,v,Istr, Iend, Jstr, Jend)

   END SUBROUTINE

!-----------------------------------------------------------------------
   SUBROUTINE LAGRANGIAN_init_main(tile)

      INTEGER :: tile
# include "compute_tile_bounds.h"
      CALL LAGRANGIAN_init(Istr,Iend,Jstr,Jend)
   END SUBROUTINE

#ifdef MPI
!-----------------------------------------------------------------------
   SUBROUTINE LAGRANGIAN_MPI_SETUP_main(tile)

      INTEGER :: tile
# include "compute_tile_bounds.h"
   !init mpi params
   !----------------
      CALL MPI_SETUP_LAG(i_X,j_E)

   END SUBROUTINE
#endif

END MODULE plug_LAGRANGIAN

#else

MODULE plug_LAGRANGIAN_empty
END MODULE plug_LAGRANGIAN_empty

#endif /* LAGRANGIAN */
