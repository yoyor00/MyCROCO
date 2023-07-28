#include "cppdefs.h"

#if defined OBSTRUCTION

      module plug_OBSTRUCTIONS
        ! interface between croco and obstruction module

      USE module_OBSTRUCTIONS

      USE initOBSTRUCTIONS, ONLY : OBSTRUCTIONS_init, OBSTRUCTIONS_init_dimension
      USE OBSTRUCTIONS, ONLY : OBSTRUCTIONS_update

      IMPLICIT NONE

      PRIVATE

      PUBLIC   OBSTRUCTIONS_update_main
      PUBLIC   OBSTRUCTIONS_init_main

CONTAINS
!
!-----------------------------------------------------------------------
!
      subroutine OBSTRUCTIONS_update_main (tile)

      integer :: tile

# include "ocean2d.h"
# include "compute_tile_bounds.h"

      CALL OBSTRUCTIONS_update(Istr, Iend, Jstr, Jend, & 
                 cm0, h, zob,     &
                 Zt_avg1,         &
                 u(:,:,:,nstp),   &
                 v(:,:,:,nstp),   &
                 Hz               &
                 )

      end subroutine
!
!-----------------------------------------------------------------------
!
      subroutine OBSTRUCTIONS_init_main (tile)

      integer :: tile
      INTEGER :: imin, imax, jmin, jmax ! compute from GLOBAL_2DARRAY definition
      real :: h0fond
# include "ocean2d.h"
# include "compute_tile_bounds.h"
      
#ifdef THREE_GHOST_POINTS
        imin = -2
        imax = Lm+3+padd_X
        jmin = -2
        jmax = Mm+3+padd_E
# ifndef MPI
#  ifdef EW_PERIODIC
#   ifndef NS_PERIODIC
            jmin = 0
            jmax = Mm+1+padd_E
#   endif
#  else
#   ifdef NS_PERIODIC
            imin = 0
            imax = Lm+1+padd_X
#   else
            imin = 0
            imax = Lm+1+padd_X
            jmin = 0
            jmax = Mm+1+padd_E
#   endif
#  endif
# endif
#else
        imin = -1
        imax = Lm+2+padd_X
        jmin = -1
        jmax = Mm+2+padd_E
# ifndef MPI
#  ifdef EW_PERIODIC
#   ifndef NS_PERIODIC
            jmin = 0
            jmax = Mm+1+padd_E
#   endif
#  else
#   ifdef NS_PERIODIC
            imin = 0
            imax = Lm+1+padd_X
#   else
            imin = 0
            imax = Lm+1+padd_X
            jmin = 0
            jmax = Mm+1+padd_E
#   endif
#  endif
# endif
#endif

# ifdef WET_DRY
            h0fond = D_wetdry
# else
            h0fond = 0.
# endif

      CALL OBSTRUCTIONS_init_dimension(imin, imax, jmin, jmax, N, stdout)
      CALL OBSTRUCTIONS_init(zob, h0fond)
      end subroutine
!-----------------------------------------------------------------------
!
      end module plug_OBSTRUCTIONS

#else

      module plug_OBSTRUCTIONS_empty
      end module plug_OBSTRUCTIONS_empty
      
#endif /* OBSTRUCTIONS */
