! Copyright (C) 2022-2026 IFREMER
! License: CeCILL-C
! See LICENSES/LICENSE_MUSTANG.txt

#include "cppdefs.h"

#if defined MUSTANG

      module plug_MUSTANG_CROCO

      USE module_substance
      USE croco_namelist, only: dt
# ifdef WET_DRY
      USE croco_namelist, only: D_wetdry
# endif
      USE initMUSTANG, ONLY : mustang_init
      USE initMUSTANG, ONLY : mustang_morphoinit
# ifdef MORPHODYN
      USE initMUSTANG, ONLY : mustang_morphoinit_dh
# endif
      USE sed_MUSTANG, ONLY : mustang_update
      USE sed_MUSTANG, ONLY : mustang_deposition
# ifdef MORPHODYN
      USE sed_MUSTANG, ONLY : mustang_morpho
# endif

      IMPLICIT NONE

      PRIVATE

      PUBLIC   mustang_update_main
      PUBLIC   mustang_deposition_main
      PUBLIC   mustang_init_main
# ifdef MORPHODYN
      PUBLIC   mustang_morpho_main
# endif

CONTAINS
!
!-----------------------------------------------------------------------
!
      subroutine mustang_update_main (tile)

      REAL    :: TEMPREF_LIN
      REAL    :: SALREF_LIN
      INTEGER :: tile
# include "compute_tile_bounds.h"

      TEMPREF_LIN = 10.0 
      SALREF_LIN  = 35.0
      CALL mustang_update (Istr, Iend, Jstr, Jend,  & 
                   t, zob,                          &
                   ubar, vbar,                      &
                   SALREF_LIN, TEMPREF_LIN, dt)
      end subroutine
!
!-----------------------------------------------------------------------
!
      subroutine mustang_deposition_main (tile)

      integer :: tile
# include "compute_tile_bounds.h"
      CALL mustang_deposition (Istr, Iend, Jstr, Jend, t)
      end subroutine
!
!-----------------------------------------------------------------------
!
      subroutine mustang_init_main (tile)

      REAL :: h0fond
      integer :: tile
# include "compute_tile_bounds.h"

# ifdef WET_DRY
            h0fond = D_wetdry
# else
            h0fond = 0.
# endif

      CALL mustang_init (Istr, Iend, Jstr, Jend,  &
                    h0fond, zob, t)

      ! call even if no morpho, initialization of hsed in all cases,
      ! even if initfromfile
      CALL mustang_morphoinit (Istr, Iend, Jstr, Jend)

# ifdef MORPHODYN
      CALL mustang_morphoinit_dh (Istr, Iend, Jstr, Jend, dh)
# endif
      end subroutine
!
!-----------------------------------------------------------------------
# ifdef MORPHODYN
      subroutine mustang_morpho_main (tile)

      integer :: tile
# include "compute_tile_bounds.h"

      CALL mustang_morpho (Istr, Iend, Jstr, Jend, dh)

      end subroutine
# endif
!-----------------------------------------------------------------------
!
      end module plug_mustang_CROCO

#else

      module plug_mustang_CROCO
      end module plug_mustang_CROCO
      
#endif /* MUSTANG */
