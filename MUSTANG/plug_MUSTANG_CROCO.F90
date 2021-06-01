#include "cppdefs.h"

#if defined MUSTANG && defined key_CROCO

      module plug_MUSTANG_CROCO

      USE module_MUSTANG
      USE sed_MUSTANG, ONLY : MUSTANG_update, MUSTANG_deposition
      USE sed_MUSTANG, ONLY : MUSTANG_morpho
      USE initMUSTANG, ONLY : MUSTANG_init_sediment

# include "coupler_define_MUSTANG.h"

      IMPLICIT NONE

      PRIVATE

      PUBLIC   mustang_update_main
      PUBLIC   mustang_deposition_main
      PUBLIC   mustang_init_sediment_main
      PUBLIC   mustang_morpho_main
      PUBLIC   mustang_update_dh

CONTAINS

!
!-----------------------------------------------------------------------
!
      subroutine mustang_update_main (tile)

      INTEGER :: tile
#include "ocean2d.h"
#include "compute_tile_bounds.h"
      CALL MUSTANG_update (Istr,Iend,Jstr,Jend,  & 
                   RESIDUAL_THICKNESS_WAT,       &
                   WATER_CONCENTRATION,Z0HYDRO,  &
                   WATER_ELEVATION,              &
#  if defined key_MUSTANG_lateralerosion || defined key_MUSTANG_bedload
                   BAROTROP_VELOCITY_U,          &
                   BAROTROP_VELOCITY_V,          &
#  endif
                   RHOREF,SALREF_LIN,TEMPREF_LIN,&
                   TRANSPORT_TIME_STEP)
      end subroutine

!
!-----------------------------------------------------------------------
!
      subroutine mustang_deposition_main (tile)
      integer :: tile
#include "ocean2d.h"
#include "compute_tile_bounds.h"
      CALL MUSTANG_deposition (Istr,Iend,Jstr,Jend, &
                   WATER_ELEVATION,                 &
                   RESIDUAL_THICKNESS_WAT,          &
                   WATER_CONCENTRATION)
      end subroutine

!
!-----------------------------------------------------------------------
!
      subroutine mustang_init_sediment_main (tile)
      integer :: tile
#include "ocean2d.h"
#include "compute_tile_bounds.h"
      CALL MUSTANG_init_sediment (Istr,Iend,Jstr,Jend,   &
                   0,WATER_ELEVATION,                    &
!                  0,BATHY_H0,WATER_ELEVATION,           &
# if (defined key_oasis && defined key_oasis_croco_ww3) || defined MORPHODYN
                   DHSED,                                &
# endif
!#if defined key_MUSTANG_debug && defined SPHERICAL
!     &     LATITUDE,LONGITUDE,
!#endif  
!#if defined key_MUSTANG_V2
!#if defined key_MUSTANG_bedload
!     &     CELL_DX,CELL_DY,
!#endif 
!#endif 
                   vname,indxT,ntrc_salt,                &
                   RESIDUAL_THICKNESS_WAT,Z0HYDRO,       &
                   WATER_CONCENTRATION)
      end subroutine

!
!-----------------------------------------------------------------------
!
      subroutine mustang_morpho_main (tile)
      integer :: tile
#include "ocean2d.h"
#include "compute_tile_bounds.h"

      CALL MUSTANG_morpho (Istr,Iend,Jstr,Jend,     &
                   WATER_ELEVATION,                 &
#      if (defined key_oasis && defined key_oasis_mars_ww3) || defined MORPHODYN
                   DHSED,                           &
#      endif                                     
                   RESIDUAL_THICKNESS_WAT)
      end subroutine

!
!-----------------------------------------------------------------------
!
      subroutine mustang_update_dh (tile)
      integer :: tile, i, j
#include "compute_tile_bounds.h"

#     ifdef MORPHODYN
       do j=J_EXT_RANGE
         do i=I_EXT_RANGE
           DHSED(i,j)=0.0_rsh
         enddo
       enddo 
#      undef I_EXT_RANGE
#      undef J_EXT_RANGE
#     endif
      end subroutine

!
!-----------------------------------------------------------------------
!


      end module plug_MUSTANG_CROCO
#else

      module plug_MUSTANG_CROCO_empty
      end module plug_MUSTANG_CROCO_empty
      
#endif /* MUSTANG */

