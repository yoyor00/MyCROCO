#include "cppdefs.h"


!======================================================================

# if defined DIAGNOSTICS_TRACER_ISO 

      subroutine set_diags_tracer_iso (tile)
      implicit none
      integer tile
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      call set_diags_tracer_iso_tile (istr,iend,jstr,jend,tile)
      return
      end

      subroutine set_diags_tracer_iso_tile (istr,iend,jstr,jend,tile)
!
      implicit none
# include "param.h"
      integer istr,iend,jstr,jend, i,j, ilc, iflux,
     &     imin,imax,jmin,jmax,tile
# ifdef SOLVE3D
     &                       , itrc, k
# endif

# include "scalars.h"
# include "ncscrum.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "diagnostics.h"

# include "compute_auxiliary_bounds.h"


! DEBUG: save fluxes for all terms

      do k=1,N
        do j=JstrR,JendR
          do i=IstrR,IendR

            ! we only save one tracer for now.
            itrc = NTA+1

            ! DEBUG: save fluxes
            TF(i,j,k,1) =  TF_xHmix(i,j,k,itrc)
            TF(i,j,k,2) =  TF_yHmix(i,j,k,itrc)
            TF(i,j,k,3) =  TF_zHmix(i,j,k,itrc)
            TF(i,j,k,4) =  TF_zVmix(i,j,k,itrc)
            TF(i,j,k,5) =  TF_Xadv(i,j,k,itrc)
            TF(i,j,k,6) =  TF_Yadv(i,j,k,itrc)
            TF(i,j,k,7) =  TF_Vadv(i,j,k,itrc)
            TF(i,j,k,8) =  0.

          enddo
        enddo
      enddo



! make the isopycnal projection for fluxes
!# include "compute_isopycnal_tracer_fluxes.h"



      return
      end
#else /* DIAGNOSTICS_TRACER_ISO  */
      subroutine set_diags_tracer_iso_empty
      end
#endif /* DIAGNOSTICS_TRACER_ISO  */










