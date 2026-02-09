#include "cppdefs.h"

      subroutine das_zeta_smooth (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call das_zeta_smooth_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_zeta_smooth_tile (Istr,Iend,Jstr,Jend)
!
!--------------------------------------------------------------------
!  Set initial conditions for momentum and tracer variables using
!  analytical expressions.
!--------------------------------------------------------------------
!
      implicit none

# include "param.h"
# include "das_param.h"
# include "scalars.h"
# include "grid.h"
# include "das_ocean.h"
# include "das_covar.h"
# include "das_ocean_smooth.h"

      integer Istr,Iend,Jstr,Jend, i,j
      integer is,ie,js,je,i0,j0
      real rr,rc,ro
      integer count,num,n0
!
      integer srad,sdim
      parameter (srad=sz_rad_len, sdim=(2*srad+1)*(2*srad+1))
      real zz(sdim),dis(sdim)
!
!
# include "compute_auxiliary_bounds.h"
!
      ro = (sz_rad*sz_rad)*( 
     *           (lonr(Lm/2,Mm/2)-lonr(Lm/2+1,Mm/2+1))
     &          *(lonr(Lm/2,Mm/2)-lonr(Lm/2+1,Mm/2+1))
     &          +(latr(Lm/2,Mm/2)-latr(Lm/2+1,Mm/2+1))
     &          *(latr(Lm/2,Mm/2)-latr(Lm/2+1,Mm/2+1)) )
      
!   
      do j=JstrR,JendR
        do i=IstrR,IendR
          count=0
          js=max(j-srad, 0)
          je=min(j+srad, Mm+1)
          is=max(i-srad,0)
          ie=min(i+srad, Lm+1)
          do j0=js,je
            do i0=is,ie
              rr=      (lonr(i,j)-lonr(i0,j0))
     &                *(lonr(i,j)-lonr(i0,j0))
     &                +(latr(i,j)-latr(i0,j0))
     &                *(latr(i,j)-latr(i0,j0)) 
              count=count+1
              zz(count)=zeta_s(i0,j0)
              dis(count)=exp(-rr/ro)
            enddo
          enddo
          rc=0.0
          zeta_sm(i,j)=0.0
          do n0=1,count
            rc=rc+dis(n0)
            zeta_sm(i,j)=zeta_sm(i,j)+zz(n0)*dis(n0)
          enddo
          zeta_sm(i,j)=zeta_sm(i,j)/rc
        enddo
      enddo
!
      return
      end
