#include "cppdefs.h"

      subroutine das_pressure_smooth (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call das_pressure_smooth_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_pressure_smooth_tile (Istr,Iend,Jstr,Jend)
!
!--------------------------------------------------------------------
! Purpose: 1 Fill the boundary grids, which are not determined by
!            DA
!          2 Smooth the field to prevent the small scale noice due to
!            the early termination of the minimization
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

      integer Istr,Iend,Jstr,Jend, i,j,k
      integer is,ie,js,je,i0,j0
      real rr,rc,ro
      integer count,num,n0
!
      integer srad,sdim
      parameter (srad=sz_rad_len, sdim=(2*srad+1)*(2*srad+1))
      real pp(sdim),dis(sdim)
!
# include "compute_auxiliary_bounds.h"
!
      ro = (sz_rad*sz_rad)* ( 
     &           (lonr(Lm/2,Mm/2)-lonr(Lm/2+1,Mm/2+1))
     &          *(lonr(Lm/2,Mm/2)-lonr(Lm/2+1,Mm/2+1))
     &          +(latr(Lm/2,Mm/2)-latr(Lm/2+1,Mm/2+1))
     &          *(latr(Lm/2,Mm/2)-latr(Lm/2+1,Mm/2+1)) )
!   
      do k=1,NDAS
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
     &                  *(lonr(i,j)-lonr(i0,j0))
     &                  +(latr(i,j)-latr(i0,j0))
     &                  *(latr(i,j)-latr(i0,j0)) 
                count=count+1
                pp(count)=p_s(i0,j0,k)
                dis(count)=exp(-rr/ro)
              enddo
            enddo
            rc=0.0
            p_sm(i,j,k)=0.0
            do n0=1,count
              rc=rc+dis(n0)
              p_sm(i,j,k)=p_sm(i,j,k)+pp(n0)*dis(n0)
            enddo
            p_sm(i,j,k)=p_sm(i,j,k)/rc
          enddo
        enddo
!
      enddo    !k

!
      return
      end
