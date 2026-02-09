#include "cppdefs.h"

      subroutine das_innov_uv_smooth (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call das_innov_uv_smooth_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_innov_uv_smooth_tile (Istr,Iend,Jstr,Jend)
!
!--------------------------------------------------------------------
! Purpose: 1 Fill the boundary grids, which are not determined by
!            DA
!          2 Smmoth the field to prevent the small scal noices due to
!            the early termination of the minimization 
!--------------------------------------------------------------------
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k
!
# include "param.h"
# include "das_param.h"
# include "scalars.h"
# include "grid.h"
# include "das_ocean.h"
# include "das_covar.h"
# include "das_ocean_smooth.h"
!
      integer is,ie,js,je,i0,j0
      real rr,rc,ro,ppm
      integer count,num,n0
!

      integer srad, sdim
      parameter (srad=6, sdim=(2*srad+1)*(2*srad+1))
      
      real dis(sdim),su(sdim),sv(sdim)
!
# include "compute_auxiliary_bounds.h"
!
!
# ifdef EW_PERIODIC
#  define IR_RANGE Istr,Iend
#  define IU_RANGE Istr,Iend
# else
#  define IR_RANGE IstrR,IendR
#  define IU_RANGE  Istr,IendR
# endif

# ifdef NS_PERIODIC
#  define JR_RANGE Jstr,Jend
#  define JV_RANGE Jstr,Jend
# else
#  define JR_RANGE JstrR,JendR
#  define JV_RANGE  Jstr,JendR
# endif

!
! NOTE: This smoothing is used to suppress noices
!       that may be introduced by psi/chi to u/v
!
      ro = (1.*1.)* (
     &           (lonr(Lm/2,Mm/2)-lonr(Lm/2+1,Mm/2+1))
     &          *(lonr(Lm/2,Mm/2)-lonr(Lm/2+1,Mm/2+1))
     &          +(latr(Lm/2,Mm/2)-latr(Lm/2+1,Mm/2+1))
     &          *(latr(Lm/2,Mm/2)-latr(Lm/2+1,Mm/2+1)))
!   
      do k=1,NDAS
        do j=JR_RANGE
          do i=IU_RANGE
            count=0
            js=max(j-srad, 1)
            je=min(j+srad, Mm)
            is=max(i-srad,1)
            ie=min(i+srad, Lm)
            do j0=js,je
              do i0=is,ie
                rr=      (lonr(i,j)-lonr(i0,j0))
     &                  *(lonr(i,j)-lonr(i0,j0))
     &                  +(latr(i,j)-latr(i0,j0))
     &                  *(latr(i,j)-latr(i0,j0)) 
                count=count+1
                su(count)=u_s(i0,j0,k)
                dis(count)=exp(-rr/ro)
              enddo
            enddo
            rc=0.0
            u_sm(i,j,k)=0.0
            do n0=1,count
              rc=rc+dis(n0)
              u_sm(i,j,k)=u_sm(i,j,k)+su(n0)*dis(n0)
            enddo
            u_sm(i,j,k)=u_sm(i,j,k)/rc
          enddo
        enddo
!
!v
!
        do j=JV_RANGE
          do i=IR_RANGE
            count=0
            js=max(j-srad, 1)
            je=min(j+srad, Mm)
            is=max(i-srad,1)
            ie=min(i+srad, Lm)
            do j0=js,je
              do i0=is,ie
                rr=      (lonr(i,j)-lonr(i0,j0))
     &                  *(lonr(i,j)-lonr(i0,j0))
     &                  +(latr(i,j)-latr(i0,j0))
     &                  *(latr(i,j)-latr(i0,j0)) 
                count=count+1
                sv(count)=v_s(i0,j0,k)
                dis(count)=exp(-rr/ro)
              enddo
            enddo
            rc=0.0
            v_sm(i,j,k)=0.0
            do n0=1,count
              rc=rc+dis(n0)                 
              v_sm(i,j,k)=v_sm(i,j,k)+sv(n0)*dis(n0)
            enddo
            v_sm(i,j,k)=v_sm(i,j,k)/rc
          enddo
        enddo
!
      enddo    !k
!
# undef IR_RANGE
# undef IU_RANGE
# undef JR_RANGE
# undef JV_RANGE
      return
      end
