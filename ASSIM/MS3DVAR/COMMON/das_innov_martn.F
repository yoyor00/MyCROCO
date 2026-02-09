#include "cppdefs.h"

      subroutine das_innov_martn (tile)
      implicit none
# include "param.h"
# include "das_param.h"
# include "scalars.h"
# include "grid.h"
# include "das_ocean.h"
# include "das_covar.h"
# include "das_innov.h"

      integer i,j,ii,jj,k,tile,kdas
      real srad,srad0,rmin,rr,crs,crt,cfs,cft
!
!
      if (tile .gt. prf_num_martn) goto 99
!
! 1  map to model grid
!
      srad=0.0185*0.0185
      srad0=0.008*0.008
! obs at > srad, not used
      rmin=srad
      ii=-999
      jj=-999
!
      do j=0,Mm+1
        do i=0,Lm+1
          rr = (lonr(i,j)-lon_martn(tile))
     &        *(lonr(i,j)-lon_martn(tile))
     &        +(latr(i,j)-lat_martn(tile))
     &        *(latr(i,j)-lat_martn(tile))
          if (rr .le. rmin) then
            rmin=rr
            ii=i
            jj=j
          endif
        enddo
      enddo
      Imartn(tile)=ii
      Jmartn(tile)=jj
!
! mask
!
      if (ii .ge. 0 ) then
        do k=1,max_prf
          if ( t_martn(tile,k) .gt. 0.0 .and.
     &        s_martn(tile,k) .gt. 0.0 ) then
            mask_martn(tile,k)=1.0
          else
            mask_martn(tile,k)=0.0
          endif
        enddo
      else
        do k=1,max_prf
          mask_martn(tile,k)=0.0
        enddo
        Imartn(tile)=0
        Jmartn(tile)=0
      endif
!
! rmask_das
!
! note: t_martn k=1:max_prf, k=1 => ndas
!
      do k=1,max_prf
        kdas=ndas-k+1
        mask_martn(tile,k)=mask_martn(tile,k)
     &        *rmask_das(Imartn(tile),Jmartn(tile),kdas)
      enddo
!
! implement simple QC and compute innovations
!
      crt=3.0
      crs=0.35
      do k=1,max_prf
        kdas=ndas-k+1
        cft=t_martn(tile,k)
     &       - t_das(Imartn(tile),Jmartn(tile),kdas,itemp)
        if (abs(cft) .gt. crt ) mask_martn(tile,k)=0.0
        t_martn(tile,k)=cft*mask_martn(tile,k)
!
        cfs=s_martn(tile,k)
     &       - t_das(Imartn(tile),Jmartn(tile),kdas,isalt)
        if (abs(cfs) .gt. crs ) mask_martn(tile,k)=0.0
        s_martn(tile,k)=cfs*mask_martn(tile,k)
      enddo
!
99    continue
      return
      end
