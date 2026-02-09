#include "cppdefs.h"
#ifdef DAS_TPSSH

      subroutine das_innov_tp (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call das_innov_tp_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_innov_tp_tile (Istr,Iend,Jstr,Jend)
!
!--------------------------------------------------------------------
!  Set initial conditions for momentum and tracer variables using
!  analytical expressions.
!--------------------------------------------------------------------
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k, itrc
      real cft,crt
      real srad,srad0,rr,rc
      real tt(1200),dis(1200)
      integer count,num,n0
# include "param.h"
# include "das_param.h"
# include "scalars.h"
# include "grid.h"
# include "das_ocean.h"
# include "das_covar.h"
# include "das_innov.h"
# include "das_ocean_smooth.h"
!
!
# include "compute_auxiliary_bounds.h"
!
# ifdef DAS_TPSSH
      do j=JstrR,JendR
        do i=IstrR,IendR
          tp_oin(i,j)=tp_ossh
        enddo
      enddo
# endif
!
!     ... form superobs
!
      srad=1.5*0.0185*0.0185
      srad0=0.08*0.08
!
!     ... form superobs
!
      do j=JstrR,JendR
        do i=IstrR,IendR
          if (rmask_das(i,j,NDAS) .gt. 0.5) then
            count=0
            do num=1,num_tp
              rr= (lonr(i,j)-lon_tp(num))
     &           *(lonr(i,j)-lon_tp(num))
     &           +(latr(i,j)-lat_tp(num))
     &           *(latr(i,j)-lat_tp(num))
              if (rr .lt. srad ) then
                count=count+1
                if (count>1200) write(*,*) 'Issue tp innov'
                tt(count)=ssh_tp_raw(num)
                dis(count)=max(rr,srad0)
              endif
            enddo
            if (count .ge. 1 ) then
              rc=0.0
              ssh_tp(i,j)=0.0
              do n0=1,count
                rc=rc+1./dis(n0)
                ssh_tp(i,j)=ssh_tp(i,j)+tt(n0)/dis(n0)
              enddo
              ssh_tp(i,j)=ssh_tp(i,j)/rc
              tp_mask(i,j)=1.0
            else
              ssh_tp(i,j)=0.0
              tp_mask(i,j)=0.0
            endif
          else
            ssh_tp(i,j)=0.0
            tp_mask(i,j)=0.0
          endif
        enddo
      enddo

      crt=0.25
!
!     ... QC
!
      do j=JstrR,JendR
        do i=IstrR,IendR
          cft=ssh_tp(i,j)-(zeta_sm(i,j)-ssh_ref(i,j))
          if (abs(cft) .gt. crt ) tp_mask(i,j)=0.0
          tp_mask(i,j)=tp_mask(i,j)*ref_mask(i,j)
          ssh_tp(i,j)=cft * tp_mask(i,j)
        enddo
      enddo
!
#else
      subroutine das_innov_tp_empty
      return
      end
#endif
