#include "cppdefs.h"
! TODO-MS in da_l0_lr the subroutine isn't protected with DAS_JASONSSH
#ifdef DAS_JASONSSH
      subroutine das_innov_js1 (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call das_innov_js1_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_innov_js1_tile (Istr,Iend,Jstr,Jend)
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
# ifdef DAS_JASONSSH
      do j=JstrR,JendR
        do i=IstrR,IendR
          js1_oin(i,j)=js1_ossh
        enddo
      enddo
# endif
!
! TODO-MS BLXD removed it's TP
!# ifdef DAS_TPSSH
!      do j=JstrR,JendR
!        do i=IstrR,IendR
!          tp_oin(i,j)=js1_ossh
!        enddo
!      enddo
!# endif
!
!     ... form superobs
!
      srad=1.5*0.0185*0.0185
      srad0=0.008*0.008
!
!     ... form superobs
!
      do j=JstrR,JendR
        do i=IstrR,IendR
          if (rmask_das(i,j,NDAS) .gt. 0.5) then
            count=0
            do num=1,num_js1
              rr= (lonr(i,j)-lon_js1(num))
     &           *(lonr(i,j)-lon_js1(num))
     &           +(latr(i,j)-lat_js1(num))
     &           *(latr(i,j)-lat_js1(num))
              if (rr .lt. srad ) then
                count=count+1
                if (count>1200) write(*,*) 'WARNING : count Issue js1 innov'
                tt(count)=ssh_js1_raw(num)
                dis(count)=max(rr,srad0)
              endif
            enddo
            if (count .ge. 1 ) then
              rc=0.0
              ssh_js1(i,j)=0.0
              do n0=1,count
                rc=rc+1./dis(n0)
                ssh_js1(i,j)=ssh_js1(i,j)+tt(n0)/dis(n0)
              enddo
              ssh_js1(i,j)=ssh_js1(i,j)/rc
              js1_mask(i,j)=1.0
            else
              ssh_js1(i,j)=0.0
              js1_mask(i,j)=0.0
            endif
          else
            ssh_js1(i,j)=0.0
            js1_mask(i,j)=0.0
          endif
        enddo
      enddo

      crt=0.20
!
!     ... QC
!
      do j=JstrR,JendR
        do i=IstrR,IendR
          cft=ssh_js1(i,j)-(zeta_das(i,j)-ssh_ref(i,j))
          if (abs(cft) .gt. crt ) js1_mask(i,j)=0.0
          js1_mask(i,j)=js1_mask(i,j)*ref_mask(i,j)
          ssh_js1(i,j)=cft * js1_mask(i,j)
        enddo
      enddo
!
      return
      end
!
#else
      subroutine das_innov_js1_empty
      return
      end
#endif
