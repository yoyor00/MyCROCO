#include "cppdefs.h"

      subroutine das_innov_calauv (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call das_innov_calauv_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_innov_calauv_tile (Istr,Iend,Jstr,Jend)
!
!--------------------------------------------------------------------
!  Set initial conditions for momentum and tracer variables using
!  analytical expressions.
!--------------------------------------------------------------------
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k, itrc
      real cft, crt,cfs, crs
      real srad,srad0,rr,rc
      real tt(1200),ss(1200),dis(1200)
      common /temp_tt_ss_dis/ tt,ss,dis
      integer count,num,n0,kdas
# include "param.h"
# include "das_param.h"
# include "scalars.h"
# include "grid.h"
# include "das_ocean.h"
# include "das_covar.h"
# include "das_innov.h"
!
!
# include "compute_auxiliary_bounds.h"
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
      srad=2.0*0.0185*0.0185
      srad0=0.008*0.008
!
      do k=1,max_prf_cal
        kdas=NDAS-k+1
        do j=JR_RANGE
          do i=IR_RANGE
            if (rmask_das(i,j,kdas) .gt. 0.5) then
              count=0
              do num=1,num_cal(k)
                rr= (lonr(i,j)-lon_cal(num,k))
     &             *(lonr(i,j)-lon_cal(num,k))
     &             +(latr(i,j)-lat_cal(num,k))
     &             *(latr(i,j)-lat_cal(num,k))
                if (rr .lt. srad ) then
                  count=count+1
                  tt(count)=t_cal_raw(num,k)
                  ss(count)=s_cal_raw(num,k)
                  dis(count)= max(rr,srad0)
                endif
              enddo
              if (count .ge. 1 ) then
                rc=0.0
                t_cal(i,j,k)=0.0
                s_cal(i,j,k)=0.0
                do n0=1,count
                  rc=rc+1./dis(n0)
                  t_cal(i,j,k)=t_cal(i,j,k)+tt(n0)/dis(n0)
                  s_cal(i,j,k)=s_cal(i,j,k)+ss(n0)/dis(n0)
                enddo
                t_cal(i,j,k)=t_cal(i,j,k)/rc
                s_cal(i,j,k)=s_cal(i,j,k)/rc
                cal_mask(i,j,k)=1.0
              else
                t_cal(i,j,k)=25.0
                s_cal(i,j,k)=35.0
                cal_mask(i,j,k)=0.0
              endif
            else
              t_cal(i,j,k)=25.0
              s_cal(i,j,k)=35.0
              cal_mask(i,j,k)=0.0
            endif
          enddo
        enddo
!
        crt=2.5
        crs=0.3
!
        do j=JR_RANGE
          do i=IR_RANGE
            cft=t_cal(i,j,k)-t_das(i,j,kdas,itemp)
            cfs=s_cal(i,j,k)-t_das(i,j,kdas,isalt)
            if (abs(cft) .gt. crt ) cal_mask(i,j,k)=0.0
            if (abs(cfs) .gt. crs ) cal_mask(i,j,k)=0.0
            t_cal(i,j,k)=cft * cal_mask(i,j,k)
            s_cal(i,j,k)=cfs * cal_mask(i,j,k)
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
