#include "cppdefs.h"

      subroutine das_innov_flight (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call das_innov_flight_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_innov_flight_tile (Istr,Iend,Jstr,Jend)
!
!--------------------------------------------------------------------
!  Set initial conditions for momentum and tracer variables using
!  analytical expressions.
!--------------------------------------------------------------------
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k, itrc
      real otmi,cft, crt
      real srad,srad0,rr,rc
      real tt,dis
      integer count,num,n0
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
      srad=1.5*0.0185*0.0185
      srad0=0.008*0.008
!
      do j=JR_RANGE
        do i=IR_RANGE
          if (rmask_das(i,j,NDAS) .gt. 0.5) then
            count=0
            tt=0.
            dis=0.
            rc=0.
            do num=1,num_flight
              rr= (lonr(i,j)-lon_fl(num))
     &           *(lonr(i,j)-lon_fl(num))
     &           +(latr(i,j)-lat_fl(num))
     &           *(latr(i,j)-lat_fl(num))
              if (rr .lt. srad ) then
                count=count+1
                dis=max(rr,srad0)
                tt=tt+sst_fl_raw(num)/dis
                rc=rc+1.0/dis
              endif
            enddo
            if (count .ge. 1 ) then
              sst_fl(i,j)=tt/rc
              fl_mask(i,j)=1.0
            else
              sst_fl(i,j)=25.0
              fl_mask(i,j)=0.0
            endif
          else
            sst_fl(i,j)=25.0
            fl_mask(i,j)=0.0
          endif
        enddo
      enddo

      crt=6.0
!
      do j=JR_RANGE
        do i=IR_RANGE
          cft=sst_fl(i,j)-t_das(i,j,NDAS,itemp)
          if (abs(cft) .gt. crt ) fl_mask(i,j)=0.0
          sst_fl(i,j)=cft * fl_mask(i,j)
        enddo
      enddo
!
# undef IR_RANGE
# undef IU_RANGE
# undef JR_RANGE
# undef JV_RANGE
      return
      end
