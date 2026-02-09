#include "cppdefs.h"

      subroutine das_ts_s_sm (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call das_ts_s_sm_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_ts_s_sm_tile (Istr,Iend,Jstr,Jend)
!
!--------------------------------------------------------------------
! Purpose: 1 Fill the boundary grids, which are not determined by
!            DA
!          2 Smooth the field to prevent the small scale noice due to
!            the early termination of the minimization
!--------------------------------------------------------------------
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k
      integer is,ie,js,je,i0,j0
      real rr,rc,ro
      integer count,num,n0
!
# include "param.h"
# include "das_param.h"
# include "scalars.h"
# include "grid.h"
# include "das_ocean.h"
# include "das_covar.h"
# include "das_ocean_smooth.h"
!
      integer srad,sdim
      parameter (srad=sts_s_rad_len, sdim=(2*srad+1)*(2*srad+1))
      real tt(sdim),ss(sdim),dis(sdim)
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
      ro = (sts_s_rad*sts_s_rad)* ( 
     &           (lonr(Lm/2,Mm/2)-lonr(Lm/2+1,Mm/2+1))
     &          *(lonr(Lm/2,Mm/2)-lonr(Lm/2+1,Mm/2+1))
     &          +(latr(Lm/2,Mm/2)-latr(Lm/2+1,Mm/2+1))
     &          *(latr(Lm/2,Mm/2)-latr(Lm/2+1,Mm/2+1)) )
!   
      do k=1,NDAS
        do j=JR_RANGE
          do i=IR_RANGE
            count=0
            js=max(j-srad, 0)
            je=min(j+srad, Mm+1)
            is=max(i-srad,0)
            ie=min(i+srad, Lm+1)
            do j0=js,je
              do i0=is,ie
                !BLXD added condition on t/S inc
                !if (rmask_das(i0,j0,k) .gt. 0.5) then
                rr=      (lonr(i,j)-lonr(i0,j0))
     &                  *(lonr(i,j)-lonr(i0,j0))
     &                  +(latr(i,j)-latr(i0,j0))
     &                  *(latr(i,j)-latr(i0,j0)) 
                count=count+1
                tt(count)=t_s(i0,j0,k,itemp)
                ss(count)=t_s(i0,j0,k,isalt)
                dis(count)=exp(-rr/ro)
                !endif
              enddo
            enddo
            rc=0.0
            t_sm(i,j,k,itemp)=0.0
            t_sm(i,j,k,isalt)=0.0
            !BLXD added condition on count 
            !if ( count .gt. 0) then

            do n0=1,count
              rc=rc+dis(n0)
              t_sm(i,j,k,itemp)=t_sm(i,j,k,itemp)+tt(n0)*dis(n0)
              t_sm(i,j,k,isalt)=t_sm(i,j,k,isalt)+ss(n0)*dis(n0)
            enddo

            if ( rc .gt. 0.) then
            t_sm(i,j,k,itemp)=t_sm(i,j,k,itemp)/rc
            t_sm(i,j,k,isalt)=t_sm(i,j,k,isalt)/rc
            else
             write(*,*)'COUNT ISSUE with t/sinc smoothing ',i,j,k
             t_sm(i,j,k,itemp)=t_s(i,j,k,itemp)
             t_sm(i,j,k,isalt)=t_s(i,j,k,isalt)
            endif

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
