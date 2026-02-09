#include "cppdefs.h"

      subroutine das_adj_ts_s_sm (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call das_adj_ts_s_sm_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_adj_ts_s_sm_tile (Istr,Iend,Jstr,Jend)
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
      real tt_adj(sdim),ss_adj(sdim)
!
      real :: temp_s_adj(GLOBAL_2D_ARRAY)
      real :: salt_s_adj(GLOBAL_2D_ARRAY)
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

        do j=1,sdim
          tt_adj(j)=0.0
          ss_adj(j)=0.0
        enddo
        do j=0,Mm+1
          do i=0,Lm+1
            temp_s_adj(i,j)=0.0
            salt_s_adj(i,j)=0.0
          enddo
        enddo

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
                dis(count)=exp(-rr/ro)
                !endif
              enddo
            enddo
            rc=0.0
            !BLXD added condition on count 
            !if ( count .gt. 0) then

            do n0=1,count
              rc=rc+dis(n0)
            enddo

           if (rc .gt. 0.) then
           t_sm_adj(i,j,k,itemp)=t_sm_adj(i,j,k,itemp)/rc
           t_sm_adj(i,j,k,isalt)=t_sm_adj(i,j,k,isalt)/rc

           do n0=1,count
             tt_adj(n0)=tt_adj(n0) + t_sm_adj(i,j,k,itemp)*dis(n0)
             ss_adj(n0)=ss_adj(n0) + t_sm_adj(i,j,k,isalt)*dis(n0)
           enddo

           t_sm_adj(i,j,k,itemp)=0.0
           t_sm_adj(i,j,k,isalt)=0.0

            count=0
            do j0=js,je
              do i0=is,ie
                count=count+1
                temp_s_adj(i0,j0)=temp_s_adj(i0,j0)+tt_adj(count)
                salt_s_adj(i0,j0)=salt_s_adj(i0,j0)+ss_adj(count)
                tt_adj(count)=0.0
                ss_adj(count)=0.0
              enddo
            enddo
            !BLXD added condition on count NO FILTERING
            else


            temp_s_adj(i,j)=t_sm_adj(i,j,k,itemp)
            salt_s_adj(i,j)=t_sm_adj(i,j,k,isalt)

            t_sm_adj(i,j,k,itemp)=0.0
            t_sm_adj(i,j,k,isalt)=0.0

            endif
!
          enddo
        enddo
!
C$OMP CRITICAL (adj_ts_s_cr)
        do j=0,Mm+1
          do i=0,Lm+1
            t_s_adj(i,j,k,itemp)=t_s_adj(i,j,k,itemp)+temp_s_adj(i,j)
            t_s_adj(i,j,k,isalt)=t_s_adj(i,j,k,isalt)+salt_s_adj(i,j)
          enddo
        enddo
C$OMP END CRITICAL (adj_ts_s_cr)
      enddo    !k
!
# undef IR_RANGE
# undef IU_RANGE
# undef JR_RANGE
# undef JV_RANGE
      return
      end
