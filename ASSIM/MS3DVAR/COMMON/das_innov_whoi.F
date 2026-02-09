#include "cppdefs.h"

      subroutine das_innov_whoi (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call das_innov_whoi_tile (Istr,Iend,Jstr,Jend)
      return
      end 


      subroutine das_innov_whoi_tile (Istr,Iend,Jstr,Jend)
      implicit none
      integer Istr,Iend,Jstr,Jend
# include "param.h"
# include "das_param.h"
# include "scalars.h"
# include "grid.h"
# include "das_ocean.h"
# include "das_covar.h"
# include "das_innov.h"

      integer i,j,k,kdas,num
      real srad,srad0,rr,crs,crt,cfs,cft,dis,ro
      real tt(max_prf), ss(max_prf), count_t(max_prf),
     &     count_s(max_prf),rc_t(max_prf),rc_s(max_prf)
!
# include "compute_auxiliary_bounds.h"
!
!
!
      ro =      (lonr(Lm/2,Mm/2)-lonr(Lm/2+1,Mm/2+1))
     &          *(lonr(Lm/2,Mm/2)-lonr(Lm/2+1,Mm/2+1))
     &          +(latr(Lm/2,Mm/2)-latr(Lm/2+1,Mm/2+1))
     &          *(latr(Lm/2,Mm/2)-latr(Lm/2+1,Mm/2+1))

      srad=9.0*ro
      srad0=0.5*ro     ! larger, smoother, 0.7*0.7
!
!WHOI
!
      do j=JstrR,JendR
        do i=IstrR,IendR
          do k=1,max_prf
            count_t(k)=0
            count_s(k)=0
            tt(k)=0.
            ss(k)=0
            rc_t(k)=0.
            rc_s(k)=0.
          enddo        !k
          dis=0.

          do num=1,prf_num_whoi
            rr= (lonr(i,j)-lon_whoi(num))                                                                  
     &         *(lonr(i,j)-lon_whoi(num))                                                                  
     &         +(latr(i,j)-lat_whoi(num))                                                                  
     &         *(latr(i,j)-lat_whoi(num))                                                                  
            if (rr .lt. srad ) then                                                                         
              dis=max(rr,srad0)
              do k=1,max_prf
                                        ! bad value -9999
                if (t_whoi_raw(num,k) .gt. 1.0 ) then
                  count_t(k)=count_t(k)+1.0
                  rc_t(k)=rc_t(k)+1.0/dis
                  tt(k)=tt(k)+t_whoi_raw(num,k)/dis
                endif
                if (s_whoi_raw(num,k) .gt. 15.0 ) then
                  count_s(k)=count_s(k)+1.0
                  rc_s(k)=rc_s(k)+1.0/dis
                  ss(k)=ss(k)+s_whoi_raw(num,k)/dis 
                endif                     
              enddo    !k
            endif      !rr               
          enddo        !num                                                                                     

          do k=1,max_prf
            if(count_t(k) .gt. 0.5 ) then
              t_whoi(i,j,k)=tt(k)/rc_t(k)
              whoi_t_mask(i,j,k)=1.
            else
              t_whoi(i,j,k)=0.0
              whoi_t_mask(i,j,k)=0.0
            endif
            if(count_s(k) .gt. 0.5 ) then
              s_whoi(i,j,k)=ss(k)/rc_s(k)
              whoi_s_mask(i,j,k)=1.
            else
              s_whoi(i,j,k)=0.0
              whoi_s_mask(i,j,k)=0.0
            endif
          enddo   !k
        enddo     !i                                                                                              
      enddo       !j
!
! implement simple QC and compute innovations
!
      crt=3.5
      crs=0.35

      do j=JstrR,JendR
        do i=IstrR,IendR
          do k=1,max_prf
            kdas=NDAS-k+1
            cft=t_whoi(i,j,k) - t_das(i,j,kdas,itemp)    
            whoi_t_mask(i,j,k)=whoi_t_mask(i,j,k)*rmask_das(i,j,kdas)
            if (abs(cft) .gt. crt ) whoi_t_mask(i,j,k)=0.0 
!
            cfs=s_whoi(i,j,k) - t_das(i,j,kdas,isalt)     
            whoi_s_mask(i,j,k)=whoi_s_mask(i,j,k)*rmask_das(i,j,kdas)
            if (abs(cfs) .gt. crs ) whoi_s_mask(i,j,k)=0.0 
!                                                        
            t_whoi(i,j,k)=cft*whoi_t_mask(i,j,k)          
            s_whoi(i,j,k)=cfs*whoi_s_mask(i,j,k)         
          enddo
        enddo
      enddo
!
      return
      end
