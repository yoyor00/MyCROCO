#include "cppdefs.h"

      subroutine das_innov_map (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call das_innov_map_tile (Istr,Iend,Jstr,Jend)
      return
      end 


      subroutine das_innov_map_tile (Istr,Iend,Jstr,Jend)
      implicit none
      integer Istr,Iend,Jstr,Jend
# include "param.h"
# include "das_param.h"
# include "scalars.h"
# include "grid.h"
# include "das_ocean.h"
# include "das_covar.h"
# include "das_innov.h"
#include "das_ocean_smooth.h"

      integer i,j,k,kdas,num
      real srad,srad0,rr,crs,crt,cfs,cft,dis, ro
      real tt(max_prf), ss(max_prf), count_t(max_prf),
     &     count_s(max_prf),rc_t(max_prf),rc_s(max_prf)
!
# include "compute_auxiliary_bounds.h"
!
!
      ro =  (lonr(Lm/2,Mm/2)-lonr(Lm/2+1,Mm/2+1))
     &     *(lonr(Lm/2,Mm/2)-lonr(Lm/2+1,Mm/2+1))
     &     +(latr(Lm/2,Mm/2)-latr(Lm/2+1,Mm/2+1))
     &     *(latr(Lm/2,Mm/2)-latr(Lm/2+1,Mm/2+1)) 

      srad=1.0*ro
      srad0=0.1*ro           ! srad, srad0 should be adjusted
                             ! larger, smoother, srad/srad0 < 2.5
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

          do num=1,prf_num_map
            rr= (lonr(i,j)-lon_map(num))                         
     &         *(lonr(i,j)-lon_map(num))                          
     &         +(latr(i,j)-lat_map(num))                        
     &         *(latr(i,j)-lat_map(num))                       
            if (rr .lt. srad ) then                            
              dis=exp(-rr/srad0)
              do k=1,max_prf
                                        ! bad value -9999
                if (t_map_raw(num,k) .gt. 2.0 ) then
                  count_t(k)=count_t(k)+1.0
                  rc_t(k)=rc_t(k)+dis
                  tt(k)=tt(k)+t_map_raw(num,k)*dis
                endif
                if (s_map_raw(num,k) .gt. 15.0 ) then
                  count_s(k)=count_s(k)+1.0
                  rc_s(k)=rc_s(k)+dis
                  ss(k)=ss(k)+s_map_raw(num,k)*dis           
                endif                                        
              enddo    !k
            endif      !rr                                                                                        
          enddo        !num                                                                                     

          do k=1,max_prf
            if(count_t(k) .gt. 0.5 ) then
              t_map(i,j,k)=tt(k)/rc_t(k)
              map_t_mask(i,j,k)=1.
            else
              t_map(i,j,k)=0.0
              map_t_mask(i,j,k)=0.0
            endif
            if(count_s(k) .gt. 0.5 ) then
              s_map(i,j,k)=ss(k)/rc_s(k)
              map_s_mask(i,j,k)=1.
            else
              s_map(i,j,k)=0.0
              map_s_mask(i,j,k)=0.0
            endif
          enddo   !k
        enddo     !i                                                                                              
      enddo       !j
!
! implement simple QC and compute innovations
!
      crt=6.0
      crs=0.8
!      crt=1.8
!      crs=0.18

      do j=JstrR,JendR
        do i=IstrR,IendR
          do k=1,max_prf
            kdas=NDAS-k+1
            cft=t_map(i,j,k) - t_sm(i,j,kdas,itemp)    
            map_t_mask(i,j,k)=map_t_mask(i,j,k)*rmask_das(i,j,kdas)
            if (abs(cft) .gt. crt ) map_t_mask(i,j,k)=0.0             
!                                                                                                                                           
            cfs=s_map(i,j,k) - t_sm(i,j,kdas,isalt)                 
            map_s_mask(i,j,k)=map_s_mask(i,j,k)*rmask_das(i,j,kdas)
            if (abs(cfs) .gt. crs ) map_s_mask(i,j,k)=0.0            
!                                                                                                                                           
            t_map(i,j,k)=cft*map_t_mask(i,j,k)                     
            s_map(i,j,k)=cfs*map_s_mask(i,j,k)                    
          enddo
        enddo
      enddo
 
!
99    continue
      return
      end
