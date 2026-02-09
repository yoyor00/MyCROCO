#include "cppdefs.h"

      subroutine das_innov_TS_obs (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call das_innov_TS_obs_tile (Istr,Iend,Jstr,Jend)
      return
      end 


      subroutine das_innov_TS_obs_tile (Istr,Iend,Jstr,Jend)
      implicit none
      integer Istr,Iend,Jstr,Jend
# include "param.h"
# include "das_param.h"
# include "scalars.h"
# include "grid.h"
# include "das_ocean.h"
# include "das_covar.h"
# include "das_innov.h"

      integer i,j,k,num
      real crs,crt,cfs,cft
      real teq, seq, t_aup,t_adn,s_aup,s_adn
!
# include "compute_auxiliary_bounds.h"
!
!
!
901   format(3I8,7F11.4)

      do j=JstrR,JendR
        do i=IstrR,IendR
          do k=1,NDAS-1
!
            obs_t_mask(i,j,k)=1.0
            obs_s_mask(i,j,k)=1.0
!
            t_obs(i,j,k)=-9999.0
            s_obs(i,j,k)=-9999.0
            dep_tobs(i,j,k)=z_das(k)
            dep_sobs(i,j,k)=z_das(k)
            err_tobs(i,j,k)=0.5;
            err_sobs(i,j,k)=0.06;
!
            do num=1,prf_num_TS_obs
              if (I_obs(num) .eq. i .and. J_obs(num) .eq. j) then
#if defined DAS_VERBOSE_INSITU
                  write(*,*) 'OK ',i,j
                  write(*,*) num, I_obs(num), J_obs(num) 
#endif
                 t_obs(i,j,k)=t_obs_raw(num,k)
                 s_obs(i,j,k)=s_obs_raw(num,k)

                 dep_tobs(i,j,k)=dep_tobs_raw(num,k)
                 dep_sobs(i,j,k)=dep_sobs_raw(num,k)

                 err_tobs(i,j,k)=err_tobs_raw(num,k)
                 err_sobs(i,j,k)=err_sobs_raw(num,k)
              endif                     
            enddo    ! prf

          enddo    !k
        enddo      !j
      enddo        !i
!
! deal with missing value in the proile with data
!
      do j=JstrR,JendR
        do i=IstrR,IendR
          do k=1,NDAS-1
              if ( (t_obs(i,j,k) .lt. -990) .or.
     &             (err_tobs(i,j,k) .lt. -990) ) then
                 t_obs(i,j,k)=0.0
                 obs_t_mask(i,j,k)=0.0
                 err_tobs(i,j,k)=0.5
              endif   
              if ( (s_obs(i,j,k) .lt. -990)  .or.
     &             (err_sobs(i,j,k) .lt. -990) ) then
                 s_obs(i,j,k)=0.0
                 obs_s_mask(i,j,k)=0.0
                 err_sobs(i,j,k)=0.06;
              endif
          enddo    !k
        enddo      !j
      enddo        !i


!
! implement simple QC and compute innovations
!
      crt=3.0
      crs=0.25
!
! up and down(dn) ---- k+1
!                   . dep obs
!                 ____ k
!

      do j=JstrR,JendR
        do i=IstrR,IendR
          do k=1,NDAS-1
            t_aup=(dep_tobs(i,j,k)-z_das(k))/(z_das(k+1)-z_das(k))
            t_adn=(z_das(k+1)-dep_tobs(i,j,k))/(z_das(k+1)-z_das(k))
            s_aup=(dep_sobs(i,j,k)-z_das(k))/(z_das(k+1)-z_das(k))
            s_adn=(z_das(k+1)-dep_sobs(i,j,k))/(z_das(k+1)-z_das(k))
!
! linear interpolation
!
            teq=t_adn*t_das(i,j,k,itemp)+t_aup*t_das(i,j,k+1,itemp)
            seq=s_adn*t_das(i,j,k,isalt)+s_aup*t_das(i,j,k+1,isalt)

            cft=t_obs(i,j,k) - teq
            obs_t_mask(i,j,k)=obs_t_mask(i,j,k)*rmask_das(i,j,k)
            if (abs(cft) .gt. crt ) then
#if defined DAS_VERBOSE_INSITU
                if ( obs_t_mask(i,j,k) .gt. 0.5 ) then
                    write(*,901),i,j,k,
     &               t_das(i,j,k,itemp),t_adn,
     &               t_das(i,j,k+1,itemp),t_aup,
     &               teq, cft, err_tobs(i,j,k)
                endif
#endif
                obs_t_mask(i,j,k)=0.0
            endif 
!
            cfs=s_obs(i,j,k) - seq
            obs_s_mask(i,j,k)=obs_s_mask(i,j,k)*rmask_das(i,j,k)
            if (abs(cfs) .gt. crs ) then
#if defined DAS_VERBOSE_INSITU
                if ( obs_s_mask(i,j,k) .gt. 0.5 ) then
                    write(*,901),i,j,k,
     &               t_das(i,j,k,isalt),s_adn,
     &               t_das(i,j,k+1,isalt),s_aup,
     &               seq, cfs, err_sobs(i,j,k)
                endif
#endif
                obs_s_mask(i,j,k)=0.0 
            endif 

!                                                        
            t_obs(i,j,k)=cft*obs_t_mask(i,j,k)          
            s_obs(i,j,k)=cfs*obs_s_mask(i,j,k)         
          enddo
        enddo
      enddo
!
      return
      end
