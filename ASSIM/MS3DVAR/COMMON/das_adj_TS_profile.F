#include "cppdefs.h"

      subroutine das_adj_TS_obs (tile)
      implicit none
      integer tile
#include "param.h"
#include "das_param.h"
#include "compute_tile_bounds.h"
      call das_adj_TS_obs_tile (Istr,Iend,Jstr,Jend) 
      return
      end

      subroutine das_adj_TS_obs_tile (Istr,Iend,Jstr,Jend) 
      implicit none
      integer Istr,Iend,Jstr,Jend,    i,j,k
      real cft, cfs, teq, seq
      real t_aup, t_adn, s_aup, s_adn, obs_ot, obs_os
      real teq_adj, seq_adj
#include "param.h"
# include "das_param.h"
#include "grid.h"
#include "scalars.h"
# include "das_ocean.h"
# include "das_covar.h"
# include "das_lbgfs.h"
# include "das_innov.h"
!
# include "compute_auxiliary_bounds.h"

      teq_adj=0.0
      seq_adj=0.0
!
      do j=JstrR,JendR
        do i=IstrR,IendR
          do k=1,NDAS-1

            t_aup=(dep_tobs(i,j,k)-z_das(k))/(z_das(k+1)-z_das(k))
            t_adn=(z_das(k+1)-dep_tobs(i,j,k))/(z_das(k+1)-z_das(k))
            obs_ot=1.0/(err_tobs(i,j,k)*err_tobs(i,j,k))

            s_aup=(dep_sobs(i,j,k)-z_das(k))/(z_das(k+1)-z_das(k))
            s_adn=(z_das(k+1)-dep_sobs(i,j,k))/(z_das(k+1)-z_das(k))
            obs_os=1.0/(err_sobs(i,j,k)*err_sobs(i,j,k))

            teq=t_adn*t_s(i,j,k,itemp)+t_aup*t_s(i,j,k+1,itemp)
            seq=s_adn*t_s(i,j,k,isalt)+s_aup*t_s(i,j,k+1,isalt)

            cft=teq-t_obs(i,j,k)
            cfs=seq-s_obs(i,j,k)
           
            teq_adj=teq_adj + cft*obs_t_mask(i,j,k)*obs_ot
            seq_adj=teq_adj + cfs*obs_s_mask(i,j,k)*obs_os
   
            t_s_adj(i,j,k,itemp)=t_s_adj(i,j,k,itemp)    +t_adn*teq_adj
            t_s_adj(i,j,k+1,itemp)=t_s_adj(i,j,k+1,itemp)+t_aup*teq_adj
            teq_adj=0.0

            t_s_adj(i,j,k,isalt)=t_s_adj(i,j,k,isalt)    +s_adn*seq_adj
            t_s_adj(i,j,k+1,isalt)=t_s_adj(i,j,k+1,isalt)+s_aup*seq_adj
            seq_adj=0.0

          enddo
        enddo
      enddo
!
!
      return
      end
