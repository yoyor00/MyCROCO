#include "cppdefs.h"

      subroutine das_adj_pressure_sm (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call das_adj_pressure_sm_tile(Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_adj_pressure_sm_tile 
     &                   (Istr,Iend,Jstr,Jend)
!
!--------------------------------------------------------------------
!  Set initial conditions for momentum and tracer variables using
!  analytical expressions.
!--------------------------------------------------------------------
!
      implicit none

# include "param.h"
# include "das_param.h"
# include "scalars.h"
# include "grid.h"
# include "das_ocean.h"
# include "das_covar.h"
# include "das_ocean_smooth.h"


      integer Istr,Iend,Jstr,Jend
      integer is,ie,js,je,i0,j0,i,j,k
      real rr,rc,ro
      integer count,num,n0
!
      integer srad,sdim
      parameter (srad=sp_rad_len, sdim=(2*srad+1)*(2*srad+1))
      real pp(sdim),dis(sdim)
      real pp_adj(sdim)
!
      real temp_p_adj(0:Lm+1, 0:Mm+1)
!
# include "compute_auxiliary_bounds.h"
!
!
      ro = (sp_rad*sp_rad)*( 
     *           (lonr(Lm/2,Mm/2)-lonr(Lm/2+1,Mm/2+1))
     &          *(lonr(Lm/2,Mm/2)-lonr(Lm/2+1,Mm/2+1))
     &          +(latr(Lm/2,Mm/2)-latr(Lm/2+1,Mm/2+1))
     &          *(latr(Lm/2,Mm/2)-latr(Lm/2+1,Mm/2+1)) )

      do k=1,NDAS
        do j=1,sdim
          pp_adj(j)=0.0
        enddo
        do j=0,Mm+1
          do i=0,Lm+1
            temp_p_adj(i,j)=0.0
          enddo
        enddo
!   
        do j=JstrR,JendR
          do i=IstrR,IendR
            count=0
            js=max(j-srad, 0)
            je=min(j+srad, Mm+1)
            is=max(i-srad,0)
            ie=min(i+srad, Lm+1)
            do j0=js,je
              do i0=is,ie
                rr=      (lonr(i,j)-lonr(i0,j0))
     &                  *(lonr(i,j)-lonr(i0,j0))
     &                  +(latr(i,j)-latr(i0,j0))
     &                  *(latr(i,j)-latr(i0,j0)) 
                count=count+1
                dis(count)=exp(-rr/ro)
              enddo
            enddo
            rc=0.0
            do n0=1,count
              rc=rc+dis(n0)
            enddo
            p_sm_adj(i,j,k)=p_sm_adj(i,j,k)/rc
            do n0=1,count  
              pp_adj(n0)=pp_adj(n0) + p_sm_adj(i,j,k)*dis(n0)
            enddo
            p_sm_adj(i,j,k)=0.0

            count=0 
            do j0=js,je  
              do i0=is,ie
                count=count+1
                temp_p_adj(i0,j0)=temp_p_adj(i0,j0)+pp_adj(count)
                pp_adj(count)=0.0
              enddo 
            enddo  
!
          enddo
        enddo
C$OMP CRITICAL (adj_psm_cr)
        do j=0,Mm+1                   
          do i=0,Lm+1 
            p_s_adj(i,j,k)=p_s_adj(i,j,k)+temp_p_adj(i,j)
          enddo      
        enddo       
C$OMP END CRITICAL (adj_psm_cr)

       enddo !k
!
      return
      end
