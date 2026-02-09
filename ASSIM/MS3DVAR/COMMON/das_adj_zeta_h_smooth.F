#include "cppdefs.h"
#ifdef DAS_PWTHZETATOT

      subroutine das_adj_zeta_h_sm (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call das_adj_zeta_h_sm_tile(Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_adj_zeta_h_sm_tile 
     &                   (Istr,Iend,Jstr,Jend)
!
!--------------------------------------------------------------------
!  Set initial conditions for momentum and tracer variables using
!  analytical expressions.
!--------------------------------------------------------------------
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j
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
      parameter (srad=sz_rad_len, sdim=(2*srad+1)*(2*srad+1))
      real zz(sdim),dis(sdim)
      real zz_adj(sdim)
!
      real temp_h_adj(0:Lm+1, 0:Mm+1)
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
      ro = (sz_rad*sz_rad)*( 
     *           (lonr(Lm/2,Mm/2)-lonr(Lm/2+1,Mm/2+1))
     &          *(lonr(Lm/2,Mm/2)-lonr(Lm/2+1,Mm/2+1))
     &          +(latr(Lm/2,Mm/2)-latr(Lm/2+1,Mm/2+1))
     &          *(latr(Lm/2,Mm/2)-latr(Lm/2+1,Mm/2+1)) )

      do j=1,sdim
        zz_adj(j)=0.0
      enddo
      do j=0,Mm+1
        do i=0,Lm+1
          temp_h_adj(i,j)=0.0
        enddo
      enddo
!   
      do j=JR_RANGE
        do i=IR_RANGE
          count=0
          js=max(j-srad, 0)
          je=min(j+srad, Mm+1)
          is=max(i-srad,0)
          ie=min(i+srad, Lm+1)
          do j0=js,je
            do i0=is,ie
              rr=      (lonr(i,j)-lonr(i0,j0))
     &                *(lonr(i,j)-lonr(i0,j0))
     &                +(latr(i,j)-latr(i0,j0))
     &                *(latr(i,j)-latr(i0,j0)) 
              count=count+1
              dis(count)=exp(-rr/ro)
            enddo
          enddo
          rc=0.0
          do n0=1,count
            rc=rc+dis(n0)
          enddo
          zeta_h_sm_adj(i,j)=zeta_h_sm_adj(i,j)/rc
          do n0=1,count  
            zz_adj(n0)=zz_adj(n0) + zeta_h_sm_adj(i,j)*dis(n0)
          enddo
          zeta_h_sm_adj(i,j)=0.0

          count=0 
          do j0=js,je  
            do i0=is,ie
              count=count+1
              temp_h_adj(i0,j0)=temp_h_adj(i0,j0)+zz_adj(count)
              zz_adj(count)=0.0
            enddo 
          enddo  
!
        enddo
      enddo
C$OMP CRITICAL (adj_zeta_h_cr)
      do j=0,Mm+1                   
        do i=0,Lm+1 
          zeta_h_adj(i,j)=zeta_h_adj(i,j)+temp_h_adj(i,j)
        enddo      
      enddo       
C$OMP END CRITICAL (adj_zeta_h_cr)
!
# undef IR_RANGE
# undef IU_RANGE
# undef JR_RANGE
# undef JV_RANGE
      return
      end
#endif /* DAS_PWTHZETATOT */
