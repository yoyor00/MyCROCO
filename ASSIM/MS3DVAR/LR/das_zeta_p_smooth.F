#include "cppdefs.h"

      subroutine das_zeta_p_smooth (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call das_zeta_p_smooth_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_zeta_p_smooth_tile (Istr,Iend,Jstr,Jend)
!
!--------------------------------------------------------------------
!  Set initial conditions for momentum and tracer variables using
!  analytical expressions.
!--------------------------------------------------------------------
!
      implicit none

!
# include "param.h"
# include "das_param.h"
# include "scalars.h"
# include "grid.h"
# include "das_ocean.h"
# include "das_covar.h"
# include "das_ocean_smooth.h"
!
      integer Istr,Iend,Jstr,Jend, i,j,k
      integer is,ie,js,je,i0,j0
      real rr,rc,ro
      integer count,num,n0
!
      integer srad,sdim
      parameter (srad=sm_rad,sdim=(3*srad+1)*(3*srad+1))
      real pp(sdim),zz(sdim),dis(sdim)
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
#define DISTANCES
#ifdef DISTANCES
      ro =  ((lonr(Lm/2,Mm/2)-lonr(Lm/2+1,Mm/2+1))
     &          *(lonr(Lm/2,Mm/2)-lonr(Lm/2+1,Mm/2+1))
     &          +(latr(Lm/2,Mm/2)-latr(Lm/2+1,Mm/2+1))
     &          *(latr(Lm/2,Mm/2)-latr(Lm/2+1,Mm/2+1)))
#endif
!
      do k=1,NDAS
        do j=JR_RANGE
          do i=IR_RANGE
            if (rmask_das(i,j,k) .gt. 0.5) then
              count=0
              js=max(j-srad, 0)
              je=min(j+srad, Mm+1)
              is=max(i-srad,0)
              ie=min(i+srad, Lm+1)
              do j0=js,je
                do i0=is,ie
                  if (rmask_das(i0,j0,k) .gt. 0.5) then
#ifdef DISTANCES
                    rr=( (lonr(i,j)-lonr(i0,j0))
     &                      *(lonr(i,j)-lonr(i0,j0))
     &                      +(latr(i,j)-latr(i0,j0))
     &                      *(latr(i,j)-latr(i0,j0)) )
#endif
                    count=count+1
                    pp(count)=p_s(i0,j0,k)
                    if (k .eq. NDAS) then
                      zz(count)=zeta_s(i0,j0)
                    endif
#ifdef DISTANCES
                    dis(count)=max(ro,rr)
#endif
                  endif
                enddo
              enddo
              rc=0.0
              p_sm(i,j,k)=0.0
              if (k .eq. NDAS) then
                zeta_sm(i,j)=0.0
              endif
              do n0=1,count
#ifdef DISTANCES
                rc=rc+1./dis(n0)
                p_sm(i,j,k)=p_sm(i,j,k)+pp(n0)/dis(n0)
                if (k .eq. NDAS) then
                  zeta_sm(i,j)=zeta_sm(i,j)+zz(n0)/dis(n0)
                endif
#else
                rc=rc+1.
                p_sm(i,j,k)=p_sm(i,j,k)+pp(n0)
                if (k .eq. NDAS) then
                  zeta_sm(i,j)=zeta_sm(i,j)+zz(n0)
                endif
#endif
              enddo
              p_sm(i,j,k)=p_sm(i,j,k)/rc
             if (k .eq. NDAS) then
               zeta_sm(i,j)=zeta_sm(i,j)/rc
             endif
            else
              p_sm(i,j,k)=0.
              if (k .eq. NDAS) then
                zeta_sm(i,j)=0.0
              endif
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
