#include "cppdefs.h"

      subroutine das_extrapolate9 (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call das_extrapolate9_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_extrapolate9_tile (Istr,Iend,Jstr,Jend)
!
!--------------------------------------------------------------------
! Purpose:  Fill the land grids for DA with the flat bottom
!--------------------------------------------------------------------
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k
      integer is,ie,js,je,i0,j0
      real rr,rc,ro
      integer count,num,n0
!
      integer srad,sdim
      parameter (srad=11, sdim=(2*srad+1)*(2*srad+1))
      real zz(sdim),tt(sdim),ss(sdim),dis(sdim)
!
# include "param.h"
# include "das_param.h"
# include "scalars.h"
# include "grid.h"
# include "das_ocean.h"
# include "das_covar.h"
# include "das_ocean_smooth.h"
#include "das_ocean9.h"
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
      ro = (5.0*5.0)* ( 
     &           (lonr(Lm/2,Mm/2)-lonr(Lm/2+1,Mm/2+1))
     &          *(lonr(Lm/2,Mm/2)-lonr(Lm/2+1,Mm/2+1))
     &          +(latr(Lm/2,Mm/2)-latr(Lm/2+1,Mm/2+1))
     &          *(latr(Lm/2,Mm/2)-latr(Lm/2+1,Mm/2+1)) )
!
!zeta
!
      do j=JR_RANGE
        do i=IR_RANGE
          if (rmask_das(i,j,NDAS) .lt. 0.5) then
            count=0
            js=max(j-srad, 0)
            je=min(j+srad, Mm+1)
            is=max(i-srad,0)
            ie=min(i+srad, Lm+1)
            do j0=js,je
              do i0=is,ie
                if (rmask_das(i0,j0,NDAS) .gt. 0.5) then
                  rr=      (lonr(i,j)-lonr(i0,j0))
     &                    *(lonr(i,j)-lonr(i0,j0))
     &                    +(latr(i,j)-latr(i0,j0))
     &                    *(latr(i,j)-latr(i0,j0))
                  count=count+1
                  zz(count)=zeta_das(i0,j0)
                  dis(count)=exp(-rr/ro)
                endif
              enddo
            enddo
            rc=0.0
            zeta_das9(i,j)=0.0
            do n0=1,count
              rc=rc+dis(n0)
              zeta_das9(i,j)=zeta_das9(i,j)+zz(n0)*dis(n0)
            enddo
            zeta_das9(i,j)=(zeta_das(i,j)+zeta_das9(i,j))/(1.0+rc)
          else
            zeta_das9(i,j)=zeta_das(i,j)
          endif    
        enddo     
      enddo
!   
! T/S
!
      do k=1,NDAS
        do j=JR_RANGE
          do i=IR_RANGE
            if (rmask_das(i,j,k) .lt. 0.5) then
              count=0
              js=max(j-srad, 0)
              je=min(j+srad, Mm+1)
              is=max(i-srad,0)
              ie=min(i+srad, Lm+1)
              do j0=js,je
                do i0=is,ie
                  if (rmask_das(i0,j0,k) .gt. 0.5) then
                    rr=      (lonr(i,j)-lonr(i0,j0))
     &                      *(lonr(i,j)-lonr(i0,j0))
     &                      +(latr(i,j)-latr(i0,j0))
     &                      *(latr(i,j)-latr(i0,j0)) 
                    count=count+1
                    tt(count)=t_das(i0,j0,k,itemp)
                    ss(count)=t_das(i0,j0,k,isalt)
                    dis(count)=exp(-rr/ro)
                  endif
                enddo
              enddo
              rc=0.0
              t_das9(i,j,k,itemp)=0.0
              t_das9(i,j,k,isalt)=0.0
              do n0=1,count
                rc=rc+dis(n0)
                t_das9(i,j,k,itemp)=t_das9(i,j,k,itemp)+tt(n0)*dis(n0)
                t_das9(i,j,k,isalt)=t_das9(i,j,k,isalt)+ss(n0)*dis(n0)
              enddo
              t_das9(i,j,k,itemp)=(t_das(i,j,k,itemp)+t_das9(i,j,k,itemp))
     &                                               /(rc+1.0)
              t_das9(i,j,k,isalt)=(t_das(i,j,k,isalt)+t_das9(i,j,k,isalt))
     &                                               /(rc+1.0)
            else
              t_das9(i,j,k,itemp)=t_das(i,j,k,itemp)
              t_das9(i,j,k,isalt)=t_das(i,j,k,isalt)
            endif      !rmask
          enddo
        enddo
      enddo    !k
!
! How about vertical?
!

!
# undef IR_RANGE
# undef IU_RANGE
# undef JR_RANGE
# undef JV_RANGE
      return
      end
