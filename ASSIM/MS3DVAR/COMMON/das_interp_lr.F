#include "cppdefs.h"

      subroutine das_interp_lr (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call das_interp_lr_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_interp_lr_tile (Istr,Iend,Jstr,Jend)
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

!
# include "param.h"
# include "das_param.h"
# include "das_param_lr.h"
# include "scalars.h"
# include "grid.h"
# include "das_ocean.h"
# include "das_ocean_lr.h"
!
      integer is,ie,js,je,i0,j0,ir,jr
      real rr,rc,ro, ppm
      integer count,num,n0
!
      integer srad,sdim
      parameter (srad=5, sdim=(2*srad+1)*(2*srad+1))

      real tt(sdim),ss(sdim),dis(sdim)
      real su(sdim),sv(sdim),zz(sdim)

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
      ro = (nratio*nratio)* ( 
     &           (lonr(Lm/2,Mm/2)-lonr(Lm/2+1,Mm/2+1))
     &          *(lonr(Lm/2,Mm/2)-lonr(Lm/2+1,Mm/2+1))
     &          +(latr(Lm/2,Mm/2)-latr(Lm/2+1,Mm/2+1))
     &          *(latr(Lm/2,Mm/2)-latr(Lm/2+1,Mm/2+1)) )
!   
      do k=1,NDAS

        do j=JR_RANGE
          do i=IR_RANGE
            count=0
            ir=(i+nhalf)/nratio
            jr=(j+nhalf)/nratio
            js=max(jr-srad, 1)
            je=min(jr+srad, Mm_lr)
            is=max(ir-srad,1)
            ie=min(ir+srad, Lm_lr)
            do j0=js,je
              do i0=is,ie
                  rr=      (lonr(i,j)-lon_lr(i0,j0))
     &                    *(lonr(i,j)-lon_lr(i0,j0))
     &                    +(latr(i,j)-lat_lr(i0,j0))
     &                    *(latr(i,j)-lat_lr(i0,j0)) 
                  count=count+1
                  tt(count)=t_lr(i0,j0,k,itemp)
                  ss(count)=t_lr(i0,j0,k,isalt)
                  dis(count)=exp(-rr/ro)
              enddo
            enddo
            rc=0.0
            t_lr_s(i,j,k,itemp)=0.0
            t_lr_s(i,j,k,isalt)=0.0
            if ( count .gt. 0) then
              do n0=1,count
                rc=rc+dis(n0)
                t_lr_s(i,j,k,itemp)=t_lr_s(i,j,k,itemp)+tt(n0)*dis(n0)
                t_lr_s(i,j,k,isalt)=t_lr_s(i,j,k,isalt)+ss(n0)*dis(n0)
              enddo
              t_lr_s(i,j,k,itemp)=t_lr_s(i,j,k,itemp)/rc
              t_lr_s(i,j,k,isalt)=t_lr_s(i,j,k,isalt)/rc
            endif 
          enddo
        enddo
!
! u
!
        do j=JR_RANGE
          do i=IU_RANGE
            count=0
            ir=(i+nhalf-1)/nratio
            jr=(j+nhalf)/nratio
            ir=i/nratio
            jr=j/nratio
            js=max(jr-srad, 1)
            je=min(jr+srad, Mm_lr)
            is=max(ir-srad,1)
            ie=min(ir+srad, Lm_lr)
            do j0=js,je
              do i0=is,ie
                rr=      (lonr(i,j)-lon_lr(i0,j0))
     &                  *(lonr(i,j)-lon_lr(i0,j0))
     &                  +(latr(i,j)-lat_lr(i0,j0))
     &                  *(latr(i,j)-lat_lr(i0,j0))
                count=count+1
                su(count)=u_lr(i0,j0,k)
                dis(count)=exp(-rr/ro)
              enddo
            enddo
            rc=0.0
            u_lr_s(i,j,k)=0.0
            if ( count .gt. 0) then
              do n0=1,count 
                rc=rc+dis(n0)
                u_lr_s(i,j,k)=u_lr_s(i,j,k)+su(n0)*dis(n0)
              enddo
              u_lr_s(i,j,k)=u_lr_s(i,j,k)/rc
            endif
          enddo
        enddo
!                                                                                                                                
!v                                                                                                                               
!                                                                                                                                
        do j=JV_RANGE                                                                                                            
          do i=IR_RANGE                                                                                                          
            count=0                                                                                                            
            ir=(i+nhalf)/nratio
            jr=(j+nhalf-1)/nratio
            js=max(jr-srad, 1)
            je=min(jr+srad, Mm_lr)
            is=max(ir-srad,1)
            ie=min(ir+srad, Lm_lr)
            do j0=js,je
              do i0=is,ie
                rr=      (lonr(i,j)-lon_lr(i0,j0))
     &                  *(lonr(i,j)-lon_lr(i0,j0))
     &                  +(latr(i,j)-lat_lr(i0,j0))
     &                  *(latr(i,j)-lat_lr(i0,j0))
                count=count+1
                sv(count)=v_lr(i0,j0,k)
                dis(count)=exp(-rr/ro)
              enddo
            enddo
            rc=0.0
            v_lr_s(i,j,k)=0.0
            if ( count .gt. 0) then
              do n0=1,count
                rc=rc+dis(n0)
                v_lr_s(i,j,k)=v_lr_s(i,j,k)+sv(n0)*dis(n0)
              enddo
              v_lr_s(i,j,k)=v_lr_s(i,j,k)/rc
            endif
          enddo
        enddo
!
      enddo    !k

!
! zeta
!
      do j=JR_RANGE
        do i=IR_RANGE
          count=0
          ir=(i+nhalf)/nratio
          jr=(j+nhalf)/nratio
          js=max(jr-srad, 1)
          je=min(jr+srad, Mm_lr)
          is=max(ir-srad,1)
          ie=min(ir+srad, Lm_lr)
          do j0=js,je
            do i0=is,ie
              rr=      (lonr(i,j)-lon_lr(i0,j0))
     &                *(lonr(i,j)-lon_lr(i0,j0))
     &                +(latr(i,j)-lat_lr(i0,j0))
     &                *(latr(i,j)-lat_lr(i0,j0))
              count=count+1
              zz(count)=zeta_lr(i0,j0)
              dis(count)=exp(-rr/ro)
            enddo
          enddo
          rc=0.0
          zeta_lr_s(i,j)=0.0
          do n0=1,count
            rc=rc+dis(n0)
            zeta_lr_s(i,j)=zeta_lr_s(i,j)+zz(n0)*dis(n0)
          enddo
          zeta_lr_s(i,j)=zeta_lr_s(i,j)/rc
        enddo
      enddo
!
# undef IR_RANGE
# undef IU_RANGE
# undef JR_RANGE
# undef JV_RANGE
      return
      end
