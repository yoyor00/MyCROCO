#include "cppdefs.h"

      subroutine das_filter (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call das_filter_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_filter_tile (Istr,Iend,Jstr,Jend)
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
      real rr,rc,ro, ppm, Lcorr
      integer count,num,n0
!
      integer srad,sdim
      parameter (srad=11, sdim=(2*srad+1)*(2*srad+1))

      real tt(sdim),ss(sdim),dis(sdim)
      real su(sdim),sv(sdim),zz(sdim)
!
# include "param.h"
# include "das_param.h"
# include "scalars.h"
# include "grid.h"
# include "das_ocean.h"
# include "das_covar.h"
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
! srad/Lcorr should not be lagere than 20, dis iand thus rc may close to zero
! 
!
      Lcorr=3.5     
      ro = (Lcorr*Lcorr)* ( 
     &           (lonr(Lm/2,Mm/2)-lonr(Lm/2+1,Mm/2+1))
     &          *(lonr(Lm/2,Mm/2)-lonr(Lm/2+1,Mm/2+1))
     &          +(latr(Lm/2,Mm/2)-latr(Lm/2+1,Mm/2+1))
     &          *(latr(Lm/2,Mm/2)-latr(Lm/2+1,Mm/2+1)) )
!   
      do k=1,NDAS

        do j=JR_RANGE
          do i=IR_RANGE
            if (rmask_das(i,j,NDAS) .gt. 0.5) then
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
              t_s(i,j,k,itemp)=0.0
              t_s(i,j,k,isalt)=0.0
              if ( count .gt. 0) then
                do n0=1,count
                  rc=rc+dis(n0)
                  t_s(i,j,k,itemp)=t_s(i,j,k,itemp)+tt(n0)*dis(n0)
                  t_s(i,j,k,isalt)=t_s(i,j,k,isalt)+ss(n0)*dis(n0)
                enddo
                t_s(i,j,k,itemp)=t_s(i,j,k,itemp)/rc
                t_s(i,j,k,isalt)=t_s(i,j,k,isalt)/rc
               else
                t_s(i,j,k,itemp)=0.0
                t_s(i,j,k,isalt)=0.0
              endif 
!            else
!              t_s(i,j,k,itemp)=0.0
!              t_s(i,j,k,isalt)=0.0
            endif      !rmask
          enddo
        enddo
!
! u
!
        do j=JR_RANGE
          do i=IU_RANGE
            if (umask_das(i,j,NDAS) .gt. 0.5) then
              count=0
              js=max(j-srad, 1)
              je=min(j+srad, Mm)
              is=max(i-srad,1)
              ie=min(i+srad, Lm)
              do j0=js,je
                do i0=is,ie
!                  ppm=pmask_das(i,j,k)*pmask_das(i,j+1,k)
                  ppm=umask_das(i0,j0,k)
                  if (ppm .gt. 0.5) then
                    rr=      (lonr(i,j)-lonr(i0,j0))
     &                      *(lonr(i,j)-lonr(i0,j0))
     &                      +(latr(i,j)-latr(i0,j0))
     &                      *(latr(i,j)-latr(i0,j0))
                    count=count+1
                    su(count)=u_das(i0,j0,k)
                    dis(count)=exp(-rr/ro)
                  endif
                enddo
              enddo
              rc=0.0
              u_s(i,j,k)=0.0
              if ( count .gt. 0) then
                do n0=1,count 
                  rc=rc+dis(n0)
                  u_s(i,j,k)=u_s(i,j,k)+su(n0)*dis(n0)
                enddo
                u_s(i,j,k)=u_s(i,j,k)/rc
              endif
            else
              u_s(i,j,k)=0.0
            endif
          enddo
        enddo
!                                                                                                                                
!v                                                                                                                               
!                                                                                                                                
        do j=JV_RANGE                                                                                                            
          do i=IR_RANGE                                                                                                          
            if (vmask_das(i,j,NDAS) .gt. 0.5) then                                                                               
              count=0                                                                                                            
              js=max(j-srad, 1)                                                                                                  
              je=min(j+srad, Mm)
              is=max(i-srad,1)
              ie=min(i+srad, Lm)
              do j0=js,je
                do i0=is,ie
!                  ppm=pmask_das(i,j,k)*pmask_das(i+1,j,k)
                  ppm=vmask_das(i0,j0,k)
                  if ( ppm .gt. 0.5) then
                    rr=      (lonr(i,j)-lonr(i0,j0))
     &                      *(lonr(i,j)-lonr(i0,j0))
     &                      +(latr(i,j)-latr(i0,j0))
     &                      *(latr(i,j)-latr(i0,j0))
                    count=count+1
                    sv(count)=v_das(i0,j0,k)
                    dis(count)=exp(-rr/ro)
                  endif
                enddo
              enddo
              rc=0.0
              v_s(i,j,k)=0.0
              if ( count .gt. 0) then
                do n0=1,count
                  rc=rc+dis(n0)
                  v_s(i,j,k)=v_s(i,j,k)+sv(n0)*dis(n0)
                enddo
                v_s(i,j,k)=v_s(i,j,k)/rc
              endif
            else
              v_s(i,j,k)=0.0
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
          if (rmask_das(i,j,NDAS) .gt. 0.5) then
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
            zeta_s(i,j)=0.0
            do n0=1,count
              rc=rc+dis(n0)
              zeta_s(i,j)=zeta_s(i,j)+zz(n0)*dis(n0)
            enddo
            zeta_s(i,j)=zeta_s(i,j)/rc
          else
            zeta_s(i,j)=0.0
          endif
        enddo
      enddo


!
# undef IR_RANGE
# undef IU_RANGE
# undef JR_RANGE
# undef JV_RANGE
      return
      end
