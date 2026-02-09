#include "cppdefs.h"

      subroutine das_calc_corr (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call das_calc_corr_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_calc_corr_tile (Istr,Iend,Jstr,Jend)
!
!--------------------------------------------------------------------
! Purpose:
!--------------------------------------------------------------------
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k
      integer is,ie,js,je,i0,j0
      real rr,rc,ro, Lcorr
      integer count
!
!     integer srad
!     parameter (srad=11)
      integer srad,sdim
      parameter (srad=11, sdim=(2*srad+1)*(2*srad+1))
! BLXD TODO eliminate count and dis array since not useful here
      real dis(sdim)

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
!      ro = (Lcorr*Lcorr)* ( 
!     &           (lonr(Lm/2,Mm/2)-lonr(Lm/2+1,Mm/2+1))
!     &          *(lonr(Lm/2,Mm/2)-lonr(Lm/2+1,Mm/2+1))
!     &          +(latr(Lm/2,Mm/2)-latr(Lm/2+1,Mm/2+1))
!     &          *(latr(Lm/2,Mm/2)-latr(Lm/2+1,Mm/2+1)) )
!   

      write(*,*) 'in das calc corr'

! T/S corr(Xi,Xi')
      Lcorr=3.5     
      ro = (Lcorr*Lcorr)*   
     &           (lonr(Lm/2,Mm/2)-lonr(Lm/2+1,Mm/2+1))
     &          *(lonr(Lm/2,Mm/2)-lonr(Lm/2+1,Mm/2+1))
!
      !ctX_das(:,:,:)=0.0
      do i=IR_RANGE
        rc=0.0
        count=0
        !if (rmask_das(i,j,NDAS) .gt. 0.5) then
        is=max(i-srad,0)
        ie=min(i+srad, Lm+1)
        do i0=is,ie
          count=count+1
          rr=     (lonr(i,Mm/2)-lonr(i0,Mm/2))
     &           *(lonr(i,Mm/2)-lonr(i0,Mm/2))
          dis(count)=exp(-rr/ro)
          rc=rc+dis(count)
          ctX_das(i,i0,itemp)=ctX_das(i,i0,itemp)+dis(count)
          ctX_das(i,i0,isalt)=ctX_das(i,i0,isalt)+dis(count)
        enddo
        ctX_das(i,i0,itemp)=ctX_das(i,i0,itemp)/rc
        ctX_das(i,i0,isalt)=ctX_das(i,i0,isalt)/rc
      enddo

! Chi corr(Xi,Xi')
      Lcorr=3.5     
      ro = (Lcorr*Lcorr)*
     &           (lonr(Lm/2,Mm/2)-lonr(Lm/2+1,Mm/2+1))
     &          *(lonr(Lm/2,Mm/2)-lonr(Lm/2+1,Mm/2+1))
!
      !cchiX_das(:,:)=0.0
      do i=IR_RANGE
        rc=0.0
        count=0
        !if (rmask_das(i,j,NDAS) .gt. 0.5) then
        is=max(i-srad,0)
        ie=min(i+srad, Lm+1)
        do i0=is,ie
          count=count+1
          rr=     (lonr(i,Mm/2)-lonr(i0,Mm/2))
     &           *(lonr(i,Mm/2)-lonr(i0,Mm/2))
          dis(count)=exp(-rr/ro)
          rc=rc+dis(count)
          cchiX_das(i,i0)=cchiX_das(i,i0)+dis(count)
        enddo
        cchiX_das(i,i0)=cchiX_das(i,i0)/rc
      enddo

! Psi corr(Xi,Xi')
      Lcorr=3.5     
      ro = (Lcorr*Lcorr)*
     &           (lonr(Lm/2,Mm/2)-lonr(Lm/2+1,Mm/2+1))
     &          *(lonr(Lm/2,Mm/2)-lonr(Lm/2+1,Mm/2+1))
!
      !cpsiX_das(:,:)=0.0
      do i=IR_RANGE
        rc=0.0
        count=0
        !if (rmask_das(i,j,NDAS) .gt. 0.5) then
        is=max(i-srad,0)
        ie=min(i+srad, Lm+1)
        do i0=is,ie
          count=count+1
          rr=     (lonr(i,Mm/2)-lonr(i0,Mm/2))
     &           *(lonr(i,Mm/2)-lonr(i0,Mm/2))
          dis(count)=exp(-rr/ro)
          rc=rc+dis(count)
          cpsiX_das(i,i0)=cpsiX_das(i,i0)+dis(count)
        enddo
        cpsiX_das(i,i0)=cpsiX_das(i,i0)/rc
      enddo


! zeta corr(Xi,Xi')
      Lcorr=3.5     
      ro = (Lcorr*Lcorr)*
     &           (lonr(Lm/2,Mm/2)-lonr(Lm/2+1,Mm/2+1))
     &          *(lonr(Lm/2,Mm/2)-lonr(Lm/2+1,Mm/2+1))
!
      !czX_das(:,:)=0.0
      do i=IR_RANGE
        rc=0.0
        count=0
        !if (rmask_das(i,j,NDAS) .gt. 0.5) then
        is=max(i-srad,0)
        ie=min(i+srad, Lm+1)
        do i0=is,ie
          count=count+1
          rr=     (lonr(i,Mm/2)-lonr(i0,Mm/2))
     &           *(lonr(i,Mm/2)-lonr(i0,Mm/2))
          dis(count)=exp(-rr/ro)
          rc=rc+dis(count)
          czX_das(i,i0)=czX_das(i,i0)+dis(count)
        enddo
        czX_das(i,i0)=czX_das(i,i0)/rc
      enddo


! T/S corr(Eta,Eta')
      Lcorr=3.5     
      ro = (Lcorr*Lcorr)* ( 
     &          +(latr(Lm/2,Mm/2)-latr(Lm/2+1,Mm/2+1))
     &          *(latr(Lm/2,Mm/2)-latr(Lm/2+1,Mm/2+1)) )
!
      !ctE_das(:,:,:)=0.0
      do j=JR_RANGE
        rc=0.0
        count=0
        js=max(j-srad, 0)
        je=min(j+srad, Mm+1)
        do j0=js,je
          count=count+1
          rr =    (latr(Lm/2,j)-latr(Lm/2,j0))
     &           *(latr(Lm/2,j)-latr(Lm/2,j0)) 
          dis(count)=exp(-rr/ro)
          rc=rc+dis(count)
          ctE_das(i,i0,itemp)=ctE_das(i,i0,itemp)+dis(count)
          ctE_das(i,i0,isalt)=ctE_das(i,i0,isalt)+dis(count)
        enddo
        ctE_das(i,i0,itemp)=ctE_das(i,i0,itemp)/rc
        ctE_das(i,i0,isalt)=ctE_das(i,i0,isalt)/rc
      enddo

! Chi corr(Eta,Eta')
      Lcorr=3.5     
      ro = (Lcorr*Lcorr)* ( 
     &          +(latr(Lm/2,Mm/2)-latr(Lm/2+1,Mm/2+1))
     &          *(latr(Lm/2,Mm/2)-latr(Lm/2+1,Mm/2+1)) )
!
      !cchiE_das(:,:)=0.0
      do j=JR_RANGE
        rc=0.0
        count=0
        js=max(j-srad, 0)
        je=min(j+srad, Mm+1)
        do j0=js,je
          count=count+1
          rr =    (latr(Lm/2,j)-latr(Lm/2,j0))
     &           *(latr(Lm/2,j)-latr(Lm/2,j0)) 
          dis(count)=exp(-rr/ro)
          rc=rc+dis(count)
          cchiE_das(i,i0)=cchiE_das(i,i0)+dis(count)
        enddo
        cchiE_das(i,i0)=cchiE_das(i,i0)/rc
      enddo

! Psi corr(Eta,Eta')
      Lcorr=3.5     
      ro = (Lcorr*Lcorr)* ( 
     &          +(latr(Lm/2,Mm/2)-latr(Lm/2+1,Mm/2+1))
     &          *(latr(Lm/2,Mm/2)-latr(Lm/2+1,Mm/2+1)) )
!
      !cpsiE_das(:,:)=0.0
      do j=JR_RANGE
        rc=0.0
        count=0
        js=max(j-srad, 0)
        je=min(j+srad, Mm+1)
        do j0=js,je
          count=count+1
          rr =    (latr(Lm/2,j)-latr(Lm/2,j0))
     &           *(latr(Lm/2,j)-latr(Lm/2,j0)) 
          dis(count)=exp(-rr/ro)
          rc=rc+dis(count)
          cpsiE_das(i,i0)=cpsiE_das(i,i0)+dis(count)
        enddo
        cpsiE_das(i,i0)=cpsiE_das(i,i0)/rc
      enddo

! zeta corr(Eta,Eta')
      Lcorr=3.5     
      ro = (Lcorr*Lcorr)* ( 
     &          +(latr(Lm/2,Mm/2)-latr(Lm/2+1,Mm/2+1))
     &          *(latr(Lm/2,Mm/2)-latr(Lm/2+1,Mm/2+1)) )
!
      !czE_das(:,:)=0.0
      do j=JR_RANGE
        rc=0.0
        count=0
        js=max(j-srad, 0)
        je=min(j+srad, Mm+1)
        do j0=js,je
          count=count+1
          rr =    (latr(Lm/2,j)-latr(Lm/2,j0))
     &           *(latr(Lm/2,j)-latr(Lm/2,j0))
          dis(count)=exp(-rr/ro)
          rc=rc+dis(count)
          czE_das(i,i0)=czE_das(i,i0)+dis(count)
        enddo
        czE_das(i,i0)=czE_das(i,i0)/rc
      enddo

! T/S corr(Z,Z'')
      Lcorr=200     
      ro = (Lcorr*Lcorr)
      !ctZ_das(:,:,itemp)=0.
!     Parametrization using mixed layer depth !
      do k=1,NDAS
        ctZ_das(k,k,itemp)=1.
      enddo

      Lcorr=200     
      ro = (Lcorr*Lcorr)
      ctZ_das(:,:,isalt)=0.
!     Parametrization using mixed layer depth !
      do k=1,NDAS
        ctZ_das(k,k,isalt)=1.
      enddo

! Chi corr(Z,Z'')
      Lcorr=200     
      ro = (Lcorr*Lcorr)
      !cchiZ_das(:,:)=0.
!     Parametrization using mixed layer depth !
      do k=1,NDAS
        cchiZ_das(k,k)=1.
      enddo

! Psi corr(Z,Z'')
      Lcorr=200     
      ro = (Lcorr*Lcorr)
      !cpsiZ_das(:,:)=0.
!     Parametrization using mixed layer depth !
      do k=1,NDAS
        cpsiZ_das(k,k)=1.
      enddo

! u/v corr(Z,Z'')
      Lcorr=200     
      ro = (Lcorr*Lcorr)
      !cuZ_das(:,:)=0.
      !cvZ_das(:,:)=0.
!     Parametrization using mixed layer depth !
      do k=1,NDAS
        cuZ_das(k,k)=1.
        cvZ_das(k,k)=1.
      enddo

      write(*,*) 'out das calc corr'

!
# undef IR_RANGE
# undef IU_RANGE
# undef JR_RANGE
# undef JV_RANGE
      return
      end
