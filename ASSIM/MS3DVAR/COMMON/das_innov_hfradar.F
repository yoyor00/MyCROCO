#include "cppdefs.h"

      subroutine das_innov_hfradar (tile)
      implicit none
      integer tile,trd, omp_get_thread_num
# include "param.h"
# include "das_param.h"
# include "das_private_scratch.h"
# include "compute_tile_bounds.h"
      trd=omp_get_thread_num()
      call das_innov_hfradar_tile (Istr,Iend,Jstr,Jend,
     &        A2d(1,1,trd),A2d(1,2,trd),A2d(1,3,trd),A2d(1,4,trd))
      return
      end

      subroutine das_innov_hfradar_tile 
     &               (Istr,Iend,Jstr,Jend,ur,vr,eerr,rmaskr)
!
!--------------------------------------------------------------------
!  Set initial conditions for momentum and tracer variables using
!  analytical expressions.
!--------------------------------------------------------------------
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k, itrc
      real cfu,cfv,crt,hfumask, hfvmask
      real srad,srad0,rr,rc, lonc,latc
      real uu,vv,err,dis
      integer count,num,n0
# include "param.h"
# include "das_param.h"
# include "scalars.h"
# include "grid.h"
# include "das_ocean.h"
# include "das_covar.h"
# include "das_innov.h"
!
      real ur(PRIVATE_2D_SCRATCH_ARRAY),                                                                                             
     &     vr(PRIVATE_2D_SCRATCH_ARRAY),
     &     eerr(PRIVATE_2D_SCRATCH_ARRAY),
     &     rmaskr(PRIVATE_2D_SCRATCH_ARRAY)

      logical :: has_non_null_superobs
      integer :: count_non_null_superobs

!
!
# include "compute_auxiliary_bounds.h"
!
      srad=2.0*0.015*0.015
      srad0=0.006*0.006
      crt=0.25
      has_non_null_superobs=.false.
      count_non_null_superobs=0
#if defined DAS_HUGEHRADAR
!     Parameter in das_param.h fhuge = 1000.
      crt = fhuge * crt
!     fvmaskhuge = - fhuge * 900. 
#endif
#if defined DAS_SYNTHRADAR
!   pseudo-obs sampled on the Croco Grid decimation 1/5
!   angular size of horizontal grid cell ~ 
!   size of the horiz. grid cell ~ 5 to 10km
!   With original srad value 
!   sqrt(srad) ~ 1.4 * 0.015 which is 2.1km and a too small radius
!   within which we will not have any pseudo-obs data to contribute 
!   to the super obs 
!   Specific issue arouse from the recalculation of u_hf obs from mass
!   point to u point, with the associated rmaskr product
!   If in i,j we have a non null superobs value, with such decimation + original
!   srad, in the adjacent cell i-1,j, we will have a null super obs and 
!   the rmaskr product will be therefore null 
!   As a result, non of the non-null superobs will be finally assimilated 
      srad  = 15.*srad 
      srad0 = 15.*srad0
#endif

!
!     ... form superobs
!
! u
!
      do j=Jstr,Jend   
        do i=Istr,IendR
          u_hf(i,j)=0.0
          hf_umask(i,j)=0.0
        enddo
      enddo
!  
! v                  
      do j=Jstr,JendR 
        do i=Istr,Iend 
          v_hf(i,j)=0.0 
          hf_vmask(i,j) = 0.0
        enddo
      enddo   
      if (num_hfradar .le. 1 ) goto 100
!
      do j=Jstr-1,JendR
        do i=Istr-1,IendR
          if (rmask_das(i,j,NDAS) .gt. 0.5) then
            count=0
            uu=0.
            vv=0.
            rc=0.0
            do num=1,num_hfradar
              rr= (lonr(i,j)-lon_hf(num))
     &           *(lonr(i,j)-lon_hf(num))
     &           +(latr(i,j)-lat_hf(num))
     &           *(latr(i,j)-lat_hf(num))
              if (rr .lt. srad ) then
                count=count+1
                dis=max(rr,srad0)
                uu=uu+u_hf_raw(num)/dis
                vv=vv+v_hf_raw(num)/dis
                rc=rc+1.0/dis
              endif
            enddo
            if (count .ge. 1 ) then

#if !defined DAS_SYNTHRADAR
!
! rotate to ROMS curvelinear
! u_r = u.*cosa + v.*sina;
! v_r = v.*cosa - u.*sina;
!
              ur(i,j)=(uu*cos(angler(i,j))+vv*sin(angler(i,j)))/rc
              vr(i,j)=(vv*cos(angler(i,j))-uu*sin(angler(i,j)))/rc
#else
              ur(i,j)=uu
              vr(i,j)=vv
#endif
              rmaskr(i,j)=1.0
            else
              ur(i,j)=0.0
              vr(i,j)=0.0
              rmaskr(i,j)=0.0
            endif
          else
            ur(i,j)=0.0
            vr(i,j)=0.0
            rmaskr(i,j)=0.0
          endif
        enddo
      enddo
!
!     ... QC and INNOVATION
!
!u
      do j=Jstr,Jend
        do i=Istr,IendR
          cfu=0.5*(ur(i-1,j)+ur(i,j))-u_das(i,j,NDAS)
!         rmaskr product eliminates zero hf vel minus bkg model issue
!         abs(cfu) = abs(- u_das(i,j,NDAS)) .gt. crt 
          hf_umask(i,j) = rmaskr(i-1,j)*rmaskr(i,j)
     &                  * pmask_das(i,j,NDAS)*pmask_das(i,j+1,NDAS)  ! exclude u adjacent to coast
! BLXD extend counting to real hfradar
#if defined DAS_SYNTHRADAR
          if (abs(cfu) .gt. crt ) then
            hf_umask(i,j)=0.0
          else
            count_non_null_superobs=count_non_null_superobs+1
          endif
#else
          if (abs(cfu) .gt. crt ) hf_umask(i,j)=0.0
#endif
          u_hf(i,j)=cfu * hf_umask(i,j)
        enddo
      enddo

!
!
! v
      do j=Jstr,JendR
        do i=Istr,Iend
          cfv=0.5*(vr(i,j-1)+vr(i,j))-v_das(i,j,NDAS) 
          hf_vmask(i,j) = rmaskr(i,j-1)*rmaskr(i,j)
     &                  * pmask_das(i,j,NDAS)*pmask_das(i+1,j,NDAS)  ! exclude v adjacent to coast
! BLXD extend counting to real hfradar
#if defined DAS_SYNTHRADAR
          if ( abs(cfv) .gt. crt) then
            hf_vmask(i,j)=0.0
          else
            count_non_null_superobs=count_non_null_superobs+1
          endif
#else
          if ( abs(cfv) .gt. crt) hf_vmask(i,j)=0.0
#endif
          v_hf(i,j)=cfv * hf_vmask(i,j)
        enddo
      enddo

! BLXD extend counting to real hfradar
#if defined DAS_SYNTHRADAR
! Test if non null superobs    
      if (count_non_null_superobs .gt. 0) then
        has_non_null_superobs=.true.
      endif
#endif

  100 continue

! BLXD extend counting to real hfradar
#if defined DAS_SYNTHRADAR
      if (has_non_null_superobs) then
       write(*,'(6x,A,I6, (/8x,A))')
     &   'OK: HAS SUPER OBS HFRADAR : ', count_non_null_superobs
      else
       write(*,'(6x,A,I6, (/8x,A))')
     &   'OOPS: NO SUPER OBS HFRADAR : ', count_non_null_superobs
      endif
#endif

!
!6KM  6KM --------------------
!
      srad=2.0*0.025*0.025
      srad0=0.01*0.01
      crt=0.25
!
!     ... form superobs
!
! u
      do j=Jstr,Jend
        do i=Istr,IendR
          u_hf6(i,j)=0.0
          hf_uierr6(i,j)=0.0
          hf_umask6(i,j)=0.0
        enddo
      enddo
!
! v
      do j=Jstr,JendR
        do i=Istr,Iend
          v_hf6(i,j)=0.0
          hf_vierr6(i,j) = 0.0
          hf_vmask6(i,j) = 0.0
        enddo
      enddo
      if (num_hfradar6 .le. 1 ) goto 200
!
      do j=Jstr-1,JendR
        do i=Istr-1,IendR
          if (rmask_das(i,j,NDAS) .gt. 0.5) then
            count=0
            uu=0.
            vv=0.
            err=0.
            rc=0.0
            do num=1,num_hfradar6
              rr= (lonr(i,j)-lon_hf6(num))
     &           *(lonr(i,j)-lon_hf6(num))
     &           +(latr(i,j)-lat_hf6(num))
     &           *(latr(i,j)-lat_hf6(num))
              if (rr .lt. srad ) then
                count=count+1
                dis=max(rr,srad0)
                uu=uu+u_hf_raw6(num)/dis
                vv=vv+v_hf_raw6(num)/dis
                err=err+amperr_hf6(num)/dis
                rc=rc+1.0/dis
              endif
            enddo
            if (count .ge. 1 ) then
#if !defined DAS_SYNTHRADAR
!
! rotate to ROMS curvelinear
! u_r = u.*cosa + v.*sina;
! v_r = v.*cosa - u.*sina;
!
              ur(i,j)=(uu*cos(angler(i,j))+vv*sin(angler(i,j)))/rc
              vr(i,j)=(vv*cos(angler(i,j))-uu*sin(angler(i,j)))/rc
#else
              ur(i,j)=uu
              vr(i,j)=vv
#endif
              eerr(i,j)=err/rc
              rmaskr(i,j)=1.0
            else
              ur(i,j)=0.0
              vr(i,j)=0.0
              eerr(i,j)=0.0
              rmaskr(i,j)=0.0
            endif
          else
            ur(i,j)=0.0
            vr(i,j)=0.0
            eerr(i,j)=0.0
            rmaskr(i,j)=0.0
          endif
        enddo
      enddo
!
!     ... QC and INNOVATION
!     
!u      
      do j=Jstr,Jend
        do i=Istr,IendR
          cfu=0.5*(ur(i-1,j)+ur(i,j))-u_das(i,j,NDAS)
          hfumask = rmaskr(i-1,j)*rmaskr(i,j)
     &                  * pmask_das(i,j,NDAS)*pmask_das(i,j+1,NDAS)  ! exclude u adjacent to coast
          if (abs(cfu) .gt. crt ) hfumask=0.0
          u_hf6(i,j)=cfu * hfumask
          err=0.5*(eerr(i-1,j)+eerr(i,j))   !!!!!NOTE: depend on amplification error
          err=max(1.0,err/0.35)
          hf_uierr6(i,j)=hfumask * hf_ouv6/(err*err)
          hf_umask6(i,j)= hfumask
        enddo
      enddo 
!           
!           
! v         
      do j=Jstr,JendR
        do i=Istr,Iend
          cfv=0.5*(vr(i,j-1)+vr(i,j))-v_das(i,j,NDAS)
          hfvmask = rmaskr(i,j-1)*rmaskr(i,j)
     &                  * pmask_das(i,j,NDAS)*pmask_das(i+1,j,NDAS)  ! exclude v adjacent to coast
          if ( abs(cfv) .gt. crt) hfvmask=0.0
          v_hf6(i,j)= cfv * hfvmask
          err=0.5*(eerr(i,j-1)+eerr(i,j))
          err=max(1.0,err/0.35)
          hf_vierr6(i,j)=hfvmask * hf_ouv6/(err*err)
          hf_vmask6(i,j) = hfvmask
        enddo   
      enddo     
                
  200 continue
!
      return
      end
