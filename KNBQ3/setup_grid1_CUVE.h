! ! ***********************************
! ! ***********************************
! ! setup_grid1_CUVE.h
! ! ***********************************
! ! ***********************************
! !
! ! ***********************************
! ! Flags
! ! ***********************************
! !
        ifl_bathy(:)=0
        ifl_bathy(1)=0
        ifl_bathy(2)=1
        ifl_bathy(3)=0
        ifl_bathy(4)=0
        ifl_bathy(5)=0
        ifl_bathy(6)=4
! !
! ! ***********************************
! ! Coordinates of specific points
! ! ***********************************
! !
        xv1(1)=2.205;yv1(1)=0.265
        xv1(2)=3.565;yv1(2)=1.185
        xv1(3)=xv1(2)-xv1(1);yv1(3)=yv1(2)-yv1(1)
        
        xv2(1)=-0.055;yv2(1)=-3.275
        xv2(2)=0.805;yv2(2)=-1.405
        xv2(3)=xv2(2)-xv2(1);yv2(3)=yv2(2)-yv2(1)

        xv3(1)=-4.805;yv3(1)=-1.415
        xv3(2)=-8.265;yv3(2)=1.765
        xv3(3)=xv3(2)-xv3(1);yv3(3)=yv3(2)-yv3(1)
        
        xv4(1)=-4.725;yv4(1)=8.195
        xv4(2)=-3.745;yv4(2)=8.195
        xv4(3)=xv4(2)-xv4(1);yv4(3)=yv4(2)-yv4(1)
        
        xv5(1)=0.595;yv5(1)=-1.545;yv5(1)=-1.405
        
        yv6(1)=0.
        
        hd2=h
! !
! ! ***********************************
! ! Copy bathy
! ! ***********************************
! !
   !     do j=JstrR,JendR
   !     do i=IstrR,IendR
          !h(i,j)=max(h(i,j),1.e-2)         ! Hmin
          ! hd2(i,j)=h(i,j)                  ! Copy to compute differences
#  ifdef MASKING
           hd3=rmask           ! Can be used to change mask...
           hd3u=umask            ! Can be used to change mask...
           hd3v=vmask             ! Can be used to change mask...
#  else
           hd3 =1.                     ! Can be used to change mask...
           hd3u=1.                     ! Can be used to change mask...
           hd3v=1.                     ! Can be used to change mask...
#  endif
  !      enddo
  !      enddo
  !      
! !
! ! ***********************************
! ! Computes Diffusion coef.
! ! ***********************************
! !            
        hd4=0.
        
        if (ifl_bathy(6)==3.or.ifl_bathy(6)==4) then
        hd=h      ! Dummy variable 
          do j=Jstr-2,Jend+1
          do i=Istr-2,Iend+1
         if (hd3(i,j).ne.0) then
             hd4(i,j)=         
     &  ( hd3u(i  ,j)*(hd(i  ,j)-hd(i-1,j))**2
     &   +hd3u(i+1,j)*(hd(i+1,j)-hd(i  ,j))**2)
     &     /max(1.,hd3u(i,j)+hd3u(i+1,j))
          endif
        enddo
        enddo
          
        call exchange_u2d_tile (Istr,Iend,Jstr,Jend,hd4(-2,-2))
        
          do j=Jstr-2,Jend+1
          do i=Istr-2,Iend+1
         if (hd3(i,j).ne.0) then
             hd4(i,j)=hd4(i,j)   
     & +( hd3v(i,j  )*(hd(i,j  )-hd(i,j-1))**2
     &   +hd3v(i,j+1)*(hd(i,j+1)-hd(i,j  ))**2)
     &     /max(1.,hd3v(i,j)+hd3v(i,j+1))
          endif
        enddo
        enddo
        call exchange_v2d_tile (Istr,Iend,Jstr,Jend,hd4(-2,-2))
        
        endif
        
        call exchange_r2d_3pts_tile (Istr,Iend,Jstr,Jend,h(-2,-2))
! !
! ! ***********************************
! ! Main iterative loop
! ! ***********************************
! !            
        do it=1,1000 ! Begins loop
! !
! ! -----------------------------------
! ! "Diffuses" bathy (coef=kappa*h)
! ! -----------------------------------
! !            
        if (ifl_bathy(1).ne.0) then
            
        hd=h      ! Dummy variable 

        do j=JstrR,JendR                         
        do i=IstrR,IendR
         if ( hd3(i,j).ne.0 ) then

! Diffusion coef.
          if (ifl_bathy(1)==1) then
             val_bathy(1)=1./5.
          elseif (ifl_bathy(1)==2) then
             val_bathy(1)=hd(i,j)/2. 
          elseif (ifl_bathy(1)==3.or.ifl_bathy(1)==4) then
             val_bathy(1)=100.*sqrt(           
     &  ( hd3u(i  ,j)*(hd(i  ,j)-hd(i-1,j))**2*pm_u(i  ,j)**2
     &   +hd3u(i+1,j)*(hd(i+1,j)-hd(i  ,j))**2*pm_u(i+1,j)**2)**2
     &     /max(1.,hd3u(i,j)+hd3u(i+1,j))
     & +( hd3v(i,j  )*(hd(i,j  )-hd(i,j-1))**2*pn_v(i,j  )**2
     &   +hd3v(i,j+1)*(hd(i,j+1)-hd(i,j  ))**2*pn_v(i,j+1)**2)**2
     &     /max(1.,hd3v(i,j)+hd3v(i,j+1))
     &                        )
             if (ifl_bathy(1)==4) val_bathy(1)=val_bathy(1)*hd(i,j)
          endif
          
! Laplacian operator (with mask):
          h(i,j)=h(i,j)+val_bathy(1)*(

     &     hd(i+1,j)*hd3(i+1,j)
     &    +hd(i,j)*(1.-hd3(i+1,j))
     &    -2.*hd(i,j)
     &     +hd(i-1,j)*hd3(i-1,j)
     &     +hd(i,j)*(1.-hd3(i-1,j))
     
     &    +hd(i,j+1)*hd3(i,j+1)
     &    +hd(i,j)*(1.-hd3(i,j+1))
     &    -2.*hd(i,j)
     &    +hd(i,j-1)*hd3(i,j-1)
     &    +hd(i,j)*(1.-hd3(i,j-1))
     &                            )

! Bi-laplacian operator:
!         h(i,j)=h(i,j)+1./val_bathy(1)*(
!     &    -hd(i+2,j)+4.*h(i+1,j)-6.*hd(i,j)+4.*hd(i-1,j)-hd(i-2,j)
!     &    -hd(i,j+2)+4.*h(i,j+1)-6.*hd(i,j)+4.*hd(i,j-1)-hd(i,j-2)
!     &                            )
        endif
        enddo
        enddo

        call exchange_r2d_3pts_tile (Istr,Iend,Jstr,Jend,h(-2,-2))
        endif

! !
! ! -----------------------------------
! ! "Diffuses" bathy (coef=kappa*grad h)
! ! -----------------------------------
! !               
        if (ifl_bathy(6).ne.0) then      
         hd4=0.     ! Diffusion coef.
        
         if (ifl_bathy(6)==3.or.ifl_bathy(6)==4) then
          hd=h      ! Dummy variable 
 
!....Computes grad-x
          do j=Jstr-2,Jend+1
          do i=Istr-2,Iend+1
           if (hd3(i,j).ne.0) then
             hd4(i,j)=         
     &        ( hd3u(i  ,j)*(hd(i  ,j)-hd(i-1,j))**2
     &         +hd3u(i+1,j)*(hd(i+1,j)-hd(i  ,j))**2)
     &        /max(1.,hd3u(i,j)+hd3u(i+1,j))
           endif
          enddo
          enddo
          call exchange_u2d_tile (Istr,Iend,Jstr,Jend,hd4(-2,-2))

!....Computes grad-y    
          do j=Jstr-2,Jend+1
          do i=Istr-2,Iend+1
           if (hd3(i,j).ne.0) then
               hd4(i,j)=hd4(i,j)   
     &          +( hd3v(i,j  )*(hd(i,j  )-hd(i,j-1))**2
     &            +hd3v(i,j+1)*(hd(i,j+1)-hd(i,j  ))**2)
     &         /max(1.,hd3v(i,j)+hd3v(i,j+1))
           endif
          enddo
          enddo
          call exchange_v2d_tile (Istr,Iend,Jstr,Jend,hd4(-2,-2))
        
       endif  !ifl_bathy(6)
        
       call exchange_r2d_3pts_tile (Istr,Iend,Jstr,Jend,h(-2,-2))
 
!....Computes diffusion coef. & diffuses      
       hd=h      ! Dummy variable 
        
       do j=Jstr,Jend
       do i=Istr,Iend
        if (hd3(i,j).ne.0) then     
          if (ifl_bathy(6)==1) then
             cff=1./2.
          elseif (ifl_bathy(6)==2) then
             cff=hd(i,j)/2. 
          elseif (ifl_bathy(6)==3) then
             cff=1./5000.+1.*sqrt(hd4(i,j))
     &             /max(1.e-20,hd(i,j)) 
          elseif (ifl_bathy(6)==4) then
             cff=1./5000.+2.*sqrt(hd4(i,j))
     &          /max(1.e-20,hd(i,j))
     &          *(exp(-(hd(i,j)-5e-3)**2/0.1**2)     
!    &           +exp(-(hd(i,j)-0.3)**2/0.025**2)    
     &           +2.*exp(-(hd(i,j)-0.5)**2/0.2**2))
     &        /2.
          endif
            
! Laplacian operator (with mask):
          h(i,j)=h(i,j)+cff*(
     &     hd(i+1,j)*hd3(i+1,j)+hd(i,j)*(1.-hd3(i+1,j))
     &    -2.*hd(i,j)
     &    +hd(i-1,j)*hd3(i-1,j)+hd(i,j)*(1.-hd3(i-1,j))
     &    +hd(i,j+1)*hd3(i,j+1)+hd(i,j)*(1.-hd3(i,j+1))
     &    -2.*hd(i,j)
     &    +hd(i,j-1)*hd3(i,j-1)+hd(i,j)*(1.-hd3(i,j-1))
     &                            )

! Bi-laplacian operator:
!         h(i,j)=h(i,j)+1./val_bathy(6)*(
!     &    -hd(i+2,j)+4.*h(i+1,j)-6.*hd(i,j)+4.*hd(i-1,j)-hd(i-2,j)
!     &    -hd(i,j+2)+4.*h(i,j+1)-6.*hd(i,j)+4.*hd(i,j-1)-hd(i,j-2)
!     &                            )
         endif
        enddo
        enddo

        call exchange_r2d_tile (Istr,Iend,Jstr,Jend,h(-2,-2))
        endif
        
        
! ! -----------------------------------
! ! "Averages" bathy (first local pass)
! ! -----------------------------------
! !        
        if (ifl_bathy(2)==1) then
        hd=h  
        do it2=1,1   
        do j=JstrR,JendR                         
        do i=IstrR,IendR
            if  (hd3(i,j).eq.1.and.(
!    Eastern shelf:
     &       (xr(i,j).ge.1.9.or.yr(i,j).le.-0.9.and.
     &        ((yv1(3)*(xr(i,j)-xv1(1))-xv1(3)*(yr(i,j)-yv1(1)).gt.0)
     &        .or.
     &         (yv2(3)*(xr(i,j)-xv2(1))-xv2(3)*(yr(i,j)-yv2(1)).gt.0.)
     &         )) 
!    Western shelf:
     &        .or.(xr(i,j).le.xv3(1)+0.5.and.
     &        (yv3(3)*(xr(i,j)-xv3(1))-xv3(3)*(yr(i,j)-yv3(1)).gt.0.))
     &        .or.(xr(i,j).le.xv4(1).and.
     &        (yv3(4)*(xr(i,j)-xv4(1))-xv4(3)*(yr(i,j)-yv4(1)).gt.0.))
     &         )) then
               h(i,j)=( hd(i,j)*hd3(i,j)
     &          +hd(i+1,j)*hd3(i+1,j)+hd(i-1,j)*hd3(i-1,j)
     &          +hd(i,j+1)*hd3(i,j+1)+hd(i,j-1)*hd3(i,j-1) )
     &          /max(1.e-20,hd3(i,j)
     &          +hd3(i+1,j)+hd3(i-1,j)
     &          +hd3(i,j+1)+hd3(i,j-1) )
        endif
        enddo
        enddo   
      
        call exchange_r2d_tile (Istr,Iend,Jstr,Jend, h(-2,-2))
        enddo
        
        endif

! !
! ! -----------------------------------
! ! "Averages" bathy (second pass)
! ! -----------------------------------
! !        
        if (ifl_bathy(3)==1) then
        hd=h
        do j=JstrR,JendR                         
        do i=IstrR,IendR
         if (hd3(i,j).ne.0.and.h(i,j).ge.0e-2) then
          h(i,j)=(
     &       hd(i  ,j  )*hd3(i  ,j  )
     &      +hd(i+1,j  )*hd3(i+1,j  ) 
     &      +hd(i  ,j+1)*hd3(i  ,j+1)
     &      +hd(i-1,j  )*hd3(i-1,j  ) 
     &      +hd(i  ,j-1)*hd3(i  ,j-1)
     &      +hd3(i+1,j+1)*hd(i+1,j+1)*hd3(i+1,j  )*hd3(i,j+1)
     &      +hd3(i-1,j+1)*hd(i-1,j+1)*hd3(i-1,j  )*hd3(i,j+1)
     &      +hd3(i-1,j-1)*hd(i-1,j-1)*hd3(i  ,j-1)*hd3(i-1,j)
     &      +hd3(i+1,j-1)*hd(i+1,j-1)*hd3(i  ,j-1)*hd3(i+1,j) ) /
     &      ( hd3(i  ,j  )  +hd3(i+1,j  )+hd3(i  ,j+1)+hd3(i-1,j  )
     &       +hd3(i  ,j-1)+hd3(i+1,j+1)*hd3(i+1,j  )*hd3(i,j+1)
     &       +hd3(i-1,j+1)*hd3(i-1,j  )*hd3(i  ,j+1)
     &       +hd3(i-1,j-1)*hd3(i  ,j-1)*hd3(i-1,j  )
     &       +hd3(i+1,j-1)*hd3(i  ,j-1)*hd3(i+1,j  )  )
        endif
        enddo
        enddo

        call exchange_r2d_3pts_tile (Istr,Iend,Jstr,Jend,h(-2,-2))

        endif
        
        enddo  ! Ends it-loop
        
        call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                        h(-2,-2))
  
! !
! ! ***********************************
! ! Computes differences & updates
! ! ***********************************
! !        
        hd4 = hd2
     !   h   = hd2      ! To be used for testing.
        hd2 = hd2-h
    !    h=hd4
        
        
