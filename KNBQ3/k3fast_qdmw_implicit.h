! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! K3FAST_qdmw_implicit.h (begin)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !********************************
! !********************************
! ! Tridiag (implicit) system RHS:
! ! Gaussian Elimination
! !********************************
! !
# if defined CVTK_DEBUG_ADV1 && defined KNBQ
      call check_tab3d(qdmw_nbq,'3d_fast bef impl qdmw_nbq',
     &  'r',ondevice=.TRUE.)
      call check_tab3d(soundspeed2_nbq,'3d_fast bef impl soundspeed',
     &  'r',ondevice=.TRUE.)
      call check_tab3d(visc2v_nbq,'3d_fast bef impl visc2v_nbq',
     &  'r',ondevice=.TRUE.)
      call check_tab3d(Hzw_nbq_inv,'3d_fast bef impl Hzw_nbq_inv',
     &  'r',ondevice=.TRUE.)
      call check_tab3d(Hzr_nbq_inv,'3d_fast bef impl Hzr_nbq_inv',
     &  'r',ondevice=.TRUE.)
# endif    
!$acc parallel loop collapse(2) default(present) async(1)
!$acc& private(k, cff, cff1, cff2)
      do j=Jstr,Jend
        do i=Istr,Iend
!
! !********************************
! !..........Bottom BC
! !********************************
! !
#   ifndef K3FAST_SEDLAYERS
        k=1
#   else
        k=-N_sl+1
#   endif
!! Calcul des coeffs
!        cff1 = dtfast*(thetaimp_nbq*thetaimp_nbq*soundspeed2_nbq(i,j,1)
!     &      *dtfast+visc2v_nbq(i,j,1))
!        cff2 = 1.D0 / cff1
!     
!           !-------------------------
!           ! \beta (normalization)
!           !-------------------------      
#   ifndef K3FAST_CSVISC2K
            cff1=1./(dtfast*(thetaimp_nbq**2*soundspeed2_nbq
     &              *dtfast+visc2v_nbq)) 
#   else
            cff1=1./(dtfast*(thetaimp_nbq**2*soundspeed2_nbq(i,j,k)
     &              *dtfast+visc2v_nbq(i,j,k))) 
            ! \beta(k+1)/\beta(k) 
            cff2=dtfast*(thetaimp_nbq**2*soundspeed2_nbq(i,j,k+1)
     &                  *dtfast+visc2v_nbq(i,j,k+1))
     &                *cff1
#   endif
!           !-------------------------
!           ! b*(k)=b(k)-a(k)*c(k-1)
!           !-------------------------        
            cff=cff1+Hzw_nbq_inv(i,j,k)       ! b(k)
     &              *(
     &                 Hzr_nbq_inv(i,j,k)
     &                +Hzr_nbq_inv(i,j,k+1)  
#   ifdef K3FAST_CSVISC2K
     &                *cff2  
#   endif
     &              )
!           !-------------------------
!           ! c*(k)=c(k)/b*(k)
!           !-------------------------
            FC3D(i,j,k)=-Hzw_nbq_inv(i,j,k+1)
     &              *Hzr_nbq_inv(i,j,k+1)/cff
#   ifdef K3FAST_CSVISC2K
     &                *cff2  
#   endif
!           !-------------------------
!           ! Y(k)=(f(k)-a(k)*Y(k-1))/b*(k)
!           !-------------------------
            DC3D(i,j,k)=qdmw_nbq(i,j,k)*cff1/cff  ! RHS=f(k)
     &              +Hzw_nbq_inv(i,j,k-1)         ! -ak(k)
     &              *Hzr_nbq_inv(i,j,k)
     &              *qdmw_nbq(i,j,k-1)/cff
              
! !
! !********************************
! !..........Inner layers
! !********************************
! !
! !--------------------------------
! ! Sediment Layer (implicit)
! !--------------------------------
! !
#   ifdef K3FAST_SEDLAYERS
        do k=-N_sl+2,0
!           !-------------------------
!           ! \beta(k) (normalization)
!           !-------------------------  
#   ifndef K3FAST_CSVISC2K
            cff1=1./(dtfast*(thetaimp_nbq**2*soundspeed2_nbq
     &              *dtfast+visc2v_nbq)) 
#   else
            cff1=1./(dtfast*(thetaimp_nbq**2*soundspeed2_nbq(i,j,k)
     &              *dtfast+visc2v_nbq(i,j,k))) 
            ! \beta(k+1)/\beta(k) 
            cff2=dtfast*(thetaimp_nbq**2*soundspeed2_nbq(i,j,k+1)
     &                  *dtfast+visc2v_nbq(i,j,k+1))
     &                *cff1
#   endif
!           !-------------------------
!           ! b*(k)=b(k)-a(k)*c*(k-1)
!           !-------------------------        
            cff=cff1+Hzw_nbq_inv(i,j,k)       ! b(k)
     &              *(
     &                Hzr_nbq_inv(i,j,k)
     &                +Hzr_nbq_inv(i,j,k+1)
#   ifdef K3FAST_CSVISC2K
     &                *cff2       
#   endif
     &              )
     &              + Hzw_nbq_inv(i,j,k-1)    ! -a(k)
     &               *Hzr_nbq_inv(i,j,k)
     &               *FC3D(i,j,k-1)           !c*(k-1)
!           !-------------------------
!           ! c*(k)=c(k)/b*(k)
!           !-------------------------
            FC3D(i,j,k)=-Hzw_nbq_inv(i,j,k+1)
     &              *Hzr_nbq_inv(i,j,k+1)/cff
#   ifdef K3FAST_CSVISC2K
     &              *cff2
#   endif
!           !-------------------------
!           ! Y(k)=(f(k)-a(k)*Y(k-1))/b*(k)
!           !-------------------------
            DC3D(i,j,k)=(qdmw_nbq(i,j,k)*cff1        ! RHS=f(k)
     &              +Hzw_nbq_inv(i,j,k-1)            ! -ak(k)
     &              *Hzr_nbq_inv(i,j,k)
!    &              *DC(i,k-1)) /cff           
     &              *DC3D(i,j,k-1)) /cff           
        enddo
#   endif
! !
! !--------------------------------
! ! Ocean inner Layer (implicit)
! !--------------------------------
! ! 
#   ifndef K3FAST_SEDLAYERS
        do k=2,N-1
#   else
        do k=1,N-1
#   endif
!           !-------------------------
!           ! \beta(k) (normalization)
!           !-------------------------  
#   ifndef K3FAST_CSVISC2K
            cff1=1./(dtfast*(thetaimp_nbq**2*soundspeed2_nbq
     &              *dtfast+visc2v_nbq)) 
#   else
            cff1=1./(dtfast*(thetaimp_nbq**2*soundspeed2_nbq(i,j,k)
     &              *dtfast+visc2v_nbq(i,j,k))) 
            ! \beta(k+1)/\beta(k) 
            cff2=dtfast*(thetaimp_nbq**2*soundspeed2_nbq(i,j,k+1)
     &                  *dtfast+visc2v_nbq(i,j,k+1))
     &                *cff1
#   endif
!           !-------------------------
!           ! b*(k)=b(k)-a(k)*c(k-1)
!           !-------------------------        
            cff=cff1+Hzw_nbq_inv(i,j,k)       ! b(k)
     &              *(
     &                Hzr_nbq_inv(i,j,k)
     &                +Hzr_nbq_inv(i,j,k+1)
#   ifdef K3FAST_CSVISC2K
     &                  *cff2       
#   endif
     &              )
     &              + Hzw_nbq_inv(i,j,k-1)    ! -a(k)
     &               *Hzr_nbq_inv(i,j,k)
     &               *FC3D(i,j,k-1)           !c(k-1)
!           !-------------------------
!           ! c*(k)=c(k)/b*(k)
!           !-------------------------
            FC3D(i,j,k)=-Hzw_nbq_inv(i,j,k+1)
     &              *Hzr_nbq_inv(i,j,k+1)/cff
#   ifdef K3FAST_CSVISC2K
     &              *cff2
#   endif
!           !-------------------------
!           ! Y(k)=(f(k)-a(k)*Y(k-1))/b*(k)
!           !-------------------------
            DC3D(i,j,k)=(qdmw_nbq(i,j,k)*cff1        ! RHS=f(k)
     &              +Hzw_nbq_inv(i,j,k-1)            ! -ak(k)
     &              *Hzr_nbq_inv(i,j,k)
!    &              *DC(i,k-1)) /cff           
     &              *DC3D(i,j,k-1)) /cff           
        enddo
! !
! !--------------------------------
! !..........Surface BC
! !--------------------------------
! !
        k=N
!           !-------------------------
!           ! \beta (normalization)
!           !-------------------------      
#   ifndef K3FAST_CSVISC2K
          cff1=1./(dtfast*(thetaimp_nbq**2*soundspeed2_nbq
     &                    *dtfast+visc2v_nbq)) 
#   else
          cff1=1./(dtfast*(thetaimp_nbq**2*soundspeed2_nbq(i,j,N)
     &                    *dtfast+visc2v_nbq(i,j,k))) 
#   endif
!           !-------------------------
!           ! b*(k)=b(k)-a(k)*c(k-1)
!           !-------------------------        
          cff=cff1+Hzw_nbq_inv(i,j,N)
     &            *Hzr_nbq_inv(i,j,N) 
     &       +Hzw_nbq_inv(i,j,N-1)
     &             *Hzr_nbq_inv(i,j,N)
     &             *FC3D(i,j,N-1) 
!           !-------------------------
!           ! c*(k)=c(k)/b*(k)
!           !-------------------------
          FC3D(i,j,N)=0. 
!           !-------------------------
!           ! Y(k)=(f(k)-a(k)*Y(k-1))/b*(k)
!           !-------------------------
!           ! WSURF
          DC3D(i,j,N)=(qdmw_nbq(i,j,N)*cff1
     &             +Hzw_nbq_inv(i,j,N-1)
     &             *Hzr_nbq_inv(i,j,N)  
     &             *DC3D(i,j,N-1))/cff
          qdmw_nbq(i,j,N)=DC3D(i,j,N) 
! !
! !********************************
! !.....Solves tri-diag system
! ! (starts above with qdmw_nbq(N))
! !********************************
! !
#   ifndef K3FAST_SEDLAYERS
        do k=N-1,1,-1
#   else
        do k=N-1,-N_sl+1,-1
#   endif
            qdmw_nbq(i,j,k)=DC3D(i,j,k)
     &                      -FC3D(i,j,k)*qdmw_nbq(i,j,k+1)
          enddo            
                        
          enddo                        
        enddo 
!$acc end parallel loop

! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! K3FAST_qdmw_implicit.h (end)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
