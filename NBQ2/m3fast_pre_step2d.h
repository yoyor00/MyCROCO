! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! m3fast_pre_step2d.h (begin)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !********************************
! ! Compute rufrc & rvfrc: 
! ! internal mode forcing for
! ! barotropic fast mode
! !
! ! During the first fast time step convert rufrc & fvfrc into forcing
! ! terms by subtracting the fast-time "rubar" and "rvbar" from them;
! ! These forcing terms are then extrapolated forward in time using
! ! optimized Adams-Bashforth weights, so that the resultant rufrc
! ! and rvfrc are centered effectively at time n+1/2. From now on,
! ! these newly computed forcing terms will remain constant during
! ! the fast time stepping and will be added to "rubar" and "rvbar"
! ! during all subsequent fast time steps.
! !********************************
! !
      if (FIRST_FAST_STEP) then
! !    
! !  Pre-step2D: extrapolate rufrc
! !    
        if (FIRST_TIME_STEP) then
          cff3=0.                        ! This version is designed
          cff2=0.                        ! for coupling during 3D
          cff1=1.                        ! predictor sub-step: here
        elseif (FIRST_TIME_STEP+1) then  ! forcing term "rufrc" is
          cff3= 0.                       ! computed as instantaneous
          cff2=-0.5                      ! value at 3D time step
          cff1= 1.5                      ! "nstp" first, and then
        else                             ! extrapolated half-step
          cff3= 0.281105                 ! forward using  AM3-like
          cff2=-0.5-2.*cff3              ! weights optimized for
          cff1= 1.5+cff3                 ! maximum stability (with
        endif                            ! special care for startup)
!$acc kernels if(compute_on_device) default(present)
# if defined M3FAST_COUPLING2D  || defined NHINT
        do j=Jstr,Jend
         do i=IstrU,Iend
           cff=rufrc(i,j)-rubar(i,j)
           rufrc(i,j)=cff1*cff + cff2*rufrc_bak(i,j,3-nstp)
     &                         + cff3*rufrc_bak(i,j,nstp)
           rufrc_bak(i,j,nstp)=cff
         enddo
        enddo
# endif
# if defined M3FAST_C3D_UVSF && defined M3FAST_COUPLING3D
        do j=Jstr,Jend
         do i=IstrU,Iend
           cff=rubar(i,j)! -rufrc(i,j)
# ifndef BODYFORCE
!#  ifdef BSTRESS_FAST
!     &        - sustr(i,j)*om_u(i,j)*on_u(i,j)
!#  else
!     &        -(sustr(i,j)-bustr(i,j))*om_u(i,j)*on_u(i,j)
!#  endif
# endif
           ru_int2d_nbq(i,j)=cff1*cff+cff2*ru_int2d_nbq_bak(i,j,3-nstp)
     &                               +cff3*ru_int2d_nbq_bak(i,j,nstp)
           ru_int2d_nbq_bak(i,j,nstp)=cff
         enddo
        enddo
# endif      
# if defined M3FAST_COUPLING2D  || defined NHINT
        do j=JstrV,Jend
         do i=Istr,Iend
           cff=rvfrc(i,j)-rvbar(i,j)
           rvfrc(i,j)=cff1*cff + cff2*rvfrc_bak(i,j,3-nstp)
     &                         + cff3*rvfrc_bak(i,j,nstp)
           rvfrc_bak(i,j,nstp)=cff
         enddo
        enddo
# endif
# if defined M3FAST_C3D_UVSF &&  defined M3FAST_COUPLING3D
        do j=JstrV,Jend
         do i=Istr,Iend
           cff=rvbar(i,j)!-rvfrc(i,j)
# ifndef BODYFORCE
!#  ifdef BSTRESS_FAST
!     &        - svstr(i,j)*om_v(i,j)*on_v(i,j)
!#  else
!     &        -(svstr(i,j)-bvstr(i,j))*om_v(i,j)*on_v(i,j)
!#  endif
# endif
           rv_int2d_nbq(i,j)=cff1*cff+cff2*rv_int2d_nbq_bak(i,j,3-nstp)
     &                               +cff3*rv_int2d_nbq_bak(i,j,nstp)
           rv_int2d_nbq_bak(i,j,nstp)=cff
         enddo
        enddo
# endif        
!$acc end kernels
#ifdef RVTK_DEBUG_ADVANCED
C$OMP BARRIER
C$OMP MASTER
        call check_tab2d(rufrc_bak(:,:,nstp),'rufrc_bak_st_fast_d','u')
         call check_tab2d(rufrc_bak(:,:,3-nstp),'rufrc_bak_d','u')
         call check_tab2d(rufrc(:,:),'rufrc st_fast_d','u')
         call check_tab2d(rvfrc(:,:),'rvfrc st_fast_d','v')
#endif
! !
! !********************************
! ! Adjust barotropic pressure force:
! ! 
! ! Since coupling requires that pressure gradient term is computed
! ! using zeta(:,:,kstp) instead of zeta_new(:,:) needed to achieve
! ! numerical stability, apply compensation to shift pressure gradient
! ! terms from "kstp" to "knew": in essense, convert the fist 2D step
! ! from Forward Euler to Forward-Backward].
! !********************************
! !   
#  define zwrk UFx
#  define rzeta  UFe
#  define rzeta2  VFe
#  define rzetaSA VFx
!$acc kernels if(compute_on_device) default(present)
        do j=JstrV-1,Jend
          do i=IstrU-1,Iend
            zwrk(i,j)=zeta(i,j,knew)-zeta(i,j,kstp)
# if defined VAR_RHO_2D && defined SOLVE3D
            rzeta(i,j)=(1.+rhoS(i,j))*zwrk(i,j)
            rzeta2(i,j)=rzeta(i,j)*(zeta(i,j,knew)+zeta(i,j,kstp))
            rzetaSA(i,j)=zwrk(i,j)*(rhoS(i,j)-rhoA(i,j))
# else
            rzeta(i,j)=zwrk(i,j)
            rzeta2(i,j)=zwrk(i,j)*(zeta(i,j,knew)+zeta(i,j,kstp))
# endif
          enddo
        enddo
        cff=0.5*g
        do j=Jstr,Jend
          do i=Istr,Iend
            rubar(i,j)=rubar(i,j) +cff*on_u(i,j)*( (h(i-1,j)+h(i,j))  
     &          *(rzeta(i-1,j)-rzeta(i,j)) +rzeta2(i-1,j)-rzeta2(i,j)
# if defined VAR_RHO_2D && defined SOLVE3D
     &              +(h(i-1,j)-h(i,j))*( rzetaSA(i-1,j)+rzetaSA(i,j)
     &                        +0.333333333333*(rhoA(i-1,j)-rhoA(i,j))
     &                                     *(zwrk(i-1,j)-zwrk(i,j)) )
# endif
     &                                                              )
!
            rvbar(i,j)=rvbar(i,j) +cff*om_v(i,j)*( (h(i,j-1)+h(i,j))
     &          *(rzeta(i,j-1)-rzeta(i,j)) +rzeta2(i,j-1)-rzeta2(i,j)

# if defined VAR_RHO_2D && defined SOLVE3D
     &              +(h(i,j-1)-h(i,j))*( rzetaSA(i,j-1)+rzetaSA(i,j)
     &                        +0.333333333333*(rhoA(i,j-1)-rhoA(i,j))
     &                                     *(zwrk(i,j-1)-zwrk(i,j)) )
# endif
     &                                                              )
          enddo
        enddo            !--> discard  zwrk, rzeta, rzeta2, rzetaSA
!$acc end kernels        
# undef rzetaSA
# undef rzeta2
# undef rzeta
# undef zwrk
      endif   !<-- FIRST_FAST_STEP    
! !
! !********************************
! !   Update internal and external forcing 
! !      terms for NBQ mode
! !
! !   Compute external forcing terms ru_ext_nbq and updated internal 
! !   forcing terms ru_int_nbq for NBQ equations
! !
! !   ru_int_nbq     : RHS (3D) ( *mask & 2D correction)
! !   ru_ext_nbq     : RHS (2D)
! !   ru_ext_nbq_old : RHS (2D) at previous time-step
! !   ru_ext_nbq_sum : time-integrated RHS (2D)
! !********************************
! !
! !--------------------------------
! !  First fast time step only
! !--------------------------------
! !
# ifdef M3FAST_C3D_UVSF
       if (FIRST_FAST_STEP) then
!$acc kernels if(compute_on_device) default(present)
        do j=Jstr,Jend
          do i=IstrU,Iend
            ru_ext_nbq_sum(i,j)=0.
            ru_ext_nbq_old(i,j)=0.
          enddo
        enddo
        do j=JstrV,Jend
          do i=Istr,Iend
            rv_ext_nbq_sum(i,j)=0.
            rv_ext_nbq_old(i,j)=0.
          enddo
        enddo 
# if defined MASKING && defined M3FAST_2D
        do k=1,N
          do j=Jstr,Jend
            do i=IstrU,Iend
              ru_int_nbq(i,j,k)=ru_int_nbq(i,j,k)*umask(i,j)
            enddo
          enddo
        enddo
        do k=1,N
          do j=JstrV,Jend
            do i=Istr,Iend
              rv_int_nbq(i,j,k)=rv_int_nbq(i,j,k)*vmask(i,j)
            enddo
          enddo
        enddo
# endif
!$acc end kernels
       endif ! FIRST_FAST_STEP
# endif
! !
! !--------------------------------
! !  All fast time steps
! !--------------------------------
! !
# ifndef M3FAST_COUPLING2D
#  define ru_ext_nbq UFx
# endif
! !       if ( FIRST_FAST_STEP) then
! ! !$acc update device( ru_int_nbq, rv_int_nbq )
! !       endif
#  if defined RVTK_DEBUG 
C$OMP BARRIER
C$OMP MASTER
!        call Check_tab3d(ru_int_nbq,'ru_int_nbq (A0)','uint')
!        call Check_tab3d(rv_int_nbq,'rv_int_nbq (A0)','vint')
       call check_tab3d(ru_int_nbq,'ru_int_nbq (A1)','uint')
       call check_tab3d(rv_int_nbq,'rv_int_nbq (A1)','vint')
C$OMP END MASTER     
#  endif  
# if defined RVTK_DEBUG && defined NBQ
C$OMP BARRIER
C$OMP MASTER
!       call check_tab3d_sedlay(Hz,'Hz','r',N_sl+1,N)
!       call check_tab2d(ru_ext_nbq_old,'ru_ext_nbq_old (A1)','u')
!       call check_tab2d(rubar,'rubar (A1)','uint')
C$OMP END MASTER     
#  endif  
  
!$acc kernels if(compute_on_device) default(present)
#ifdef M3FAST_C3D_UVSF
      do j=Jstr,Jend
        do i=IstrU,Iend
          ru_ext_nbq(i,j)=
# ifndef M3FAST_COUPLING3D
     &                     (rufrc(i,j)+rubar(i,j))     
# else
     &              (rubar(i,j)-ru_int2d_nbq(i,j)) ! 2* delta with slow 
# endif     
     &                                *pm_u(i,j)*pn_u(i,j)
     &                                /(Drhs(i,j)+Drhs(i-1,j))
#  ifdef MASKING
     &                                *umask(i,j)
#  endif
   !  &                    (rubar(i,j)+rufrc(i,j))
          ru_ext_nbq_old(i,j)=ru_ext_nbq(i,j)-ru_ext_nbq_old(i,j)  ! 2* delta (m+1 - m)
          ru_ext_nbq_sum(i,j)=ru_ext_nbq_sum(i,j)+ru_ext_nbq(i,j)
        enddo
      enddo
      do k=1,N
        do j=Jstr,Jend
          do i=IstrU,Iend
            ru_int_nbq(i,j,k)=ru_int_nbq(i,j,k)
     &                        +ru_ext_nbq_old(i,j)       ! 2* delta (m+1 - m)
     &                        *(Hz(i-1,j,k)+Hz(i,j,k))
          enddo
        enddo
      enddo
       do j=Jstr,Jend
         do i=IstrU,Iend
           ru_ext_nbq_old(i,j)=ru_ext_nbq(i,j)   ! delta with slow
         enddo
       enddo
# endif /* M3FAST_C3D_UVSF */
!
# ifndef M3FAST_COUPLING2D
#  undef  ru_ext_nbq      
#  define rv_ext_nbq UFx
# endif
!
# ifdef M3FAST_C3D_UVSF
      do j=JstrV,Jend
        do i=Istr,Iend
          rv_ext_nbq(i,j)=
# ifndef M3FAST_COUPLING3D
     &                     (rvfrc(i,j)+rvbar(i,j))
# else
     &                    (rvbar(i,j)-rv_int2d_nbq(i,j))
# endif     
     &                                *pm_v(i,j)*pn_v(i,j)
     &                                /(Drhs(i,j)+Drhs(i,j-1))
#  ifdef MASKING
     &                                *vmask(i,j)
#  endif
          rv_ext_nbq_old(i,j)=rv_ext_nbq(i,j)-rv_ext_nbq_old(i,j)
          rv_ext_nbq_sum(i,j)=rv_ext_nbq_sum(i,j)+rv_ext_nbq(i,j)
        enddo
      enddo    
      do k=1,N
        do j=JstrV,Jend
          do i=Istr,Iend
            rv_int_nbq(i,j,k)=rv_int_nbq(i,j,k)
     &                       +rv_ext_nbq_old(i,j)
     &                         *(Hz(i,j,k)+Hz(i,j-1,k))
          enddo
        enddo
      enddo
      do j=JstrV,Jend
        do i=Istr,Iend
          rv_ext_nbq_old(i,j)=rv_ext_nbq(i,j)
        enddo
      enddo
# endif /* M3FAST_C3D_UVSF */
!
# ifndef M3FAST_COUPLING2D
#  undef rv_ext_nbq  
# endif
!
!$acc end kernels
# if defined RVTK_DEBUG && defined NBQ
C$OMP BARRIER
C$OMP MASTER
       call check_tab3d(ru_int_nbq(:,:,1),'ru_int_nbq (A2)','uint'
     &    ,ondevice=.TRUE.)
       call check_tab3d(rv_int_nbq,'rv_int_nbq (A2)','vint'
     &    ,ondevice=.TRUE.)
C$OMP END MASTER     
#  endif  

#  ifdef M3FAST_C3D_UVSF
#   if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI 
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,  
     &                        ru_int_nbq(START_2D_ARRAY,1))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,  
     &                        rv_int_nbq(START_2D_ARRAY,1))
#   endif 
#  endif 

# if defined RVTK_DEBUG && defined NBQ
C$OMP BARRIER
C$OMP MASTER
       call check_tab3d(ru_int_nbq,'ru_int_nbq (A)','uint'
     &    ,ondevice=.TRUE.)
       call check_tab3d(rv_int_nbq,'rv_int_nbq (A)','vint'
     &    ,ondevice=.TRUE.)
C$OMP END MASTER     
#  endif  
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! m3fast_pre_step2d.h (end)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
