 
#ifdef AGRIF
        irhox = Agrif_Irhox()
        irhoy = Agrif_Irhoy()
        irhot = Agrif_Irhot()
#endif

      if (FIRST_2D_STEP) then         ! Meaning of temporal indices
        kbak=kstp                     ! ------- -- -------- -------
        kold=kstp                     ! m-2     m-1      m      m+1
        cff1=1.                       ! kold    kbak    kstp    knew
        cff2=0.
        cff3=0.
      elseif (FIRST_2D_STEP+1) then
        kbak=kstp-1
        if (kbak.lt.1) kbak=4
        kold=kbak
        cff1=1.                  ! Logically AB2-AM3 forward-backward 
        cff2=0.                  ! scheme with coefficients chosen for
        cff3=0.                  ! maximum stability ... (see below)
      else
        kbak=kstp-1
        if (kbak.lt.1) kbak=4
        kold=kbak-1
        if (kold.lt.1) kold=4
#ifdef M2FILTER_NONE
        cff1= 1.5+mybeta
        cff2=-2.0*mybeta-0.5
        cff3= mybeta
#else
        cff1= 1.781105
        cff2=-1.06221
        cff3= 0.281105
#endif
      endif

#ifdef NBQ_ZETAW
      if (iic==1.and.iif==1) then
        knew2=2
        kstp2=1
        kbak2=4
        kold2=3
      else
        kold2=kbak2
        kbak2=kbak2+1
        if (kbak2==5) kbak2=1
      endif

      knew2=knew
      kstp2=kstp
#endif

!#ifndef NBQ_ZETAW

#if defined M2FILTER_NONE && defined NBQ_ZETAW
!        cff4= 1.
!        cff5= 0.
!        cff6= 0.
!        cff4= 1.5+mybeta
!        cff5=-2.0*mybeta-0.5
!        cff6= mybeta
      if (FIRST_2D_STEP) then
#ifdef SOLVE3D
        cff0=0.                !---> Compute pressure-gradient
        cff1=1.                !  terms using just zeta(:,:,kstp)
#else
        cff0=1.
        cff1=0.
#endif
        cff2=0.
        cff3=0.
      elseif (FIRST_2D_STEP+1) then
        cff0= 1.0833333333333    ! Logically AB2-AM3 forward-backward
        cff1=-0.1666666666666    ! scheme with coefficients chosen for
        cff2= 0.0833333333333    ! maximum stability, while maintaining
        cff3= 0.                 ! third-accuracy; alpha_max=1.73
      else
#ifdef M2FILTER_NONE
        cff0=0.5+2.*myepsilon+mygamma+2.*myalpha
        cff1=1.-cff0-mygamma-myepsilon
        cff2=mygamma
        cff3=myepsilon
#else
        cff0=0.614
        cff1=0.285
        cff2=0.088
        cff3=0.013
#endif
      endif

      if (iif.eq.1) then
        cff8=1.D0
        cff9=0.D0
        cff10=0.D0
      elseif (iif.eq.1+1) then
        cff8=1.D0
        cff9=0.D0
        cff10=0.D0
      else
        cff8= 1.5D0+mybeta
        cff9=-2.0D0*mybeta-0.5D0
        cff10= mybeta
      endif

#endif
!#endif

!
#ifdef RVTK_DEBUG_ADVANCED
C$OMP BARRIER
C$OMP MASTER
       call check_tab2d(zeta(:,:,kstp),'zeta step2d #0','r')
# ifndef NBQ_ZETAW
       call check_tab2d(ubar(:,:,kstp),'ubar step2d #0','u')
       call check_tab2d(vbar(:,:,kstp),'vbar step2d #0','v')
# endif
C$OMP END MASTER       
#endif  
!
!-----------------------------------------------------------------------  
! Computes D,u & v for RHS
!-----------------------------------------------------------------------
!
!
#if defined EW_PERIODIC || defined MPI || !defined M2_HADV_UP3
      imin=IstrU-2
      imax=Iend+1
# define IV_EXT_RANGE Istr-1,Iend+1
#else
      imin=max(IstrU-3,0)
      imax=min(Iend+2,Lm+1)
# define IV_EXT_RANGE max(Istr-2,0),min(Iend+2,Lm+1)
#endif
#if defined NS_PERIODIC || defined MPI || !defined M2_HADV_UP3
      jmin=JstrV-2
      jmax=Jend+1
# define JU_EXT_RANGE Jstr-1,Jend+1
#else
      jmin=max(JstrV-3,0)
      jmax=min(Jend+2,Mm+1)
# define JU_EXT_RANGE max(Jstr-2,0),min(Jend+2,Mm+1)
#endif
!
#ifndef NBQ
      do j=jmin,jmax
        do i=imin,imax
          Drhs(i,j)= cff1*(zeta(i,j,kstp)+h(i,j))
     &              +cff2*(zeta(i,j,kbak)+h(i,j))
     &              +cff3*(zeta(i,j,kold)+h(i,j))
        enddo
      enddo
#else /* NBQ */
      do j=jmin,jmax
        do i=imin,imax
# ifdef NBQ_MASS
#  ifndef NBQ_ZETAW
          Drhs(i,j)= cff1*(zeta(i,j,kstp)+h(i,j))*rhobar_nbq(i,j,kstp)
     &              +cff2*(zeta(i,j,kbak)+h(i,j))*rhobar_nbq(i,j,kbak)
     &              +cff3*(zeta(i,j,kold)+h(i,j))*rhobar_nbq(i,j,kold)
#  else /*  ! NBQ_ZETAW */
#   ifndef NBQ_AB3
          Drhs(i,j)=(zeta(i,j,kstp2)+h(i,j)) *rhobar_nbq(i,j,kstp2)
#   else
          Drhs(i,j)= cff8 *(zeta(i,j,kstp2)+h(i,j))
     &                    *rhobar_nbq(i,j,kstp2)
     &              +cff9 *(zeta(i,j,kbak2)+h(i,j))
     &                    *rhobar_nbq(i,j,kbak2)
     &              +cff10*(zeta(i,j,kold2)+h(i,j))
     &                    *rhobar_nbq(i,j,kold2)
#   endif
#  endif /* NBQ_ZETAW */
# else /* ! NBQ_MASS */
#  ifndef NBQ_ZETAW
          Drhs(i,j)= cff1*(zeta(i,j,kstp)+h(i,j))
     &              +cff2*(zeta(i,j,kbak)+h(i,j))
     &              +cff3*(zeta(i,j,kold)+h(i,j))
#  else /*  ! NBQ_ZETAW */
#   ifndef NBQ_AB3
          Drhs(i,j)=(zeta(i,j,kstp2)+h(i,j)) 
#   else
          Drhs(i,j)= cff8 *(zeta(i,j,kstp2)+h(i,j))
     &              +cff9 *(zeta(i,j,kbak2)+h(i,j))
     &              +cff10*(zeta(i,j,kold2)+h(i,j))
#   endif
#  endif /* NBQ_ZETAW */
# endif /* NBQ_MASS */
        enddo
      enddo
#endif /* NBQ */

#ifdef RVTK_DEBUG_ADVANCED
C$OMP BARRIER
C$OMP MASTER
!       call check_tab2d(Drhs(:,:),'Drhs step2d #00','r')
!       call check_tab2d(urhs(:,:),'urhs step2d #1','u')
C$OMP END MASTER       
#endif 

# if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
#  ifdef NBQ_ZETAW
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                   ubar(START_2D_ARRAY,kstp2))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                   vbar(START_2D_ARRAY,kstp2))
#  elif defined NBQ
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                   ubar(START_2D_ARRAY,kstp))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                   vbar(START_2D_ARRAY,kstp))
#  endif
#endif

      do j=JU_EXT_RANGE
        do i=imin+1,imax
#ifndef NBQ_ZETAW
          urhs(i,j)=cff1*ubar(i,j,kstp) +cff2*ubar(i,j,kbak)
     &                                  +cff3*ubar(i,j,kold)
#else
!         urhs(i,j)=cff4*ubar(i,j,kstp2) +cff5*ubar(i,j,kbak2)
!    &                                   +cff6*ubar(i,j,kold2)
        !   urhs(i,j)=DU_nbq(i,j)    !ubar(i,j,kstp2)
#ifndef NBQ_AB3
          urhs(i,j)=ubar(i,j,kstp2)
#else
          urhs(i,j)=cff8*ubar(i,j,kstp2) +cff9*ubar(i,j,kbak2)
     &                                  +cff10*ubar(i,j,kold2)
#endif
#endif
#if defined MRL_WCI && defined MASKING
          urhs(i,j)=urhs(i,j)*umask(i,j)+ust2d(i,j)*(umask(i,j)-1.0)
#endif
          DUon(i,j)=0.5*(Drhs(i,j)+Drhs(i-1,j))*on_u(i,j)*( urhs(i,j)
#ifdef MRL_WCI
     &                                                   + ust2d(i,j)
#endif
     &                                                              )
        enddo
      enddo
#undef JU_EXT_RANGE

      do j=jmin+1,jmax
        do i=IV_EXT_RANGE
#ifndef NBQ_ZETAW
          vrhs(i,j)=cff1*vbar(i,j,kstp) +cff2*vbar(i,j,kbak)
     &                                  +cff3*vbar(i,j,kold)
#else
!         vrhs(i,j)=cff4*vbar(i,j,kstp2) +cff5*vbar(i,j,kbak2)
!    &                                   +cff6*vbar(i,j,kold2)

#ifndef NBQ_AB3
         !  vrhs(i,j)=DV_nbq(i,j)   !vbar(i,j,kstp2)
            vrhs(i,j)=vbar(i,j,kstp2)
#else
          vrhs(i,j)=cff8*vbar(i,j,kstp2) +cff9*vbar(i,j,kbak2)
     &                                  +cff10*vbar(i,j,kold2)
#endif
#endif
#if defined MRL_WCI && defined MASKING
          vrhs(i,j)=vrhs(i,j)*vmask(i,j)+vst2d(i,j)*(vmask(i,j)-1.0)
#endif
          DVom(i,j)=0.5*(Drhs(i,j)+Drhs(i,j-1))*om_v(i,j)*( vrhs(i,j)
#ifdef MRL_WCI
     &                                                   + vst2d(i,j)
#endif
     &                                                              )
        enddo
      enddo
#undef IV_EXT_RANGE

#ifdef M2_HADV_UP3
# if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                   urhs(START_2D_ARRAY))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                   vrhs(START_2D_ARRAY))
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                   Duon (START_2D_ARRAY))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                   Dvom(START_2D_ARRAY))
# endif
#endif

#ifdef OBC_VOLCONS
      call set_DUV_bc_tile (Istr,Iend,Jstr,Jend, Drhs, DUon,DVom)
#endif

#ifdef RVTK_DEBUG_ADVANCED
C$OMP BARRIER
C$OMP MASTER
       call check_tab2d(zeta(:,:,kstp),'zeta step2d #1','r')
# ifndef NBQ_ZETAW
       call check_tab2d(ubar(:,:,kstp),'ubar step2d #1','u')
       call check_tab2d(vbar(:,:,kstp),'vbar step2d #1','v')
# endif
C$OMP END MASTER       
#endif 
