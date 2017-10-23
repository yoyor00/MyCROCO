
 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
#ifdef NBQ_MASS
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
!# ifdef NBQ_ZETAW
!          zetaw_nbq(i,j)= zetaw_nbq(i,j) SWITCH rmask(i,j)
!# endif

!***********************************************************************
# ifndef NBQ_ZETAW
!***********************************************************************
          Dnew(i,j)=zeta_new(i,j)
          rhobar_nbq(i,j,knew) = 2.*rhobar_nbq(i,j,kstp)-rhobar_nbq(i,j,kbak)
!***********************************************************************
# else 
!***********************************************************************
       !   rhobar_nbq(i,j,knew2) = 2.*rhobar_nbq(i,j,kstp2)-rhobar_nbq(i,j,kbak2)
       !   zwrk(i,j)= zeta(i,j,kstp2)/rhobar_nbq(i,j,kstp2)-h(i,j)
           zwrk(i,j)=zeta(i,j,kstp2)
!          zwrk(i,j)= 
!     &              +cff4*zeta(i,j,kstp2)/rhobar_nbq(i,j,kstp2)
!     &              +cff5*zeta(i,j,kbak2)/rhobar_nbq(i,j,kbak2)
!     &              +cff6*zeta(i,j,kold2)/rhobar_nbq(i,j,kold2)
!     &              -h(i,j)

#  ifdef MASKING
          Dnew(i,j)=(zeta(i,j,kstp2)*rmask(i,j)+h(i,j))*rhobar_nbq(i,j,kstp2)  ! CAUTION: rhobar_nbq must 1 on Mask !!
          zwrk(i,j)=zwrk(i,j)*rmask(i,j)
#  else
          Dnew(i,j)=(zeta(i,j,kstp2)+h(i,j))*rhobar_nbq(i,j,kstp2)
#  endif

!***********************************************************************
# endif /* ! NBQ_ZETAW */
!***********************************************************************

# if defined VAR_RHO_2D && defined SOLVE3D
          rzeta(i,j)=(1.+rhoS(i,j))*zwrk(i,j)
          rzeta2(i,j)=rzeta(i,j)*zwrk(i,j)
          rzetaSA(i,j)=zwrk(i,j)*(rhoS(i,j)-rhoA(i,j))
# else
          rzeta(i,j)=zwrk(i,j)
          rzeta2(i,j)=zwrk(i,j)*zwrk(i,j)
# endif
        enddo
      enddo
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
#else /* ! NBQ_MASS */
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
# ifndef NBQ
          zeta_new(i,j)=zeta_new(i,j) SWITCH rmask(i,j)
# endif

!***********************************************************************
# ifndef NBQ_ZETAW
!***********************************************************************
          Dnew(i,j)=zeta_new(i,j)+h(i,j)
          zwrk(i,j)=cff0*zeta_new(i,j) +cff1*zeta(i,j,kstp)
     &             +cff2*zeta(i,j,kbak)+cff3*zeta(i,j,kold)  
!*********************************************************************** 
# else
!***********************************************************************
          zetaw_nbq(i,j,knew2)= zetaw_nbq(i,j,knew2)
     &            SWITCH rmask(i,j)   ! Utilité ?????
!         Dnew(i,j)=(zeta(i,j,kstp2)+h(i,j))*rhobar_nbq(i,j,kstp2)
          Dnew(i,j)=(zeta(i,j,kstp2)+h(i,j))
          zwrk(i,j)=zeta(i,j,kstp2) SWITCH rmask(i,j)  
!         zwrk(i,j)=zeta_new(i,j) SWITCH rmask(i,j)  
!          zwrk(i,j)= cff4*zeta(i,j,kstp2)
!     &              +cff5*zeta(i,j,kbak2)+cff6*zeta(i,j,kold2)  
!          zwrk(i,j)= cff4*zeta(i,j,knew2)+cff5*zeta(i,j,kstp2)
!     &              +cff6*zeta(i,j,kbak2)+cff7*zeta(i,j,kold2)  
!***********************************************************************
# endif /* ! NBQ_ZETAW*/
!***********************************************************************

# if defined VAR_RHO_2D && defined SOLVE3D
          rzeta(i,j)=(1.+rhoS(i,j))*zwrk(i,j)
          rzeta2(i,j)=rzeta(i,j)*zwrk(i,j)
          rzetaSA(i,j)=zwrk(i,j)*(rhoS(i,j)-rhoA(i,j))
# else
          rzeta(i,j)=zwrk(i,j)
          rzeta2(i,j)=zwrk(i,j)*zwrk(i,j)
# endif
        enddo
      enddo

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
#endif /* NBQ_MASS */
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------


# ifdef NBQ_DTDRHO
        if (iic==1.and.iif==1) then
         z_nbq(:,:,:,kstp)=z_nbq(:,:,:,knew)
         hrho_nbq(:,:,:,kstp)=hrho_nbq(:,:,:,knew)
        endif
        do k=1,N
          do j=jstrq_nh-1,jendq_nh+1
            do i=istrq_nh-1,iendq_nh+1
              hrho_nbq(i,j,k,knew)=Hzr_half_nbq(i,j,k)*rho(i,j,k)/rho0
            enddo  
          enddo  
        enddo
        do k=0,N
          do j=jstrq_nh-1,jendq_nh+1
            do i=istrq_nh-1,iendq_nh+1
              z_nbq(i,j,k,knew)=zw_half_nbq(i,j,k)
            enddo
          enddo  
        enddo
# endif /* NBQ_DTDRHO */
      

