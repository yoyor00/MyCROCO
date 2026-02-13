!
!===============================================================
!
! Compute 5th order horizontal advection
!
!===============================================================
!
#ifdef NS_PERIODIC
      jmin=1
      jmax=LOCALMM+1
#else
# ifdef MPI
      if (SOUTH_INTER) then
        jmin=1
      else
        jmin=3
      endif
      if (NORTH_INTER) then
        jmax=Mmmpi+1
      else
        jmax=Mmmpi-1
      endif
# else
      jmin=3
      jmax=Mm-1
# endif
#endif
#ifdef EW_PERIODIC
      imin=1
      imax=LOCALLM+1
#else
# ifdef MPI
      if (WEST_INTER) then
        imin=1
      else
        imin=3
      endif
      if (EAST_INTER) then
        imax=Lmmpi+1
      else
        imax=Lmmpi-1
      endif
# else
      imin=3
      imax=Lm-1
# endif
#endif
!
!--------------------------------------------------------------------
!  j loop: FE
!--------------------------------------------------------------------
!
      DO j = Jstr,Jend+1
                                                !
        IF ( j.ge.jmin .and. j.le.jmax ) THEN   ! full stencil
                                                !
          DO i = Istr,Iend+1
# ifdef MRL_CEW
            cev=0.5*om_v(i,j)*((cge(i,j)+cge(i,j-1))+vwave(i,j,wmod))
# else
            cev=0.5*om_v(i,j)*( cge(i,j)+cge(i,j-1) )
# endif
            flx5=cev*FLUX5(VAR(i,j-3), VAR(i,j-2),
     &                     VAR(i,j-1), VAR(i,j  ),
     &                     VAR(i,j+1), VAR(i,j+2), cev)
#ifdef MASKING
            flx3=cev*FLUX3(VAR(i,j-2), VAR(i,j-1),
     &                     VAR(i,j  ), VAR(i,j+1), cev)
            flx2=cev*FLUX2(VAR(i,j-1), VAR(i,j  ), cev)
            mask1=rmask(i,j-2)*rmask(i,j+1)
            mask2=rmask(i,j-3)*rmask(i,j+2)
            mask0=mask1*mask2
            FE(i,j)=mask0*flx5+(1-mask0)*mask1*flx3+
     &                         (1-mask0)*(1-mask1)*flx2
#else
            FE(i,j)=flx5
#endif
          ENDDO

        ELSE IF ( (j.eq.jmin-1 .and. jmax.ge.jmin)       ! 3rd order flux 
     &                           .or. j.eq.jmax+1) THEN  ! 2 points in from
                                                         ! south/north borders
          DO i = Istr,Iend+1
# ifdef MRL_CEW
            cev=0.5*om_v(i,j)*((cge(i,j)+cge(i,j-1))+vwave(i,j,wmod))
# else
            cev=0.5*om_v(i,j)*( cge(i,j)+cge(i,j-1) )
# endif
            flx3=cev*FLUX3(VAR(i,j-2),VAR(i,j-1),
     &                     VAR(i,j  ),VAR(i,j+1), cev)
# ifdef MASKING
            flx2=cev*FLUX2(VAR(i,j-1),VAR(i,j  ), cev)
            mask1=rmask(i,j-2)*rmask(i,j+1)
            FE(i,j)=mask1*flx3+(1-mask1)*flx2
# else
            FE(i,j)=flx3
# endif
          ENDDO

        ELSE IF ( j.eq.jmin-2 .or. j.eq.jmax+2) THEN  ! 1rst order flux next 
                                                      ! to south/north borders
          DO i = Istr,Iend+1
# ifdef MRL_CEW
            cev=0.5*om_v(i,j)*((cge(i,j)+cge(i,j-1))+vwave(i,j,wmod))
# else
            cev=0.5*om_v(i,j)*( cge(i,j)+cge(i,j-1) )
# endif
            FE(i,j)=cev*FLUX2(VAR(i,j-1),VAR(i,j), cev)
          ENDDO

        ENDIF
      ENDDO   ! j loop
!
!--------------------------------------------------------------------
!  i loop: FX
!--------------------------------------------------------------------
!
      DO i = Istr,Iend+1
                                              !
        IF ( i.ge.imin .and. i.le.imax ) THEN ! use full stencil
                                              !
          DO j = Jstr,Jend+1
# ifdef MRL_CEW
            cxu=0.5*on_u(i,j)*((cgx(i,j)+cgx(i-1,j))+uwave(i,j,wmod))
# else
            cxu=0.5*on_u(i,j)*( cgx(i,j)+cgx(i-1,j) )
# endif
            flx5=cxu*FLUX5(VAR(i-3,j), VAR(i-2,j),
     &                     VAR(i-1,j), VAR(i  ,j),
     &                     VAR(i+1,j), VAR(i+2,j), cxu)
# ifdef MASKING
            flx3=cxu*FLUX3(VAR(i-2,j), VAR(i-1,j),
     &                     VAR(i  ,j), VAR(i+1,j), cxu)
            flx2=cxu*FLUX2(VAR(i-1,j), VAR(i  ,j), cxu)
            mask1=rmask(i-2,j)*rmask(i+1,j)
            mask2=rmask(i-3,j)*rmask(i+2,j)
            mask0=mask1*mask2
            FX(i,j)=mask0*flx5+(1-mask0)*mask1*flx3+
     &                         (1-mask0)*(1-mask1)*flx2
# else
            FX(i,j)=flx5
# endif
          ENDDO

        ELSE IF ( (i.eq.imin-1 .and. imax.ge.imin)      ! 3rd order flux
     &                           .or. i.eq.imax+1) THEN ! 2 points in from
                                                        ! west/east borders
          DO j = Jstr,Jend+1
# ifdef MRL_CEW
            cxu=0.5*on_u(i,j)*((cgx(i,j)+cgx(i-1,j))+uwave(i,j,wmod))
# else
            cxu=0.5*on_u(i,j)*( cgx(i,j)+cgx(i-1,j) )
# endif
            flx3=cxu*FLUX3(VAR(i-2,j),VAR(i-1,j),
     &                     VAR(i  ,j),VAR(i+1,j), cxu)
# ifdef MASKING
            flx2=cxu*FLUX2(VAR(i-1,j),VAR(i  ,j), cxu)
            mask1=rmask(i-2,j)*rmask(i+1,j)
            FX(i,j)=mask1*flx3+(1-mask1)*flx2
# else
            FX(i,j)=flx3
# endif
          ENDDO

        ELSE IF ( i.eq.imin-2 .or. i.eq.imax+2 ) THEN  ! 1rst order flux next
                                                       ! to east/west borders

          DO j = Jstr,Jend+1
# ifdef MRL_CEW
            cxu=0.5*on_u(i,j)*((cgx(i,j)+cgx(i-1,j))+uwave(i,j,wmod))
# else
            cxu=0.5*on_u(i,j)*( cgx(i,j)+cgx(i-1,j) )
# endif
            FX(i,j)=cxu*FLUX2(VAR(i-1,j),VAR(i,j), cxu)
          ENDDO

        ENDIF
      ENDDO ! i loop

