        DO k=1,N
#  ifndef EW_PERIODIC
          if (WESTERN_EDGE) then
            do j=J_EXT_RANGE
              WORK(Istr-1,j,k)=WORK(Istr,j,k)
            enddo
          endif
          if (EASTERN_EDGE) then
            do j=J_EXT_RANGE
              WORK(Iend+1,j,k)=WORK(Iend,j,k)
            enddo
          endif
#  endif
#  ifndef NS_PERIODIC
          if (SOUTHERN_EDGE) then
            do i=I_EXT_RANGE
              WORK(i,Jstr-1,k)=WORK(i,Jstr,k)
            enddo
          endif
          if (NORTHERN_EDGE) then
            do i=I_EXT_RANGE
              WORK(i,Jend+1,k)=WORK(i,Jend,k)
            enddo
          endif
#   ifndef EW_PERIODIC
          if (WESTERN_EDGE.and.SOUTHERN_EDGE) then
            WORK(Istr-1,Jstr-1,k)=WORK(Istr,Jstr,k)
          endif
          if (WESTERN_EDGE.and.NORTHERN_EDGE) then
            WORK(Istr-1,Jend+1,k)=WORK(Istr,Jend,k)
          endif
          if (EASTERN_EDGE.and.SOUTHERN_EDGE) then
            WORK(Iend+1,Jstr-1,k)=WORK(Iend,Jstr,k)
          endif
          if (EASTERN_EDGE.and.NORTHERN_EDGE) then
            WORK(Iend+1,Jend+1,k)=WORK(Iend,Jend,k)
          endif
#   endif
#  endif
          DO j=Jstr-1,Jend+1
            DO i=Istr,Iend+1
              FX(i,j)=( WORK(i  ,j,k)
     &               -  WORK(i-1,j,k) ) SWITCH umask(i,j)
            ENDDO
          ENDDO
          DO j=Jstr,Jend+1
            DO i=Istr-1,Iend+1
              FE1(i,j)=( WORK(i,j  ,k)
     &                 - WORK(i,j-1,k) ) SWITCH vmask(i,j)
            ENDDO
            DO i=Istr,Iend
              FE(i,j)=FE1(i,j)
     &              + smoo_a*( FX(i+1,j)+FX(i  ,j-1)
     &                        -FX(i  ,j)-FX(i+1,j-1))
            ENDDO
          ENDDO
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
              FX(i,j)=FX(i,j)
     &              + smoo_a*( FE1(i,j+1)+FE1(i-1,j  )
     &                        -FE1(i,j  )-FE1(i-1,j+1))
            ENDDO
            DO i=Istr,Iend
              WORK(i,j,k)=WORK(i,j,k)
     &                   + smoo_b*( FX(i+1,j)-FX(i,j)
     &                             +FE(i,j+1)-FE(i,j) )
     &                                SWITCH rmask(i,j)
#  ifdef PERMEABILITY
     &                              *(1.-pena_r(i,j,k))
#  endif
            ENDDO
          ENDDO      !--> discard FE1
        ENDDO        !--> k loop

