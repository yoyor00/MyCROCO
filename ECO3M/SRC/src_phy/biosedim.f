C
C
C **********************************************************************
C **********************************************************************
       subroutine BIOSEDIM
C
C Appele par biologie.f a chaque passage: calcul des flux de sedimenta-
C tion
C
C ======================================================================
C
        Use comrunmod
        Use comdynmod
        Use mod_varphy_coupl

C        write(*,*) "Sed speed ", vitsed
C
       do jtr=1,jptract  ! barrier.n nbrprono -> jptract
         SEDFLUX(1,jtr)=-VITSED(jtr)*TENEUR(1,jtr)/DE3T(1)
c         do jk=2,nzt
         do jk=2,nzt-2
           SEDFLUX(jk,jtr)=VITSED(jtr)*(TENEUR(jk-1,jtr)-TENEUR(jk,jtr))
     &            /DE3T(jk)
         enddo
           SEDFLUX(nzt-1,jtr)=VITSED(jtr)*(TENEUR(nzt-2,jtr))
     &            /DE3T(nzt-1)
           SEDFLUX(nzt,jtr)=0.0
       enddo
       
C
       do jtr=1,jptract  ! barrier.n nbrprono -> jptract
         do jk=1,nzt
           SEDBIO(jk,jtr)=SEDBIO(jk,jtr)+SEDFLUX(jk,jtr)
         enddo
       enddo
C
       do jtr=1,jptract  ! barrier.n nbrprono -> jptract
         do jk=1,nzt-1
           TENEUR(jk,jtr)=TENEUR(jk,jtr)+dts*SEDFLUX(jk,jtr)
         enddo
       enddo
C
       return
       end
