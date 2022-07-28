C
C
C **********************************************************************
C **********************************************************************
       subroutine BIOINVERSE
C
C Appele par BIODYN: inversion d'une matrive triangulaire (resolution
C    implicate de la diffusion verticale)
C
C ======================================================================
C
       Use comrunmod
       Use comdynmod
       USE mod_varphy_coupl
       Implicit none
       Integer:: jk

C
C Initialisation pour le traceur
C
       do jk=1,nzt-1
         DZW(jk)=-dts*DAVT(jk)/DE3T(jk)/DE3W(jk)
         DZY(jk)=-dts*DAVT(jk+1)/DE3T(jk)/DE3W(jk+1)
         DZX(jk)=dun-DZW(jk)-DZY(jk)
       enddo
C 
C Traitement de surface: flux deja pris en compte
C
       DZW(1)=0.
       DZX(1)=dun-DZY(1)
C
C Traitement du fond: pas de flux
C
       DZY(nzt-1)=0.
       DZX(nzt-1)=dun-DZW(nzt-1)
C
C Inversion matrice tridiagonale
C
       DZR(1)=DZX(1)
       do jk=2,nzt-1
         DZR(jk)=DZX(jk)-DZW(jk)*DZY(jk-1)/DZR(jk-1)
       enddo

       DZX(1)=WORK1(1)
       do jk=2,nzt-1
         DZX(jk)=WORK1(jk)-DZW(jk)/DZR(jk-1)*DZX(jk-1)
       enddo
       DZW(nzt-1)=DZX(nzt-1)/DZR(nzt-1)
       do jk=nzt-2,1,-1
         DZW(jk)=(DZX(jk)-DZY(jk)*DZW(jk+1))/DZR(jk)
       enddo
C
C Update
C
       do jk=1,nzt-1
!MB         WORK(jk)=SNGL(DZW(jk))
         WORK(jk)=DZW(jk)
       enddo
C
       return
       end
