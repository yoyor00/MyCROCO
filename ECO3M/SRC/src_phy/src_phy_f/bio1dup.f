C     
C
C***********************************************************************
C***********************************************************************
       subroutine BIO1DUP
C
C ------------------------------------------------
C Moyenne et stockage des champs 1D de la biologie
C ------------------------------------------------
C nind = indice temporel
C
C       include 'comrun'
C       include 'comdyn'
        Use comrunmod
        Use comdynmod
        Use mod_varphy_coupl
C       include 'combio'
C

       if (nind0.ne.0) then
         DEPTHZE(nind)=DEPTHZE(nind)/nind0
         XNEBUL(nind)=XNEBUL(nind)/nind0
         XPAR0PLUS(nind)=XPAR0PLUS(nind)/nind0
       endif
       do jtr=1,jptract  ! barrier.n nbrprono -> jptract
         FLUXAIR(nind,jtr)=FLUXAIR(nind,jtr)/nsave
       enddo
C       VPCO2(nind)=VPCO2(nind)/nsave
C       VITTRANS(nind)=100.*hour*VITTRANS(nind)/nsave
C
       return
       end
