C
C ----------
C biologie.f
C ----------
C 27/02/97: Version 1 provenant de biotke.f (lm)
C
C Appele a chaque pas de temps par main.f.
C
C Programme en quatre parties:
C   1. Calcul des profils optiques PAR et PUR, ainsi que des moyennes
C           dans couche melangee, et couche euphotique
C            (subroutines dans optique.f)
C           PARPLUS = calcul de PAR0+ par ciel clair
C           PURPROF = calculs des profils de PAR et de PUR
C   2. Recherche et calage temporel des profils verticaux de rappels pour
C          les variables biogeochimiques (dissoutes) = BIORAPPEL
C   3. Calcul des flux a l'interface air - mer = BIOSURFACE
C   4. Estimation des termes non conservatifs biogeochimiques =
C                BIOMASSE
C
C ======================================================================

       subroutine BIOLOGIE
       Use comrunmod
       Use comdynmod
       Use mod_eco3m
       USE mod_varphy_coupl

C Calcul du PAR0+ et PAR0-
C       call OPTPARPLUS
C	write(*,*) "biologie xpar0p ",xpar0p

C Calcul du profil vertical de PAR et de PUR    
C       npurday0 = npurday0 + 1
C       if (xpar0p.gt.0.) then
C         call OPTPURPROF
C         DEPTHZE(nind)=DEPTHZE(nind)+seuph
C         XNEBUL(nind)=XNEBUL(nind)+trpar
C         XPAR0PLUS(nind)=XPAR0PLUS(nind)+xpar0p
C         nind0=nind0+1
C         npurd0=npurd0+1
C         fzede0=fzede0+seuph
C         fmlde0=fmlde0+xdml
C         fmlpur0=fmlpur0+fpurml
C        fzepur0=fzepur0+fpurze
C       endif

C Flux a la surface
       call BIOSURFACE

C Calcul des flux biogeochimiques
C       call BIOMASSE
       call eco3m_conc_update

       
       if (npurday0.eq.npurday) then
         npurday0=0
         npurd0=0
       endif

       return
       end
