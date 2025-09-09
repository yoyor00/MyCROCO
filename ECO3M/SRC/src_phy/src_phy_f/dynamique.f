C
C -----------
C dynamique.f
C -----------
C 27/02/97: Version 1 provenant de biotke.f (lm)
C
C Appele a chaque pas de temps par main.f
C
C Programme en quatre parties:
C   1. Diagnostique sur stratification, longueur de melange et diffusion
C              verticale.
C   2. Equation pronostique de TKE
C   3. Equation pronostique de U et V (couche d'Ekman)
C   4. Equation pronostique pour T et S
C
C ======================================================================
C
       subroutine DYNAMIQUE
       Use comrunmod
       Use comdynmod
C
C 1111111111111111111111111111111111111111
C
C Calcul de la diffusion verticale
C
       call DYNKZ
C
C Calcul de la profondeur de la couche melangee
C
       call DYNMIXLAY
       DEPTHML(nind)=DEPTHML(nind)+xdml
C
C 1111111111111111111111111111111111111111
C
C 2222222222222222222222222222222222222222
C
C Lecture des forcages atmospheriques
C
       call DYNFLUXSUR
       FLUXVENTX(nind)=FLUXVENTX(nind)+windx
       FLUXVENTY(nind)=FLUXVENTY(nind)+windy
       FLUXVENT(nind)=FLUXVENT(nind)+wind
C
C Integration de TKE
C
       call DYNTKE
C
C 2222222222222222222222222222222222222222
C
C 3333333333333333333333333333333333333333
C
C Integration des vitesses
C
       call DYNVIT
C
C 3333333333333333333333333333333333333333
C
C 4444444444444444444444444444444444444444
C
C Integration temperature, salinite
C
       call DYNTS
C
C 4444444444444444444444444444444444444444
       
       return
       end
