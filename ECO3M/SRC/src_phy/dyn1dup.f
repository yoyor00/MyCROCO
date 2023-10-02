C
C
C***********************************************************************
C***********************************************************************
       subroutine DYN1DUP
C
C -------------------------------------------------
C Moyenne et stockage des champs 1D de la dynamique
C -------------------------------------------------
C nind = indice temporel
C
       Use comrunmod
       Use comdynmod
C       include 'combio'
C
       DEPTHML(nind)=DEPTHML(nind)/nsave
       FLUXSOL(nind)=FLUXSOL(nind)/nsave
       FLUXNSOL(nind)=FLUXNSOL(nind)/nsave
       FLUXADVT(nind)=FLUXADVT(nind)/nsave
       FLUXSAL(nind)=FLUXSAL(nind)/nsave
       FLUXADVS(nind)=FLUXADVS(nind)/nsave
       FLUXVENTX(nind)=FLUXVENTX(nind)/nsave
       FLUXVENTY(nind)=FLUXVENTY(nind)/nsave
       FLUXVENT(nind)=FLUXVENT(nind)/nsave
       FLUXSURT(nind)=FLUXSURT(nind)/nsave
       FLUXSURS(nind)=FLUXSURS(nind)/nsave
       FLUXSURTKE(nind)=FLUXSURTKE(nind)/nsave
c
       return
       end