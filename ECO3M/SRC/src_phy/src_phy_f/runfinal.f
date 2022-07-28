C
C
C***********************************************************************
C***********************************************************************
       subroutine RUNFINAL 
C
C Appele par main.f a la fin de la simulation: sauvegarde des champs 1d
C   d'evolution temporelle et ecruture finale de statistiques sur les
C   valeurs sous seuil des variables biologiques.
C
C ======================================================================
C
C ----------------------------------------------
C Sauvegarde et ecriture finales a la fin du run
C ----------------------------------------------
C
       Use comrunmod
       Use comdynmod
C       include 'combio'
C
C Sauvegarde des champs 2D dynamiques sur nundyn1d
C
       write(nundyn1d,1001) (float(jt),jt=1,nind-1)
       write(nundyn1d,1001) (-DEPTHML(jt),jt=1,nind-1)
       write(nundyn1d,1002) (FLUXVENTX(jt),jt=1,nind-1)
       write(nundyn1d,1002) (FLUXVENTY(jt),jt=1,nind-1)
       write(nundyn1d,1002) (FLUXVENT(jt),jt=1,nind-1)
       write(nundyn1d,1003) (FLUXSOL(jt),jt=1,nind-1)
       write(nundyn1d,1003) (FLUXNSOL(jt),jt=1,nind-1)
       write(nundyn1d,1003) (FLUXADVT(jt),jt=1,nind-1)
       write(nundyn1d,1003) (FLUXSURT(jt),jt=1,nind-1)
       write(nundyn1d,1003) (FLUXSAL(jt),jt=1,nind-1)
       write(nundyn1d,1003) (FLUXADVS(jt),jt=1,nind-1)
       write(nundyn1d,1003) (FLUXSURS(jt),jt=1,nind-1)
       write(nundyn1d,1003) (FLUXSURTKE(jt),jt=1,nind-1)
 1001  format(0p,10f12.1)
 1002  format(0p,10f12.3)
 1003  format(1p,10e12.4)
C
C Sauvegarde des champs 2D biologiques sur nunbio1d
C
C Lignes commentees par MB
C MB       write(nunbio1d,1001) (float(jt),jt=1,nind-1)
C MB       write(nunbio1d,1001) (-DEPTHZE(jt),jt=1,nind-1)
C MB       write(nunbio1d,1002) (XNEBUL(jt),jt=1,nind-1)
C MB       write(nunbio1d,1003) (XPAR0PLUS(jt),jt=1,nind-1)
C       do jtr=1,nbrbio
C         if (MSTOCK(jtr).eq.1) then
C           write(89,1003) (FLUXAIR(jt,jtr),jt=1,nind-1)
C         endif
C       enddo
C
C Ecriture des statistiques finales sur les parametres biologiques < 0
C durant l'integration
C
       call WRUNFIN
C
       return
       end
