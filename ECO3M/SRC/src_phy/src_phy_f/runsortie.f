C
C
C***********************************************************************
C***********************************************************************
       subroutine RUNSORTIE
C
C Appele par main.f a chaque pas de temps: Moyenne et sauvegarde des
C   champs dynamiques et biologiques.
C   N'est effectif que tous les nsave pas de temps.
C
C ======================================================================
C
C ----------
C Sauvegarde
C ----------
C
C s'efectue tous les nsave pas de temps, ie nsa = nsave
C
       Use comrunmod
       Use comdynmod
       USE mod_varphy_coupl
C       include 'combio'
C
       nsa=nsa+1
       if (nsa.lt.nsave) goto 101
C
C Champs 1D
C
       call DYN1DUP
       call BIO1DUP
C
       nind=nind+1
C
C Champs 3D traceurs
C
       call DYN3DSAVE
C MB   call BIO3DSAVE
C
C Ecriture tous les nwrite*nsave pas de temps
C
       nwrite0=nwrite0+1
       if (mflagwrite.eq.1.and.nwrite0.eq.nwrite) then
         call WDYNPRINT
         call WBIOPRINT
         nwrite0=0
       endif
C
C Initialisation champs 3D traceurs apres ecriture
C
       call DYN3DUP
       call BIO3DUP
C
C Champs 3D tendances traceurs
C
       call DYNTEND
C MB       call BIOTEND
C
C Initialisation des indices apres sauvegarde
C
       nsa=0
       nind0=0
C
 101   continue
       return
       end
