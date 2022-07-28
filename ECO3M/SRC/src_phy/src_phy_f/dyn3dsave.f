C
C
C***********************************************************************
C***********************************************************************
       subroutine DYN3DSAVE
C
C -------------------------------------------
C Sauvegarde sur nundyn des champs dynamiques
C -------------------------------------------
C
       Use comrunmod
       Use comdynmod
C       include 'combio'
C
C Temperature
C
       do jk=1,nzt
         TMOY(jk)=TMOY(jk)/nsave
       enddo
       if (mt.eq.1) write(nundyn) (TMOY(jk),jk=nzt,1,-1)
C
C Salinite - Densite
C
       do jk=1,nzt
         SMOY(jk)=SMOY(jk)/nsave
         DENSITE(jk)=DENSITE(jk)/nsave
       enddo
       if (ms.eq.1) then
         write(nundyn) (SMOY(jk),jk=nzt,1,-1)
         write(nundyn) (DENSITE(jk),jk=nzt,1,-1)
       endif
C
C TKE          
C
       do jk=1,nzt
         TKEMOY(jk)=ALOG10(AMAX1(TKEMOY(jk)/nsave,emin))
       enddo
       if (mtke.eq.1) write(nundyn) (TKEMOY(jk),jk=nzt,1,-1)
C
C Diffusion traceur
C
       do jk=1,nzt
         DIFFT(jk)=ALOG10(AMAX1(DIFFT(jk)/nsave,avtb))
       enddo
       if (mdifft.eq.1) write(nundyn) (DIFFT(jk),jk=nzt,1,-1)
C
       return
       end
