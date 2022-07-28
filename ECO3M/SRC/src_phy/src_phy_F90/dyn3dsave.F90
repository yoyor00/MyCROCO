!
!
!***********************************************************************
!***********************************************************************
       subroutine DYN3DSAVE
!
! -------------------------------------------
! Sauvegarde sur nundyn des champs dynamiques
! -------------------------------------------
!
       Use comrunmod
       Use comdynmod
!       include 'combio'
!
! Temperature
!
       do jk=1,nzt
         TMOY(jk)=TMOY(jk)/nsave
       enddo
       if (mt.eq.1) write(nundyn) (TMOY(jk),jk=nzt,1,-1)
!
! Salinite - Densite
!
       do jk=1,nzt
         SMOY(jk)=SMOY(jk)/nsave
         DENSITE(jk)=DENSITE(jk)/nsave
       enddo
       if (ms.eq.1) then
         write(nundyn) (SMOY(jk),jk=nzt,1,-1)
         write(nundyn) (DENSITE(jk),jk=nzt,1,-1)
       endif
!
! TKE          
!
       do jk=1,nzt
         TKEMOY(jk)=ALOG10(AMAX1(TKEMOY(jk)/nsave,emin))
       enddo
       if (mtke.eq.1) write(nundyn) (TKEMOY(jk),jk=nzt,1,-1)
!
! Diffusion traceur
!
       do jk=1,nzt
         DIFFT(jk)=ALOG10(AMAX1(DIFFT(jk)/nsave,avtb))
       enddo
       if (mdifft.eq.1) write(nundyn) (DIFFT(jk),jk=nzt,1,-1)
!
       return
       end
