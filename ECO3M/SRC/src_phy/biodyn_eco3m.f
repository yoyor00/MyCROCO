C
C
C **********************************************************************
C **********************************************************************
       subroutine BIODYN 
C
C Appele par biologie.f a chaque passage: calcul des tendances dynamiques
C des variables biologiques et update des variables
C
C ======================================================================
C
       Use comrunmod
       Use comdynmod
       USE mod_varphy_coupl


       do jtr=1,jptract  ! barrier.n -> replacement of nbrprono by jptract
         workb=BIOAIR(jtr)*dts
         do jz=1,nzt
           WORK(jz)=TENEUR(jz,jtr)
           WORK1(jz)=TENEUR(jz,jtr)
         enddo
         WORK(1)=WORK(1)+workb
         WORK1(1)=WORK1(1)+workb


         call BIOINVERSE

         do jz=1,nzt
           TENEUR(jz,jtr)=WORK(jz)
         enddo
         FLUXAIR(nind,jtr)=FLUXAIR(nind,jtr)+BIOAIR(jtr)
         do jz=1,nzt
           DIFFBIO(jz,jtr)=DIFFBIO(jz,jtr)+(WORK(jz)-WORK1(jz))/dts
         enddo
       enddo

C
C Diagnostique sur rapport C/Chl
C
C MB       if (nchlapro.eq.1) then
C MB         do jz=1,nzt
C MB           if (TENEUR(jz,niparchl).gt.teneurmin) then
C MB             TENEUR(jz,nicchl)=carbmass*REDCN(niphy)*TENEUR(jz,niphy)/
C MB     &            TENEUR(jz,niparchl)
C MB           else
C MB             TENEUR(jz,nicchl)=cchlmin
C MB           endif
C MB         enddo
C MB       endif
C
C Stockage en carbone
C
C       VPCO2(nind)=VPCO2(nind)+xpco2
C       VITTRANS(nind)=VITTRANS(nind)+vitgaz
C
       return
       end
