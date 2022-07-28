C
C
C***********************************************************************
C***********************************************************************
       subroutine BIO3DSAVE
C
C ------------------------------------------------
C Sauvegarde sur nunbio des champs biogeochimiques
C ------------------------------------------------
C
C       include 'comrun'
C       include 'comdyn'
C       include 'combio'
       Use comrunmod
       Use comdynmod
       Use mod_varphy_coupl
C
C Moyenne
C
       do jtr=1,nbrbio
         do jk=1,nzt
           BIOAVER(jk,jtr)=BIOAVER(jk,jtr)/nsave
         enddo
       enddo
CCCC!!!!!!!!!!!!!!!!!!!!!!!!
C PUR et PAR: moyenne sur periode eclairee de la journee
C
        do jk=1,nzt
          BIOAVER(jk,13)=BIOAVER(jk,13)*nsave/nind0
          BIOAVER(jk,14)=BIOAVER(jk,14)*nsave/nind0
        enddo
CCCC!!!!!!!!!!!!!!!!!!!!!!!!
C
C Sauvegarde sur nunbio
C
       do jtr=1,nbrbio
         if (MSTOCK(jtr).eq.1) then
           write(nunbio,*) (BIOAVER(jk,jtr),jk=nzt,1,-1)
         endif
       enddo
C
       return
       end
