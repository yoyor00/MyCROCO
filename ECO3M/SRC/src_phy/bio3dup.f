C
C
C***********************************************************************
C***********************************************************************
       subroutine BIO3DUP 
C
C ------------------------------------------------------
C Initialisation des moyennes apres sauvegarde effectuee
C ------------------------------------------------------
C Traceurs biologiques
C
C       include 'comrun'
C       include 'comdyn'
C       include 'combio'
       Use comrunmod
       Use comdynmod
       Use mod_varphy_coupl

C
       do jtr=1,nbrbio
         do jk=1,nzt
           BIOAVER(jk,jtr)=0.
         enddo
       enddo
C
       return
       end
