C
C
C***********************************************************************
C***********************************************************************
       subroutine DYN3DUP 
C
C ------------------------------------------------------
C Initialisation des moyennes apres sauvegarde effectuee
C ------------------------------------------------------
C Traceurs dynamiques
C
       Use comrunmod
       Use comdynmod
C       include 'combio'
C
       do jk=1,nzt
         TMOY(jk)=0.
       enddo
       do jk=1,nzt
         SMOY(jk)=0.
         DENSITE(jk)=0.
       enddo
       do jk=1,nzt
         TKEMOY(jk)=0.
       enddo
       do jk=1,nzt
         DIFFT(jk)=0.
       enddo
C
       return
       end
