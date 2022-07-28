C
C
C***********************************************************************
C***********************************************************************
       subroutine BIOSUM
C       include 'comrun'
C       include 'comdyn'
C       include 'combio'
       Use comrunmod
       Use comdynmod
       Use mod_varphy_coupl
C
       do jk=1,nzt
         do jtr=1,nbrbio
           BIOAVER(jk,jtr)=BIOAVER(jk,jtr)+TENEUR(jk,jtr)
         enddo
       enddo
C
       return
       end
