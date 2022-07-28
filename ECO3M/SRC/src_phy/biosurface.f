C
C
C***********************************************************************
C***********************************************************************
       subroutine BIOSURFACE
       Use comrunmod
       Use comdynmod
       USE mod_varphy_coupl
C
       do jtr=1,nbrbio
         BIOAIR(jtr)=0.
       enddo
C
C       call BIOFLUXPCO2
C
       return
       end
