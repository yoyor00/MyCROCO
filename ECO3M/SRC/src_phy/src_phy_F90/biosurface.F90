!
!
!***********************************************************************
!***********************************************************************
       subroutine BIOSURFACE
       Use comrunmod
       Use comdynmod
       USE mod_varphy_coupl
!
       do jtr=1,nbrbio
         BIOAIR(jtr)=0.
       enddo
!
!       call BIOFLUXPCO2
!
       return
       end
