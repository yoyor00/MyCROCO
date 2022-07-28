C
C
C***********************************************************************
C***********************************************************************
       subroutine DYNMIXLAY    
       Use comrunmod
       Use comdynmod
C
C       do jk=2,nzt-1
C         if (AVT(jk).lt.AVT(jk-1)) goto 11
C       enddo
C       xdml=DEPW(2)
C       return
CC
C 11    jk0=jk-1
C       do jk=jk0,nzt-1
C         if (AVT(jk).le.xlmkz) goto 13
C       enddo
C 13    if (jk.eq.jk0) then
C         xdml=DEPW(2)
C       else
C         xdml=DEPW(jk)
C       endif
C
        do jk=2,nzt-1
          if (EN(jk).lt.2.*emin) goto 11
        enddo
 11     xdml=DEPW(jk)
C
       return
       end
