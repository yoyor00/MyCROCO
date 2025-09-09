C
C
C***********************************************************************
C***********************************************************************
       subroutine DYNRAPPEL
       Use comrunmod
       Use comdynmod
C
       if (dayflux0.ge.dayts0.and.dayflux0.lt.dayts1) return
 102   continue
       read(95,*,end=101) dayts0,dayts1
       do jk=1,nzt
         read(95,*) XTS(jk),TRP(jk),SRP(jk)
	TRP(jk) = TRP(jk) 
       enddo
       return
C
 101   rewind(95)
       goto 102
C
       end
