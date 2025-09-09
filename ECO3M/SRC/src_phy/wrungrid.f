C     
C=======================================================================
C     
      subroutine WRUNGRID
C     
      Use comrunmod
      Use comdynmod
C     
      do jk=1,nzt-1,2
         write(99,1005) jk,DEPW(jk),DEPT(jk),E3W(jk),E3T(jk)
      enddo
c-- Ajout MB
c      do jk=1,nzt-1
c         write(199,1006) jk,DEPW(jk),DEPT(jk),E3W(jk),E3T(jk),DE3T(jk)
c      enddo
 1005 format(i5,4f8.1)  
c-- Ajout MB
c 1006 format(i5,2X,4(f10.5,1X),E16.10)
C     
      return
      end
