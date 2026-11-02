  FUNCTION tool_julien(jou,moi,ia)

   !&E---------------------------------------------------------------------
   !&E                 ***  FUNCTION tool_julien  ***
   !&E
   !&E ** Purpose : compute the Julian day number from the date entered as 
   !&E              function arguments in the form of day, month and year
   !&E
   !&E ** Description :
   !&E
   !&E ** Called by : paramet
   !&E
   !&E ** External calls : 
   !&E
   !&E ** Used ij-arrays : 
   !&E
   !&E ** Modified variables : tool_julien
   !&E
   !&E ** Reference :
   !&E
   !&E ** History :
   !&E       !  2004-08-18
   !&E
   !&E---------------------------------------------------------------------
   !! * Modules used
   !USE parameters

   IMPLICIT NONE
   INTEGER,PARAMETER     :: rsh=8

   !! * Declarations function
   INTEGER              :: tool_julien

   !! * Arguments
   INTEGER, INTENT(in)  :: ia,jou,moi

   !! * Local declarations
   INTEGER,DIMENSION(12,2),PARAMETER :: jo=                            &
      RESHAPE((/0,31,59,90,120,151,181,212,243,273,304,334             &
               ,0,31,60,91,121,152,182,213,244,274,305,335/),(/12,2/))
   INTEGER                           :: i,ibs,m,iy,j,mois,jour
   REAL(kind=rsh)                    :: a,b

   !!----------------------------------------------------------------------
   !! * Executable part

   mois=moi
   jour=jou
   IF(mois == 1) THEN
     ibs=1
     IF (ia < 1582 .OR. (ia == 1582 .AND. jou < 283)) GOTO 10
     ibs=MOD(ia,4)+2
     IF(ibs /= 2) ibs=1
10 CONTINUE

     DO i=1,12
       mois=12-i+1
       jour=jou-jo(mois,ibs)
       IF(jour > 0) GOTO 20
     ENDDO
   ENDIF

20 CONTINUE

   a=REAL(ia,rsh)+REAL(mois,rsh)/100.0_rsh+REAL(jour,rsh)/10000.0_rsh
   b=0.0_rsh
   iy=ia
   m=mois
   IF(mois <= 2) THEN
     iy=iy-1
     m=m+12
   ENDIF

   IF(a >= 1582.1015_rsh) THEN
     j=iy/100
     b=2.0_rsh-REAL(j,rsh)+REAL(j,rsh)/4.0_rsh
   ENDIF

   tool_julien=INT(365.25_rsh*iy)+INT(30.6001_rsh*(m+1))+jour+1720995+INT(b)

  END FUNCTION tool_julien
