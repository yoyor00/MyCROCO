!======================================================================
! CROCO is a branch of ROMS developped at IRD, INRIA, 
! Ifremer, CNRS and Univ. Toulouse III  in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
!
! CROCO website : http://www.croco-ocean.org
!======================================================================
!
       subroutine biology_tile(Istr,Iend,Jstr,Jend)
!------------------------------------------------------------------
!
!   ROUTINE biology_pisces : PISCES MODEL
!   *************************************
!
!
!     PURPOSE.
!     --------
!          *ROMS_PISCES ECOSYSTEM MODEL FOR THE WHOLE OCEAN
!                       THIS ROUTINE COMPUTES INTERACTIONS
!                       BETWEEN THE DIFFERENT COMPARTMENTS OF THE
!                       MODEL
!----------------------------------------------------------------
       USE sms_pisces
       USE trcsms_pisces

C       implicit none
       INTEGER Istr,Iend,Jstr,Jend
       INTEGER i, j, jk, k, jn, ilc

#include "ocean2pisces.h90"

      CALL ocean_2_pisces( Istr,Iend,Jstr,Jend) 

      DO jn = 1, jptra
         DO jk = KRANGE
            DO j = JRANGE
               DO i = IRANGE          ! masked grid volume
                  trb(i,j,K,jn) 
     &           = MAX( 0., trb(i,j,K,jn) * 1.e-6 )
               ENDDO
            ENDDO
         ENDDO
      END DO
      DO jk = KRANGE
         DO j = JRANGE
            DO i = IRANGE          ! masked grid volume
                trb(i,j,K,jpno3) 
     &        = trb(i,j,K,jpno3) / rno3
# if ! defined key_pisces_light
                trb(i,j,K,jpnh4) 
     &        = trb(i,j,K,jpnh4) / rno3
                trb(i,j,K,jppo4) 
     &        = trb(i,j,K,jppo4) / po4r
#endif
            ENDDO
         ENDDO
       ENDDO

!       ilc=1+iic-ntstart   ! number of time step since restart
!       ! To be conform with clasic (no xios) croco file first
!       ! time step must be writen
!       if (ilc == 0) ilc=1
!       CALL iom_setkt( ilc )
       CALL trc_sms_pisces( iic )

      DO jn = 1, jptra
         DO jk = KRANGE
            DO j = JRANGE
               DO i = IRANGE          ! masked grid volume
                  trb(i,j,K,jn) 
     &          = trb(i,j,K,jn) * 1.e6
               ENDDO
            ENDDO
         ENDDO
      END DO
      DO jk = KRANGE
         DO j = JRANGE
            DO i = IRANGE          ! masked grid volume
               trb(i,j,K,jpno3)
     &       = trb(i,j,K,jpno3) * rno3
# if ! defined key_pisces_light
               trb(i,j,K,jpnh4) 
     &       = trb(i,j,K,jpnh4) * rno3
               trb(i,j,K,jppo4) 
     &       = trb(i,j,K,jppo4) * po4r
# endif
            ENDDO
         ENDDO
      ENDDO


      RETURN       
      END
