!======================================================================
! CROCO is derived from the ROMS-AGRIF branch of ROMS.
! ROMS-AGRIF was developed by IRD and Inria. CROCO also inherits
! from the UCLA branch (Shchepetkin et al.) and the Rutgers
! University branch (Arango et al.), both under MIT/X style license.
! Copyright (C) 2005-2026 CROCO Development Team
! License: CeCILL-2.1 - see LICENSE.txt
!
! CROCO website : https://www.croco-ocean.org
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
       USE oce_trc
       USE trc

# include "ocean2pisces.h90"

       implicit none
       INTEGER, intent(in)  :: Istr,Iend,Jstr,Jend
       INTEGER :: i, j, k, jn
       INTEGER :: Kbb, Kmm
!----------------------------------------------------------------
!----------------------------------------------------------------


      Kbb = nnew     ;      Kmm = nnew     ;       Krhs = nrhs

      CALL ocean_2_pisces( iic, Istr,Iend,Jstr,Jend) 

      DO jn = 1, jptra
         DO k= 1, jpk ; DO j=Jstrp,Jendp ; DO i=Istrp,Iendp
           tr(i,j,k,jn,Kbb) 
     &  = MAX( 0., tr(i,j,k,jn,Kbb) * 1.e-6 )
         END DO   ;   END DO   ;   END DO         
      END DO

      DO k= 1, jpk ; DO j=Jstrp,Jendp ; DO i=Istrp,Iendp
         tr(i,j,k,jpno3,Kbb) 
     &   = tr(i,j,k,jpno3,Kbb) / rno3
# if ! defined key_pisces_npzd
         tr(i,j,k,jpnh4,Kbb) 
     &   = tr(i,j,k,jpnh4,Kbb) / rno3
         tr(i,j,k,jppo4,Kbb)
     &   = tr(i,j,k,jppo4,Kbb) / po4r
#endif
      END DO   ;   END DO   ;   END DO         

      DO jn = 1, jptra
         DO k= 1, jpk ; DO j=Jstrp,Jendp ; DO i=Istrp,Iendp
             tra(i,j,k,jn)
     &       = tr(i,j,k,jn,Krhs)
         END DO   ;   END DO   ;   END DO         
         !
         DO k= 1, jpk ; DO j=Jstrp,Jendp ; DO i=Istrp,Iendp
            tr(i,j,k,jn,Krhs) = 0.
         END DO   ;   END DO   ;   END DO         
       ENDDO

      CALL trc_sms_pisces( iic, Kbb, Kmm, Krhs )


      DO jn = 1, jptra
         DO k= 1, jpk ; DO j=Jstrp,Jendp ; DO i=Istrp,Iendp
           tr(i,j,k,jn,Kbb) 
     &   = tr(i,j,k,jn,Kbb) * 1.e6 
         END DO   ;   END DO   ;   END DO         
      END DO

      DO k= 1, jpk ; DO j=Jstrp,Jendp ; DO i=Istrp,Iendp
         tr(i,j,k,jpno3,Kbb)
     &   = tr(i,j,k,jpno3,Kbb) * rno3
# if ! defined key_pisces_npzd
         tr(i,j,k,jpnh4,Kbb) 
     &  = tr(i,j,k,jpnh4,Kbb) * rno3
         tr(i,j,k,jppo4,Kbb) 
     &  = tr(i,j,k,jppo4,Kbb) * po4r
#endif
      END DO   ;   END DO   ;   END DO         
      !
      DO jn = 1, jptra
         DO k= 1, jpk ; DO j=Jstrp,Jendp ; DO i=Istrp,Iendp
            tr(i,j,k,jn,Krhs)
     &    = tra(i,j,k,jn)
         END DO   ;   END DO   ;   END DO         
      ENDDO

      RETURN       
      END
