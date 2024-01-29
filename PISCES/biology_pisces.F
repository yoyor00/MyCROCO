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
       USE oce_trc
       USE trc

# include "ocean2pisces.h90"
# include "do_loop_substitute.h90"

       implicit none
       INTEGER :: Istr,Iend,Jstr,Jend
       INTEGER :: ji, jk, jn
       INTEGER :: Kbb, Kmm
!----------------------------------------------------------------
!----------------------------------------------------------------


      Kbb = nnew     ;      Kmm = nnew     ;       Krhs = nrhs

      CALL ocean_2_pisces( Istr,Iend,Jstr,Jend) 

      DO jn = 1, jptra
         DO_3D( 0, 0, 0, 0, 1, jpk)
           tr(ji,jj,jk,jn,Kbb) 
     &  = MAX( 0., tr(ji,jj,jk,jn,Kbb) * 1.e-6 )
         END_3D
      END DO

      DO_3D( 0, 0, 0, 0, 1, jpk)
         tr(ji,jj,jk,jpno3,Kbb) 
     &   = tr(ji,jj,jk,jpno3,Kbb) / rno3
# if ! defined key_pisces_npzd
         tr(ji,jj,jk,jpnh4,Kbb) 
     &   = tr(ji,jj,jk,jpnh4,Kbb) / rno3
         tr(ji,jj,jk,jppo4,Kbb)
     &   = tr(ji,jj,jk,jppo4,Kbb) / po4r
#endif
      END_3D

!       ilc=1+iic-ntstart   ! number of time step since restart
!       ! To be conform with clasic (no xios) croco file first
!       ! time step must be writen
!       if (ilc == 0) ilc=1
!       CALL iom_setkt( ilc )
      DO jn = 1, jptra
         DO_3D( 0, 0, 0, 0, 1, jpk)
            tra(ji,jj,jk,jn)
     &    = tr(ji,jj,jk,jn,Krhs)
         END_3D
         !
         DO_3D( 0, 0, 0, 0, 1, jpk)
            tr(ji,jj,jk,jn,Krhs) = 0.
         END_3D
      ENDDO

      CALL trc_sms_pisces( iic, Kbb, Kmm, Krhs )


      DO jn = 1, jptra
         DO_3D( 0, 0, 0, 0, 1, jpk)
           tr(ji,jj,jk,jn,Kbb) 
     &   = tr(ji,jj,jk,jn,Kbb) * 1.e6 
         END_3D
      END DO

      DO_3D( 0, 0, 0, 0, 1, jpk)
         tr(ji,jj,jk,jpno3,Kbb)
     &   = tr(ji,jj,jk,jpno3,Kbb) * rno3
# if ! defined key_pisces_npzd
         tr(ji,jj,jk,jpnh4,Kbb) 
     &  = tr(ji,jj,jk,jpnh4,Kbb) * rno3
         tr(ji,jj,jk,jppo4,Kbb) 
     &  = tr(ji,jj,jk,jppo4,Kbb) * po4r
#endif
      END_3D
      !
      DO jn = 1, jptra
         DO_3D( 0, 0, 0, 0, 1, jpk)
            tr(ji,jj,jk,jn,Krhs)
     &    = tra(ji,jj,jk,jn)
        END_3D
      ENDDO

      RETURN       
      END
