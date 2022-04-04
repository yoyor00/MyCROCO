! $Id: wetdry.F 1408 2013-12-20 12:41:14Z marchesiello $
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
! 
! CROCO website : http://www.croco-ocean.org
!
! This routine was adapted from John Warner's code
! by P. Marchesiello, IRD 2013
!======================================================================
!
#include "cppdefs.h"
#ifdef WAVE_DRY
!***********************************************************************
      SUBROUTINE wavedry_tile (Istr,Iend,Jstr,Jend)
!***********************************************************************
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k
# include "param.h"
# include "grid.h"
# include "scalars.h"
# include "ocean2d.h"
# include "coupling.h"
!
      real   cff, eps, cff1,cff2,cff3
      parameter (eps=1.e-10)
# ifdef MASKING
#  define SWITCH *
# else
#  define SWITCH !
# endif
!
#include "compute_auxiliary_bounds.h"
!
!-----------------------------------------------------------------------
! If wet/drying, compute new masks for cells with depth < Dcrit.
!-----------------------------------------------------------------------
!
# ifdef SOLVE3D
      IF (iif.le.nfast) THEN
# endif
!
!  Wet/dry mask at RHO-points.
!
        DO j=Jstr-1,JendR
          DO i=Istr-1,IendR
            wave_wetdry(i,j)=1.
# ifdef MASKING
            wave_wetdry(i,j)=wave_wetdry(i,j)*rmask(i,j)
# endif
            IF ((zeta(i,j,knew)+h(i,j)).le.Dcrit_wave(i,j)+eps) THEN
              wave_wetdry(i,j)=0.
            END IF
          END DO
        END DO

        DO j=JstrR,JendR
          DO i=IstrR,IendR
            rmask_wavewet(i,j)=wave_wetdry(i,j)
          END DO
        END DO
!
!  Wet/dry mask at PSI-points.
!
        DO j=Jstr,JendR
          DO i=Istr,IendR
            pmask_wavewet(i,j)=
     &             wave_wetdry(i-1,j  )*wave_wetdry(i  ,j  )*  
     &             wave_wetdry(i-1,j-1)*wave_wetdry(i  ,j-1)
          END DO
        END DO
!
!  Wet/dry mask at U-points.
!
        DO j=JstrR,JendR
          DO i=Istr,IendR
            umask_wavewet(i,j)=wave_wetdry(i-1,j)+wave_wetdry(i,j)
            IF (umask_wavewet(i,j).eq.1.) THEN
              umask_wavewet(i,j)=wave_wetdry(i-1,j)-wave_wetdry(i,j)
            END IF
          END DO
        END DO
!
!  Wet/dry mask at V-points.
!
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            vmask_wavewet(i,j)=wave_wetdry(i,j-1)+wave_wetdry(i,j)
            IF (vmask_wavewet(i,j).eq.1.) THEN
              vmask_wavewet(i,j)=wave_wetdry(i,j-1)-wave_wetdry(i,j)
            END IF
          END DO
        END DO
!
# ifdef SOLVE3D
      END IF
# endif
!
!--------------------------------------------------------------------
!  Averaged wet/dry mask for 3D slow mode if M2 filter
!--------------------------------------------------------------------
!
# if defined SOLVE3D && !defined M2FILTER_NONE
!
!  Wet/dry mask at RHO-points, averaged over all fast time-steps.
!
      IF (iif.le.nfast) THEN
        IF (FIRST_2D_STEP) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              rmask_wavewet_avg(i,j)=wave_wetdry(i,j)
            END DO
          END DO
        ELSE
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              rmask_wavewet_avg(i,j)=rmask_wavewet_avg(i,j)
     &                               +wave_wetdry(i,j)
            END DO
          END DO
        END IF
      END IF
!
!  If done fast time-stepping, scale mask by 2 nfast.
!
      IF (iif.eq.nfast+1) THEN
        cff=1./REAL(nfast)
        DO j=Jstr-1,JendR
          DO i=Istr-1,IendR
            wave_wetdry(i,j)=AINT(rmask_wavewet_avg(i,j)*cff)
          END DO
        END DO
!
!  Wet/dry mask at RHO-points, averaged over all fast time-steps.
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            rmask_wavewet(i,j)=wave_wetdry(i,j)
          END DO
        END DO
!
!  Wet/dry mask at PSI-points, averaged over all fast time-steps.
!
        DO j=Jstr,JendR
          DO i=Istr,IendR
            pmask_wavewet(i,j)=wave_wetdry(i-1,j  )
     &                        *wave_wetdry(i  ,j  )* 
     &                     wave_wetdry(i-1,j-1)*wave_wetdry(i  ,j-1)
          END DO
        END DO
!
!  Wet/dry mask at U-points, averaged over all fast time-steps.
!
        DO j=JstrR,JendR
          DO i=Istr,IendR
            cff1=wave_wetdry(i-1,j)+wave_wetdry(i,j)
            IF (cff1.eq.1.) THEN
              cff1=wave_wetdry(i-1,j)-wave_wetdry(i,j)
            END IF
            cff2=ABS(ABS(cff1)-1.)
            cff3=0.5+SIGN(0.5,DU_avg1(i,j,nnew))*cff1
            umask_wavewet(i,j)=0.5*cff1*cff2+cff3*(1.-cff2)
!                                                    catch lone ponds
            IF (DU_avg1(i,j,nnew).eq.0. .and.
     &          wave_wetdry(i-1,j)+wave_wetdry(i,j).le.1.) THEN
              umask_wavewet(i,j)=0.
            END IF
          END DO
        END DO
!
!  Wet/dry mask at V-points, averaged over all fast time-steps.
!
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            cff1=wave_wetdry(i,j-1)+wave_wetdry(i,j)
            IF (cff1.eq.1.) THEN
              cff1=wave_wetdry(i,j-1)-wave_wetdry(i,j)
            END IF
            cff2=ABS(ABS(cff1)-1.)
            cff3=0.5+SIGN(0.5,DV_avg1(i,j,nnew))*cff1
            vmask_wavewet(i,j)=0.5*cff1*cff2+cff3*(1.-cff2)
!                                                    catch lone ponds
            IF (DV_avg1(i,j,nnew).eq.0. .and.
     &          wave_wetdry(i,j-1)+wave_wetdry(i,j).le.1.) THEN
              vmask_wavewet(i,j)=0.
            END IF
          END DO
        END DO
!
      END IF
# endif
!
!  Exchange boundary data.
!
# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,  rmask_wavewet)
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,  rmask_wavewet_avg)
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,  umask_wavewet)
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,  vmask_wavewet)
      call exchange_p2d_tile (Istr,Iend,Jstr,Jend,  pmask_wavewet)
# endif
!
      RETURN
      END SUBROUTINE wavedry_tile
!     
#endif

