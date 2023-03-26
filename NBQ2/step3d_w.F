! $Id:$
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
! 
! CROCO website : http://www.croco-ocean.org
!======================================================================
!
#include "cppdefs.h"
#if defined SOLVE3D && defined M3SLOW_W

# define NBQ_WDIF

      SUBROUTINE step3d_w (tile)

!======================================================================
!                       ***  Subroutine STEP3D_W  ***
! NBQ mode :  Advance the vertical velocity wz to time n+1 
!======================================================================
! History :   2016-11  (F. Lemarié)  Original code
!----------------------------------------------------------------------
      
      IMPLICIT NONE
      integer tile, trd, omp_get_thread_num
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      trd=omp_get_thread_num()

      CALL step3d_w_tile (Istr,Iend,Jstr,Jend, A3d(1,3,trd), !<-- A3d(1,3,trd) contains rw
     &                           A2d(1,1,trd), A2d(1,2,trd),
     &                           A2d(1,3,trd), A2d(1,4,trd)  )

      return
      END
!
!----------------------------------------------------------------------        
      SUBROUTINE step3d_w_tile (Istr,Iend,Jstr,Jend, rw, BC,CF,FC,DC)
!----------------------------------------------------------------------
!     
      IMPLICIT NONE
# include "param.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "coupling.h"
# include "forces.h"
# include "mixing.h"
# include "scalars.h"
# include "sources.h"
# if defined M3NUDGING && defined M3CLIMATOLOGY
#  include "climat.h"
# endif
# ifdef DIAGNOSTICS_UV
#  include "diagnostics.h"
# endif
# ifdef M3FAST_W
#  include "nbq.h"
# endif
      INTEGER Istr,Iend,Jstr,Jend, i,j,k,k1
# ifdef PSOURCE
     &       ,is
# endif
      REAL  ::  BC   (PRIVATE_1D_SCRATCH_ARRAY,0:N  )
      REAL  ::  CF   (PRIVATE_1D_SCRATCH_ARRAY,0:N  )
      REAL  ::  FC   (PRIVATE_1D_SCRATCH_ARRAY,0:N+1)
      REAL  ::  DC   (PRIVATE_1D_SCRATCH_ARRAY,0:N  )
      REAL  ::  rw   (PRIVATE_2D_SCRATCH_ARRAY,0:N  )       
      REAL  ::  grad (PRIVATE_2D_SCRATCH_ARRAY      )
      REAL  ::  dpth,cff,cff1,cff2
      REAL  ::  my_cadv_max   
      INTEGER :: jzi,jzj,jzk 
      !
# include "compute_auxiliary_bounds.h"
      !
# ifdef MASKING
#  define SWITCH *
# else
#  define SWITCH !
# endif

!$acc kernels if(compute_on_device) default(present)
!
!---------------------------------------------------------------------- 
!  Apply right-hand-side
!---------------------------------------------------------------------- 
!
      DO j=Jstr,Jend
        DO k=1,N
          DO i=Istr,Iend
            wz(i,j,k,nnew)=wz(i,j,k,nnew)
     &                +dt*pn(i,j)*pm(i,j)*(rw(i,j,k)
# ifdef M3FAST_C3D_WFS
     &                +rw_nbq_avg2(i,j,k)
# endif   
     &                                    )
          ENDDO
        ENDDO
      ENDDO  

! ! Bottom boundary condition:
      DO j=Jstr,Jend
        DO i=Istr,Iend
# ifdef NBQ_FREESLIP 
#  ifndef M3FAST_SEDLAYERS
          wz(i,j,0,nnew)=qdmw0_nbq(i,j)  
#  else
! !  Another BC is possible based on Slow-Mode Kinematic relation
           wz(i,j,0,nnew)= 0.25*(    
     &                 u(i,j,1,nnew) * pm_u(i,j) ! Caution: u is in m/s
     &                * ( z_w(i  ,j,0)-z_w(i-1,j,0) )
     &                +u(i+1,j,1,nnew) * pm_u(i+1,j) 
     &                * ( z_w(i+1,j,0)-z_w(i,j,0) )
     &                + v(i,j,1,nnew) * pm_v(i,j) 
     &                * ( z_w(i,j  ,0)-z_w(i,j-1,0) )
     &                + v(i,j+1,1,nnew) * pm_v(i,j+1) 
     &                * ( z_w(i,j+1,0)-z_w(i,j,0) ))
     &                    *Hz(i,j,1)
! !
! ! Association with fast mode
! !
! !        wz(i,j,0,nnew)=qdmw_nbq(i,j,0)
#  endif
# else
           wz(i,j,0,nnew)=0.
! !
! ! Association with fast mode
! !
! !        wz(i,j,0,nnew)=qdmw_nbq(i,j,0)
# endif          
        ENDDO
      ENDDO  
  
# ifdef NBQ_WDIF
!
!====================================================================     
!== Turbulent vertical diffusion for wz
!== At this point wz contains rho.Hz.wz with units kg.m-1.s-1  
!== The implicit integration contains a division by rho.Hz
!==================================================================== 
!
      DO j=Jstr,Jend
      
        ! Off-diagonal terms
        DO k=1,N
          DO i=Istr,Iend
            FC(i,k) = -0.5*dt*Hz(i,j,k)*( Akv(i,j,k)+Akv(i,j,k-1) )
     &                             /    ( z_w(i,j,k)-z_w(i,j,k-1) )**2
          ENDDO
        ENDDO
        DO i=Istr,Iend
          FC(i,0  )=0.
          FC(i,N+1)=0.
        ENDDO
         
        ! Rhs and diagonal term         
        DO k=1,N-1
          DO i=Istr,Iend  
            DC(i,k)=wz(i,j,k,nnew)                                    !<-- rhs
            BC(i,k)=0.5*(Hz(i,j,k)+Hz(i,j,k+1)) - FC(i,k+1) - FC(i,k) !<-- diagonal term                
          ENDDO
        ENDDO  
        DO i=Istr,Iend            
          BC(i,N)=0.5*Hz(i,j,N) - FC(i,N+1) - FC(i,N)    !<-- Hz contains rho.Hz
          DC(i,N)=wz(i,j,N,nnew)
          BC(i,0)=1.  !0.5*Hz(i,j,1)-FC(i,1)-FC(i,0)
#  ifdef NBQ_FREESLIP
          DC(i,0)=wz(i,j,0,nnew)  
#  else    
          DC(i,0)=0.
#  endif                 
        ENDDO 
        
        ! Gaussian elimination
        DO i=Istr,Iend
          cff     = 1.                      !<-- 1./b(0)
          CF(i,0) = 0.                      !<-- q(0) = c(0) / b(0)
          DC(i,0) = cff*DC(i,0)             !<-- f(0) = f(0) / b(0)
        ENDDO
#ifdef M3FAST_BOTH
        
        DO k=1,N-1
          DO i=Istr,Iend 
            cff     = 1./(BC(i,k)-FC(i,k)*CF(i,k-1))
            CF(i,k) = cff* FC(i,k+1) 
            DC(i,k) = cff*(DC(i,k)-FC(i,k)*DC(i,k-1))
          ENDDO
        ENDDO
        k=N
        DO i=Istr,Iend 
            cff     = 1./(BC(i,k)-FC(i,k)*CF(i,k-1))
            CF(i,k) = cff* FC(i,k+1) 
            DC(i,k) = DC(i,k)*2./Hz(i,j,k)
          ENDDO
#else
        DO k=1,N
          DO i=Istr,Iend 
            cff     = 1./(BC(i,k)-FC(i,k)*CF(i,k-1))
            CF(i,k) = cff* FC(i,k+1) 
            DC(i,k) = cff*(DC(i,k)-FC(i,k)*DC(i,k-1))
          ENDDO
        ENDDO
#endif
        DO k=N-1,1,-1
          DO i=Istr,Iend
            DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1) 
          ENDDO
        ENDDO
        DO k=1,N         
          DO i=Istr,Iend
            wz(i,j,k,nnew)=DC(i,k)
          ENDDO
        ENDDO
#  if defined NBQ_FREESLIP
        DO i=Istr,Iend  
          wz(i,j,0,nnew)=2.*wz(i,j,0,nnew)/Hz(i,j,1)
        ENDDO
#  endif

      ENDDO      ! <-- j

# else

      DO j=Jstr,Jend
        DO k=1,N-1
          DO i=Istr,Iend
            wz(i,j,k,nnew)=2.*wz(i,j,k,nnew)/(Hz(i,j,k)+Hz(i,j,k+1))
          ENDDO
        ENDDO
        DO i=Istr,Iend
          wz(i,j,N,nnew)=2.*wz(i,j,N,nnew)/Hz(i,j,N)
#  ifdef NBQ_FREESLIP
         wz(i,j,0,nnew)=2.*wz(i,j,0,nnew)/Hz(i,j,1)
#   else
         wz(i,j,0,nnew)=0.
#  endif
        ENDDO
      ENDDO

# endif /* NBQ_WDIF */
!$acc end kernels
!
!--------------------------------------------------------------------
! Set PHYSICAL lateral boundary conditions.
!--------------------------------------------------------------------
!
      call w3dbc_tile (Istr,Iend,Jstr,Jend, grad)  
!
!--------------------------------------------------------------------
! Exchange periodic boundaries and computational margins.
!--------------------------------------------------------------------
!
# if defined BD_EXCHANGE
#  ifdef THREE_GHOST_POINTS_UV
      call exchange_w3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                        wz(START_2D_ARRAY,0,nnew))
#  else
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                        wz(START_2D_ARRAY,0,nnew))
#  endif
# endif

      
      return
      end
#else
      subroutine step3d_w_empty 
      return
      end
#endif /* SOLVE3D && M3SLOW_W */
