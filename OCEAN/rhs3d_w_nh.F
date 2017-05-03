#include "cppdefs.h"

#   define W_HADV
#   define W_VADV
#   define WADV_SPLINES
#   define NBQ_WNOSLIP

#if defined SOLVE3D && defined NBQ

      SUBROUTINE rhs3d_w_nh(tile)

      !!======================================================================
      !!                       ***  Subroutine RHS3D  ***
      !! NBQ mode :  compute right-hand-side for the vertical velocity wz 
      !!======================================================================
      !! History :   2016-11  (F. Lemarié)  Original code
      !!----------------------------------------------------------------------

      implicit none
      integer tile, trd, omp_get_thread_num
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      trd=omp_get_thread_num()
      CALL rhs3d_w_tile (Istr,Iend,Jstr,Jend,  A3d(1,3,trd),    !<-- A3d(1,1,trd), A3d(1,2,trd) contain ru,rv 
     &             A2d(1,1,trd), A2d(1,2,trd), A2d(1,3,trd),    !<-- A3d(1,3,trd) contains rw 
     &             A2d(1,4,trd), A2d(1,5,trd), A2d(1,6,trd), 
     &             A2d(1,7,trd), A3d(1,4,trd), A3d(1,5,trd)     !<-- A3d(1,3,trd) and A3d(1,4,trd) should 
     &                                                     )    !<-- be available at this point ...
      return
      END
      
      !!----------------------------------------------------------------------        
      
      SUBROUTINE rhs3d_w_tile (istr, iend,   jstr,   jend, rw ,
     &                           CF,   FC,  We_r ,    lap, Hzw, 
     &                          WFx,  WFe, HUon_w, HVom_w       )      
      IMPLICIT NONE
      !!
      INTEGER  :: Istr,Iend,Jstr,Jend
      INTEGER  :: i,j,k,kp
      INTEGER  :: imin,imax,jmin,jmax
      !!
# include "param.h"
      REAL     :: rw     (PRIVATE_2D_SCRATCH_ARRAY,0:N  )
      REAL     :: CF     (PRIVATE_1D_SCRATCH_ARRAY,0:N  )
      REAL     :: FC     (PRIVATE_1D_SCRATCH_ARRAY,0:N+1)  
      REAL     :: We_r   (PRIVATE_1D_SCRATCH_ARRAY,0:N+1)        
      REAL     :: lap    (PRIVATE_2D_SCRATCH_ARRAY      )
      REAL     :: Hzw    (PRIVATE_1D_SCRATCH_ARRAY,0:N  )      
      REAL     :: WFx    (PRIVATE_2D_SCRATCH_ARRAY      ) 
      REAL     :: WFe    (PRIVATE_2D_SCRATCH_ARRAY      )
      REAL     :: HUon_w (PRIVATE_2D_SCRATCH_ARRAY,0:N  )
      REAL     :: HVom_w (PRIVATE_2D_SCRATCH_ARRAY,0:N  )
      !!    
      REAL     :: gamma, epsil, cff, cff1, cff2, cff3, cff4, Omeg_r
      REAL     :: wz_bot,wz_sfc
      PARAMETER (gamma=0.25  )
      PARAMETER (epsil=1.E-16)
      !!
# include "grid.h"
# include "ocean3d.h"
# include "coupling.h"
# include "forces.h"
# include "scalars.h"  
# include "nbq.h"
# include "compute_auxiliary_bounds.h"
      !!
# ifdef EW_PERIODIC 
#   define IU_RANGE Istr,Iend
#   define IV_RANGE Istr,Iend
# else
#   define IU_RANGE Istr,IendR
#   define IV_RANGE IstrR,IendR
# endif

# ifdef NS_PERIODIC 
#   define JU_RANGE Jstr,Jend
#   define JV_RANGE Jstr,Jend
# else
#   define JU_RANGE JstrR,JendR
#   define JV_RANGE Jstr,JendR
# endif

      !-----------------------------------------------------------------------
      DO j = jstr,jend
         DO i = istr,iend+1
            HUon_w(i,j,0)=0.5* HUon(i ,j ,1 )                
            HUon_w(i,j,N)=0.5* HUon(i ,j ,N ) 
         ENDDO
      ENDDO
      !
      DO j = jstr,jend+1
         DO i = istr,iend
            HVom_w(i,j,0)=0.5* HVom(i ,j ,1 )             
            HVom_w(i,j,N)=0.5* HVom(i ,j ,N ) 
         ENDDO
      ENDDO 
      !
      DO k = 1,N-1
         DO j = jstr,jend
            DO i = istr,iend+1
               HUon_w(i,j,k)=0.5*( HUon(i,j,k  )+HUon(i,j,k+1) )               
            ENDDO
         ENDDO
         
         DO j = jstr,jend+1
            DO i = istr,iend
               HVom_w(i,j,k)=0.5*( HVom(i,j,k  )+HVom(i,j,k+1) )               
            ENDDO
         ENDDO         
      ENDDO
      !-----------------------------------------------------------------------


# if defined UV_COR
      !-----------------------------------------------------------------------
      !  Add in Coriolis  FWcor = +e (U cos(angle) - V sin(angle) ) 
      !  with e = 2 Omega cos Phi
      !-----------------------------------------------------------------------
      DO k = 0, N
         DO j = jstr,jend   
            DO i = istr,iend        
              rw(i,j,k) =   0.5*eomn(i,j)*(
     &            cosa(i,j)*pn(i,j)*(HUon_w(i,j,k)+HUon_w(i+1,j,k)) 
     &          - sina(i,j)*pm(i,j)*(HVom_w(i,j,k)+HVom_w(i,j+1,k))
     &                                    )             
            ENDDO
         ENDDO
      ENDDO      
      !!----------------------------------------------------------------------         
# else
      DO k = 0, N
         DO j = jstr,jend   
            DO i = istr,iend        
              rw(i,j,k) =   0.            
            ENDDO
         ENDDO
      ENDDO 
# endif

# if defined UV_ADV && defined W_HADV
      !-----------------------------------------------------------------------      
      !  Horizontal advection
      !----------------------------------------------------------------------- 
      !!
#ifndef EW_PERIODIC
      IF (WESTERN_EDGE) THEN       
         imin=istr                  
      ELSE                        
         imin=istr-1             
      ENDIF                      
      IF (EASTERN_EDGE) THEN      
         imax=iend                 
      ELSE                  
         imax=iend+1           
      ENDIF
#else
      imin=istr-1
      imax=iend+1
#endif
      !!
#ifndef NS_PERIODIC
      IF (SOUTHERN_EDGE) THEN
         jmin=jstr
      ELSE
         jmin=jstr-1
      ENDIF
      IF (NORTHERN_EDGE) then
         jmax=jend
      ELSE
         jmax=jend+1
      ENDIF
#else
      jmin=jstr-1
      jmax=jend+1
#endif
      !!      
      DO k=0,N             
         !
         DO j=jstr,jend
            DO i=imin,imax+1
              WFx(i,j) = ( wz(i,j,k,nrhs) - wz(i-1,j,k,nrhs) )
#ifdef MASKING
     &                                             *umask(i,j)
#endif
            ENDDO
         ENDDO  
         !    
#ifndef EW_PERIODIC
         IF (WESTERN_EDGE) THEN
            DO j=jstr,jend
              WFx(istr-1,j) = WFx(istr  ,j)
            ENDDO
         ENDIF
         IF (EASTERN_EDGE) then
            DO j=jstr,jend
              WFx(iend+2,j) = WFx(iend+1,j)
            ENDDO
         ENDIF
#endif
         !
         DO j=Jstr,Jend
            DO i=Istr-1,Iend+1          
               lap(i,j) = WFx(i+1,j) - WFx(i,j)
            ENDDO
         ENDDO  
         !
         DO j=Jstr,Jend
            DO i=Istr,Iend+1
               IF ( Huon_w(i,j,k).gt. 0. ) THEN
                  cff = lap( i-1, j )
               ELSE
                  cff = lap( i  , j )
               ENDIF            
               WFx(i,j)=0.5*( wz(i,j,k,nrhs)+wz(i-1,j,k,nrhs)
     &                              -0.25*cff )*Huon_w(i,j,k)
            ENDDO    
         ENDDO
         
         
         DO j=jmin,jmax+1
            DO i=istr,iend
               WFe(i,j) = ( Wz(i,j,k,nrhs)-Wz(i,j-1,k,nrhs) ) 
#  ifdef MASKING
     &                                         *vmask(i,j)
#  endif
            ENDDO 
         ENDDO         
#ifndef NS_PERIODIC
         IF (SOUTHERN_EDGE) then
            DO i=istr,iend
              WFe(i,jstr-1) = WFe(i,jstr)
            ENDDO
         ENDIF
         IF (NORTHERN_EDGE) then
            DO i=istr,iend
              WFe(i,jend+2) = WFe(i,jend+1)
            ENDDO
         ENDIF
#endif         
         !
         DO j=jstr-1,jend+1
            DO i=istr,iend
              lap(i,j) = WFe(i,j+1) - WFe(i,j  )
            ENDDO
         ENDDO  
         !
         DO j=jstr,jend+1
            DO i=istr,iend
               IF ( HVom_w(i,j,k).gt. 0. ) THEN
                  cff = lap( i , j-1 )
               ELSE
                  cff = lap( i , j   )
               ENDIF 
               WFe(i,j)=0.5*( Wz(i,j,k,nrhs)+Wz(i,j-1,k,nrhs)
     &                              -0.25*cff )*Hvom_w(i,j,k)
            ENDDO
         ENDDO 
         !           
         DO j=Jstr,Jend
            DO i=Istr,Iend
              rw(i,j,k)= rw(i,j,k) - WFx(i+1,j) + WFx(i,j)
     &                             - WFe(i,j+1) + WFe(i,j)    
            ENDDO
         ENDDO
         ! 
       ENDDO          !<-- outer loop k
# endif /* W_HADV */




# if defined UV_ADV && defined W_VADV

# define udH_x WFx
# define vdH_y WFe  
# define DUon(i,j) Huon_w(i,j,0)  
# define DVom(i,j) Hvom_w(i,j,0)
# define udZ_x(i,j) Huon_w(i,j,1)  
# define vdZ_y(i,j) Hvom_w(i,j,1)
       !-----------------------------------------------------------------------      
       !  Bottom boundary condition for wz
       !----------------------------------------------------------------------- 
       IF(nrhs.eq.3) then    !<-- corrector
          DO j=Jstr,Jend 
             DO i=Istr,Iend+1 
# ifndef NBQ_WNOSLIP
                udH_x(i,j) = - u(i,j,1,nrhs) * pm_u(i,j) 
     &                                     * ( z_w(i  ,j,0)-z_w(i-1,j,0) )  
# else
                udH_x(i,j) = 0. 
# endif
                DUon(i,j) = 2.*DU_avg2(i,j) 
     &             / ( rhobar_nbq_avg1(i,j) + rhobar_nbq_avg1(i-1,j) ) 
                udZ_x(i,j) = u(i,j,N,nrhs) * pm_u(i,j) 
     &                                     * ( z_w(i  ,j,N)-z_w(i-1,j,N) ) 
             ENDDO
          ENDDO

          DO j=Jstr,Jend+1 
             DO i=Istr,Iend 
# ifndef NBQ_WNOSLIP
                vdH_y(i,j) = - v(i,j,1,nrhs) * pn_v(i,j) 
     &                                     * ( z_w(i,j,0)-z_w(i,j-1,0) )
# else
                vdH_y(i,j) = 0. 
# endif
                DVom(i,j)  = 2.*DV_avg2(i,j)
     &           / ( rhobar_nbq_avg1(i,j) + rhobar_nbq_avg1(i,j-1) )
                vdZ_y(i,j) = v(i,j,N,nrhs) * pn_v(i,j) 
     &                                     * ( z_w(i,j,N)-z_w(i,j-1,N) ) 
             ENDDO
         ENDDO
       ELSE                  !<-- predictor
          DO j=Jstr,Jend 
             DO i=Istr,Iend+1 
# ifndef NBQ_WNOSLIP
                udH_x(i,j) = - u(i,j,1,nstp) * pm_u(i,j) 
     &                                     * ( z_w(i  ,j,0)-z_w(i-1,j,0) )   
# else
                udH_x(i,j) = 0.
# endif
                DUon(i,j) = 2.*DU_avg1(i,j,nstp) 
     &             / ( rhobar_nbq_avg1(i,j) + rhobar_nbq_avg1(i-1,j) )  
                udZ_x(i,j) = u(i,j,N,nstp) * pm_u(i,j) 
     &                                     * ( z_w(i  ,j,N)-z_w(i-1,j,N) ) 
             ENDDO
          ENDDO

          DO j=Jstr,Jend+1 
             DO i=Istr,Iend 
# ifndef NBQ_WNOSLIP
                vdH_y(i,j) = - v(i,j,1,nstp) * pn_v(i,j) 
     &                                     * ( z_w(i,j,0)-z_w(i,j-1,0) )
# else
                vdH_y(i,j) = 0.
# endif
                DVom(i,j) = 2.*DV_avg1(i,j,nstp)
     &           / ( rhobar_nbq_avg1(i,j) + rhobar_nbq_avg1(i,j-1) )
                vdZ_y(i,j) = v(i,j,N,nstp) * pn_v(i,j) 
     &                                     * ( z_w(i,j,N)-z_w(i,j-1,N) )      
             ENDDO
          ENDDO
       ENDIF    

      !-----------------------------------------------------------------------      
      !  Vertical advection
      !----------------------------------------------------------------------- 
       DO j=Jstr,Jend 
 
          DO i=Istr,Iend          
             Hzw (i,0  ) = 0.5 * HZR(i,j,1) 
             Hzw (i,N  ) = 0.5 * HZR(i,j,N)  
             We_r(i,0  ) =        We(i,j,0)
             We_r(i,N  ) = 0.5 * (We(i,j,N)+We(i,j,N-1))             
             We_r(i,N+1) =        We(i,j,N)        
          ENDDO

 
          DO k=1,N-1
             DO i=Istr,Iend
                Hzw (i,k) = 0.5d0 * ( HZR(i,j,k) + HZR(i,j,k+1) )
                We_r(i,k) = 0.5d0 * ( We (i,j,k) + We (i,j,k-1) ) 
             ENDDO        
          ENDDO
         
#undef WADV_SPLINES 
#  ifdef WADV_SPLINES          
          DO i=Istr,Iend      
             wz_bot  = - 0.5d0 * ( udH_x(i,j) + udH_x(i+1,j)
     &                         +   vdH_y(i,j) + vdH_y(i,j+1) )  
#     ifdef MOVING_BATHY
             wz_bot  =  wz_bot - dh(i,j)     
#     endif  

             cff     =   0.5 / ( Hzw (i,1) + Hzw (i,0) )         !<--- 1 / b(1)
             CF(i,1) =   cff * Hzw(i,0)                          !<--- q(1) = c(1) / b(1)
             FC(i,1) =   cff * (  3.*Hzw(i,0)*wz(i,j,1,nrhs)     
     &                          + 3.*Hzw(i,1)*wz(i,j,0,nrhs)
     &                          -    Hzw(i,1)*wz_bot             !<--- f(1) / b(1) 
     &                         )     
          ENDDO          
          

          DO k=2,N-1
             DO i=Istr,Iend          
                cff     = 1./( 2.*Hzw(i,k-1)+Hzw(i,k)*(2.-CF(i,k-1)) )   !<--- p = 1 / ( b(k)+a(k)*q(k-1) )
                CF(i,k) = cff* Hzw(i,k-1)                                !<--- c(k) * p
                FC(i,k) = cff*( 3.*( Hzw(i,k-1)*wz(i,j,k  ,nrhs)         !<--- f(k)=( f(k)-a(k)*f(k-1) )*p
     &                             + Hzw(i,k  )*wz(i,j,k-1,nrhs) )
     &                             - Hzw(i,k  )*FC(i,k-1)  )
             ENDDO
          ENDDO                   
          
          
          DO i=Istr,Iend           
             wz_sfc  = pn(i,j)*pm(i,j)*( DUon(i,j)-DUon(i+1,j  )
     &                               +   DVom(i,j)-DVom(i  ,j+1) )    !
     &               + 0.5d0 * ( udZ_x(i,j) + udZ_x(i+1,j)
     &                       +   vdZ_y(i,j) + vdZ_y(i,j+1) )

             cff     = 1./( 2.*Hzw(i,N-1)+Hzw(i,N)*(2.-CF(i,N-1)) )   !<--- p = 1 / ( b(k)+a(k)*q(k-1) )
             FC(i,N) = cff*( 3.*( Hzw(i,N-1)*wz(i,j,N  ,nrhs)         !<--- f(k)=( f(k)-a(k)*f(k-1) )*p
     &                          + Hzw(i,N  )*wz(i,j,N-1,nrhs) )
     &                          - Hzw(i,N-1)*wz_sfc
     &                          - Hzw(i,N  )*FC(i,N-1)  )    
          ENDDO          

          DO k=N-1,1,-1
             DO i=Istr,Iend 
                FC(i,k  )=FC(i,k)-CF(i,k)*FC(i,k+1)
                FC(i,k+1)=FC(i,k+1)*We_r(i,k+1)  
             ENDDO
          ENDDO

          DO i=Istr,Iend 
             FC(i,1  )=FC(i,1  )*We_r(i,1)
             FC(i,0  )=0. 
             FC(i,N+1)=0.
          ENDDO
#  else      
          DO k=2,N-1
             DO i=Istr,Iend
                FC(i,k) = 0.5*(wz(i,j,k,nrhs) + wz(i,j,k-1,nrhs))*We_r(i,k)
             ENDDO
          ENDDO

          DO i=Istr,Iend        
  
!             wz_bot  = - 0.5d0 * ( udH_x(i,j) + udH_x(i+1,j)
!     &                         +   vdH_y(i,j) + vdH_y(i,j+1) )  
!#     ifdef MOVING_BATHY
!             wz_bot  =  wz_bot - dh(i,j)     
!#     endif 
           ! FC(i,1) = 0.5*( wz(i,j,1,nrhs) + wz_bot  )*We_r(i,1)
             FC(i,1) = 0.5*( wz(i,j,1,nrhs) + wz(i,j,0,nrhs)  )*We_r(i,1)

!             wz_sfc  = pn(i,j)*pm(i,j)*( DUon(i,j)-DUon(i+1,j  )
!     &                               +   DVom(i,j)-DVom(i  ,j+1) )    !
!     &               + 0.5d0 * ( udZ_x(i,j) + udZ_x(i+1,j)
!     &                       +   vdZ_y(i,j) + vdZ_y(i,j+1) )
!            FC(i,N) = 0.5*( wz_sfc + wz(i,j,N-1,nrhs))*We_r(i,N)
             FC(i,N) = 0.5*( wz(i,j,N,nrhs) + wz(i,j,N-1,nrhs))*We_r(i,N)
          ENDDO          

          DO i=Istr,Iend
             FC(i,0  )=0.
             FC(i,N+1)=0.
          ENDDO          
#  endif
# undef udH_x
# undef vdH_y   
# undef DUon 
# undef DVom 
# undef udZ_x  
# undef vdZ_y
      
          DO k=1,N
             DO i=Istr,Iend
                rw(i,j,k) = rw(i,j,k)  - FC(i,k+1) + FC(i,k)
#  ifdef MASKING              
                rw(i,j,k) = rw(i,j,k) * rmask(i,j)
#  endif             
             ENDDO
          ENDDO
!          DO i=Istr,Iend
!             rw(i,j,0) = rw(i,j,0)  + FC(i,0)
!#  ifdef MASKING             
!             rw(i,j,0) = rw(i,j,0) * rmask(i,j)
!#  endif             
!          ENDDO

       ENDDO
                     
# endif /* UV_ADV && W_VADV */
      return
      end  
#else
      subroutine rhs3d_w_empty
      end
#endif /* SOLVE3D && NBQ */       
