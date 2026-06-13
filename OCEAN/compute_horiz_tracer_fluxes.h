#if defined TS_HADV_UP5 || defined TS_HADV_C6 || defined TS_HADV_WENO5
!----------------------------------------------------------
# ifdef PREDICTOR
!
!----------------------------------------------------------
! Sixth order advection scheme [PREDICTOR]
!  ... except for WENO5
!----------------------------------------------------------
!
#  ifdef TS_HADV_WENO5
#   define FLUX5 flux5_weno
#   define FLUX3 flux3_weno
#   define FLUX2 flux1
#   define UP5_MASKING
#  else
#   define FLUX5 flux6
#   define FLUX3 flux4
#   define FLUX2 flux2
#   undef  UP5_MASKING
#  endif
          cdif=1.
#  include "t3dadv_order5.h"
#  undef FLUX5
#  undef FLUX3
#  undef FLUX2
#  undef UP5_MASKING
# else
!----------------------------------------------------------
! 5th or 6th order or WENO5 advection schemes [CORRECTOR]
!----------------------------------------------------------
!
#  ifdef TS_HADV_C6
#   define FLUX5 flux6
#   define FLUX3 flux4
#   define FLUX2 flux2
#   undef  UP5_MASKING
#  elif defined TS_HADV_WENO5
#   define FLUX5 flux5_weno
#   define FLUX3 flux3_weno
#   define FLUX2 flux1
#   define UP5_MASKING
#  else
#   define FLUX5 flux5
#   define FLUX3 flux3
#   define FLUX2 flux1
#   define UP5_MASKING
#  endif
          cdif=1.
#  include "t3dadv_order5.h"
#  undef FLUX5
#  undef FLUX3
#  undef FLUX2
#  undef UP5_MASKING
# endif
!----------------------------------------------------------
#elif defined TS_HADV_C2
!
!----------------------------------------------------------
! Second order advection scheme
!----------------------------------------------------------
!
          do j=Jstr-2,Jend+1
            do i=Istr-2,Iend+1
              zdxt(i,j)=(t(i+1,j,k,nadv,itrc)-t(i,j,k,nadv,itrc))
#  ifdef MASKING
     &                                               *umask(i+1,j)
#  endif
              zdyt(i,j)=(t(i,j+1,k,nadv,itrc)-t(i,j,k,nadv,itrc))
#  ifdef MASKING
     &                                               *vmask(i,j+1)
#  endif
            enddo
          enddo

# if defined PREDICTOR
          rdtstp = (1.0-gamma)*dt
# else
          rdtstp = dt
# endif

          do j=Jstr-1,Jend+1
            do i=Istr-1,Iend+1
              zzslpx = ( zdxt(i,j) + zdxt(i-1,j) )
     &          * ( 0.25 + SIGN( 0.25, zdxt(i,j) * zdxt (i-1,j) ) )
              zzslpy = ( zdyt(i,j) + zdyt(i,j-1) )
     &          * ( 0.25 + SIGN( 0.25, zdyt(i,j) * zdyt (i,j-1) ) )
              zslpx(i,j) = SIGN( 1.0, zzslpx) * MIN ( ABS( zzslpx),
     &                     2. * ABS(zdxt(i-1,j)), 2.*ABS(zdxt(i,j)) )
              zslpy(i,j) = SIGN( 1.0, zzslpy) * MIN ( ABS( zzslpy),
     &                     2. * ABS(zdyt(i,j-1)), 2.*ABS(zdyt(i,j)) )
             ENDDO
          ENDDO

          DO j=Jstr,Jend
             DO i=Istr-1,Iend
               z0u = SIGN( 0.5, Huon(i+1,j,k) )
               zalpha = 0.5 - z0u
               zu  = z0u-0.5*Huon(i+1,j,k)*rdtstp*pn_u(i+1,j)
     &               *pm_u(i+1,j)/(0.5*(Hz(i,j,k)+Hz(i+1,j,k)))
               zzwx = t(i+1,j,k,nadv,itrc) + zu * zslpx(i+1,j)
               zzwy = t(i,j,k,nadv,itrc) + zu * zslpx(i,j)
               FX(i+1,j) = Huon(i+1,j,k)*(zalpha*zzwx+(1.-zalpha)*zzwy )
             END DO
          END DO

          DO j=Jstr-1,Jend
             DO i=Istr,Iend
               !
               !
               z0v = SIGN( 0.5, Hvom(i,j+1,k) )
               zalpha = 0.5 - z0v
               zv  = z0v-0.5*Hvom(i,j+1,k)*rdtstp*pn_v(i,j+1)
     &               *pm_v(i,j+1)/(0.5*(Hz(i,j,k)+Hz(i,j+1,k)))
               zzwx = t(i,j+1,k,nadv,itrc)+zv * zslpy(i,j+1)
               zzwy = t(i,j,k,nadv,itrc)+zv * zslpy(i,j)
               FE(i,j+1) = Hvom(i,j+1,k)*(zalpha*zzwx+(1.-zalpha)*zzwy)
            enddo
        enddo


#else   /*  --> UP3 (default) or C4 */
!
!----------------------------------------------------------
! Fourth or Third order advection scheme
!----------------------------------------------------------


!----------------------------------------------------------
# ifdef PREDICTOR
!----------------------------------------------------------
! PREDICTOR [Fourth order advection scheme]
!
#  define grad WORK
!----------------------------------------------------------
# else
!----------------------------------------------------------
!  CORRECTOR
!
#  if ! defined TS_HADV_C4
#   define TS_HADV_UP3
#  endif
!
#  ifdef TS_HADV_UP3
#   define curv WORK
#  else
#   define grad WORK
#  endif
!------------------------------------------------------------
# endif
!------------------------------------------------------------


!------------------------------------------------------------
# ifdef EW_PERIODIC
#  define I_EXT_RANGE Istr-1,Iend+2
# else
#  ifdef MPI
          if (WEST_INTER) then
            imin=Istr-1
          else
            imin=max(Istr-1,1)
          endif
          if (EAST_INTER) then
            imax=Iend+2
          else
            imax=min(Iend+2,Lmmpi+1)
          endif
#   define I_EXT_RANGE imin,imax
#  else
#   define I_EXT_RANGE max(Istr-1,1),min(Iend+2,Lm+1)
#  endif
# endif
# ifdef NS_PERIODIC
#  define J_EXT_RANGE Jstr-1,Jend+2
# else
#  ifdef MPI
          if (SOUTH_INTER) then
            jmin=Jstr-1
          else
            jmin=max(Jstr-1,1)
          endif
          if (NORTH_INTER) then
            jmax=Jend+2
          else
            jmax=min(Jend+2,Mmmpi+1)
          endif
#   define J_EXT_RANGE jmin,jmax
#  else
#   define J_EXT_RANGE max(Jstr-1,1),min(Jend+2,Mm+1)
#  endif
# endif
!-------------------------------------------------------------------
!$acc loop collapse(2)
          do j=Jstr,Jend
            do i=I_EXT_RANGE
              FX(i,j)=(t(i,j,k,nadv,itrc)-t(i-1,j,k,nadv,itrc))
# ifdef MASKING
     &                                               *umask(i,j)
# endif
            enddo
          enddo

# undef I_EXT_RANGE
# ifndef EW_PERIODIC
          if (WESTERN_EDGE) then
            do j=Jstr,Jend
              FX(0,j)=FX(1,j)
            enddo
          endif
          if (EASTERN_EDGE) then
#  ifdef MPI
            do j=Jstr,Jend
              FX(Lmmpi+2,j)=FX(Lmmpi+1,j)
            enddo
#  else
             do j=Jstr,Jend
              FX(Lm+2,j)=FX(Lm+1,j)
            enddo
#  endif
          endif
# endif
!---------------------------------------------------------------------
!$acc loop collapse(2)
          do j=Jstr,Jend 
            do i=Istr-1,Iend+1
# if (defined TS_HADV_C4 || defined PREDICTOR)
              grad(i,j)=0.5*(FX(i+1,j)+FX(i,j))
# elif defined TS_HADV_UP3
              curv(i,j)=FX(i+1,j)-FX(i,j)
# endif
            enddo
          enddo             !--> discard FX
!$acc loop collapse(2)
          do j=Jstr,Jend
            do i=Istr,Iend+1
# if (defined TS_HADV_UP3 && !defined PREDICTOR)
              if (Huon(i,j,k) .gt. 0.) then
                cff=curv(i-1,j)
              else
                cff=curv(i,j)
              endif
              FX(i,j)=0.5*( t(i,j,k,nadv,itrc)+t(i-1,j,k,nadv,itrc)
     &                           -0.333333333333*cff )*Huon(i,j,k)
# else
              FX(i,j)=0.5*( t(i,j,k,nadv,itrc)+t(i-1,j,k,nadv,itrc)
     &                     -0.333333333333*(grad(i,j)-grad(i-1,j))
     &                                                )*Huon(i,j,k)
# endif
            enddo
          enddo            !--> discard grad
!---------------------------------------------------------------------
!$acc loop collapse(2)
          do j=J_EXT_RANGE
            do i=Istr,Iend
              FE(i,j)=(t(i,j,k,nadv,itrc)-t(i,j-1,k,nadv,itrc))
# ifdef MASKING
     &                                               *vmask(i,j)
# endif
            enddo
          enddo
# undef J_EXT_RANGE
# ifndef NS_PERIODIC
          if (SOUTHERN_EDGE) then
            do i=Istr,Iend
              FE(i,0)=FE(i,1)
            enddo
          endif
          if (NORTHERN_EDGE) then
#  ifdef MPI
            do i=Istr,Iend
              FE(i,Mmmpi+2)=FE(i,Mmmpi+1)
            enddo
#  else
            do i=Istr,Iend
              FE(i,Mm+2)=FE(i,Mm+1)
            enddo
#  endif
          endif
# endif
!---------------------------------------------------------------------
!$acc loop collapse(2)
          do j=Jstr-1,Jend+1           !<-- C4 [only for pred]
            do i=Istr,Iend
# if (defined TS_HADV_C4 || defined PREDICTOR)
              grad(i,j)=0.5*(FE(i,j+1)+FE(i,j))
# elif defined TS_HADV_UP3
              curv(i,j)=FE(i,j+1)-FE(i,j)
# endif
            enddo
          enddo            !--> discard FE
          do j=Jstr,Jend+1
            do i=Istr,Iend
# if (defined TS_HADV_UP3 && !defined PREDICTOR)
              if (Hvom(i,j,k) .gt. 0.) then
                cff=curv(i,j-1)
              else
                cff=curv(i,j)
              endif
              FE(i,j)=0.5*( t(i,j,k,nadv,itrc)+t(i,j-1,k,nadv,itrc)
     &                          -0.333333333333*cff )*Hvom(i,j,k)
#  undef curv
# else
              FE(i,j)=0.5*( t(i,j,k,nadv,itrc)+t(i,j-1,k,nadv,itrc)
     &                     -0.333333333333*(grad(i,j)-grad(i,j-1))
     &                                               )*Hvom(i,j,k)
#  undef grad
# endif
            enddo
          enddo            !--> discard grad
!---------------------------------------------------------------------
#endif /* TS_HADV_UP5 */


