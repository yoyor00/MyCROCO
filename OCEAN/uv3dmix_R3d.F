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
#include "cppdefs.h"
#ifndef CHILD_SPG
      subroutine uv3dmix (tile)
      implicit none
      integer tile, trd, omp_get_thread_num
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      trd=omp_get_thread_num()
!
# ifdef MASKING
#  define SWITCH *
# else
#  define SWITCH !
# endif
!
# ifdef AGRIF
       if (AGRIF_Root()) then
        call uv3dmix_tile (Istr,Iend,Jstr,Jend,
     &        A2d(1,1,trd), A2d(1,2,trd),A2d(1,3,trd), A2d(1,4,trd)
#  ifdef NBQ
     &                       ,A3d(1,3,trd), A3d(1,4,trd)
#  endif
     &                                                 )
       else
        call uv3dmix_child_tile (Istr,Iend,Jstr,Jend,
     &        A2d(1,1,trd), A2d(1,2,trd),A2d(1,3,trd), A2d(1,4,trd)
#  ifdef NBQ
     &                       ,A3d(1,3,trd), A3d(1,4,trd)
#  endif
     &                                                 )
       endif
      return
      end
# else
      call uv3dmix_tile (Istr,Iend,Jstr,Jend,
     &      A2d(1,1,trd), A2d(1,2,trd),A2d(1,3,trd), A2d(1,4,trd)
#  ifdef NBQ
     &                       ,A3d(1,3,trd), A3d(1,4,trd)
#  endif
     &                                                 )
      return
      end
# endif   /* AGRIF */
!
!---------------------------------------------------------------------
!*********************************************************************
!---------------------------------------------------------------------
!
!PARENT
!
      subroutine uv3dmix_tile (Istr,Iend,Jstr,Jend,UFx,UFe,VFx,VFe
# ifdef NBQ
     &                                      ,UFz,VFz
# endif
     &                                                 )
# undef CLIMAT_UV_MIXH_FINE
!
#else
!
! CHILD
!
      subroutine uv3dmix_child_tile (Istr,Iend,Jstr,Jend,UFx,UFe,VFx,VFe
# ifdef NBQ
     &                                        ,UFz,VFz
# endif
     &                                                 )
! Diffusion always applied on U-UCLM in fine grid
# if !defined UV_HADV_RSUP3
#  undef CLIMAT_UV_MIXH_FINE
# endif
!
#endif /* CHILD_SPG */
!
!--------------------------------------------------------------------
! (1) Compute 3D viscous terms for momentum, along constant S-surfaces
! as 3D divergence of the full Reynolds stress tensor.
! Components of the stress tensor are:
!                       du                  dv                  dw
!         s_xx = 2 Akv ---- , s_yy = 2 Akv ---- , s_zz = 2 Akv ----
!                       dx                  dy                  dz
!
!                               du     dv
!         s_xy =  s_yx = Akv ( ---- + ---- )
!                               dy     dx
!
!                               du     dw
!         s_xz =  s_zx = Akv ( ---- + ---- )
!                               dz     dx
!
!                               dv     dw
!         s_yz =  s_zy = Akv ( ---- + ---- )
!                               dz     dy
!
! The pieces du/dz and dv/dz in s_xz and s_yz are handled implicitly
!     in the u and v equation and explicitly in the w equation
!
! (2) Compute 3D TKE production term by shear
!
! Future developments :
!   - Energy consistent discretization of the Sprod3d term following
!                        Burchard (2002) and Marsaleix et al. (2008)
!
!---------------------------------------------------------------------
!---------------------------------------------------------------------
! ******************************Common Code***************************
!---------------------------------------------------------------------
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k, indx
#include "param.h"
      real UFe(PRIVATE_2D_SCRATCH_ARRAY),
     &     UFx(PRIVATE_2D_SCRATCH_ARRAY),
     &     VFe(PRIVATE_2D_SCRATCH_ARRAY),
     &     VFx(PRIVATE_2D_SCRATCH_ARRAY)
#ifdef NBQ
      real UFz(PRIVATE_2D_SCRATCH_ARRAY,0:N),
     &     VFz(PRIVATE_2D_SCRATCH_ARRAY,0:N)
#endif
#include "scalars.h"
#include "grid.h"
#include "ocean3d.h"
#include "coupling.h"
#include "mixing.h"
#ifdef CLIMAT_UV_MIXH
# include "climat.h"
#endif
#ifdef AGRIF
# include "zoom.h"
#endif
      real     Akv_loc,Hz_loc
      real     cff,cff1,cff2,cff3,cff4
      real     divE,divX,divZ
      integer  kp,km
!
#include "compute_auxiliary_bounds.h"
!
#ifdef CHILD_SPG
# define UCLM usponge
# define VCLM vsponge
#else
# define UCLM uclm
# define VCLM vclm
#endif
!
      indx=3-nstp    !--> time index for target arrays;
!$acc kernels if(compute_on_device) default(present)
!
!-------------------------------------------------------------------------------
!  Compute flux-components of the Reynolds stress tensor (m4/s2)
!-------------------------------------------------------------------------------
!
#ifdef NBQ
      DO k = 0,N
        kp=MIN(k+1,N); km=MAX(k  ,1)
        DO j=JstrV-1,Jend+1
          DO i=IstrU-1,Iend+1
            Akv_loc = 0.25*( Akv(i  ,j,km) + Akv(i  ,j,kp)
     &                     + Akv(i-1,j,km) + Akv(i-1,j,kp) )  ! Akv at i-1/2,j,k+1/2 locations
# ifdef SPONGE_VIS2
            Akv_loc = MAX( Akv_loc, 0.5*( visc2_sponge_r(i  ,j)
     &                                 +  visc2_sponge_r(i-1,j) ) )
# endif
            Hz_loc  = 0.25*( Hz(i  ,j,km) + Hz(i  ,j,kp)
     &                     + Hz(i-1,j,km) + Hz(i-1,j,kp) )   ! Hz at i-1/2,j,k+1/2 locations
            !
            !-- UFz = s_xz = Akv (Hz/pn) d_x( wz/Hz )    (the vertical piece is missing as it is treated implicitly in step3D_uv)
            !
            UFz(i,j,k) = Akv_loc*Hz_loc*on_u(i,j)*(
     &            2.*wz(i  ,j,k,nstp)/(Hz(i  ,j,kp)+Hz(i  ,j,km))
     &         -  2.*wz(i-1,j,k,nstp)/(Hz(i-1,j,kp)+Hz(i-1,j,km)) )  ! m4/s2

            Akv_loc = 0.25*( Akv(i,j  ,km) + Akv(i,j  ,kp)
     &                     + Akv(i,j-1,km) + Akv(i,j-1,kp) )  ! Akv at i,j-1/2,k+1/2 locations
# ifdef SPONGE_VIS2
            Akv_loc = MAX( Akv_loc,0.5*( visc2_sponge_r(i,j  )
     &                                +  visc2_sponge_r(i,j-1) ) )
# endif
            Hz_loc  = 0.25*( Hz(i,j  ,km) + Hz(i,j  ,kp)
     &                     + Hz(i,j-1,km) + Hz(i,j-1,kp) )  ! Hz at i,j-1/2,k+1/2 locations
            !
            !-- VFz = s_yz = Akv (Hz/pm) d_y( wz/Hz )   (the vertical piece is missing as it is treated implicitly in step3D_uv)
            !
            VFz(i,j,k) = Akv_loc*Hz_loc*om_v(i,j)*(
     &           2.*wz(i,j  ,k,nstp)/(Hz(i,j  ,kp)+Hz(i,j  ,km))
     &        -  2.*wz(i,j-1,k,nstp)/(Hz(i,j-1,kp)+Hz(i,j-1,km)) )  ! m4/s2
          ENDDO
        ENDDO
      ENDDO
#endif
      !--
      DO k=1,N
        !=======================================================================
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
#ifdef SPONGE_VIS2
            Akv_loc  = MAX(Akv(i,j,k),visc2_sponge_r(i,j))
#else
            Akv_loc  = Akv(i,j,k)
#endif
            cff      = 2.*Akv_loc*Hz(i,j,k)
            !
            !-- UFx = s_xx = 2 Akv (Hz pm / pn) d_x( u )
            !
            UFx(i,j) = cff*pmon_r(i,j)*(u(i+1,j,k,nstp)-u(i,j,k,nstp)) ! m4/s2
            !
            !-- VFe = s_yy = 2 Akv (Hz pn / pm) d_y( v )
            !
            VFe(i,j) = cff*pnom_r(i,j)*(v(i,j+1,k,nstp)-v(i,j,k,nstp)) ! m4/s2
          ENDDO
        ENDDO
        !=======================================================================
        DO j=Jstr,Jend+1
          DO i=Istr,Iend+1
            Hz_loc  = 0.25*( Hz(i  ,j,k) + Hz(i  ,j-1,k)
     &                      +Hz(i-1,j,k) + Hz(i-1,j-1,k) )      ! Hz at i-1/2,j-1/2,k locations
            Akv_loc = 0.25*( Akv(i  ,j-1,k) + Akv(i  ,j  ,k)
     &                     + Akv(i-1,j  ,k) + Akv(i-1,j-1,k) )  ! Akv at i-1/2,j-1/2,k locations
#ifdef SPONGE_VIS2
            Akv_loc = MAX(Akv_loc,visc2_sponge_p(i,j))
#endif
            !
            !-- UFe = s_xy = Akv Hz ( d_x(pn v)/pn + d_y(pm u) (pn/pm^2) )
            !
            UFe(i,j) = Akv_loc*Hz_loc*(
     &               on_p(i,j)*(   pn_v(i  ,j)*v(i  ,j,k,nstp)
     &                           - pn_v(i-1,j)*v(i-1,j,k,nstp) )
     &  +pnom_p(i,j)*om_p(i,j)*(   pm_u(i,j  )*u(i,j  ,k,nstp)
     &                           - pm_u(i,j-1)*u(i,j-1,k,nstp) )
     &                                                           )  ! m4/s2
            !
            !-- VFx = s_yx = Akv Hz ( d_x(pn v)*(pm/pn^2) + d_y(pm u)/pm )
            !
            VFx(i,j) = Akv_loc*Hz_loc*(
     &   pmon_p(i,j)*on_p(i,j)*(   pn_v(i  ,j)*v(i  ,j,k,nstp)
     &                           - pn_v(i-1,j)*v(i-1,j,k,nstp) )
     &  +            om_p(i,j)*(   pm_u(i,j  )*u(i,j  ,k,nstp)
     &                           - pm_u(i,j-1)*u(i,j-1,k,nstp) )
     &                                                           )  ! m4/s2
          ENDDO
        ENDDO
        !=======================================================================
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            divX = UFx(i,j  )-UFx(i-1,j); divE = UFe(i,j+1)-UFe(i  ,j)
#ifdef NBQ
            divZ = UFz(i,j,k)-UFz(i,j,k-1)
#else
            divZ = 0.
#endif
            cff  = divX + divE + divZ
            cff1 = pm_u(i,j)*pn_u(i,j)*cff
            ! add the contribution of the viscous stress terms to the barotropic forcing
            rufrc(i,j)=rufrc(i,j) + cff
            ! WARNING : at this stage u(i,j,k,indx) contains (Hz u)(i,j,k,indx)
            u(i,j,k,indx)=u(i,j,k,indx) + dt*cff1  ! Apply flux divergence
          ENDDO
        ENDDO
        !=======================================================================
        DO j=JstrV,Jend
          DO i=Istr,Iend
            divX = VFx(i+1,j)-VFx(i,j  ); divE = VFe(i  ,j)-VFe(i,j-1)
#ifdef NBQ
            divZ = VFz(i,j,k)-VFz(i,j,k-1)
#else
            divZ = 0.
#endif
            cff  = divX + divE + divZ
            cff1 = pm_v(i,j)*pn_v(i,j)*cff
            ! add the contribution of the viscous stress terms to the barotropic forcing
            rvfrc(i,j)=rvfrc(i,j) + cff
            ! WARNING : at this stage v(i,j,k,indx) contains (Hz v)(i,j,k,indx)
            v(i,j,k,indx)=v(i,j,k,indx) + dt*cff1
          ENDDO
        ENDDO
        !=======================================================================
      ENDDO
      !--
#ifdef NBQ
#  define WFx UFx
#  define WFe VFe
      DO k = 0,N
        kp=MIN(k+1,N); km=MAX(k  ,1)
        !=======================================================================
        DO j=jstr,jend+1
          DO i=istr,iend+1
            Akv_loc = 0.25*( Akv(i  ,j,km) + Akv(i  ,j,kp)
     &                     + Akv(i-1,j,km) + Akv(i-1,j,kp) )  ! Akv at i-1/2,j,k+1/2 locations
# ifdef SPONGE_VIS2
            Akv_loc = MAX( Akv_loc, 0.5*( visc2_sponge_r(i  ,j)
     &                                 +  visc2_sponge_r(i-1,j) ) )
# endif
            Hz_loc  = 0.25*( Hz(i  ,j,km) + Hz(i  ,j,kp)
     &                     + Hz(i-1,j,km) + Hz(i-1,j,kp) )  ! Hz at i-1/2,j,k+1/2 locations
            !
            !-- WFx = s_zx = Akv( d_z(u)/pn + (pm Hz^2/pn) d_x(w/Hz)  )
            !
            WFx(i,j) = pm_u(i,j)*Hz_loc*UFz(i,j,k)+
     &                 Akv_loc*on_u(i,j)*(u(i,j,kp,nstp)-u(i,j,km,nstp))
            !--
            Akv_loc = 0.25*( Akv(i,j  ,km) + Akv(i,j  ,kp)
     &                     + Akv(i,j-1,km) + Akv(i,j-1,kp) )   ! Akv at i,j-1/2,k+1/2 locations
# ifdef SPONGE_VIS2
            Akv_loc = MAX( Akv_loc,0.5*( visc2_sponge_r(i,j  )
     &                                +  visc2_sponge_r(i,j-1) ) )
# endif
            Hz_loc  = 0.25*( Hz(i,j  ,km) + Hz(i,j  ,kp)
     &                     + Hz(i,j-1,km) + Hz(i,j-1,kp) )  ! Hz at i,j-1/2,k+1/2 locations
            !
            !-- WFe = s_zy = Akv( d_z(v)/pm + (pn Hz^2/pm) d_y(w/Hz)  )
            !
            WFe(i,j) = pn_v(i,j)*Hz_loc*VFz(i,j,k)+
     &                 Akv_loc*om_v(i,j)*(v(i,j,kp,nstp)-v(i,j,km,nstp))
          ENDDO
        ENDDO
        !=======================================================================
        DO j=jstr,jend
          DO i=istr,iend
            cff = pn(i,j)*pm(i,j)
            wz(i,j,k,indx)=wz(i,j,k,indx)+dt*cff
     &            *(WFx(i+1,j)-WFx(i,j)+WFe(i,j+1)-WFe(i,j))
          ENDDO
        ENDDO
        !=======================================================================
      ENDDO
#  undef WFx
#  undef WFe
#endif

#ifdef TKE3D_MIXING
!
!-------------------------------------------------------------------------------
!  Compute 3D TKE shear production term (naive way for the moment)
!-------------------------------------------------------------------------------
!
      DO k=1,N
        kp = k+1; IF( k==N ) kp = k
        km = k-1; IF( k==1 ) km = k
        DO j=jstr,jend
          DO i=istr,iend
            ! Contribution from diagonal terms
            divX    = pm(i,j)*(u(i+1,j,k,nstp)-u(i,j,k,nstp))
            divE    = pn(i,j)*(v(i,j+1,k,nstp)-v(i,j,k,nstp))
            divZ    = (wz(i,j,k,nstp)-wz(i,j,k-1,nstp)) / Hz(i,j,k)
            !
            cff     = 2.*( divX*divX+divE*divE+divZ*divZ ) ! diagonal terms
            ! (dudy+dvdx) term (SUM the contributions from psi points to rho points)
            cff1    = 0.5*( pn_u(i,j+1)+pn_u(i  ,j  ) )*
     &                            ( u(i,j+1,k,nstp) - u(i  ,j  ,k,nstp) )
     &             +  0.5*( pm_v(i,j+1)+pm_v(i-1,j+1) )*
     &                            ( v(i,j+1,k,nstp) - v(i-1,j+1,k,nstp) ) !-> (i-1/2,j+1/2,k)
            cff2    = 0.5*( pn_u(i+1,j+1)+pn_u(i+1,j  ) )*
     &                            ( u(i+1,j+1,k,nstp) - u(i+1,j  ,k,nstp) )
     &             +  0.5*( pm_v(i+1,j+1)+pm_v(i,j+1) )*
     &                            ( v(i+1,j+1,k,nstp) - v(i,j+1,k,nstp) ) !-> (i+1/2,j+1/2,k)
            cff3    = 0.5*( pn_u(i,j)+pn_u(i  ,j-1) )*
     &                            ( u(i,j,k,nstp) - u(i  ,j-1,k,nstp) )
     &             +  0.5*( pm_v(i,j)+pm_v(i-1,j) )*
     &                            ( v(i,j,k,nstp) - v(i-1,j,k,nstp) )     !-> (i-1/2,j-1/2,k)
            cff4    = 0.5*( pn_u(i+1,j)+pn_u(i+1,j-1) )*
     &                            ( u(i+1,j,k,nstp) - u(i+1,j-1,k,nstp) )
     &             +  0.5*( pm_v(i+1,j)+pm_v(i,j) )*
     &                            ( v(i+1,j,k,nstp) - v(i,j,k,nstp) )     !-> (i+1/2,j-1/2,k)
      !
            cff = cff + 0.25*Akv_loc*( cff1*cff1+cff2*cff2+cff3*cff3+cff4*cff4 )
            ! (dudz+dwdx) term
            Hz_loc = 0.25*(Hz(i,j,k)+Hz(i+1,j,k)+Hz(i,j,kp)+Hz(i+1,j,kp))
            cff1 = pm_u(i+1,j)*(wz(i+1,j,k  ,nstp) - wz(i  ,j,k  ,nstp))
     &                     +   ( u(i+1,j,kp,nstp) -  u(i+1,j,k  ,nstp))/Hz_loc  !-> (i+1/2,j,k+1/2)
            Hz_loc = 0.25*(Hz(i-1,j,k)+Hz(i,j,k)+Hz(i-1,j,kp)+Hz(i,j,kp))
            cff2 = pm_u(i  ,j)*(wz(i  ,j,k  ,nstp) - wz(i-1,j,k  ,nstp))
     &                     +   ( u(i  ,j,kp,nstp) -  u(i  ,j,k  ,nstp))/Hz_loc  !-> (i-1/2,j,k+1/2)
            Hz_loc = 0.25*(Hz(i,j,km)+Hz(i+1,j,km)+Hz(i,j,k)+Hz(i+1,j,k))
            cff3 = pm_u(i+1,j)*(wz(i+1,j,km,nstp) - wz(i  ,j,km,nstp))
     &                     +   ( u(i+1,j,k  ,nstp) -  u(i+1,j,km,nstp))/Hz_loc  !-> (i+1/2,j,k-1/2)
            Hz_loc = 0.25*(Hz(i-1,j,km)+Hz(i,j,km)+Hz(i-1,j,k)+Hz(i,j,k))
            cff4 = pm_u(i  ,j)*(wz(i  ,j,km,nstp) - wz(i-1,j,km,nstp))
     &                     +   ( u(i  ,j,k  ,nstp) -  u(i  ,j,km,nstp))/Hz_loc  !-> (i-1/2,j,k-1/2)
            !
            cff = cff + 0.25*( cff1*cff1+cff2*cff2+cff3*cff3+cff4*cff4 )
            ! (dvdz+dwdy) term
            Hz_loc = 0.25*(Hz(i,j,k)+Hz(i,j+1,k)+Hz(i,j,kp)+Hz(i,j+1,kp))
            cff1 = pn_v(i,j+1)*(wz(i,j+1,k  ,nstp) - wz(i  ,j,k  ,nstp))
     &                     +   ( v(i,j+1,kp,nstp) -  v(i,j+1,k  ,nstp))/Hz_loc  !-> (i,j+1/2,k+1/2)
            Hz_loc = 0.25*(Hz(i,j-1,k)+Hz(i,j,k)+Hz(i,j-1,kp)+Hz(i,j,kp))
            cff2 = pn_v(i  ,j)*(wz(i  ,j,k  ,nstp) - wz(i,j-1,k  ,nstp))
     &                     +   ( v(i  ,j,kp,nstp) -  v(i  ,j,k  ,nstp))/Hz_loc  !-> (i,j-1/2,k+1/2)
            Hz_loc = 0.25*(Hz(i,j,km)+Hz(i,j+1,km)+Hz(i,j,k)+Hz(i,j+1,k))
            cff3 = pn_v(i,j+1)*(wz(i,j+1,km,nstp) - wz(i  ,j,km,nstp))
     &                     +   ( v(i,j+1,k  ,nstp) -  v(i,j+1,km,nstp))/Hz_loc  !-> (i,j+1/2,k-1/2)
            Hz_loc = 0.25*(Hz(i,j-1,km)+Hz(i,j,km)+Hz(i,km,k)+Hz(i,j,k))
            cff4 = pn_v(i  ,j)*(wz(i  ,j,km,nstp) - wz(i,j-1,km,nstp))
     &                     +   ( v(i  ,j,k  ,nstp) -  v(i  ,j,km,nstp))/Hz_loc  !-> (i,j-1/2,k-1/2)
            !
            cff = cff + 0.25*( cff1*cff1+cff2*cff2+cff3*cff3+cff4*cff4 )
            !
            Sprod3d(i,j,k) = Akv(i,j,k)*cff
            !
          ENDDO
        ENDDO
      ENDDO
      !
      !------------------------------------------------
      !  Apply lateral boundary conditions
      !------------------------------------------------
      !
# ifndef EW_PERIODIC
      IF (WESTERN_EDGE) THEN
        DO j=jstr,jend; DO k=1,N
          Sprod3d(istr-1,j,k)=Sprod3d(istr,j,k)
        ENDDO; ENDDO
      ENDIF
      IF (EASTERN_EDGE) THEN
        DO j=jstr,jend; DO k=1,N
          Sprod3d(iend+1,j,k)=Sprod3d(iend,j,k)
        ENDDO; ENDDO
      ENDIF
# endif
# ifndef NS_PERIODIC
      IF (SOUTHERN_EDGE) THEN
        DO i=istr,iend; DO k=1,N
          Sprod3d(i,jstr-1,k)=Sprod3d(i,jstr,k)
        ENDDO; ENDDO
      ENDIF
      IF (NORTHERN_EDGE) THEN
        DO i=istr,iend; DO k=1,N
          Sprod3d(i,jend+1,k)=Sprod3d(i,jend,k)
        ENDDO; ENDDO
      ENDIF
#  ifndef EW_PERIODIC
      IF (WESTERN_EDGE .and. SOUTHERN_EDGE) THEN
        DO k=1,N
          Sprod3d(istr-1,jstr-1,k)=Sprod3d(istr,jstr,k)
        ENDDO
      ENDIF
      IF (WESTERN_EDGE .and. NORTHERN_EDGE) THEN
        DO k=1,N
          Sprod3d(istr-1,jend+1,k)=Sprod3d(istr,jend,k)
        ENDDO
      ENDIF
      IF (EASTERN_EDGE .and. SOUTHERN_EDGE) THEN
        DO k=1,N
          Sprod3d(iend+1,jstr-1,k)=Sprod3d(iend,jstr,k)
        ENDDO
      ENDIF
      IF (EASTERN_EDGE .and. NORTHERN_EDGE) THEN
        DO k=1,N
          Sprod3d(iend+1,jend+1,k)=Sprod3d(iend,jend,k)
        ENDDO
      ENDIF
#  endif
# endif
!
! Exchange computational margines and/or periodic boundaries:
!--------- ------------- -------- ------ -------- -----------
!
# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      CALL exchange_r3d_tile (Istr,Iend,Jstr,Jend,Sprod3d)
# endif
#endif
!$acc end kernels
!
      return
      end
!
#ifndef CHILD_SPG
# undef UCLM
# undef VCLM
# define CHILD_SPG
# ifdef AGRIF
#  include "uv3dmix_R3d.F"
# endif
# undef CHILD_SPG
#endif  /* !CHILD_SPG */
!
